/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PREPROCESS_SESAME      } from '../subworkflows/local/preprocess_sesame'
include { PREPROCESS_MINFI       } from '../subworkflows/local/preprocess_minfi'
include { SAMPLE_QC              } from '../subworkflows/local/sample_qc'
include { PROBE_FILTER           } from '../subworkflows/local/probe_filter'
include { CELL_COMPOSITION       } from '../subworkflows/local/cell_composition'
include { CLOCKS                 } from '../subworkflows/local/clocks'
include { EXPLORATION            } from '../subworkflows/local/exploration'
include { EWAS                   } from '../subworkflows/local/ewas'
include { DMR                    } from '../subworkflows/local/dmr'
include { ENRICHMENT             } from '../subworkflows/local/enrichment'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_methylarray_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METHYLARRAY {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_phenotypes  // channel: phenotype table, or []
    ch_models      // channel: one map per validated model (empty without --models)
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // Arrays excluded in the samplesheet are left out of every step
    //
    def ch_arrays = ch_samplesheet.filter { meta, _idat_red, _idat_grn -> !meta.exclude }

    //
    // SUBWORKFLOW: Preprocessing
    // sesame: SeSAMe only. minfi, enmix: minfi import. combined: the minfi
    // import with ENmix correction, masked with the SeSAMe result.
    //
    def ch_sesame_se = channel.value([])
    def ch_preprocessing_qc = channel.empty()
    def ch_route_se = channel.empty()
    if (params.preprocessing in ['sesame', 'combined']) {
        PREPROCESS_SESAME(ch_arrays, params.sesame_cache, params.platform)
        ch_sesame_se = PREPROCESS_SESAME.out.se.map { _meta, se -> se }
        if (params.preprocessing == 'sesame') {
            ch_preprocessing_qc = PREPROCESS_SESAME.out.qc.map { _meta, qc -> qc }
            ch_route_se = PREPROCESS_SESAME.out.se
            ch_multiqc_files = ch_multiqc_files.mix(PREPROCESS_SESAME.out.multiqc.map { _meta, table -> table })
        }
    }
    if (params.preprocessing in ['minfi', 'enmix', 'combined']) {
        PREPROCESS_MINFI(ch_arrays, params.preprocessing, ch_sesame_se)
        ch_preprocessing_qc = PREPROCESS_MINFI.out.qc.map { _meta, qc -> qc }
        ch_route_se = PREPROCESS_MINFI.out.se
        ch_multiqc_files = ch_multiqc_files.mix(PREPROCESS_MINFI.out.multiqc.map { _meta, table -> table })
    }

    //
    // SUBWORKFLOW: Sample QC
    //
    SAMPLE_QC(ch_samplesheet, ch_preprocessing_qc, params.skip_ewastools ? null : params.ewastools_source)
    ch_multiqc_files = ch_multiqc_files.mix(SAMPLE_QC.out.multiqc.map { _meta, files -> files }.flatten())

    //
    // SUBWORKFLOW: Probe annotation and filtering
    //
    PROBE_FILTER(ch_route_se, SAMPLE_QC.out.retained, params.probe_annotation, params.probe_manifest_source, params.probe_mask_source, params.platform)
    ch_multiqc_files = ch_multiqc_files.mix(PROBE_FILTER.out.multiqc.map { _meta, table -> table })

    //
    // SUBWORKFLOW: Blood cell-type composition
    //
    def ch_cell_fractions = channel.value([])
    if (params.tissue == 'blood') {
        CELL_COMPOSITION(PROBE_FILTER.out.se, PROBE_FILTER.out.annotation)
        ch_cell_fractions = CELL_COMPOSITION.out.fractions.map { _meta, fractions -> fractions }
        ch_multiqc_files = ch_multiqc_files.mix(CELL_COMPOSITION.out.multiqc.map { _meta, table -> table })
    }

    //
    // SUBWORKFLOW: Epigenetic clocks
    //
    if (!params.skip_clocks) {
        CLOCKS(ch_arrays, ch_route_se, SAMPLE_QC.out.retained, PROBE_FILTER.out.annotation, ch_cell_fractions,
            params.dnamethyage_source, params.dunedinpace_source)
        ch_multiqc_files = ch_multiqc_files.mix(CLOCKS.out.multiqc.map { _meta, table -> table })
    }

    //
    // SUBWORKFLOW: Exploratory analysis (control-probe PCs, PCA associations)
    //
    EXPLORATION(ch_arrays, PROBE_FILTER.out.se, SAMPLE_QC.out.qc, SAMPLE_QC.out.controls, ch_cell_fractions, ch_phenotypes)
    ch_multiqc_files = ch_multiqc_files.mix(EXPLORATION.out.multiqc)

    //
    // SUBWORKFLOW: Epigenome-wide association, one model per entry of --models
    //
    EWAS(ch_models, ch_arrays, PROBE_FILTER.out.se, SAMPLE_QC.out.qc, ch_phenotypes, ch_cell_fractions,
        EXPLORATION.out.ctrl_pcs.map { _meta, ctrl_pcs -> ctrl_pcs })
    ch_multiqc_files = ch_multiqc_files.mix(EWAS.out.multiqc)

    //
    // SUBWORKFLOW: Differentially methylated regions of every model
    //
    if (!params.skip_dmr) {
        DMR(EWAS.out.stats, EWAS.out.resolved, PROBE_FILTER.out.se, params.dmrff_source)
    }

    //
    // SUBWORKFLOW: Gene set enrichment of every model
    //
    if (!params.skip_enrichment) {
        ENRICHMENT(EWAS.out.stats, EWAS.out.resolved, params.gene_set_gmts)
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'methylarray_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'methylarray'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
