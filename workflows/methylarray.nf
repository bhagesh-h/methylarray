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
include { DNAM_SCORE             } from '../subworkflows/local/dnam_scores'
include { EXPLORATION            } from '../subworkflows/local/exploration'
include { EWAS                   } from '../subworkflows/local/ewas'
include { DMR                    } from '../subworkflows/local/dmr'
include { ENRICHMENT             } from '../subworkflows/local/enrichment'
include { EXPORT_GEO             } from '../modules/local/export/geo'
include { ANCESTRY_PCS           } from '../modules/local/ancestry/pcs'
include { MEQTL_PREPARE          } from '../modules/local/meqtl/prepare'
include { BIOMARKERS_ASSOCIATE   } from '../modules/local/biomarkers/associate'
include { REPORT_QUARTO          } from '../modules/local/report/quarto'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_methylarray_pipeline'
include { sampleTable            } from '../subworkflows/local/utils_nfcore_methylarray_pipeline'

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
    def ch_report_files = channel.empty()

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
    ch_report_files = ch_report_files.mix(SAMPLE_QC.out.qc.map { _meta, qc -> qc }, SAMPLE_QC.out.snp)
    ch_multiqc_files = ch_multiqc_files.mix(SAMPLE_QC.out.snp_mqc)

    //
    // SUBWORKFLOW: Probe annotation and filtering
    //
    PROBE_FILTER(ch_route_se, SAMPLE_QC.out.retained, params.probe_annotation, params.probe_manifest_source, params.probe_mask_source, params.platform,
        params.illumina_flagged_probes, params.illumina_mapping_inaccuracy)
    ch_multiqc_files = ch_multiqc_files.mix(PROBE_FILTER.out.multiqc.map { _meta, table -> table })
    ch_report_files = ch_report_files.mix(PROBE_FILTER.out.attrition.map { _meta, table -> table })

    //
    // SUBWORKFLOW: Blood cell-type composition
    //
    def ch_cell_fractions = channel.value([])
    if (params.tissue == 'blood') {
        CELL_COMPOSITION(PROBE_FILTER.out.se, PROBE_FILTER.out.annotation, params.skip_idol ? null : params.idol_reference)
        ch_cell_fractions = CELL_COMPOSITION.out.fractions.map { _meta, fractions -> fractions }
        ch_multiqc_files = ch_multiqc_files.mix(CELL_COMPOSITION.out.multiqc.map { _meta, table -> table })
        ch_report_files = ch_report_files.mix(CELL_COMPOSITION.out.fractions.map { _meta, f -> f }, CELL_COMPOSITION.out.coverage.map { _meta, f -> f },
            CELL_COMPOSITION.out.idol.map { _meta, f -> f })
    }

    //
    // MODULE: CpGs with a cis-meQTL, from a user-supplied meQTL table (--meqtl)
    //
    def ch_meqtl_cpgs = channel.value([])
    if (params.meqtl) {
        MEQTL_PREPARE(
            PROBE_FILTER.out.annotation.map { annotation -> [ [ id: 'meqtl' ], file(params.meqtl, checkIfExists: true), annotation ] }
        )
        ch_meqtl_cpgs = MEQTL_PREPARE.out.cpgs.map { _meta, f -> f }
        ch_report_files = ch_report_files.mix(MEQTL_PREPARE.out.summary.map { _meta, f -> f })
        ch_multiqc_files = ch_multiqc_files.mix(MEQTL_PREPARE.out.multiqc.map { _meta, f -> f })
    }

    //
    // SUBWORKFLOW: Epigenetic clocks
    //
    def ch_clocks = channel.value([])
    def ch_clocks_robust = channel.value([])
    if (!params.skip_clocks) {
        CLOCKS(ch_arrays, ch_route_se, SAMPLE_QC.out.retained, PROBE_FILTER.out.annotation, ch_cell_fractions, ch_meqtl_cpgs,
            params.dnamethyage_source, params.dunedinpace_source)
        ch_clocks = CLOCKS.out.clocks.map { _meta, f -> f }
        if (params.meqtl) {
            ch_report_files = ch_report_files.mix(CLOCKS.out.robust.map { _meta, f -> f })
            // optional output: an empty channel would hold back every step that combines it
            ch_clocks_robust = CLOCKS.out.robust.map { _meta, f -> f }.ifEmpty([])
        }
        ch_multiqc_files = ch_multiqc_files.mix(CLOCKS.out.multiqc.map { _meta, table -> table })
        ch_report_files = ch_report_files.mix(CLOCKS.out.clocks.map { _meta, f -> f }, CLOCKS.out.coverage.map { _meta, f -> f },
            CLOCKS.out.weights.map { _meta, f -> f })
    }

    //
    // SUBWORKFLOW: DNA methylation scores (smoking, alcohol, BMI, IL-6)
    //
    def ch_dnam_scores = channel.value([])
    if (!params.skip_dnam_scores) {
        DNAM_SCORE(ch_route_se, SAMPLE_QC.out.retained, PROBE_FILTER.out.annotation)
        ch_dnam_scores = DNAM_SCORE.out.scores.map { _meta, scores -> scores }
        ch_multiqc_files = ch_multiqc_files.mix(DNAM_SCORE.out.multiqc.map { _meta, table -> table })
        ch_report_files = ch_report_files.mix(DNAM_SCORE.out.scores.map { _meta, f -> f }, DNAM_SCORE.out.coverage.map { _meta, f -> f },
            DNAM_SCORE.out.weights.map { _meta, f -> f })
    }

    //
    // SUBWORKFLOW: Exploratory analysis (control-probe PCs, PCA associations)
    //
    EXPLORATION(ch_arrays, PROBE_FILTER.out.se, SAMPLE_QC.out.qc, SAMPLE_QC.out.controls, ch_cell_fractions, ch_phenotypes, params.skip_variance_partition)
    ch_multiqc_files = ch_multiqc_files.mix(EXPLORATION.out.multiqc)
    ch_report_files = ch_report_files.mix(EXPLORATION.out.variance.map { _meta, f -> f }, EXPLORATION.out.associations.map { _meta, f -> f },
        EXPLORATION.out.scores.map { _meta, f -> f }, EXPLORATION.out.outliers.map { _meta, f -> f }, EXPLORATION.out.variance_partition)

    //
    // MODULE: Genetic ancestry PCs from the SNP CpGs and rs probes (EpiAnceR+ approach)
    //
    def ch_ancestry_pcs = channel.value([])
    if (!params.skip_ancestry) {
        ANCESTRY_PCS(
            ch_route_se
                .combine(PROBE_FILTER.out.annotation)
                .combine(sampleTable(ch_arrays, false))
                .combine(SAMPLE_QC.out.qc.map { _meta, qc -> qc })
                .combine(EXPLORATION.out.ctrl_pcs.map { _meta, pcs -> pcs })
                .combine(ch_cell_fractions.map { f -> [ f ] })
                .combine(ch_phenotypes.map { f -> [ f ] })
                .map { meta, se, annotation, samples, qc, ctrl_pcs, cells, phenotypes -> [ meta, se, annotation, samples, qc, ctrl_pcs, cells, phenotypes ] }
        )
        ch_ancestry_pcs = ANCESTRY_PCS.out.pcs.map { _meta, f -> f }
        ch_multiqc_files = ch_multiqc_files.mix(ANCESTRY_PCS.out.multiqc.map { _meta, files -> files }.flatten())
        ch_report_files = ch_report_files.mix(ANCESTRY_PCS.out.pcs.map { _meta, f -> f }, ANCESTRY_PCS.out.variance.map { _meta, f -> f },
            ANCESTRY_PCS.out.summary.map { _meta, f -> f }, ANCESTRY_PCS.out.associations.map { _meta, f -> f })
    }

    //
    // MODULE: clocks, scores and cell fractions against each other, the batches and the phenotypes
    //
    if (!params.skip_clocks || !params.skip_dnam_scores) {
        BIOMARKERS_ASSOCIATE(
            sampleTable(ch_arrays, false)
                .combine(SAMPLE_QC.out.qc.map { _meta, qc -> qc })
                .combine(ch_clocks.map { f -> [ f ] })
                .combine(ch_cell_fractions.map { f -> [ f ] })
                .combine(ch_dnam_scores.map { f -> [ f ] })
                .combine(ch_phenotypes.map { f -> [ f ] })
                .combine(ch_ancestry_pcs.map { f -> [ f ] })
                .combine(ch_clocks_robust.map { f -> [ f ] })
                .map { samples, qc, clocks, cells, scores, phenotypes, ancestry, robust ->
                    [ [ id: 'biomarkers' ], samples, qc, clocks, cells, scores, phenotypes, ancestry, robust ]
                }
        )
        ch_multiqc_files = ch_multiqc_files.mix(BIOMARKERS_ASSOCIATE.out.multiqc.map { _meta, files -> files }.flatten())
        def ch_biomarker_tables = BIOMARKERS_ASSOCIATE.out.agreement
            .mix(BIOMARKERS_ASSOCIATE.out.cells_measures, BIOMARKERS_ASSOCIATE.out.cells_phenotypes, BIOMARKERS_ASSOCIATE.out.batch,
                BIOMARKERS_ASSOCIATE.out.phenotypes, BIOMARKERS_ASSOCIATE.out.quartiles, BIOMARKERS_ASSOCIATE.out.contamination,
                BIOMARKERS_ASSOCIATE.out.ancestry, BIOMARKERS_ASSOCIATE.out.robust)
            .map { _meta, f -> f }
        ch_report_files = ch_report_files.mix(ch_biomarker_tables)
    }

    //
    // SUBWORKFLOW: Epigenome-wide association, one model per entry of --models
    //
    EWAS(ch_models, ch_arrays, PROBE_FILTER.out.se, SAMPLE_QC.out.qc, ch_phenotypes, ch_cell_fractions,
        EXPLORATION.out.ctrl_pcs.map { _meta, ctrl_pcs -> ctrl_pcs }, ch_ancestry_pcs, ch_dnam_scores, ch_meqtl_cpgs)
    ch_multiqc_files = ch_multiqc_files.mix(EWAS.out.multiqc)
    ch_report_files = ch_report_files.mix(EWAS.out.summary, EWAS.out.stats.map { _meta, f -> f }, EWAS.out.top.map { _meta, f -> f },
        EWAS.out.design_check.map { _meta, f -> f }, EWAS.out.cells.map { _meta, f -> f }, EWAS.out.cells_summary, EWAS.out.meqtl,
        EWAS.out.characteristics.map { _meta, f -> f }, EWAS.out.comparison)

    //
    // SUBWORKFLOW: Differentially methylated regions of every model
    //
    if (!params.skip_dmr) {
        DMR(EWAS.out.stats, EWAS.out.resolved, PROBE_FILTER.out.se, params.dmrff_source)
        ch_report_files = ch_report_files.mix(DMR.out.consensus.map { _meta, f -> f }, DMR.out.summary.map { _meta, f -> f })
        ch_multiqc_files = ch_multiqc_files.mix(DMR.out.multiqc.map { _meta, f -> f })
    }

    //
    // SUBWORKFLOW: Gene set enrichment of every model
    //
    if (!params.skip_enrichment) {
        ENRICHMENT(EWAS.out.stats, EWAS.out.resolved, params.gene_set_gmts, PROBE_FILTER.out.annotation, params.skip_kycg, params.kycg_cache, params.platform)
        ch_report_files = ch_report_files.mix(ENRICHMENT.out.enrichment.map { _meta, f -> f }, ENRICHMENT.out.kycg.map { _meta, f -> f })
        ch_multiqc_files = ch_multiqc_files.mix(ENRICHMENT.out.multiqc.map { _meta, f -> f }, ENRICHMENT.out.kycg_mqc.map { _meta, f -> f })
    }

    //
    // MODULE: GEO submission files (opt-in)
    //
    if (params.export_geo) {
        EXPORT_GEO(
            ch_route_se
                .combine(sampleTable(ch_arrays, true))
                .combine(ch_arrays.map { _meta, idat_red, idat_grn -> [ idat_red, idat_grn ] }.collect().map { idats -> [ idats ] })
                .combine(SAMPLE_QC.out.retained.map { _meta, retained -> retained })
                .map { meta, se, samples, idats, retained -> [ meta, samples, idats, se, retained ] }
        )
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
    //
    // MODULE: HTML report (Quarto) from the result tables, versions, methods and parameters
    //
    if (!params.skip_report) {
        def ch_params_json = channel.value(groovy.json.JsonOutput.toJson(params.collectEntries { name, value -> [ name, value?.toString() ] })).collectFile(name: 'params.json')
        REPORT_QUARTO(
            file("${projectDir}/assets/report/methylarray_report.qmd", checkIfExists: true),
            ch_report_files
                .mix(ch_collated_versions, ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'), ch_params_json)
                .collect()
        )
    }

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
