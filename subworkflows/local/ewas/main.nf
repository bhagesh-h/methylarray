//
// Epigenome-wide association: one limma task per model of --models, then one
// summary table across the models
//

include { EWAS_LIMMA   } from '../../../modules/local/ewas/limma'
include { EWAS_SUMMARY } from '../../../modules/local/ewas/summary'
include { sampleTable  } from '../utils_nfcore_methylarray_pipeline'

workflow EWAS {

    take:
    ch_models         // channel: one map per validated model
    ch_arrays         // channel: [ meta, idat_red, idat_grn ] of the processed arrays
    ch_se             // channel: [ meta, filtered SummarizedExperiment, one row per CpG ]
    ch_sample_qc      // channel: [ meta, sample QC table ]
    ch_phenotypes     // channel: phenotype table, or []
    ch_cell_fractions // channel: cell fraction table, or []
    ch_ctrl_pcs       // channel: control-probe PCs table
    ch_ancestry_pcs   // channel: ancestry PCs table, or []
    ch_dnam_scores    // channel: DNA methylation scores table, or []
    ch_meqtl_cpgs     // channel: CpGs with a cis-meQTL, or []

    main:
    def ch_model_files = ch_models
        .collectFile { model -> [ "${model.name}.json", groovy.json.JsonOutput.toJson(model) ] }
        .map { json -> [ [ id: json.baseName ], json ] }

    EWAS_LIMMA(
        ch_model_files
            .combine(ch_se.map { _meta, se -> se })
            .combine(sampleTable(ch_arrays, false))
            .combine(ch_sample_qc.map { _meta, qc -> qc })
            .combine(ch_phenotypes.map { phenotypes -> [ phenotypes ] })
            .combine(ch_cell_fractions.map { fractions -> [ fractions ] })
            .combine(ch_ctrl_pcs)
            .combine(ch_ancestry_pcs.map { pcs -> [ pcs ] })
            .combine(ch_dnam_scores.map { scores -> [ scores ] })
            .map { meta, model, se, samples, sample_qc, phenotypes, fractions, ctrl_pcs, ancestry_pcs, scores ->
                [ meta, model, se, samples, sample_qc, phenotypes, fractions, ctrl_pcs, ancestry_pcs, scores ]
            }
    )

    EWAS_SUMMARY(
        EWAS_LIMMA.out.summary.map { _meta, summary -> summary }.collect(),
        EWAS_LIMMA.out.stats.map { _meta, stats -> stats }.collect(),
        ch_meqtl_cpgs
    )

    emit:
    stats    = EWAS_LIMMA.out.stats    // channel: [ meta, full statistics per model ]
    resolved = EWAS_LIMMA.out.resolved // channel: [ meta, resolved model ]
    top     = EWAS_LIMMA.out.top      // channel: [ meta, top CpGs per model ]
    design_check = EWAS_LIMMA.out.design_check // channel: [ meta, design checks per model ]
    characteristics = EWAS_LIMMA.out.characteristics // channel: [ meta, baseline characteristics per model ]
    missing  = EWAS_LIMMA.out.missing  // channel: [ meta, arrays per combination of missing model variables ]
    cells    = EWAS_LIMMA.out.cells    // channel: [ meta, cell-type-specific statistics per model ]
    summary = EWAS_SUMMARY.out.summary // channel: summary table across models
    multiqc = EWAS_SUMMARY.out.multiqc // channel: MultiQC table
    meqtl   = EWAS_SUMMARY.out.meqtl   // channel: hits with a cis-meQTL per model (header only without --meqtl)
    comparison = EWAS_SUMMARY.out.comparison // channel: overlap and agreement of every pair of model contrasts
    cells_summary = EWAS_SUMMARY.out.cells // channel: cell-type-specific hits per model and lineage
}
