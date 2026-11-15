//
// Exploratory analysis: control-probe PCs and the association of the main
// principal components with technical and biological variables
//

include { EXPLORATION_PCA } from '../../../modules/local/exploration/pca'
include { sampleTable     } from '../utils_nfcore_methylarray_pipeline'

workflow EXPLORATION {

    take:
    ch_arrays         // channel: [ meta, idat_red, idat_grn ] of the processed arrays
    ch_se             // channel: [ meta, filtered SummarizedExperiment, one row per CpG ]
    ch_sample_qc      // channel: [ meta, sample QC table ]
    ch_controls       // channel: control-probe intensities, or [] without ewastools
    ch_cell_fractions // channel: cell fraction table, or [] without cell composition
    ch_phenotypes     // channel: phenotype table, or []

    main:
    EXPLORATION_PCA(
        ch_se
            .combine(sampleTable(ch_arrays, false))
            .combine(ch_sample_qc.map { _meta, qc -> qc })
            .combine(ch_controls.map { controls -> [ controls ] })
            .combine(ch_cell_fractions.map { fractions -> [ fractions ] })
            .combine(ch_phenotypes.map { phenotypes -> [ phenotypes ] })
            .map { meta, se, samples, sample_qc, controls, fractions, phenotypes ->
                [ meta, se, samples, sample_qc, controls, fractions, phenotypes ]
            }
    )

    emit:
    ctrl_pcs     = EXPLORATION_PCA.out.ctrl_pcs                                                // channel: [ meta, control-probe PCs ]
    variance     = EXPLORATION_PCA.out.variance                                                // channel: [ meta, variance per component ]
    scores       = EXPLORATION_PCA.out.scores                                                  // channel: [ meta, PC scores ]
    associations = EXPLORATION_PCA.out.associations                                            // channel: [ meta, PC associations ]
    multiqc      = EXPLORATION_PCA.out.multiqc.map { _meta, files -> files }.flatten()         // channel: MultiQC custom content
}
