//
// Exploratory analysis: control-probe PCs, the association of the main
// principal components with technical and biological variables, and the
// variance partition of the most variable CpGs
//

include { EXPLORATION_PCA    } from '../../../modules/local/exploration/pca'
include { VARIANCE_PARTITION } from '../../../modules/local/variance/partition'
include { sampleTable     } from '../utils_nfcore_methylarray_pipeline'

workflow EXPLORATION {

    take:
    ch_arrays         // channel: [ meta, idat_red, idat_grn ] of the processed arrays
    ch_se             // channel: [ meta, filtered SummarizedExperiment, one row per CpG ]
    ch_sample_qc      // channel: [ meta, sample QC table ]
    ch_controls       // channel: control-probe intensities, or [] without ewastools
    ch_cell_fractions // channel: cell fraction table, or [] without cell composition
    ch_phenotypes     // channel: phenotype table, or []
    skip_variance_partition // boolean: skip the variance partition

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

    def ch_variance_partition = channel.empty()
    def ch_variance_partition_mqc = channel.empty()
    if (!skip_variance_partition) {
        VARIANCE_PARTITION(
            ch_se
                .combine(sampleTable(ch_arrays, false))
                .combine(ch_sample_qc.map { _meta, qc -> qc })
                .combine(ch_cell_fractions.map { fractions -> [ fractions ] })
                .map { meta, se, samples, sample_qc, fractions -> [ meta, se, samples, sample_qc, fractions ] }
        )
        ch_variance_partition = VARIANCE_PARTITION.out.variance.mix(VARIANCE_PARTITION.out.summary).map { _meta, f -> f }
        ch_variance_partition_mqc = VARIANCE_PARTITION.out.multiqc.map { _meta, f -> f }
    }

    emit:
    ctrl_pcs     = EXPLORATION_PCA.out.ctrl_pcs                                                // channel: [ meta, control-probe PCs ]
    variance     = EXPLORATION_PCA.out.variance                                                // channel: [ meta, variance per component ]
    scores       = EXPLORATION_PCA.out.scores                                                  // channel: [ meta, PC scores ]
    associations = EXPLORATION_PCA.out.associations                                            // channel: [ meta, PC associations ]
    outliers     = EXPLORATION_PCA.out.outliers                                                // channel: [ meta, PCA outlier flags (report only) ]
    multiqc      = EXPLORATION_PCA.out.multiqc.map { _meta, files -> files }.flatten().mix(ch_variance_partition_mqc) // channel: MultiQC custom content
    variance_partition = ch_variance_partition                                                 // channel: variance partition table and summary (empty when skipped)
}
