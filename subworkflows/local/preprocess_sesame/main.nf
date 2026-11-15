//
// SeSAMe preprocessing: reference cache, per-array preprocessing, and one
// HDF5-backed SummarizedExperiment for all arrays
//

include { SESAME_CACHE      } from '../../../modules/local/sesame/cache'
include { SESAME_PREPROCESS } from '../../../modules/local/sesame/preprocess'
include { SESAME_BUILDSE    } from '../../../modules/local/sesame/buildse'
include { sampleTable       } from '../utils_nfcore_methylarray_pipeline'

workflow PREPROCESS_SESAME {

    take:
    ch_arrays    // channel: [ meta, idat_red, idat_grn ]
    sesame_cache // string: directory from an earlier SESAME_CACHE run, or null
    platform     // string: array platform, e.g. EPICv2

    main:
    def ch_cache = channel.empty()
    if (sesame_cache) {
        ch_cache = channel.value(file(sesame_cache, checkIfExists: true))
    } else {
        SESAME_CACHE(platform)
        ch_cache = SESAME_CACHE.out.cache
    }

    SESAME_PREPROCESS(ch_arrays, ch_cache)

    def ch_samples = sampleTable(ch_arrays, false)

    def ch_rds = SESAME_PREPROCESS.out.rds.map { _meta, rds -> rds }.collect().map { files -> [ files ] }
    def ch_qc = SESAME_PREPROCESS.out.qc.map { _meta, qc -> qc }.collect().map { files -> [ files ] }

    SESAME_BUILDSE(
        ch_samples
            .combine(ch_rds)
            .combine(ch_qc)
            .map { samples, rds, qc -> [ [ id: 'sesame' ], samples, rds, qc ] }
    )

    emit:
    se        = SESAME_BUILDSE.out.se           // channel: [ meta, se directory ]
    qc        = SESAME_BUILDSE.out.qc           // channel: [ meta, sample QC table ]
    multiqc   = SESAME_BUILDSE.out.multiqc      // channel: [ meta, MultiQC general statistics ]
    controls  = SESAME_PREPROCESS.out.controls  // channel: [ meta, control-probe intensities ]
    cache     = ch_cache                        // channel: sesame cache directory
}
