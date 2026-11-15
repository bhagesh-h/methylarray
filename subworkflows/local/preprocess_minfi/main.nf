//
// minfi-based preprocessing (routes minfi, enmix and combined): all arrays are
// read together and written to one HDF5-backed SummarizedExperiment
//

include { MINFI_PREPROCESS } from '../../../modules/local/minfi/preprocess'
include { sampleTable      } from '../utils_nfcore_methylarray_pipeline'

workflow PREPROCESS_MINFI {

    take:
    ch_arrays    // channel: [ meta, idat_red, idat_grn ]
    method       // string: minfi, enmix or combined
    ch_sesame_se // channel: SeSAMe SummarizedExperiment directory (combined route), or [] otherwise

    main:
    def ch_samples = sampleTable(ch_arrays, true)
    def ch_idats = ch_arrays
        .map { _meta, idat_red, idat_grn -> [ idat_red, idat_grn ] }
        .collect()
        .map { files -> [ files ] }

    MINFI_PREPROCESS(
        ch_samples
            .combine(ch_idats)
            .combine(ch_sesame_se.map { se -> [ se ] })
            .map { samples, idats, sesame_se -> [ [ id: method ], samples, idats, sesame_se ] }
    )

    emit:
    se      = MINFI_PREPROCESS.out.se      // channel: [ meta, se directory ]
    qc      = MINFI_PREPROCESS.out.qc      // channel: [ meta, sample QC table ]
    multiqc = MINFI_PREPROCESS.out.multiqc // channel: [ meta, MultiQC general statistics ]
}
