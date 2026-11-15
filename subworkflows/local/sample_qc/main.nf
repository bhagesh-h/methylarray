//
// Sample QC: ewastools control metrics, sex and SNP genotypes (optional), then
// one summary with a PASS, FLAG, FAIL or EXCLUDED status for every array
//

include { EWASTOOLS_QC     } from '../../../modules/local/ewastools/qc'
include { SAMPLEQC_SUMMARY } from '../../../modules/local/sampleqc/summary'
include { sampleTable      } from '../utils_nfcore_methylarray_pipeline'

workflow SAMPLE_QC {

    take:
    ch_samplesheet      // channel: [ meta, idat_red, idat_grn ] for every array, including excluded ones
    ch_preprocessing_qc // channel: QC table of the preprocessing route
    ewastools_source    // string: ewastools source tarball (URL or path), or null to skip the ewastools QC

    main:
    def ch_arrays = ch_samplesheet.filter { meta, _idat_red, _idat_grn -> !meta.exclude }

    def ch_ewastools_qc = channel.value([])
    def ch_snp_conflicts = channel.value([])
    def ch_controls = channel.value([])
    if (ewastools_source) {
        def ch_idats = ch_arrays
            .map { _meta, idat_red, idat_grn -> [ idat_red, idat_grn ] }
            .collect()
            .map { files -> [ files ] }
        EWASTOOLS_QC(
            sampleTable(ch_arrays, true)
                .combine(ch_idats)
                .map { samples, idats -> [ [ id: 'ewastools' ], samples, idats ] },
            file(ewastools_source, checkIfExists: true)
        )
        ch_ewastools_qc = EWASTOOLS_QC.out.qc.map { _meta, qc -> qc }
        ch_snp_conflicts = EWASTOOLS_QC.out.conflicts.map { _meta, conflicts -> conflicts }
        ch_controls = EWASTOOLS_QC.out.controls.map { _meta, controls -> controls }
    }

    // Every samplesheet row, so that excluded arrays are listed too
    def ch_samples = ch_samplesheet
        .map { meta, _idat_red, _idat_grn ->
            [ meta.id, meta.reported_sex ?: '', meta.replicate_group ?: '', meta.exclude ? 'true' : 'false' ].join('\t')
        }
        .collectFile(name: 'samples.tsv', seed: 'sample\treported_sex\treplicate_group\texclude', newLine: true, sort: true)

    SAMPLEQC_SUMMARY(
        ch_samples
            .combine(ch_preprocessing_qc)
            .combine(ch_ewastools_qc.map { qc -> [ qc ] })
            .combine(ch_snp_conflicts.map { conflicts -> [ conflicts ] })
            .map { samples, preprocessing_qc, ewastools_qc, snp_conflicts ->
                [ [ id: 'sample_qc' ], samples, preprocessing_qc, ewastools_qc, snp_conflicts ]
            }
    )

    emit:
    qc       = SAMPLEQC_SUMMARY.out.qc       // channel: [ meta, sample_qc.tsv ]
    retained = SAMPLEQC_SUMMARY.out.retained // channel: [ meta, sample_qc_retained.txt ]
    multiqc  = SAMPLEQC_SUMMARY.out.multiqc  // channel: [ meta, [ MultiQC custom content ] ]
    controls = ch_controls                   // channel: control-probe intensities, or [] without ewastools
}
