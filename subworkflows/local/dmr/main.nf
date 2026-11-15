//
// Differentially methylated regions (DMRcate, dmrff and their consensus) for
// every tested association model
//

include { DMR_CALL } from '../../../modules/local/dmr/call'

workflow DMR {

    take:
    ch_stats      // channel: [ meta, EWAS statistics ] per model
    ch_resolved   // channel: [ meta, resolved model ] per model
    ch_se         // channel: [ meta, filtered SummarizedExperiment, one row per CpG ]
    dmrff_source  // string: dmrff source tarball (URL or path)

    main:
    DMR_CALL(
        ch_stats
            .join(ch_resolved)
            .combine(ch_se.map { _meta, se -> se }),
        file(dmrff_source, checkIfExists: true)
    )

    emit:
    dmrcate   = DMR_CALL.out.dmrcate   // channel: [ meta, DMRcate regions ]
    dmrff     = DMR_CALL.out.dmrff     // channel: [ meta, dmrff regions ]
    consensus = DMR_CALL.out.consensus // channel: [ meta, consensus regions ]
    summary   = DMR_CALL.out.summary   // channel: [ meta, per-model summary ]
}
