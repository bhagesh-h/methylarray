//
// Blood cell-type composition with EpiDISH on the filtered CpG-level data
//

include { EPIDISH_RPC } from '../../../modules/local/epidish/rpc'

workflow CELL_COMPOSITION {

    take:
    ch_se         // channel: [ meta, filtered SummarizedExperiment, one row per CpG ]
    ch_annotation // channel: probe annotation table
    idol_reference // string: IDOL reference table (URL or path), or null to skip the IDOL cross-check

    main:
    EPIDISH_RPC(
        ch_se
            .combine(ch_annotation)
            .map { meta, se, annotation -> [ meta, se, annotation, idol_reference ? file(idol_reference, checkIfExists: true) : [] ] }
    )

    emit:
    fractions = EPIDISH_RPC.out.fractions // channel: [ meta, cell fractions per array ]
    coverage  = EPIDISH_RPC.out.coverage  // channel: [ meta, reference CpG coverage ]
    idol      = EPIDISH_RPC.out.idol_agreement // channel: [ meta, agreement of the IDOL fractions with the lineages ]
    multiqc   = EPIDISH_RPC.out.multiqc   // channel: [ meta, MultiQC custom content ]
}
