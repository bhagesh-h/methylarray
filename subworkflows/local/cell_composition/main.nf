//
// Blood cell-type composition with EpiDISH on the filtered CpG-level data
//

include { EPIDISH_RPC } from '../../../modules/local/epidish/rpc'

workflow CELL_COMPOSITION {

    take:
    ch_se         // channel: [ meta, filtered SummarizedExperiment, one row per CpG ]
    ch_annotation // channel: probe annotation table

    main:
    EPIDISH_RPC(ch_se.combine(ch_annotation))

    emit:
    fractions = EPIDISH_RPC.out.fractions // channel: [ meta, cell fractions per array ]
    coverage  = EPIDISH_RPC.out.coverage  // channel: [ meta, reference CpG coverage ]
    multiqc   = EPIDISH_RPC.out.multiqc   // channel: [ meta, MultiQC custom content ]
}
