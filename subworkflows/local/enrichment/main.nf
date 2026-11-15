//
// Gene set enrichment (missMethyl: GO, KEGG, MSigDB collections) of the CpGs
// significant in every tested association model
//

include { ENRICHMENT_MISSMETHYL } from '../../../modules/local/enrichment/missmethyl'

workflow ENRICHMENT {

    take:
    ch_stats      // channel: [ meta, EWAS statistics ] per model
    ch_resolved   // channel: [ meta, resolved model ] per model
    gene_set_gmts // string: comma-separated MSigDB GMT files (URLs or paths, Entrez IDs), or empty

    main:
    def gmts = gene_set_gmts ? gene_set_gmts.tokenize(',').collect { gmt -> file(gmt.trim(), checkIfExists: true) } : []
    ENRICHMENT_MISSMETHYL(ch_stats.join(ch_resolved), gmts)

    emit:
    enrichment = ENRICHMENT_MISSMETHYL.out.enrichment // channel: [ meta, enrichment table ]
    summary    = ENRICHMENT_MISSMETHYL.out.summary    // channel: [ meta, per-model summary ]
}
