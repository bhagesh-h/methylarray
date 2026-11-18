//
// Gene set enrichment (missMethyl: GO, KEGG, MSigDB collections) and genomic
// context enrichment (KnowYourCG) of the CpGs significant in every tested
// association model
//

include { ENRICHMENT_MISSMETHYL } from '../../../modules/local/enrichment/missmethyl'
include { KYCG_CACHE            } from '../../../modules/local/kycg/cache'
include { ENRICHMENT_KYCG       } from '../../../modules/local/enrichment/kycg'

workflow ENRICHMENT {

    take:
    ch_stats      // channel: [ meta, EWAS statistics ] per model
    ch_resolved   // channel: [ meta, resolved model ] per model
    gene_set_gmts // string: comma-separated MSigDB GMT files (URLs or paths, Entrez IDs), or empty
    ch_annotation // channel: probe annotation table
    skip_kycg     // boolean: skip the KnowYourCG enrichment
    kycg_cache    // string: directory from an earlier KYCG_CACHE run, or null
    platform      // string: array platform

    main:
    def gmts = gene_set_gmts ? gene_set_gmts.tokenize(',').collect { gmt -> file(gmt.trim(), checkIfExists: true) } : []
    ENRICHMENT_MISSMETHYL(ch_stats.join(ch_resolved), gmts)

    def ch_kycg = channel.empty()
    def ch_kycg_multiqc = channel.empty()
    if (!skip_kycg) {
        def ch_cache
        if (kycg_cache) {
            ch_cache = channel.value(file(kycg_cache, checkIfExists: true))
        } else {
            // Only download the knowledgebases when a model was tested
            KYCG_CACHE(ch_stats.first().map { _meta, _stats -> platform })
            ch_cache = KYCG_CACHE.out.cache
        }
        ENRICHMENT_KYCG(ch_stats.join(ch_resolved), ch_cache, ch_annotation)
        ch_kycg = ENRICHMENT_KYCG.out.enrichment
        ch_kycg_multiqc = ENRICHMENT_KYCG.out.multiqc
    }

    emit:
    enrichment = ENRICHMENT_MISSMETHYL.out.enrichment // channel: [ meta, enrichment table ]
    summary    = ENRICHMENT_MISSMETHYL.out.summary    // channel: [ meta, per-model summary ]
    multiqc    = ENRICHMENT_MISSMETHYL.out.multiqc    // channel: [ meta, MultiQC table of the top sets ]
    kycg       = ch_kycg                              // channel: [ meta, KnowYourCG enrichment table ]
    kycg_mqc   = ch_kycg_multiqc                      // channel: [ meta, MultiQC table of the top KnowYourCG sets ]
}
