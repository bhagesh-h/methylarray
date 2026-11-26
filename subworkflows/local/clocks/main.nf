//
// Epigenetic clocks (dnaMethyAge, DunedinPACE) on the probe-filtered betas of the
// arrays that passed sample QC, as in Qi et al. 2025 and the studies that follow it
//

include { CLOCKS_COMPUTE } from '../../../modules/local/clocks/compute'
include { sampleTable    } from '../utils_nfcore_methylarray_pipeline'

workflow CLOCKS {

    take:
    ch_arrays          // channel: [ meta, idat_red, idat_grn ] of the processed arrays
    ch_se              // channel: [ meta, SummarizedExperiment directory ] after the probe filter, collapsed to CpGs
    ch_retained        // channel: [ meta, list of the arrays that passed sample QC ]
    ch_annotation      // channel: probe annotation table
    ch_cell_fractions  // channel: cell fraction table, or [] without cell composition
    ch_meqtl_cpgs      // channel: CpGs with a cis-meQTL, or [] without --meqtl
    dnamethyage_source // string: dnaMethyAge source tarball (URL or path)
    dunedinpace_source // string: DunedinPACE source tarball (URL or path)

    main:
    CLOCKS_COMPUTE(
        ch_se
            .combine(ch_retained.map { _meta, retained -> retained })
            .combine(ch_annotation)
            .combine(sampleTable(ch_arrays, false))
            .combine(ch_cell_fractions.map { fractions -> [ fractions ] })
            .combine(ch_meqtl_cpgs.map { cpgs -> [ cpgs ] })
            .map { meta, se, retained, annotation, samples, fractions, meqtl -> [ meta, se, retained, annotation, samples, fractions, meqtl ] },
        file(dnamethyage_source, checkIfExists: true),
        file(dunedinpace_source, checkIfExists: true)
    )

    emit:
    clocks   = CLOCKS_COMPUTE.out.clocks   // channel: [ meta, clock values per array ]
    coverage = CLOCKS_COMPUTE.out.coverage // channel: [ meta, CpG coverage per clock ]
    weights  = CLOCKS_COMPUTE.out.weights  // channel: [ meta, standardised CpG weights per clock ]
    robust   = CLOCKS_COMPUTE.out.robust   // channel: [ meta, genotype-robust clocks ] (only with a meQTL list)
    multiqc  = CLOCKS_COMPUTE.out.multiqc  // channel: [ meta, MultiQC custom content ]
}
