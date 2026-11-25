//
// Probe annotation (Peters 2024 manifest and Zhou lab masks, or a saved
// annotation table) and probe filtering of the arrays that passed sample QC
//

include { PROBES_ANNOTATE } from '../../../modules/local/probes/annotate'
include { PROBES_FILTER   } from '../../../modules/local/probes/filter'

workflow PROBE_FILTER {

    take:
    ch_se             // channel: [ meta, SummarizedExperiment directory ] of the preprocessing route
    ch_retained       // channel: [ meta, list of the arrays that passed sample QC ]
    probe_annotation  // string: saved probe annotation table, or null to build it
    manifest_source   // string: Peters et al. 2024 EPIC v2 manifest .Rda (URL or path)
    probe_mask_source // string: Zhou lab EPIC v2 mask (URL or path)
    platform          // string: array platform
    illumina_flagged  // string: Illumina flagged-probe list (path), or null
    illumina_mapping  // string: Illumina mapping-inaccuracy list (path), or null
    exclude_probes    // string: list of probes or CpGs to remove (URL or path), or null

    main:
    def ch_annotation = channel.empty()
    if (probe_annotation) {
        ch_annotation = channel.value(file(probe_annotation, checkIfExists: true))
    } else {
        PROBES_ANNOTATE(channel.value([
            [ id: platform ],
            file(manifest_source, checkIfExists: true),
            file(probe_mask_source, checkIfExists: true),
            illumina_flagged ? file(illumina_flagged, checkIfExists: true) : [],
            illumina_mapping ? file(illumina_mapping, checkIfExists: true) : []
        ]))
        ch_annotation = PROBES_ANNOTATE.out.annotation.map { _meta, annotation -> annotation }
    }

    PROBES_FILTER(
        ch_se
            .combine(ch_retained.map { _meta, retained -> retained })
            .combine(ch_annotation)
            .map { meta, se, retained, annotation -> [ meta, se, retained, annotation, exclude_probes ? file(exclude_probes, checkIfExists: true) : [] ] }
    )

    emit:
    se             = PROBES_FILTER.out.se             // channel: [ meta, filtered SummarizedExperiment, collapsed to CpGs ]
    se_uncollapsed = PROBES_FILTER.out.se_uncollapsed // channel: [ meta, filtered SummarizedExperiment, one row per probe ]
    attrition      = PROBES_FILTER.out.attrition      // channel: [ meta, probe attrition table ]
    annotation     = ch_annotation                    // channel: probe annotation table
    multiqc        = PROBES_FILTER.out.multiqc        // channel: [ meta, MultiQC custom content ]
}
