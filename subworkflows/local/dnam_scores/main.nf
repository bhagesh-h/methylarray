//
// DNA methylation scores: weighted sums over published CpG sets (assets/dnam_scores)
//

include { DNAM_SCORES } from '../../../modules/local/dnam/scores'

workflow DNAM_SCORE {

    take:
    ch_se         // channel: [ meta, SummarizedExperiment directory ] of the preprocessing route
    ch_retained   // channel: [ meta, list of the arrays that passed sample QC ]
    ch_annotation // channel: probe annotation table

    main:
    def weights = file("${projectDir}/assets/dnam_scores/*.tsv").findAll { f -> f.name != 'sources.tsv' }
    DNAM_SCORES(
        ch_se
            .combine(ch_retained.map { _meta, retained -> retained })
            .combine(ch_annotation),
        weights,
        file("${projectDir}/assets/dnam_scores/sources.tsv", checkIfExists: true)
    )

    emit:
    scores   = DNAM_SCORES.out.scores   // channel: [ meta, scores per array ]
    coverage = DNAM_SCORES.out.coverage // channel: [ meta, CpG coverage per score ]
    multiqc  = DNAM_SCORES.out.multiqc  // channel: [ meta, MultiQC table ]
}
