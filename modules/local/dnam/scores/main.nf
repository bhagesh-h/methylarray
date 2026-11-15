process DNAM_SCORES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-hdf5array_bioconductor-summarizedexperiment_procps-ng_r-base_r-data.table:f08a921e65e6bc1c' :
        'community.wave.seqera.io/library/bioconductor-hdf5array_bioconductor-summarizedexperiment_procps-ng_r-base_r-data.table:460ab90a75a0c805' }"

    input:
    tuple val(meta), path(se), path(retained), path(annotation)
    path weights, stageAs: 'weights/*'
    path sources

    output:
    tuple val(meta), path("${prefix}.tsv")             , emit: scores
    tuple val(meta), path("${prefix}_coverage.tsv")    , emit: coverage
    tuple val(meta), path("${prefix}_weights.tsv")     , emit: weights
    tuple val(meta), path("${prefix}_coverage_mqc.tsv"), emit: multiqc
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "dnam_scores"
    template 'dnam_scores.R'

    stub:
    prefix = task.ext.prefix ?: "dnam_scores"
    """
    printf 'sample\\tdnam_smoking\\n' > ${prefix}.tsv
    printf 'score\\tcpgs\\tpresent\\n' > ${prefix}_coverage.tsv
    printf 'score\\tcpg\\tweight\\tsd_beta\\tstandardised\\n' > ${prefix}_weights.tsv
    printf 'Score\\tCpGs\\n' > ${prefix}_coverage_mqc.tsv
    printf '"%s":\\n    SummarizedExperiment: stub\\n    data.table: stub\\n' "${task.process}" > versions.yml
    """
}
