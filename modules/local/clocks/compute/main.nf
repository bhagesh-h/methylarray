process CLOCKS_COMPUTE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-hdf5array_bioconductor-impute_bioconductor-preprocesscore_bioconductor-summarizedexperiment_pruned:4ca057f283ed6e32' :
        'community.wave.seqera.io/library/bioconductor-hdf5array_bioconductor-impute_bioconductor-preprocesscore_bioconductor-summarizedexperiment_pruned:07d497e1c23ca6c1' }"

    input:
    tuple val(meta), path(se), path(retained), path(annotation), path(samples), path(cell_fractions), path(meqtl_cpgs)
    path dnamethyage_source, stageAs: 'dnaMethyAge.tar.gz'
    path dunedinpace_source, stageAs: 'DunedinPACE.tar.gz'

    output:
    tuple val(meta), path("${prefix}.tsv")             , emit: clocks
    tuple val(meta), path("${prefix}_coverage.tsv")    , emit: coverage
    tuple val(meta), path("${prefix}_weights.tsv")     , emit: weights
    tuple val(meta), path("${prefix}_robust.tsv")      , emit: robust, optional: true
    tuple val(meta), path("${prefix}_*_mqc.tsv")       , emit: multiqc
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "clocks"
    template 'clocks_compute.R'

    stub:
    prefix = task.ext.prefix ?: "clocks"
    """
    printf 'sample\\tage\\n' > ${prefix}.tsv
    printf 'clock\\tcpgs\\tpresent\\n' > ${prefix}_coverage.tsv
    printf 'clock\\tcpg\\tgene\\tweight\\tsd_beta\\tstandardised\\tsnp\\tmeqtl\\n' > ${prefix}_weights.tsv
    printf 'Clock\\tCpGs\\n' > ${prefix}_coverage_mqc.tsv
    ${meqtl_cpgs ? "printf 'sample\\tage\\n' > ${prefix}_robust.tsv" : ''}
    printf '"%s":\\n    dnaMethyAge: stub\\n    DunedinPACE: stub\\n    preprocessCore: stub\\n    impute: stub\\n    RPMM: stub\\n' "${task.process}" > versions.yml
    """
}
