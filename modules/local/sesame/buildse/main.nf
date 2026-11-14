process SESAME_BUILDSE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/bhagesh/methylarray-sesame:1.28.1'

    input:
    tuple val(meta), path(samples), path(rds, stageAs: 'rds/*'), path(qc, stageAs: 'qc/*')

    output:
    tuple val(meta), path("${prefix}_se")          , emit: se
    tuple val(meta), path("${prefix}_qc.tsv")      , emit: qc
    tuple val(meta), path("${prefix}_qc_mqc.tsv")  , emit: multiqc
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'sesame_buildse.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_se
    touch ${prefix}_se/assays.h5 ${prefix}_se/se.rds
    printf 'sample\\tplatform\\n' > ${prefix}_qc.tsv
    printf 'Sample\\n' > ${prefix}_qc_mqc.tsv
    printf '"%s":\\n    SummarizedExperiment: stub\\n    HDF5Array: stub\\n' "${task.process}" > versions.yml
    """
}
