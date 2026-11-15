process MINFI_PREPROCESS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // bhagesh/methylarray-minfi:1.56.0, built from environment.yml with containers/minfi/Dockerfile.
    container 'docker.io/bhagesh/methylarray-minfi@sha256:fbbd3f2e40d24e5ed33807f34367919ccba2e3a37e9b64e0edd80ced7793cbdc'

    input:
    tuple val(meta), path(samples), path(idats, stageAs: 'idat/*'), path(sesame_se)

    output:
    tuple val(meta), path("${prefix}_se")          , emit: se
    tuple val(meta), path("${prefix}_qc.tsv")      , emit: qc
    tuple val(meta), path("${prefix}_qc_mqc.tsv")  , emit: multiqc
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'minfi_preprocess.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_se
    touch ${prefix}_se/assays.h5 ${prefix}_se/se.rds
    printf 'sample\\tplatform\\n' > ${prefix}_qc.tsv
    printf 'Sample\\n' > ${prefix}_qc_mqc.tsv
    printf '"%s":\\n    minfi: stub\\n    ENmix: stub\\n    SummarizedExperiment: stub\\n    HDF5Array: stub\\n' "${task.process}" > versions.yml
    """
}
