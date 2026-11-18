process VARIANCE_PARTITION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Built with containers/conda.Dockerfile from environment.yml (dev/nfdev image): bhagesh/methylarray-variancepartition:1.0.0
    container 'docker.io/bhagesh/methylarray-variancepartition@sha256:a5d5c713c5ce40b3fde2a1b29dcc6a30278cd302e31f58508ee62cea8f446943'

    input:
    tuple val(meta), path(se), path(samples), path(sample_qc), path(cell_fractions)

    output:
    tuple val(meta), path("${prefix}.tsv.gz")     , emit: variance
    tuple val(meta), path("${prefix}_summary.tsv"), emit: summary
    tuple val(meta), path("${prefix}_mqc.tsv")    , emit: multiqc, optional: true
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "variance_partition"
    template 'variance_partition.R'

    stub:
    prefix = task.ext.prefix ?: "variance_partition"
    """
    printf 'cpg_id\\tresidual\\n' | gzip -n > ${prefix}.tsv.gz
    printf 'variable\\ttype\\tmedian\\tmean\\tq90\\tcpgs_over_10\\tarrays\\tcpgs\\tformula\\tnote\\n' > ${prefix}_summary.tsv
    printf 'variable\\tresidual\\n' > ${prefix}_mqc.tsv
    printf '"%s":\\n    variancePartition: stub\\n    lme4: stub\\n    HDF5Array: stub\\n' "${task.process}" > versions.yml
    """
}
