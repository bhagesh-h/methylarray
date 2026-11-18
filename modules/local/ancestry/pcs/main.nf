process ANCESTRY_PCS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // The image of EXPLORATION_PCA (same environment): bhagesh/methylarray-exploration:1.0.0
    container 'docker.io/bhagesh/methylarray-exploration@sha256:caf0bd980d570c61b20e9925303ad719b3dd33a83ff95047419d44b11eb3fde4'

    input:
    tuple val(meta), path(se), path(annotation), path(samples), path(sample_qc), path(ctrl_pcs), path(cell_fractions), path(phenotypes)

    output:
    tuple val(meta), path("${prefix}.tsv")             , emit: pcs
    tuple val(meta), path("${prefix}_variance.tsv")    , emit: variance
    tuple val(meta), path("${prefix}_summary.tsv")     , emit: summary
    tuple val(meta), path("${prefix}_associations.tsv"), emit: associations
    tuple val(meta), path("${prefix}_*_mqc.*")         , emit: multiqc, optional: true
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "ancestry_pcs"
    template 'ancestry_pcs.R'

    stub:
    prefix = task.ext.prefix ?: "ancestry_pcs"
    """
    printf 'sample\\tancestry_pc1\\n' > ${prefix}.tsv
    printf 'pc\\tvariance_explained\\tcumulative\\tcovariate\\n' > ${prefix}_variance.tsv
    printf 'metric\\tvalue\\n' > ${prefix}_summary.tsv
    printf 'variable\\tpc\\tn\\tr_squared\\tp_value\\n' > ${prefix}_associations.tsv
    printf 'pc\\tvariance\\n' > ${prefix}_scree_mqc.tsv
    printf '"%s":\\n    R: stub\\n    SummarizedExperiment: stub\\n    HDF5Array: stub\\n' "${task.process}" > versions.yml
    """
}
