process EXPLORATION_PCA {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Built with containers/conda.Dockerfile from environment.yml (dev/nfdev image), because
    // Wave rate-limited the Singularity build: bhagesh/methylarray-exploration:1.0.0
    container 'docker.io/bhagesh/methylarray-exploration@sha256:caf0bd980d570c61b20e9925303ad719b3dd33a83ff95047419d44b11eb3fde4'

    input:
    tuple val(meta), path(se), path(samples), path(sample_qc), path(controls), path(cell_fractions), path(phenotypes)

    output:
    tuple val(meta), path("${prefix}_ctrl_pcs.tsv")          , emit: ctrl_pcs
    tuple val(meta), path("${prefix}_ctrl_pca_variance.tsv") , emit: ctrl_variance
    tuple val(meta), path("${prefix}_scores.tsv")            , emit: scores
    tuple val(meta), path("${prefix}_variance.tsv")          , emit: variance
    tuple val(meta), path("${prefix}_associations.tsv")      , emit: associations
    tuple val(meta), path("${prefix}_*_mqc.*")               , emit: multiqc
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "pca"
    template 'exploration_pca.R'

    stub:
    prefix = task.ext.prefix ?: "pca"
    """
    printf 'sample\\n' > ${prefix}_ctrl_pcs.tsv
    printf 'pc\\tvariance_explained\\tcumulative\\n' > ${prefix}_ctrl_pca_variance.tsv
    printf 'sample\\tPC1\\n' > ${prefix}_scores.tsv
    printf 'pc\\tvariance_explained\\tcumulative\\ttested\\n' > ${prefix}_variance.tsv
    printf 'variable\\tgroup\\tpc\\ttest\\tn\\tr_squared\\tp_value\\n' > ${prefix}_associations.tsv
    printf 'PC\\tvariance\\n' > ${prefix}_scree_mqc.tsv
    printf '"%s":\\n    isva: stub\\n    SummarizedExperiment: stub\\n    HDF5Array: stub\\n' "${task.process}" > versions.yml
    """
}
