process EPIDISH_RPC {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-epidish_bioconductor-hdf5array_bioconductor-summarizedexperiment_procps-ng_pruned:95c4af4791c4592c' :
        'community.wave.seqera.io/library/bioconductor-epidish_bioconductor-hdf5array_bioconductor-summarizedexperiment_procps-ng_pruned:0716f5da89501159' }"

    input:
    tuple val(meta), path(se), path(annotation)

    output:
    tuple val(meta), path("${prefix}.tsv")                   , emit: fractions
    tuple val(meta), path("${prefix}_reference_coverage.tsv"), emit: coverage
    tuple val(meta), path("${prefix}_mqc.tsv")               , emit: multiqc
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "cell_fractions"
    template 'epidish_rpc.R'

    stub:
    prefix = task.ext.prefix ?: "cell_fractions"
    """
    printf 'sample\\tlineage_B\\n' > ${prefix}.tsv
    printf 'reference\\treference_cpgs\\tsex_chromosome\\tnot_on_array\\tremoved_by_filter\\tfound\\n' > ${prefix}_reference_coverage.tsv
    printf 'Sample\\tB\\n' > ${prefix}_mqc.tsv
    printf '"%s":\\n    EpiDISH: stub\\n    SummarizedExperiment: stub\\n    HDF5Array: stub\\n' "${task.process}" > versions.yml
    """
}
