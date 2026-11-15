process SESAME_CACHE {
    tag "${platform}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-experimenthub_bioconductor-hdf5array_bioconductor-sesame_bioconductor-sesamedata_pruned:0c0cc322759e48bc' :
        'community.wave.seqera.io/library/bioconductor-experimenthub_bioconductor-hdf5array_bioconductor-sesame_bioconductor-sesamedata_pruned:8ca814a29ce1f4ed' }"

    input:
    val platform

    output:
    path "sesame_cache"         , emit: cache
    path "sesame_resources.tsv" , emit: resources
    path "versions.yml", emit: versions, topic: versions

    script:
    template 'sesame_cache.R'

    stub:
    """
    mkdir sesame_cache
    printf 'title\\teh_id\\tfile\\tmd5\\n' > sesame_resources.tsv
    printf '"%s":\\n    sesame: stub\\n    sesameData: stub\\n    ExperimentHub: stub\\n' "${task.process}" > versions.yml
    """
}
