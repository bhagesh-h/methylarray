process SESAME_CACHE {
    tag "${platform}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/bhagesh/methylarray-sesame:1.28.1'

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
