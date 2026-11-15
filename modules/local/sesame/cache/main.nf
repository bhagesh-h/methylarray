process SESAME_CACHE {
    tag "${platform}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // bhagesh/methylarray-sesame:1.28.1, built from environment.yml with containers/sesame/Dockerfile.
    // Pinned by digest only: Singularity and Apptainer reject a reference with both tag and digest.
    container 'docker.io/bhagesh/methylarray-sesame@sha256:75cf491e12071a0c3eca3479b227ecedaa35e521d8b6266b53497c5681716ea9'

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
