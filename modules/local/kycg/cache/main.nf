process KYCG_CACHE {
    tag "${platform}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Built with containers/conda.Dockerfile from environment.yml (dev/nfdev image): bhagesh/methylarray-kycg:1.0.0
    container 'docker.io/bhagesh/methylarray-kycg@sha256:2026eb27c96964b18d3e8befae919209b102f809a0fda7afb01d35e3b22d15a9'

    input:
    val platform

    output:
    path "kycg_cache"        , emit: cache
    path "kycg_resources.tsv", emit: resources
    path "versions.yml", emit: versions, topic: versions

    script:
    template 'kycg_cache.R'

    stub:
    """
    mkdir kycg_cache
    printf 'database\\ttitle\\teh_id\\tfile\\tmd5\\n' > kycg_resources.tsv
    printf '"%s":\\n    knowYourCG: stub\\n    ExperimentHub: stub\\n' "${task.process}" > versions.yml
    """
}
