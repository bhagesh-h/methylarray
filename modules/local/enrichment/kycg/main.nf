process ENRICHMENT_KYCG {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // Built with containers/conda.Dockerfile from environment.yml (dev/nfdev image): bhagesh/methylarray-kycg:1.0.0
    container 'docker.io/bhagesh/methylarray-kycg@sha256:2026eb27c96964b18d3e8befae919209b102f809a0fda7afb01d35e3b22d15a9'

    input:
    tuple val(meta), path(stats), path(resolved)
    path cache, stageAs: 'kycg_cache'
    path annotation

    output:
    tuple val(meta), path("${prefix}.tsv")         , emit: enrichment
    tuple val(meta), path("${prefix}_summary.json"), emit: summary
    tuple val(meta), path("${prefix}_mqc.tsv")     , emit: multiqc, optional: true
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "kycg_${meta.id}"
    template 'enrichment_kycg.R'

    stub:
    prefix = task.ext.prefix ?: "kycg_${meta.id}"
    """
    printf 'contrast\\tquery\\tdatabase\\tset\\tlog2_or\\tp_value\\tfdr\\toverlap\\tquery_probes\\tset_probes\\n' > ${prefix}.tsv
    printf '{"model": "${meta.id}", "status": "stub"}\\n' > ${prefix}_summary.json
    printf '"%s":\\n    knowYourCG: stub\\n    ExperimentHub: stub\\n' "${task.process}" > versions.yml
    """
}
