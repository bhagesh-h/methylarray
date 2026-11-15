process ENRICHMENT_MISSMETHYL {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // Built with containers/conda.Dockerfile from environment.yml (dev/nfdev image): missMethyl
    // depends on every Illumina annotation package, beyond the Wave build time limit
    container 'docker.io/bhagesh/methylarray-enrichment@sha256:81293c7b16ad744c38bd8dca9cb06e2a598e7970c03cbf979662a6b7bd1a883e'

    input:
    tuple val(meta), path(stats), path(resolved)
    path gmts, stageAs: 'gmt/*'

    output:
    tuple val(meta), path("${prefix}.tsv")         , emit: enrichment
    tuple val(meta), path("${prefix}_summary.json"), emit: summary
    tuple val(meta), path("${prefix}_mqc.tsv")     , emit: multiqc, optional: true
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "enrichment_${meta.id}"
    template 'enrichment_missmethyl.R'

    stub:
    prefix = task.ext.prefix ?: "enrichment_${meta.id}"
    """
    printf 'collection\\tcontrast\\tset\\tterm\\tontology\\tn\\tde\\tp_value\\tfdr\\n' > ${prefix}.tsv
    printf '{"model": "${meta.id}", "status": "stub"}\\n' > ${prefix}_summary.json
    printf '"%s":\\n    missMethyl: stub\\n    org.Hs.eg.db: stub\\n    GO.db: stub\\n' "${task.process}" > versions.yml
    """
}
