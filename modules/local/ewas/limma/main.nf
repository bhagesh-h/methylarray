process EWAS_LIMMA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // Built with containers/conda.Dockerfile from environment.yml (dev/nfdev image), because
    // Wave rate-limited anonymous builds: bhagesh/methylarray-ewas:1.0.0
    container 'docker.io/bhagesh/methylarray-ewas@sha256:e208cdcd0107e83b7aaa316d5d67146f1e780e5e74cd0048a1d1ff277cee4fa2'

    input:
    tuple val(meta), path(model), path(se), path(samples), path(sample_qc), path(phenotypes), path(cell_fractions), path(ctrl_pcs), path(ancestry_pcs), path(dnam_scores)

    output:
    tuple val(meta), path("${prefix}_stats.tsv.gz")    , emit: stats
    tuple val(meta), path("${prefix}_top.tsv")         , emit: top
    tuple val(meta), path("${prefix}_design_check.tsv"), emit: design_check
    tuple val(meta), path("${prefix}_resolved.json")   , emit: resolved
    tuple val(meta), path("${prefix}_summary.json")    , emit: summary
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "ewas_${meta.id}"
    template 'ewas_limma.R'

    stub:
    prefix = task.ext.prefix ?: "ewas_${meta.id}"
    """
    printf '\\x1f\\x8b\\x08\\x00\\x00\\x00\\x00\\x00\\x00\\x03\\x03\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00' > ${prefix}_stats.tsv.gz
    touch ${prefix}_top.tsv ${prefix}_design_check.tsv
    printf '{}\\n' > ${prefix}_resolved.json
    printf '{"model": "${meta.id}", "status": "stub"}\\n' > ${prefix}_summary.json
    printf '"%s":\\n    limma: stub\\n    bacon: stub\\n    SummarizedExperiment: stub\\n' "${task.process}" > versions.yml
    """
}
