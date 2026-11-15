process EWAS_SUMMARY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Same image as EWAS_LIMMA (identical environment.yml)
    container 'docker.io/bhagesh/methylarray-ewas:1.0.0'

    input:
    path summaries, stageAs: 'summaries/*'

    output:
    path "${prefix}.tsv"    , emit: summary
    path "${prefix}_mqc.tsv", emit: multiqc
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "ewas_summary"
    template 'ewas_summary.R'

    stub:
    prefix = task.ext.prefix ?: "ewas_summary"
    """
    printf 'model\\tn\\tstatus\\n' > ${prefix}.tsv
    printf 'Model\\n' > ${prefix}_mqc.tsv
    printf '"%s":\\n    jsonlite: stub\\n' "${task.process}" > versions.yml
    """
}
