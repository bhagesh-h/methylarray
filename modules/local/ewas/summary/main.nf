process EWAS_SUMMARY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Same image as EWAS_LIMMA (identical environment.yml)
    container 'docker.io/bhagesh/methylarray-ewas@sha256:e208cdcd0107e83b7aaa316d5d67146f1e780e5e74cd0048a1d1ff277cee4fa2'

    input:
    path summaries, stageAs: 'summaries/*'
    path stats, stageAs: 'stats/*'
    path meqtl_cpgs

    output:
    path "${prefix}.tsv"    , emit: summary
    path "${prefix}*_mqc.tsv", emit: multiqc
    path "${prefix}_cells.tsv", emit: cells
    path "${prefix}_meqtl.tsv", emit: meqtl
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "ewas_summary"
    template 'ewas_summary.R'

    stub:
    prefix = task.ext.prefix ?: "ewas_summary"
    """
    printf 'model\\tn\\tstatus\\n' > ${prefix}.tsv
    printf 'Model\\n' > ${prefix}_mqc.tsv
    printf 'model\\tcell\\thits_fdr\\n' > ${prefix}_cells.tsv
    printf 'model\\tcontrast\\ttested\\ttested_meqtl\\thits_fdr\\thits_meqtl\\tshare_tested\\tshare_hits\\todds_ratio\\tp_value\\n' > ${prefix}_meqtl.tsv
    printf '"%s":\\n    jsonlite: stub\\n' "${task.process}" > versions.yml
    """
}
