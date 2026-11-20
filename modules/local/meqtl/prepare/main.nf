process MEQTL_PREPARE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/procps-ng_r-base_r-data.table:c9210996d5e64487' :
        'community.wave.seqera.io/library/procps-ng_r-base_r-data.table:ec708683a2a90437' }"

    input:
    tuple val(meta), path(meqtl), path(annotation)

    output:
    tuple val(meta), path("${prefix}.tsv")        , emit: cpgs
    tuple val(meta), path("${prefix}_summary.tsv"), emit: summary
    tuple val(meta), path("${prefix}_mqc.tsv")    , emit: multiqc
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "meqtl_cpgs"
    template 'meqtl_prepare.R'

    stub:
    prefix = task.ext.prefix ?: "meqtl_cpgs"
    """
    printf 'cpg_id\\tmeqtl_name\\tsnps\\ttop_snp\\ttop_p\\ttop_beta\\n' > ${prefix}.tsv
    printf 'metric\\tvalue\\n' > ${prefix}_summary.tsv
    printf "# id: 'meqtl_summary'\\n# plot_type: 'table'\\nTable\\trows_read\\nmeQTL\\t0\\n" > ${prefix}_mqc.tsv
    printf '"%s":\\n    R: stub\\n    data.table: stub\\n' "${task.process}" > versions.yml
    """
}
