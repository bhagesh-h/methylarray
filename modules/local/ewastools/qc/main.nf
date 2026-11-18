process EWASTOOLS_QC {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-illuminaio_procps-ng_r-base_r-data.table_pruned:6a3be0c701c401e0' :
        'community.wave.seqera.io/library/bioconductor-illuminaio_procps-ng_r-base_r-data.table_pruned:ef0153fbffb03c9a' }"

    input:
    tuple val(meta), path(samples), path(idats, stageAs: 'idat/*')
    path ewastools_source

    output:
    tuple val(meta), path("${prefix}_qc.tsv")           , emit: qc
    tuple val(meta), path("${prefix}_snp_conflicts.tsv"), emit: conflicts
    tuple val(meta), path("${prefix}_controls.tsv.gz")  , emit: controls
    tuple val(meta), path("${prefix}_snp_summary.tsv")  , emit: snp_summary
    tuple val(meta), path("${prefix}_snp_histograms.tsv"), emit: snp_histograms
    tuple val(meta), path("${prefix}_relatedness.tsv.gz"), emit: relatedness
    tuple val(meta), path("${prefix}_relatedness_summary.tsv"), emit: relatedness_summary
    tuple val(meta), path("${prefix}_*mqc.tsv")         , emit: multiqc
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'ewastools_qc.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tplatform\\n' > ${prefix}_qc.tsv
    printf 'group\\tsample1\\texpected_donor1\\tsample2\\texpected_donor2\\tagreement\\n' > ${prefix}_snp_conflicts.tsv
    printf '\\x1f\\x8b\\x08\\x00\\x00\\x00\\x00\\x00\\x00\\x03\\x03\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00' > ${prefix}_controls.tsv.gz
    printf 'sample\\tsex_used\\tsnp_outlier\\thet_share_autosomal\\thet_share_x\\tmax_agreement\\tsnps_autosomal\\tsnps_x_like\\n' > ${prefix}_snp_summary.tsv
    printf 'histogram\\tbin_start\\tbin_end\\tcount\\n' > ${prefix}_snp_histograms.tsv
    printf 'Sample\\n' > ${prefix}_snp_mqc.tsv
    # the gzip-compressed header line: the Singularity image has no gzip command
    printf '\\x1f\\x8b\\x08\\x00\\x00\\x00\\x00\\x00\\x02\\xff\\x1d\\x8c\\x41\\x0a\\x80\\x30\\x0c\\x04\\xcf\\xf1\\x37\\xfa\\xa1\\x10\\x9b\\x55\\x02\\x6d\\x5a\\x92\\x22\\x3e\\x5f\\xe9\\x69\\x66\\x2e\\x93\\xd2\\x46\\xc5\\x4e\\xb9\\x78\\x50\\xfa\\x48\\xb2\\x33\\xb9\\x41\\x7c\\xc9\\x23\\x41\\xb5\\x2b\\x9b\\xc2\\xa7\\x15\\xa9\\xab\\x2e\\x8b\\x9c\\xac\\xb8\\x03\\xa0\\x40\\x95\\x69\\xdd\\x09\\xef\\x40\\x99\\x50\\xfe\\x77\\x60\\xed\\xde\\x63\\xfb\\x00\\x33\\x3d\\x3b\\x74\\x62\\x00\\x00\\x00' > ${prefix}_relatedness.tsv.gz
    printf 'metric\\tvalue\\n' > ${prefix}_relatedness_summary.tsv
    printf '"%s":\\n    ewastools: stub\\n    illuminaio: stub\\n' "${task.process}" > versions.yml
    """
}
