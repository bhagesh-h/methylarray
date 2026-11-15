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
    tuple val(meta), path("${prefix}_snp_mqc.tsv")      , emit: multiqc
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
    printf '"%s":\\n    ewastools: stub\\n    illuminaio: stub\\n' "${task.process}" > versions.yml
    """
}
