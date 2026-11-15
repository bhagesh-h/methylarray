process SAMPLEQC_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/procps-ng_r-base_r-jsonlite:c5d69391e11768e2' :
        'community.wave.seqera.io/library/procps-ng_r-base_r-jsonlite:78d6964c424dbe27' }"

    input:
    tuple val(meta), path(samples), path(preprocessing_qc), path(ewastools_qc), path(snp_conflicts)

    output:
    tuple val(meta), path("${prefix}.tsv")         , emit: qc
    tuple val(meta), path("${prefix}_retained.txt"), emit: retained
    tuple val(meta), path("${prefix}*_mqc.*")      , emit: multiqc
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'sampleqc_summary.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tstatus\\treasons\\n' > ${prefix}.tsv
    touch ${prefix}_retained.txt
    printf 'Sample\\tqc_status\\n' > ${prefix}_mqc.tsv
    printf '"%s":\\n    R: stub\\n' "${task.process}" > versions.yml
    """
}
