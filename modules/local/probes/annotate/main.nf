process PROBES_ANNOTATE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/procps-ng_r-base_r-data.table:c9210996d5e64487' :
        'community.wave.seqera.io/library/procps-ng_r-base_r-data.table:ec708683a2a90437' }"

    input:
    tuple val(meta), path(manifest, stageAs: 'EPICv2manifest.Rda'), path(probe_mask)

    output:
    tuple val(meta), path("${prefix}.tsv.gz")        , emit: annotation
    tuple val(meta), path("${prefix}_summary.tsv")   , emit: summary
    tuple val(meta), path("${prefix}_resources.tsv") , emit: resources
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "probe_annotation"
    template 'probes_annotate.R'

    stub:
    prefix = task.ext.prefix ?: "probe_annotation"
    """
    printf 'probe_id\\tcpg_id\\tprobe_type\\n' > ${prefix}.tsv
    Rscript -e 'con <- gzfile("${prefix}.tsv.gz", "w"); writeLines(readLines("${prefix}.tsv"), con); close(con)'
    rm ${prefix}.tsv
    printf 'flag\\tn\\n' > ${prefix}_summary.tsv
    printf 'resource\\tfile\\tmd5\\n' > ${prefix}_resources.tsv
    printf '"%s":\\n    data.table: stub\\n' "${task.process}" > versions.yml
    """
}
