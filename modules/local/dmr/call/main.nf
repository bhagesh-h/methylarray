process DMR_CALL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Built with containers/conda.Dockerfile from environment.yml (dev/nfdev image), because
    // Wave rate-limited anonymous builds: bhagesh/methylarray-dmr:1.0.0
    container 'docker.io/bhagesh/methylarray-dmr@sha256:1edf71d7343a0078ee854b7a55d0bc4787b29923bac5a7f8fc7f50f3504bd96d'

    input:
    tuple val(meta), path(stats), path(resolved), path(se)
    path dmrff_source, stageAs: 'dmrff.tar.gz'

    output:
    tuple val(meta), path("${prefix}_dmrcate.tsv")  , emit: dmrcate
    tuple val(meta), path("${prefix}_dmrff.tsv")    , emit: dmrff
    tuple val(meta), path("${prefix}_consensus.tsv"), emit: consensus
    tuple val(meta), path("${prefix}_consensus.bed"), emit: bed
    tuple val(meta), path("${prefix}_summary.json") , emit: summary
    tuple val(meta), path("${prefix}_mqc.tsv")      , emit: multiqc
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "dmr_${meta.id}"
    template 'dmr_call.R'

    stub:
    prefix = task.ext.prefix ?: "dmr_${meta.id}"
    """
    touch ${prefix}_dmrcate.tsv ${prefix}_dmrff.tsv ${prefix}_consensus.tsv ${prefix}_consensus.bed
    printf '{"model": "${meta.id}", "status": "stub"}\\n' > ${prefix}_summary.json
    printf "# id: 'dmr_summary'\\n# plot_type: 'table'\\nModel\\tstatus\\n${meta.id}\\tstub\\n" > ${prefix}_mqc.tsv
    printf '"%s":\\n    DMRcate: stub\\n    dmrff: stub\\n' "${task.process}" > versions.yml
    """
}
