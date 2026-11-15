process PROBES_FILTER {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-hdf5array_bioconductor-summarizedexperiment_procps-ng_r-base_r-data.table:f08a921e65e6bc1c' :
        'community.wave.seqera.io/library/bioconductor-hdf5array_bioconductor-summarizedexperiment_procps-ng_r-base_r-data.table:460ab90a75a0c805' }"

    input:
    tuple val(meta), path(se), path(retained), path(annotation)

    output:
    tuple val(meta), path("${prefix}_se")               , emit: se
    tuple val(meta), path("${prefix}_uncollapsed_se")   , emit: se_uncollapsed
    tuple val(meta), path("${prefix}_attrition.tsv")    , emit: attrition
    tuple val(meta), path("${prefix}_exclusions.tsv.gz"), emit: exclusions
    tuple val(meta), path("${prefix}_attrition_mqc.tsv"), emit: multiqc
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}_filtered"
    template 'probes_filter.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_filtered"
    """
    mkdir ${prefix}_se ${prefix}_uncollapsed_se
    touch ${prefix}_se/assays.h5 ${prefix}_se/se.rds ${prefix}_uncollapsed_se/assays.h5 ${prefix}_uncollapsed_se/se.rds
    printf 'step\\tremoved\\tremaining\\tn_with_reason\\n' > ${prefix}_attrition.tsv
    Rscript -e 'con <- gzfile("${prefix}_exclusions.tsv.gz", "w"); writeLines("probe_id\\tfirst_reason\\treasons", con); close(con)'
    printf 'step\\tremoved\\n' > ${prefix}_attrition_mqc.tsv
    printf '"%s":\\n    SummarizedExperiment: stub\\n    HDF5Array: stub\\n    data.table: stub\\n' "${task.process}" > versions.yml
    """
}
