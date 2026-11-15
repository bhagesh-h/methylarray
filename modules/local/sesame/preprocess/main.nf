process SESAME_PREPROCESS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-experimenthub_bioconductor-hdf5array_bioconductor-sesame_bioconductor-sesamedata_pruned:0c0cc322759e48bc' :
        'community.wave.seqera.io/library/bioconductor-experimenthub_bioconductor-hdf5array_bioconductor-sesame_bioconductor-sesamedata_pruned:8ca814a29ce1f4ed' }"

    input:
    tuple val(meta), path(idat_red), path(idat_grn)
    path sesame_cache

    output:
    tuple val(meta), path("*.sesame.rds")       , emit: rds
    tuple val(meta), path("*.qc.tsv")           , emit: qc
    tuple val(meta), path("*.controls.tsv.gz")  , emit: controls
    path "versions.yml", emit: versions, topic: versions

    script:
    template 'sesame_preprocess.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sesame.rds
    printf 'sample\\tplatform\\n${meta.id}\\tEPICv2\\n' > ${prefix}.qc.tsv
    # an empty gzip stream: the Singularity image has no gzip command
    printf '\\x1f\\x8b\\x08\\x00\\x00\\x00\\x00\\x00\\x00\\x03\\x03\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00' > ${prefix}.controls.tsv.gz
    printf '"%s":\\n    sesame: stub\\n    sesameData: stub\\n' "${task.process}" > versions.yml
    """
}
