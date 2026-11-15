process SESAME_PREPROCESS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // bhagesh/methylarray-sesame:1.28.1, built from environment.yml with containers/sesame/Dockerfile.
    // Pinned by digest only: Singularity and Apptainer reject a reference with both tag and digest.
    container 'docker.io/bhagesh/methylarray-sesame@sha256:75cf491e12071a0c3eca3479b227ecedaa35e521d8b6266b53497c5681716ea9'

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
    echo | gzip > ${prefix}.controls.tsv.gz
    printf '"%s":\\n    sesame: stub\\n    sesameData: stub\\n' "${task.process}" > versions.yml
    """
}
