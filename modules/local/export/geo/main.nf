process EXPORT_GEO {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // Same image as MINFI_PREPROCESS (identical environment.yml; data.table comes with minfi)
    container 'docker.io/bhagesh/methylarray-minfi@sha256:fbbd3f2e40d24e5ed33807f34367919ccba2e3a37e9b64e0edd80ced7793cbdc'

    input:
    tuple val(meta), path(samples), path(idats, stageAs: 'idat/*'), path(se), path(retained)

    output:
    tuple val(meta), path("${prefix}_processed_matrix.tsv.gz")  , emit: processed
    tuple val(meta), path("${prefix}_signal_intensities.tsv.gz"), emit: intensities
    tuple val(meta), path("${prefix}_samples.tsv")              , emit: samples
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "geo"
    template 'export_geo.R'

    stub:
    prefix = task.ext.prefix ?: "geo"
    """
    printf 'ID_REF\\n' | gzip > ${prefix}_processed_matrix.tsv.gz
    printf 'ID_REF\\n' | gzip > ${prefix}_signal_intensities.tsv.gz
    printf 'sample\\n' > ${prefix}_samples.tsv
    printf '"%s":\\n    minfi: stub\\n    data.table: stub\\n' "${task.process}" > versions.yml
    """
}
