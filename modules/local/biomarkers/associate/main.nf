process BIOMARKERS_ASSOCIATE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/procps-ng_r-base_r-jsonlite:c5d69391e11768e2' :
        'community.wave.seqera.io/library/procps-ng_r-base_r-jsonlite:78d6964c424dbe27' }"

    input:
    tuple val(meta), path(samples), path(sample_qc), path(clocks), path(cell_fractions), path(dnam_scores), path(phenotypes)

    output:
    tuple val(meta), path("${prefix}_agreement.tsv")       , emit: agreement
    tuple val(meta), path("${prefix}_cells_measures.tsv")  , emit: cells_measures
    tuple val(meta), path("${prefix}_cells_phenotypes.tsv"), emit: cells_phenotypes
    tuple val(meta), path("${prefix}_batch.tsv")           , emit: batch
    tuple val(meta), path("${prefix}_phenotypes.tsv")      , emit: phenotypes
    tuple val(meta), path("${prefix}_quartiles.tsv")       , emit: quartiles
    tuple val(meta), path("${prefix}_*_mqc.json")          , emit: multiqc, optional: true
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "biomarkers"
    template 'biomarkers_associate.R'

    stub:
    prefix = task.ext.prefix ?: "biomarkers"
    """
    printf 'measure1\\tmeasure2\\trho\\tp\\tn\\n' > ${prefix}_agreement.tsv
    printf 'cell\\tmeasure\\trho\\tp\\tn\\tq\\n' > ${prefix}_cells_measures.tsv
    printf 'cell\\ttrait\\ttype\\testimate\\tp\\tn\\tq\\n' > ${prefix}_cells_phenotypes.tsv
    printf 'measure\\tbatch\\tadj_r2\\tp\\tlevels\\tn\\n' > ${prefix}_batch.tsv
    printf 'measure\\ttrait\\ttype\\tversion\\testimate\\tp\\tn\\tq\\n' > ${prefix}_phenotypes.tsv
    printf 'measure\\ttrait\\ttype\\tgroup\\tvalue\\tmean_z\\tp_trend\\tq\\n' > ${prefix}_quartiles.tsv
    printf '"%s":\\n    R: stub\\n    jsonlite: stub\\n' "${task.process}" > versions.yml
    """
}
