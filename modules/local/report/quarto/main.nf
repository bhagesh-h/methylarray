process REPORT_QUARTO {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Built with containers/conda.Dockerfile from environment.yml (dev/nfdev image), because
    // Wave rate-limited anonymous builds: bhagesh/methylarray-report:1.0.0
    container 'docker.io/bhagesh/methylarray-report@sha256:9326a68532d0034dbc2e4bf1760e2a282b703c265c4f51f1b6b803a8eb239d42'

    input:
    path qmd, stageAs: 'report_template.qmd'
    path inputs, stageAs: 'inputs/*'

    output:
    path "methylarray_report.html", emit: report
    path "figures"                , emit: figures
    // Not on the 'versions' topic: the report reads the collated versions, which
    // are only complete when every process on the topic has finished
    path "versions.yml", emit: versions

    script:
    """
    export HOME=\$PWD XDG_CACHE_HOME=\$PWD/.cache XDG_DATA_HOME=\$PWD/.local XDG_CONFIG_HOME=\$PWD/.config
    # Quarto finds its tools (deno, pandoc, sass) through conda activation scripts;
    # Apptainer and Singularity start the image without activating the environment
    if [ -z "\${QUARTO_DENO:-}" ]; then
        prefix=\$(dirname "\$(dirname "\$(command -v quarto)")")
        for script in quarto.sh deno.sh deno_dom.sh; do
            [ -f "\$prefix/etc/conda/activate.d/\$script" ] && . "\$prefix/etc/conda/activate.d/\$script"
        done
    fi
    cp -L ${qmd} methylarray_report.qmd
    quarto render methylarray_report.qmd --to html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        R: \$(Rscript -e 'cat(as.character(getRversion()))')
    END_VERSIONS
    """

    stub:
    """
    touch methylarray_report.html
    mkdir figures
    printf '"%s":\\n    quarto: stub\\n    R: stub\\n' "${task.process}" > versions.yml
    """
}
