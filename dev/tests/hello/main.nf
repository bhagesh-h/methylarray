#!/usr/bin/env nextflow

// Smoke test for the developer execution modes (dev/nfdev).
// One task runs in its own container and writes a file; the workflow
// publishes it so ownership and path handling can be checked.

params.outdir = null
params.greeting = 'hello from methylarray'

process HELLO {
    conda 'conda-forge::coreutils=9.5'
    container 'ubuntu:24.04'

    input:
    val greeting

    output:
    path 'hello.txt'

    script:
    """
    echo '${greeting}' > hello.txt
    echo "container user: \$(id -u):\$(id -g)" >> hello.txt
    echo "conda environment: \${CONDA_PREFIX:-none}" >> hello.txt
    """
}

workflow {
    main:
    ch_hello = HELLO(params.greeting)

    publish:
    hello = ch_hello
}

output {
    hello {
        path '.'
    }
}
