# bhagesh-h/methylarray: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of bhagesh-h/methylarray, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Samplesheet validation for IDAT input: matching red and green files, unique IDATs, Sentrix ID and position from the file names, manual exclusions with a reason.
- Optional phenotype table (`--phenotypes`) and model definitions (`--models`), checked against the samplesheet and each other before any task starts.
- SeSAMe preprocessing (`--preprocessing sesame`): sesameData cache (`SESAME_CACHE`, reusable offline with `--sesame_cache`), per-array `prepSesame` (`SESAME_PREPROCESS`, default `QCDPB`) and one HDF5-backed SummarizedExperiment with the sample QC statistics (`SESAME_BUILDSE`), plus a MultiQC general statistics table.

### `Fixed`

- Numeric Sentrix IDs, plates and batches in the samplesheet are accepted (#3).
- The test profiles stage the GEO test IDATs four at a time to avoid HTTP 503 responses (#7).

### `Dependencies`

- Nextflow 26.04.2 or later (#4).
- SeSAMe 1.28.1, sesameData 1.28.0, ExperimentHub 3.0.0, SummarizedExperiment 1.40.0 and HDF5Array 1.38.0 on R 4.5.3, with the reference (netlib) BLAS.

### `Deprecated`
