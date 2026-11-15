# bhagesh-h/methylarray: Developer guide

This guide covers how to develop and test the pipeline on a workstation where only Docker is installed. End users do not need any of this: they run the pipeline with their own Nextflow and a container or conda profile, as described in [usage.md](usage.md).

## Requirements

- Linux with Docker (the developer machine needs nothing else; Nextflow, nf-core tools, nf-test, prek and Apptainer run inside containers).
- About 30 GB of free memory and 50 GB of free disk for the `test` profile and the module tests.
- Network access to GitHub, NCBI GEO (test IDATs), Bioconductor ExperimentHub (SeSAMe resources) and `community.wave.seqera.io` (task images).

## Setup

All developer commands go through `dev/nfdev`, which starts the devkit image with the repository and one output folder mounted at the same absolute paths as on the host.

1. Create `.nfdev.env` in the repository root. The file is ignored by git.

   ```bash
   # Every output, work directory, test result, cache and log goes here
   NFDEV_OUT_ROOT=/absolute/path/to/methylarray-dev
   # Optional: folders mounted read-only, for example local IDAT files
   NFDEV_DATA_DIRS=/absolute/path/to/idats
   ```

   Further settings, all optional, are listed in the header of `dev/nfdev`. `NFDEV_TASK_MEMORY_GB` (default 28) caps the memory of all tasks of a run together; with the Nextflow JVM on top, a run stays under 32 GB.

2. Build the devkit images:

   ```bash
   dev/nfdev build                    # methylarray-devkit:dev (Nextflow, nf-core tools, nf-test, prek, nf-metro)
   dev/nfdev --engine apptainer build # also methylarray-apptainer-kit:dev, for the apptainer and singularity profiles
   ```

3. Check the setup with the small test pipeline in `dev/tests/hello`:

   ```bash
   dev/nfdev nextflow run "$PWD/dev/tests/hello" -profile docker
   ```

Nothing is installed on the host. The repository is never written to by a run: Nextflow launches in `NFDEV_OUT_ROOT`, and a prek hook (`block-pipeline-outdir`) refuses commits that contain pipeline output folders.

## nfdev commands

| Command                                  | What it does                                                                                                          |
| ---------------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| `dev/nfdev nextflow <args>`              | Runs Nextflow with `NFDEV_OUT_ROOT` as launch directory                                                               |
| `dev/nfdev test <args>`                  | Runs `nf-test test <args>` in the repository; test work directories go to `NFDEV_OUT_ROOT/tests/nf-test`              |
| `dev/nfdev lint`                         | Exports the staged tree and runs `nf-core pipelines lint` and `nextflow lint` on it, then prek on the repository      |
| `dev/nfdev container <module> [--write]` | Builds the module images from its `environment.yml` with Wave and prints or writes the `container` directive          |
| `dev/nfdev metro`                        | Renders `assets/metro_map.mmd` to `docs/images/metro_map.svg` and two preview PNGs                                    |
| `dev/nfdev exec <cmd>`                   | Runs any command in the devkit, for example `dev/nfdev exec nf-core pipelines schema lint`                            |
| `dev/nfdev shell`                        | Interactive shell in the devkit                                                                                       |
| `dev/nfdev clean`                        | Deletes work, results and test folders, download caches and all but the 10 newest logs; keeps images and environments |

Every command writes its full output to `NFDEV_OUT_ROOT/logs/`.

### Engines

`--engine` selects how tasks run. Use it with the matching `-profile`.

| Engine             | Profile                    | Notes                                                                                          |
| ------------------ | -------------------------- | ---------------------------------------------------------------------------------------------- |
| `docker` (default) | `docker`                   | Task containers run on the host Docker daemon                                                  |
| `apptainer`        | `apptainer`, `singularity` | Apptainer inside a privileged container; images are cached in `NFDEV_OUT_ROOT/cache/apptainer` |
| `conda`            | `conda`, `mamba`           | micromamba environments inside the devkit, cached in `NFDEV_OUT_ROOT/cache/conda`              |
| `dind`             | `docker`                   | Docker in Docker sidecar, for testing without the host daemon (`dev/nfdev dind-stop` ends it)  |

### Nextflow version

The devkit pins one Nextflow release. To run another one, as the CI matrix does, set `NXF_VER`:

```bash
NXF_VER=26.04.2 dev/nfdev test tests/default.nf.test --profile +docker
```

The release is downloaded once into `NFDEV_OUT_ROOT/cache/nextflow-dist/`. The pipeline needs Nextflow 26.04.2 or later.

## Running the pipeline

Use one shared work directory and `-resume`, so staged IDATs and finished tasks are reused between runs:

```bash
OUT=/absolute/path/to/methylarray-dev   # NFDEV_OUT_ROOT
dev/nfdev nextflow run "$PWD" -profile test,docker \
    --outdir "$OUT/results/sesame" -w "$OUT/work" -resume
```

Change the route with `--preprocessing minfi|enmix|combined` and a separate `--outdir`. Use `--engine apptainer` with `-profile test,apptainer`, or `--engine conda` with `-profile test,conda`, to check the other engines. Run `dev/nfdev clean` when a feature is finished.

## Repository layout

| Path                                                    | Content                                                                                                  |
| ------------------------------------------------------- | -------------------------------------------------------------------------------------------------------- |
| `main.nf`                                               | Entry point: initialisation, `METHYLARRAY`, completion                                                   |
| `workflows/methylarray.nf`                              | The main workflow: routes the arrays through preprocessing and QC, collects MultiQC files                |
| `subworkflows/local/<name>/main.nf`                     | Pipeline steps built from several modules (`preprocess_sesame`, `preprocess_minfi`, `sample_qc`)         |
| `subworkflows/local/utils_nfcore_methylarray_pipeline/` | Input validation (samplesheet, phenotypes, models), the shared `sampleTable` helper and the methods text |
| `modules/local/<tool>/<subtool>/`                       | One process per folder (see below)                                                                       |
| `modules/nf-core/`, `subworkflows/nf-core/`             | Installed from nf-core with `nf-core modules install`; never edited by hand                              |
| `conf/base.config`                                      | Resource labels (`process_single`, `process_low`, `process_medium`, `process_high`)                      |
| `conf/modules.config`                                   | Per-process `ext.args` and `publishDir`                                                                  |
| `conf/test.config`, `conf/test_full.config`             | The `test` and `test_full` profiles                                                                      |
| `nextflow.config`, `nextflow_schema.json`               | Parameters with defaults, and their schema (help text, types, ranges)                                    |
| `assets/schema_*.json`                                  | Schemas of the samplesheet, phenotype table and model file                                               |
| `tests/`                                                | Pipeline-level nf-test files, snapshots, `.nftignore` and the test samplesheets                          |
| `dev/`                                                  | `nfdev`, the devkit Dockerfiles, `wave_container.py` and the hello-world test pipeline                   |
| `docs/`                                                 | `usage.md`, `output.md`, this guide                                                                      |

## Writing a local module

A module is one process in `modules/local/<tool>/<subtool>/`, for example `modules/local/ewastools/qc/` for the process `EWASTOOLS_QC`.

```text
modules/local/ewastools/qc/
├── main.nf               # the process
├── environment.yml       # conda packages with exact versions
├── meta.yml              # inputs, outputs, tools, licences
├── templates/
│   └── ewastools_qc.R    # the script, run by `template`
└── tests/
    ├── main.nf.test      # nf-test cases
    ├── main.nf.test.snap # snapshots, written by nf-test
    └── samples.tsv       # small test inputs, if needed
```

### main.nf

```groovy
process TOOL_SUBTOOL {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/...' :
        'community.wave.seqera.io/library/...' }"

    input:
    tuple val(meta), path(samples), path(idats, stageAs: 'idat/*')

    output:
    tuple val(meta), path("${prefix}_qc.tsv"), emit: qc
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'tool_subtool.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tplatform\\n' > ${prefix}_qc.tsv
    printf '"%s":\\n    tool: stub\\n' "${task.process}" > versions.yml
    """
}
```

Conventions:

- Every output has an `emit` name. `versions.yml` is emitted on the `versions` topic; the main workflow collects the topic, so modules do not pass versions around.
- The stub creates every output file with the same names, so `-stub` runs test the wiring without computing anything.
- Tool options come from `task.ext.args`, set in `conf/modules.config`, never from `params` inside the module.
- Optional inputs are passed as `[]` and checked in the script (`nzchar("${input}")` in R).

### R templates

The pipeline's analysis code is R, run with `template`. Nextflow fills in `${...}` before the script runs, so:

- Use `x[["column"]]` instead of `x$column`; a bare `$` is read by Nextflow.
- Write a literal backslash twice (`"\\t"`).
- Read options with the `parse_args` function used by the existing templates: `ext.args` holds `--name value` pairs, each name must have a default, and unknown names stop the task.
- Set `Sys.setenv(HOME = getwd())` before loading packages that write caches.
- Sort samples by name, so the output order does not depend on staging order.
- Write `versions.yml` at the end of the template, one line per tool, under the key `"${task.process}"`.
- Never drop an array or a probe silently: mask values or flag samples and record the reason in an output table.

### environment.yml

```yaml
# yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json
channels:
  - conda-forge
  - bioconda
dependencies:
  - bioconda::bioconductor-sesame=1.28.1
  - conda-forge::procps-ng=4.0.6
  - conda-forge::r-base=4.5.3
```

Pin every package the script loads to an exact version, sorted by name, with its channel. Use the same `r-base` and Bioconductor release in all modules. Find versions with:

```bash
docker run --rm mambaorg/micromamba:2.10.0-ubuntu24.04 \
    micromamba repoquery search -c conda-forge -c bioconda --platform linux-64 bioconductor-sesame
```

Bioconductor packages that call preprocessCore (SeSAMe dye-bias correction) also need the reference BLAS, `conda-forge::libblas=3.11.0=*_netlib` and `conda-forge::liblapack=3.11.0=*_netlib`, because preprocessCore cannot start its threads next to OpenBLAS on recent glibc.

### Containers

Task images are built from `environment.yml` by Wave and stored in the public Seqera Containers registry. Nobody builds or pushes images by hand:

```bash
dev/nfdev container modules/local/tool/subtool --write
```

The command requests a Docker image and a Singularity image, waits for both (a few minutes for a new environment, instantly for a known one) and writes the `container` directive into `main.nf`. Run it again whenever `environment.yml` changes. Wave stops builds after about 15 minutes. The one environment that exceeds this is `MINFI_PREPROCESS`: the EPIC v2 manifest and annotation data packages alone take about 9 minutes to install. Its image is therefore built from `containers/minfi/Dockerfile`, pushed once to Docker Hub as `bhagesh/methylarray-minfi` and referenced by digest (Apptainer and Singularity convert the Docker image). Rebuild and push it only when that module's `environment.yml` changes:

```bash
docker build -f containers/minfi/Dockerfile -t bhagesh/methylarray-minfi:<minfi version> .
docker push bhagesh/methylarray-minfi:<minfi version>
docker inspect --format '{{json .RepoDigests}}' bhagesh/methylarray-minfi:<minfi version>
```

A package that is not on conda but is plain R (no compiled code) is installed by the module from a staged source file instead, as `EWASTOOLS_QC` does with `--ewastools_source`: the source tarball is a pinned input file, so it is cached, works offline from a local copy, and the image only holds the conda dependencies.

### meta.yml

Describe every input and output in the order of `main.nf`, with `type`, `description`, `pattern` and EDAM `ontologies` for files, and list each tool with homepage, documentation, DOI and licence. Copy the layout from `modules/local/minfi/preprocess/meta.yml`. Add the tool's reference to `CITATIONS.md`.

## Adding a step to the pipeline

1. Create the module as above, with `environment.yml` first, then `dev/nfdev container ... --write`.
2. Add its options to `conf/modules.config`:

   ```groovy
   withName: 'TOOL_SUBTOOL' {
       ext.args   = { "--platform ${params.platform} --threshold ${params.some_threshold}" }
       publishDir = [
           path: { "${params.outdir}/qc/tool" },
           mode: params.publish_dir_mode,
           saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
       ]
   }
   ```

3. Add new parameters to `nextflow.config` with a default, and to `nextflow_schema.json` in the matching group, with `description`, `help_text`, type and range. Booleans default to `false`. Format the schema with `dev/nfdev exec prek run prettier --files nextflow_schema.json`.
4. Call the module from a subworkflow in `subworkflows/local/<name>/main.nf` and the subworkflow from `workflows/methylarray.nf`. Inputs shared by many arrays are built with `sampleTable(ch_arrays, with_idat_prefix)` and `.collect()`. Add MultiQC files (`*_mqc.tsv`, `*_mqc.json`) to `ch_multiqc_files`.
5. Document outputs in `docs/output.md`, parameters in `docs/usage.md`, the change in `CHANGELOG.md`, and update `assets/metro_map.mmd` (then `dev/nfdev metro`). If the tool is cited in the methods text, extend `toolCitationText()` and `toolBibliographyText()` in `subworkflows/local/utils_nfcore_methylarray_pipeline/main.nf`.
6. Write the tests (next section), run lint, commit.

## Testing with nf-test

There are three levels of tests. Each `main.nf.test` has a snapshot file next to it; nf-test compares outputs against it and fails on any difference.

| Level       | Location                                | What it runs                                                                                                   |
| ----------- | --------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| Module      | `modules/local/<tool>/<subtool>/tests/` | One process on two public EPIC v2 arrays (GEO GSE286313) or small TSV fixtures, plus a stub test               |
| Subworkflow | `subworkflows/local/<name>/tests/`      | Input validation and helpers, offline                                                                          |
| Pipeline    | `tests/*.nf.test`                       | The whole pipeline with `-profile test`; only `tests/default.nf.test` uses real data, the others are stub runs |

### Writing a module test

```groovy
nextflow_process {

    name "Test Process TOOL_SUBTOOL"
    script "../main.nf"
    process "TOOL_SUBTOOL"

    tag "modules"
    tag "modules_local"
    tag "tool"
    tag "tool/subtool"

    test("two EPIC v2 arrays") {
        when {
            process {
                """
                input[0] = [
                    [ id: 'test' ],
                    file("$moduleTestDir/samples.tsv", checkIfExists: true),
                    [ file('https://ftp.ncbi.nlm.nih.gov/geo/samples/...', checkIfExists: true) ]
                ]
                """
            }
        }
        then {
            assertAll(
                { assert process.success },
                { assert snapshot(process.out.qc, process.out.versions).match() }
            )
        }
    }

    test("two EPIC v2 arrays - stub") {
        options "-stub"
        when {
            process {
                """
                input[0] = [ [ id: 'test' ], file("$moduleTestDir/samples.tsv", checkIfExists: true), [] ]
                """
            }
        }
        then {
            assertAll(
                { assert process.success },
                { assert snapshot(process.out).match() }
            )
        }
    }
}
```

Guidelines:

- Assert the biology, not only the snapshot: predicted sex of the two test arrays (one female, one male), sample order, expected flags. A snapshot alone also passes when a wrong result is stable.
- Snapshot only files with stable content. Binary or time-stamped files (HDF5, RDS, HTML) go into `tests/.nftignore` for pipeline tests, or are checked by name for module tests.
- Keep real-data tests to two arrays. New pipeline-level tests are stub runs (`options "-stub"`), so CI stays short.
- Module-specific configuration (for example `ext.args`) goes into `tests/nextflow.config` next to the test, loaded with `config "./nextflow.config"`.

### Running tests

Always pass `--profile +docker` (or `test,docker`): the `+` adds the engine to the `test` profile set in `nf-test.config`. A bare `--profile docker` drops the test resource limits and `process_high` tasks then ask for 72 GB.

```bash
# One module
dev/nfdev test modules/local/ewastools/qc/tests/main.nf.test --profile +docker

# Everything a change can affect, compared with dev
dev/nfdev test --changed-since origin/dev --profile +docker

# By tag
dev/nfdev test --tag sesame --profile +docker

# Pipeline test
dev/nfdev test tests/default.nf.test --profile +docker

# Other engines
dev/nfdev --engine apptainer test modules/local/ewastools/qc/tests/main.nf.test --profile +apptainer
dev/nfdev --engine conda test modules/local/ewastools/qc/tests/main.nf.test --profile +conda
```

After an intended change of output, rewrite the snapshots and review the diff before committing:

```bash
dev/nfdev test modules/local/ewastools/qc/tests/main.nf.test --profile +docker --update-snapshot
git diff modules/local/ewastools/qc/tests/main.nf.test.snap
```

Remove snapshots of renamed tests (`--clean-snapshot`, or delete the entry). A module change is ready when its tests pass with Docker, Apptainer and conda against the same snapshot.

### Lint

```bash
git add -A
dev/nfdev lint
```

`dev/nfdev lint` checks the staged tree: `nf-core pipelines lint` (0 failures required; accepted warnings are listed in `docs/development/lint_exceptions.md`), `nextflow lint` (no errors) and the prek hooks (prettier, whitespace, Nextflow lint, and blocks on pipeline output folders, IDAT and spreadsheet files).

## Git and pull requests

- Branch from `dev`, one branch per feature (for example a preprocessing route or a QC subworkflow). Fixes found while working on it are extra commits on the same branch.
- Commit with a sign-off (`git commit -s`), one topic per commit, with a message that says what changed and why. Never commit IDAT files, spreadsheets, study data or pipeline output.
- Push the branch and open one pull request into `dev`. CI then runs:
  - `nf-core` lint and prek (`linting.yml`);
  - nf-test (`nf-test.yml`) for the tests affected by the change, in up to 7 shards. Changes to `nextflow.config`, `nextflow_schema.json`, `conf/test.config` or the CI files run every test.
- Pull requests into `dev` test with Docker and Nextflow 26.04.2. Pull requests into `main` and releases also test conda, Singularity and the latest Nextflow.
- The maintainer reviews and merges. GitHub issues track bioinformatics steps and major features; reference them in the pull request and close them with the commit IDs after the merge.

## Cleaning up

```bash
dev/nfdev clean
```

Keeps the Apptainer images, conda environments, SeSAMe resources and Nextflow releases under `NFDEV_OUT_ROOT/cache`, and the 10 newest logs. Everything else under `NFDEV_OUT_ROOT` (work directories, results, nf-test folders) is deleted.
