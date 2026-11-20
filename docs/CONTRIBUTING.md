---
title: Contributing
markdownPlugin: checklist
---

# `bhagesh-h/methylarray`: Contributing guidelines

Hi there!
Thanks for taking an interest in improving bhagesh-h/methylarray.

This page describes the recommended nf-core way to contribute to both bhagesh-h/methylarray and nf-core pipelines in general, including:

- [General contribution guidelines](#general-contribution-guidelines): common procedures or guides across all nf-core pipelines.
- [Pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines): procedures or guides specific to the development conventions of bhagesh-h/methylarray.

## General contribution guidelines

### Contribution quick start

To contribute code to any nf-core pipeline:

- [ ] Ensure you have Nextflow (26.04.2 or later), nf-core tools, and nf-test installed. See the [nf-core/tools repository](https://github.com/nf-core/tools) for instructions, or use the Docker-only setup of the [developer guide](developer.md).
- [ ] Check whether a GitHub [issue](https://github.com/bhagesh-h/methylarray/issues) about your idea already exists. If an issue does not exist, create one so that others are aware you are working on it.
- [ ] [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [bhagesh-h/methylarray repository](https://github.com/bhagesh-h/methylarray) to your GitHub account.
- [ ] Create a branch on your forked repository and make your changes following [pipeline conventions](#pipeline-contribution-conventions) (if applicable).
- [ ] To fix major bugs, name your branch `patch` and follow the [patch release](#patch-release) process.
- [ ] Update relevant documentation within the `docs/` folder, `nextflow_schema.json`, `CHANGELOG.md` and `CITATIONS.md`.
- [ ] Run and/or update tests. See [Testing](#testing) for more information.
- [ ] [Lint](#lint-tests) your code with nf-core/tools.
- [ ] Submit a pull request (PR) against the `dev` branch and request a review.

If you are not used to this workflow with Git, see the [GitHub documentation](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or [Git resources](https://try.github.io/) for more information.

### GitHub Codespaces

You can contribute to bhagesh-h/methylarray without installing a local development environment on your machine by using [GitHub Codespaces](https://github.com/codespaces).

[GitHub Codespaces](https://github.com/codespaces) is an online developer environment that runs in your browser, complete with VS Code and a terminal.
Most nf-core repositories include a devcontainer configuration, which creates a GitHub Codespaces environment specifically for Nextflow development.
The environment includes pre-installed nf-core tools, Nextflow, and a few other helpful utilities via a Docker container.

To get started, open the repository in [Codespaces](https://github.com/bhagesh-h/methylarray/codespaces).

### Documentation website

The documentation website ([bhagesh-h.github.io/methylarray](https://bhagesh-h.github.io/methylarray/)) is built from the Markdown files in this repository: `README.md` is the home page, `docs/*.md` the guides, and `CITATIONS.md` and `CHANGELOG.md` their own pages. Edit those files, not the website. `dev/docs/build_site.py` assembles them into `.site/src` and rewrites the links between them; `mkdocs.yml` sets the navigation and theme. Preview your changes in Docker with:

```bash
dev/nfdev docs          # build into .site/html
dev/nfdev docs serve    # preview at http://localhost:8000
```

Every pull request builds the website in strict mode, so a broken link or anchor fails the check. The website is published only from `main`, so it always describes the released pipeline. A new page needs an entry in `PAGES` of `dev/docs/build_site.py` and in `nav` of `mkdocs.yml`.

### Testing

Once you have made your changes, run the tests with nf-test locally.
For additional information, use the `--verbose` flag to view the Nextflow console log output.

```bash
# the tests that need no download (stubs, fixtures, input validation)
nf-test test --tag fast --profile +docker --verbose

# the pipeline test on the six public test arrays
nf-test test tests/default.nf.test --profile +docker --verbose
```

If you have added new functionality, ensure you update the test assertions in the `.nf.test` files.
Update the snapshots with the following command, and review the diff before committing:

```bash
nf-test test tests/default.nf.test --profile +docker --verbose --update-snapshot
```

When you create a pull request with changes, GitHub Actions will run automatic tests.
Pull requests are typically reviewed when these tests are passing.

Two types of tests are typically run:

#### Lint tests

nf-core has a [set of guidelines](https://nf-co.re/docs/specifications/overview) which all pipelines must follow.
To enforce these, run linting with nf-core/tools:

```bash
nf-core pipelines lint <pipeline_directory>
```

If you encounter failures or warnings, follow the linked documentation printed to screen.
For more information about linting tests, see [nf-core/tools API documentation](https://nf-co.re/docs/nf-core-tools/api_reference/latest/pipeline_lint_tests/actions_awsfulltest).

#### Pipeline tests

Each nf-core pipeline should be set up with a minimal set of test data.
GitHub Actions runs the pipeline on this data to ensure it runs through and exits successfully.
If there are any failures then the automated tests fail.
Pull requests into `dev` run the `fast` tests with Docker and the minimum required Nextflow version. Pull requests into `main` and releases run every test with Docker, conda and Singularity, on the minimum required and the latest Nextflow version.

### Patch release

> [!WARNING]
> Only in the unlikely event of a release that contains a critical bug.

- [ ] Create a new branch `patch` on your fork based on `upstream/main` or `upstream/master`.
- [ ] Fix the bug and use nf-core/tools to bump the version to the next semantic version, for example, `1.2.3` → `1.2.4`.
- [ ] Open a Pull Request from `patch` directly to `main`/`master` with the changes.

### Pipeline contribution conventions

nf-core semi-standardises how you write code and other contributions to make the bhagesh-h/methylarray code and processing logic more understandable for new contributors and to ensure quality.

#### Add a new pipeline step

To contribute a new step to the pipeline, follow the general nf-core coding procedure.
Please also refer to the [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines):

- [ ] Define the corresponding [input channel](#channel-naming-schemes) into your new process from the expected previous process channel.
- [ ] Install a module with nf-core/tools, or write a local module (see [default processes resource requirements](#default-processes-resource-requirements)), and call it from a subworkflow in `subworkflows/local/`.
- [ ] Define the output channel if needed. Emit `versions.yml` on the `versions` topic and mix relevant files into `ch_multiqc_files`.
- [ ] Add new or updated parameters to `nextflow.config` with a [default value](#default-parameter-values).
- [ ] Add new or updated parameters and relevant help text to `nextflow_schema.json` with [nf-core/tools](#default-parameter-values).
- [ ] Add validation for relevant parameters to the `PIPELINE_INITIALISATION` workflow in `subworkflows/local/utils_nfcore_methylarray_pipeline/main.nf`.
- [ ] Perform local tests to validate that the new code works as expected.
  - [ ] If applicable, add a new test in the `tests` directory.
- [ ] Update `docs/usage.md`, `docs/output.md`, `CHANGELOG.md` and `CITATIONS.md` as appropriate.
- [ ] [Lint](#lint-tests) the code with nf-core/tools.
- [ ] Update the workflow diagram (`assets/metro_map.mmd`, rendered to `docs/images/metro_map.svg`) as necessary.
- [ ] Update MultiQC config `assets/multiqc_config.yml` so relevant suffixes, file name cleanup, and module plots are in the appropriate order.
- [ ] If applicable, create a [MultiQC](https://seqera.io/multiqc/) module.
- [ ] Add a description of the output files and, if relevant, images from the MultiQC report to `docs/output.md`.

To update the minimum required Nextflow version, see the [Nextflow version bumping](#nextflow-version-bumping) section below. For more information about pipeline contributions, see [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines).

#### Channel naming schemes

Use the following naming schemes for channels to make the channel flow easier to understand:

- Initial process channel: `ch_output_from_<process>`
- Intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

#### Default parameter values

Parameters should be initialised and defined with default values within the `params` scope in `nextflow.config`.
They should also be documented in the pipeline JSON schema.

To update `nextflow_schema.json`, run:

```bash
nf-core pipelines schema build
```

The schema builder interface that loads in your browser should automatically update the defaults in the parameter documentation.

#### Default processes resource requirements

If you write a local module, specify a default set of resource requirements for the process.

Sensible defaults for process resource requirements (CPUs, memory, time) should be defined in `conf/base.config`.
Specify these with generic `withLabel:` selectors, so they can be shared across multiple processes and steps of the pipeline.

nf-core provides a set of standard labels that you should follow where possible, as seen in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/conf/base.config).
These labels define resource defaults for single-core processes, modules that require a GPU, and different levels of multi-core configurations with increasing memory requirements.

Values assigned within these labels can be dynamically passed to a tool using the `${task.cpus}` and `${task.memory}` Nextflow variables in the `script:` block of a module (see an example in the [modules repository](https://github.com/nf-core/modules/blob/bd1b6a40f55933d94b8c9ca94ec8c1ea0eaf4b82/modules/nf-core/samtools/bam2fq/main.nf#L30)).

#### Nextflow version bumping

If you use a new feature from core Nextflow, bump the minimum required Nextflow version in the pipeline with:

```bash
nf-core pipelines bump-version --nextflow . <min_nf_version>
```

#### Images and figures guidelines

If you update images or graphics, follow the nf-core [style guidelines](https://nf-co.re/docs/community/brand/workflow-schematics).

## Pipeline specific contribution guidelines

See the [developer guide](developer.md): running and testing the pipeline in Docker with `dev/nfdev`, the module structure, building task images with Wave, and the nf-test and pull request workflow.
