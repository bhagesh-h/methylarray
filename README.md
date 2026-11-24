# methylarray

[![GitHub Actions CI Status](https://github.com/bhagesh-h/methylarray/actions/workflows/nf-test.yml/badge.svg)](https://github.com/bhagesh-h/methylarray/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/bhagesh-h/methylarray/actions/workflows/linting.yml/badge.svg)](https://github.com/bhagesh-h/methylarray/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.04.2-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.2.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.2.2)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c.svg?labelColor=000000)](https://apptainer.org/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24399752.svg)](https://doi.org/10.5281/zenodo.24399752)

Documentation website: [bhagesh-h.github.io/methylarray](https://bhagesh-h.github.io/methylarray/) (built from the files in this repository, for the released version on `main`).

## Introduction

**methylarray** is a Nextflow pipeline for Illumina Infinium DNA methylation arrays. It starts from raw IDAT files and produces:

- preprocessed beta values, M-values and detection p-values
- sample QC with a status and the reasons for every array
- probe filtering with an attrition table
- blood cell-type estimates
- epigenetic clocks and DNA methylation scores
- exploratory principal component analyses
- epigenome-wide association results, differentially methylated regions and gene set enrichment, when association models are given
- optional GEO submission files
- a MultiQC report and an HTML analysis report

The pipeline is built for the **Infinium MethylationEPIC v2.0 BeadChip** and handles the specific features of that array: replicate probes, the GRCh38 manifest, the `nv` and `rs` probe classes, renamed CpGs and the updated cross-hybridisation annotation. It is written in Nextflow DSL2, built from the nf-core template, and every step runs with pinned software through Docker, Singularity, Apptainer, Podman or conda.

![methylarray workflow](docs/images/metro_map.svg)

> **Status:** released, the latest version is 1.0.2. [CHANGELOG.md](CHANGELOG.md) lists what each version contains.

## Pipeline summary

1. Validate the samplesheet, phenotype table and model definitions against JSON schemas ([nf-schema](https://github.com/nextflow-io/nf-schema)) and against each other.
2. Preprocess the arrays with one of four routes (see [Preprocessing options](#preprocessing-options)): [minfi](https://bioconductor.org/packages/minfi/) (default), [SeSAMe](https://bioconductor.org/packages/sesame/), [ENmix](https://bioconductor.org/packages/ENmix/) or a combined route. No array or probe is removed here; failed values are masked.
3. Sample QC, with [ewastools](https://github.com/hhhh5/ewastools) for all checks except detection:
   - fraction of undetected probes, from the preprocessing route
   - Illumina BeadArray Controls Reporter metrics
   - predicted against reported sex
   - SNP-based contamination score
   - SNP genotype identity checks against `replicate_group`

   Every array gets a status (`PASS`, `FLAG`, `FAIL` or `EXCLUDED`) with its reasons. `PASS` and `FLAG` arrays go on to the analysis.

4. Annotate the probes from the Peters et al. (2024) extended EPIC v2 manifest and the Zhou lab InfiniumAnnotation masks, then remove, in order:
   - `rs`, `nv`, control and non-CG probes, and probes without an hg38 position
   - probes on chrX and chrY (kept with `--keep_sex_chr`)
   - probes failing detection or the bead count in too many arrays
   - SNP-affected probes (global or population-specific masks), the Zhou general mask and cross-hybridising probes

   Probes with manifest position mismatches (Peters et al. 2024 `MismatchPos`), poorly titrating probes (Zhou) and, from user-supplied lists, Illumina's flagged and mapping-inaccuracy probes are annotated and counted, never removed. EPIC v2 replicate probes are then averaged per CpG. An attrition table documents each step.

5. Estimate leukocyte fractions with [EpiDISH](https://bioconductor.org/packages/EpiDISH/) (12-type `cent12CT` reference, with the 7-type `centDHSbloodDMC` reference for comparison).
6. Compute epigenetic clocks with [dnaMethyAge](https://github.com/yiluyucheng/dnaMethyAge) (Horvath, Hannum, PhenoAge, skin and blood, PC clocks, DNAmTL, epiTOC2 and others) and [DunedinPACE](https://github.com/danbelsky/DunedinPACE), with age acceleration and the probe coverage of every clock.
7. Compute DNA methylation scores for smoking, alcohol and BMI (McCartney et al. 2018) and IL-6 (Stevenson et al. 2021), with their probe coverage.
8. Exploratory analysis:
   - control-probe principal components, used as model covariates
   - genetic ancestry components from the SNP probes (EpiAnceR+ approach), used as model covariates with `@ancestry`
   - with a user-supplied meQTL table (`--meqtl`): the clock weight on CpGs under genetic control, genotype-robust clocks and the association hits with a cis-meQTL
   - principal components of the most variable CpGs (number estimated with [isva](https://cran.r-project.org/package=isva)) and their associations with technical, biological and phenotype variables
   - PCA outlier flags (report only, no array is removed)
   - variance partition of the most variable CpGs with [variancePartition](https://bioconductor.org/packages/variancePartition/): the share explained by plate, BeadChip, row, sex, age and cell types
9. When a models file is given, for every model:
   - [limma](https://bioconductor.org/packages/limma/) on M-values, with effect sizes reported as differences in beta values, and a design check (confounding with plate and Sentrix row, collinearity)
   - bias and inflation control with [bacon](https://bioconductor.org/packages/bacon/)
   - a study-wide threshold across models from the effective number of independent tests (Meff)
   - differentially methylated regions with [DMRcate](https://bioconductor.org/packages/DMRcate/) and [dmrff](https://github.com/perishky/dmrff), and their consensus
   - gene set enrichment with [missMethyl](https://bioconductor.org/packages/missMethyl/): GO, KEGG (downloaded at run time) and MSigDB collections
   - genomic context enrichment with [KnowYourCG](https://bioconductor.org/packages/knowYourCG/): chromatin states, transcription factor binding sites, histone marks and other features
10. Optionally write GEO submission files (`--export_geo`).
11. Build a [MultiQC](https://seqera.io/multiqc/) report and an HTML analysis report with [Quarto](https://quarto.org/), including a methods text with references.

## Preprocessing options

Select one with `--preprocessing`.

| Option            | Route                                                                                                                           | When to use it                                                                                    |
| ----------------- | ------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------- |
| `minfi` (default) | minfi import, minfi detection p-values, Noob (Funnorm, quantile or SWAN with `--minfi_norm`)                                    | Established route; comparable with most published array studies                                   |
| `sesame`          | SeSAMe per-array preprocessing (`prepSesame`, default `QCDPB`): quality mask, channel inference, dye bias, pOOBAH, noob         | Each array is processed on its own, so results do not depend on which other arrays are in the run |
| `enmix`           | minfi import, ENmix out-of-band background with RELIC dye-bias correction, optional quantile normalisation, RCP probe-type step | ENmix background and probe-type correction                                                        |
| `combined`        | the `enmix` route, with the SeSAMe quality mask and pOOBAH detection mask applied to it                                         | ENmix correction with SeSAMe masking                                                              |

The minfi-based routes read the arrays in chunks to limit memory use.

## Requirements

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) 26.04.2 or later (earlier releases can wrongly report remote HTTP or HTTPS input files as missing).
- One of [Docker](https://docs.docker.com/engine/installation/), [Singularity](https://www.sylabs.io/guides/3.0/user-guide/), [Apptainer](https://apptainer.org/docs/user/latest/), [Podman](https://podman.io/) or [conda](https://conda.io/miniconda.html). Containers are recommended.
- Memory: per-array SeSAMe tasks need about 6 GB. Steps that read all arrays together use the `process_high` label, which requests 12 CPUs and 72 GB by default. On a smaller machine, cap the requests with `process.resourceLimits` in a config passed with `-c`.
- Internet access on the first run, to fetch containers and reference resources. Later runs can be offline, see [Offline runs](docs/usage.md#offline-runs).

## Quick start

1. Test the pipeline on the bundled test profile (six public EPIC v2 arrays):

   ```bash
   nextflow run bhagesh-h/methylarray -profile test,docker --outdir results
   ```

   Replace `docker` with `singularity`, `apptainer`, `podman` or `conda` to match your system. Institutional profiles from [nf-core/configs](https://github.com/nf-core/configs) also work.

2. Prepare a samplesheet and, for association analyses, a phenotype table and a models file (see [Input](#input)).

3. Run the pipeline:

   ```bash
   nextflow run bhagesh-h/methylarray \
       -r <version> \
       -profile docker \
       --input samplesheet.csv \
       --phenotypes phenotypes.tsv \
       --models models.yml \
       --outdir results
   ```

Pipeline parameters must be passed on the command line or in a `-params-file`. A `-c` config file should only change infrastructure settings such as resources and executors.

## Input

### Samplesheet

A CSV, TSV, JSON or YAML file with one row per array. Required columns are `sample`, `idat_red` and `idat_grn`.

```csv
sample,idat_red,idat_grn,plate,reported_sex
S001,/data/idat/200000000001/200000000001_R01C01_Red.idat,/data/idat/200000000001/200000000001_R01C01_Grn.idat,P01,F
S002,/data/idat/200000000001/200000000001_R02C01_Red.idat,/data/idat/200000000001/200000000001_R02C01_Grn.idat,P01,M
```

| Column                      | Required | Description                                                                |
| --------------------------- | -------- | -------------------------------------------------------------------------- |
| `sample`                    | yes      | Unique sample identifier                                                   |
| `idat_red`                  | yes      | Red channel IDAT (`_Red.idat` or `_Red.idat.gz`), local path or URL        |
| `idat_grn`                  | yes      | Green channel IDAT with the same prefix as the red file                    |
| `sentrix_id`                | no       | Chip barcode. Derived from the file name if omitted                        |
| `sentrix_position`          | no       | Position such as `R01C01`. Derived from the file name if omitted           |
| `plate`                     | no       | Processing plate, used as a technical covariate                            |
| `well`                      | no       | Plate well                                                                 |
| `batch`                     | no       | Any further processing batch                                               |
| `reported_sex`              | no       | `M` or `F`, compared with the sex predicted from the array                 |
| `reported_age`              | no       | Age in years, used for age acceleration and the PCA associations           |
| `replicate_group`           | no       | Shared value for technical replicates, used by the genotype identity check |
| `exclude`, `exclude_reason` | no       | Manual exclusion; an excluded array needs a reason and is never processed  |

### Phenotypes

An optional CSV or TSV file with a `sample` column and any number of analysis variables. Phenotype values are never written to logs.

### Models

An optional YAML file listing the association models to run. Each model names a trait, its type, an optional transformation and the covariates. Shortcut tokens expand to derived covariates: `@cells` (five cell lineages, granulocytes as reference), `@technical` (plate and Sentrix row), `@ctrlpcs` (control-probe components), `@ancestry` (genetic ancestry components from the SNP probes), `@sex` (predicted sex) and `@smoking` (DNA methylation smoking score).

```yaml
- name: bmi
  trait: bmi
  trait_type: continuous
  covariates: [age, "@sex", "@cells", "@technical", "@ctrlpcs", "@smoking"]
- name: crp
  trait: hscrp
  trait_type: continuous
  transform: log
  covariates: [age, "@sex", bmi, "@cells", "@technical", "@ctrlpcs"]
```

Without a models file the pipeline runs preprocessing, QC, probe filtering, cell composition, clocks, DNA methylation scores and the exploratory analysis, and no association models.

## Usage notes

- Review the sample QC (`qc/sample_qc.tsv` and the reports). To leave arrays out of a later run, set `exclude` and `exclude_reason` in the samplesheet.
- Add `-resume` to reuse completed tasks after a change or failure.
- Batch effects are handled as model covariates (`@technical`, `@ctrlpcs`). The pipeline does not apply ComBat or other batch correction to the data.
- Run `nextflow run bhagesh-h/methylarray --help` for all parameters. The full description of every option is in [docs/usage.md](docs/usage.md), and the parameter schema is [nextflow_schema.json](nextflow_schema.json).

## Outputs

```text
results/
├── preprocessing/      HDF5-backed SummarizedExperiment of the preprocessing route
├── qc/                 preprocessing QC, ewastools QC, sample QC status
├── probe_filter/       filtered data (CpG and probe level), attrition and exclusion tables
├── cell_composition/   cell-type fractions and reference coverage
├── clocks/             epigenetic clocks, age acceleration, coverage and standardised weights
├── dnam_scores/        DNA methylation scores, coverage and standardised weights
├── biomarkers/         clocks, scores and cell fractions against each other, the batches and the phenotypes
├── exploration/        control-probe PCs, PCA scores, variance, associations, outlier flags
├── ewas/               statistics, top CpGs and design checks per model, summary across models
├── dmr/                DMRcate, dmrff and consensus regions per model
├── enrichment/         gene set enrichment per model
├── geo/                GEO submission files (--export_geo)
├── references/         probe annotation and SeSAMe cache (--save_reference)
├── report/             HTML analysis report
├── multiqc/            MultiQC report
└── pipeline_info/      execution reports, parameters, software versions
```

Every file is described in [docs/output.md](docs/output.md).

## Testing

The pipeline is tested with [nf-test](https://www.nf-test.com) at the module, subworkflow and pipeline level. The `test` profile runs six public EPIC v2 blood arrays from GEO GSE286313, without phenotypes or models, within 4 CPUs and 15 GB of memory. The `test_full` profile runs a larger set of arrays from the same GEO series.

## Reproducibility

- Each module pins its software in an `environment.yml` and a matching container image.
- Packages that are not on conda (ewastools, dnaMethyAge, DunedinPACE, dmrff) are installed from source packages at fixed GitHub commits.
- Software versions are collected from every task and reported in MultiQC, the HTML report and `pipeline_info/`.
- Arrays are sorted by sample name before they are combined, so results do not depend on the order of the samplesheet or of task completion.
- The parameters of every run are written to `pipeline_info/`.

## Data privacy

The pipeline processes data where it runs and does not upload samples, phenotypes or results anywhere. When working with human data, keep IDATs, samplesheets with identifiers, phenotype tables and results out of version control.

## Credits

methylarray is developed and maintained by Bhagesh Hunakunti ([bhagesh-h](https://github.com/bhagesh-h), ORCID [0000-0002-5957-8005](https://orcid.org/0000-0002-5957-8005)).

The workflow diagram is drawn with [nf-metro](https://github.com/seqeralabs/nf-metro) from [assets/metro_map.mmd](assets/metro_map.mmd).

## Contributions and support

Bug reports and feature requests are welcome through GitHub issues. Development happens on the `dev` branch: open pull requests against `dev`, never against `main`. See [docs/CONTRIBUTING.md](docs/CONTRIBUTING.md) and the [developer guide](docs/developer.md) before opening a pull request.

## Citations

To cite methylarray, use its Zenodo DOI, [10.5281/zenodo.24399752](https://doi.org/10.5281/zenodo.24399752), which always resolves to the latest version; version 1.0.2 is [10.5281/zenodo.24516701](https://doi.org/10.5281/zenodo.24516701), version 1.0.1 is [10.5281/zenodo.24458958](https://doi.org/10.5281/zenodo.24458958) and version 1.0.0 is [10.5281/zenodo.24399753](https://doi.org/10.5281/zenodo.24399753). The full citation is in [CITATION.cff](CITATION.cff), also shown by the "Cite this repository" button on the GitHub page. Please also cite the tools, methods, reference data and frameworks the pipeline uses, listed in [CITATIONS.md](CITATIONS.md).

## License

methylarray is released under the [MIT License](LICENSE).
