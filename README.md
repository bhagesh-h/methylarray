# methylarray

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A525.10.4-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c.svg)](https://apptainer.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Introduction

**methylarray** is a Nextflow pipeline for Illumina Infinium DNA methylation arrays. It starts from raw IDAT files and produces:

- quality-controlled and normalised methylation matrices
- sample and probe QC reports
- blood cell-type estimates
- DNA methylation biomarkers (epigenetic clocks and DNAm scores)
- epigenome-wide association results
- differentially methylated regions
- pathway and chromatin enrichment
- a MultiQC report and an HTML analysis report

The pipeline is built for the **Infinium MethylationEPIC v2.0 BeadChip** and handles the specific features of that array: replicate probes, the GRCh38 manifest, the `nv` and `rs` probe classes, and the updated cross-hybridisation annotation.

It is written in Nextflow DSL2 with strict syntax and built from the nf-core template. Every step runs in a container with pinned software versions. The only host requirement is Docker.

> **Status:** under active development. Version 1.0.0 has not been released. [CHANGELOG.md](CHANGELOG.md) lists what is implemented.

## Pipeline summary

1. Validate the samplesheet, phenotype table and model definitions against JSON schemas (nf-schema).
2. Prepare reference resources once and cache them for offline runs: the SeSAMe data cache, the EPIC v2 manifest and probe masks.
3. Read and preprocess each sample independently with SeSAMe:
   - bead-count masking
   - quality masking
   - Infinium I channel inference
   - non-linear dye-bias correction
   - pOOBAH detection masking
   - noob background correction
4. Assemble all samples into an HDF5-backed `SummarizedExperiment`.
5. Run sample QC:
   - detection rate and signal intensity
   - Illumina control-probe metrics (ewastools)
   - predicted versus reported sex
   - SNP-based identity and duplicate detection
   - contamination score
   - outlier detection
   - technical-replicate concordance

   Every exclusion is recorded with a reason.
6. Filter probes, in order:
   - detection failure
   - low bead count
   - cross-hybridising, SNP-affected, `nv`, `rs` and non-CpG probes
   - sex chromosomes (optional)

   EPIC v2 replicate probes are then collapsed to one value per CpG. A probe attrition table documents each step.
7. Estimate the proportions of 12 leukocyte subtypes with EpiDISH and the cent12CT reference.
8. Compute epigenetic clocks (Horvath, Hannum, skin and blood, PhenoAge, PC clocks, DunedinPACE) and DNAm scores (CRP, IL-6, smoking, alcohol, BMI), reporting the probe coverage of each on the array.
9. Explore technical and biological variation:
   - principal components and their associations with covariates
   - control-probe principal components
   - variance partitioning
   - a design-balance check for every model
10. Run epigenome-wide association studies:
    - limma on M-values, with effect sizes reported as differences in beta values
    - bias and inflation control with bacon
    - FDR and EPIC-specific significance thresholds
    - lookup of known associations
11. Call differentially methylated regions with DMRcate and dmrff, and report regions supported by both methods.
12. Test for enrichment with missMethyl (GO, KEGG, MSigDB Hallmark) and SeSAMe KnowYourCG.
13. Build a MultiQC report and an HTML report with an automatically generated methods section.

## Requirements

- Docker Engine 24 or later. Nextflow, Java and every analysis tool run inside containers, so nothing else is installed on the host.
- Alternatively, Nextflow 25.10.4 or later with Docker, Podman, Singularity, Apptainer or conda.
- Memory:
  - per-sample tasks need about 4 GB
  - whole-cohort steps need about 32 GB for 300 EPIC v2 samples
- Internet access on the first run, to fetch containers and reference resources. Later runs can be offline.

## Quick start

1. Test the installation with the bundled test profile:

   ```bash
   scripts/nf-docker run -profile test,docker --outdir /absolute/path/to/results
   ```

   `scripts/nf-docker` starts the official Nextflow container. It mounts the Docker socket and the required paths at identical locations, then runs the pipeline. Tasks start as sibling containers on the host Docker daemon.

2. Prepare a samplesheet (see [Input](#input)).

3. Run the pipeline:

   ```bash
   scripts/nf-docker run \
       -profile docker \
       --input samplesheet.csv \
       --phenotypes phenotypes.tsv \
       --models models.yaml \
       --outdir /absolute/path/to/results
   ```

If Nextflow is already installed on your system, the equivalent command is:

```bash
nextflow run bhagesh-h/methylarray -r <version> -profile docker --input samplesheet.csv --outdir results
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

| Column | Required | Description |
| --- | --- | --- |
| `sample` | yes | Unique sample identifier |
| `idat_red` | yes | Red channel IDAT (`_Red.idat` or `_Red.idat.gz`) |
| `idat_grn` | yes | Green channel IDAT with the same prefix as the red file |
| `sentrix_id` | no | Chip barcode. Derived from the file name if omitted |
| `sentrix_position` | no | Position such as `R01C01`. Derived from the file name if omitted |
| `plate` | no | Processing plate, used as a batch variable |
| `well` | no | Plate well |
| `batch` | no | Any further processing batch |
| `reported_sex` | no | `M` or `F`, compared with the sex predicted from the array |
| `reported_age` | no | Chronological age, used for clock-based checks |
| `replicate_group` | no | Shared value for technical replicates |
| `exclude`, `exclude_reason` | no | Manual exclusion with a recorded reason |

### Phenotypes

An optional CSV or TSV file with a `sample` column and any number of analysis variables. Phenotype values are never written to logs.

### Models

An optional YAML file listing the association models to run. Each model names a trait, its type, an optional transformation and the covariates. Shortcut tokens expand to derived covariates: `@cells` (cell composition), `@technical` (plate and Sentrix row), `@ctrlpcs` (control-probe components), `@sex` (predicted sex) and `@smoking` (DNAm smoking score).

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

Without a models file the pipeline stops after QC, preprocessing, cell composition, biomarkers and exploratory analysis.

## Usage notes

- `--qc_only` stops after sample QC and reporting. Use it to review the QC report and write a sample exclusion file, then pass that file with `--sample_exclusions` in the full run.
- Add `-resume` to reuse completed tasks after a change or failure.
- `--sesame_cache` reuses a reference cache from an earlier run, which allows fully offline execution.
- Batch effects are modelled as covariates by default. ComBat is available through `--combat_batch` but only runs for models that pass the design-balance check, because ComBat on unbalanced designs can create false positives.
- Run `nextflow run bhagesh-h/methylarray --help` for all parameters. The full description of every option is in [docs/usage.md](docs/usage.md), and the parameter schema is [nextflow_schema.json](nextflow_schema.json).

## Execution modes

| Mode | When to use | How it works |
| --- | --- | --- |
| Docker outside of Docker (default) | Workstations and servers where you can use Docker | The Nextflow container drives the host Docker daemon through the mounted socket. Each task runs in its own pinned container |
| Docker in Docker | CI and isolated environments | A separate Docker daemon runs in a sidecar container. The pipeline runs against that daemon only |
| Conda in Docker | Hosts where the Docker socket cannot be shared | All tasks run inside one container, using pinned conda environments built with micromamba |
| HPC | Clusters with Singularity or Apptainer | Standard Nextflow installation with `-profile singularity` or `-profile apptainer` and an institutional config |

All paths given to the pipeline must be absolute. In the Docker-based modes, the wrapper mounts them at the same location inside the containers so that task containers resolve them correctly.

## Outputs

```text
results/
├── pipeline_info/      execution reports, software versions, validated inputs
├── qc/                 sample QC, control metrics, sex and identity checks, probe attrition
├── preprocessing/      HDF5-backed SummarizedExperiment objects and exported matrices
├── cell_composition/   cell-type proportions and plots
├── biomarkers/         epigenetic clocks and DNAm scores
├── exploratory/        PCA, variance partitioning, control-probe components, design checks
├── ewas/               summary statistics, top hits and plots per model, cross-model summary
├── dmr/                DMRcate, dmrff and consensus regions per model
├── enrichment/         gene set and CpG-centric enrichment per model
├── multiqc/            MultiQC report
└── report/             HTML analysis report with methods text
```

Every file is described in [docs/output.md](docs/output.md).

## Testing

The pipeline is tested with [nf-test](https://www.nf-test.com) at the module, subworkflow and pipeline level. The `test` profile uses a small public EPIC v2 whole-blood dataset with synthetic phenotypes and runs on a 4 CPU, 15 GB machine. The `test_full` profile runs the complete public dataset.

## Reproducibility

- Each module pins its software in an `environment.yml` and a matching container image.
- Software versions are collected from every task and reported in MultiQC and in `pipeline_info/`.
- Random seeds are fixed, and inputs are sorted before aggregation.
- The exact command, parameters, container digests and pipeline revision are written with every run.

## Data privacy

The pipeline processes data locally and does not upload samples, phenotypes or results anywhere. Container builds only use environment files. When working with human cohort data, keep IDATs, samplesheets with identifiers, phenotype tables and results out of version control.

## Credits

methylarray is developed and maintained by [bhagesh-h](https://github.com/bhagesh-h).

The design draws on the nf-core/methylarray development pipeline and its contributors, the Garvan Institute EPIC v2 workflow, and the EPIC v2 preprocessing work of the Kobor lab.

## Contributions and support

Bug reports and feature requests are welcome through GitHub issues. See [docs/CONTRIBUTING.md](docs/CONTRIBUTING.md) before opening a pull request.

## Citations

All tools, methods, reference data and frameworks used by the pipeline are listed in [CITATIONS.md](CITATIONS.md). A Zenodo DOI for the pipeline will be added at the first release.

## License

methylarray is released under the [MIT License](LICENSE).
