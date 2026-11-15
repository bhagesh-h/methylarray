# bhagesh-h/methylarray: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The pipeline takes up to three input files:

| Parameter      | Required | Content                                                          |
| -------------- | -------- | ---------------------------------------------------------------- |
| `--input`      | yes      | Samplesheet: one row per array with its red and green IDAT files |
| `--phenotypes` | no       | Phenotype table: one row per sample with traits and covariates   |
| `--models`     | no       | Model definitions: the association models to run, in YAML        |

All three are validated against JSON schemas before any task starts, and checked against each other: every sample, column and model name must resolve, otherwise the run stops with a message that lists every problem.

## Samplesheet input

```bash
--input '[path to samplesheet file]'
```

The samplesheet is a CSV, TSV, YAML or JSON file with one row per array. The IDAT files can be local paths or URLs, and may be gzipped.

```csv title="samplesheet.csv"
sample,idat_red,idat_grn,plate,reported_sex,reported_age
S001,/data/idat/207059430159/207059430159_R01C01_Red.idat,/data/idat/207059430159/207059430159_R01C01_Grn.idat,P01,F,52
S002,/data/idat/207059430159/207059430159_R02C01_Red.idat,/data/idat/207059430159/207059430159_R02C01_Grn.idat,P01,M,47
```

| Column             | Required    | Description                                                                                                      |
| ------------------ | ----------- | ---------------------------------------------------------------------------------------------------------------- |
| `sample`           | yes         | Unique sample name. Letters, digits, `.`, `_` and `-` only                                                       |
| `idat_red`         | yes         | Red channel IDAT, ending in `_Red.idat` or `_Red.idat.gz`                                                        |
| `idat_grn`         | yes         | Green channel IDAT with the same prefix as the red file                                                          |
| `sentrix_id`       | no          | Chip barcode (10 to 14 digits). Taken from the file name when omitted; an error if it contradicts the file name  |
| `sentrix_position` | no          | Array position such as `R01C01`. Taken from the file name when omitted; an error if it contradicts the file name |
| `plate`            | no          | Processing plate, used as a technical covariate                                                                  |
| `well`             | no          | Plate well, such as `A1` or `H12`                                                                                |
| `batch`            | no          | Any further processing batch                                                                                     |
| `reported_sex`     | no          | `M` or `F`. Compared with the sex predicted from the array                                                       |
| `reported_age`     | no          | Age in years                                                                                                     |
| `replicate_group`  | no          | Technical replicates share a value; used for replicate concordance                                               |
| `exclude`          | no          | `true` to exclude an array before analysis                                                                       |
| `exclude_reason`   | conditional | Required when `exclude` is `true`; recorded in the QC report                                                     |

Each IDAT file may appear in only one row. The pipeline derives `sentrix_row` (such as `R01`) and `sentrix_col` (such as `C01`) from the array position.

An [example samplesheet](../assets/samplesheet.csv) is provided with the pipeline.

## Phenotype table

```bash
--phenotypes '[path to phenotype table]'
```

A CSV or TSV file with a `sample` column that matches the samplesheet and any number of further columns. Each further column is a trait or covariate that models can name. Missing values are written as `NA` or left empty.

```tsv title="phenotypes.tsv"
sample	age	bmi	smoker
S001	52	24.1	no
S002	47	NA	yes
```

- Column names may not reuse samplesheet column names (for example `plate` or `reported_sex`) or the derived covariate names listed below.
- Samples in the phenotype table that are not in the samplesheet are ignored with a warning.
- Samples without a phenotype row are still processed, and are left out of the models that use phenotype columns.
- Only column names and sample counts are written to the logs, never phenotype values.

## Model definitions

```bash
--models '[path to models file]'
```

A YAML list with one entry per association model. Each model tests one trait against DNA methylation at every CpG, adjusted for its covariates.

```yaml title="models.yml"
- name: bmi
  trait: bmi
  trait_type: continuous
  transform: log
  covariates: [age, "@sex", "@cells", "@technical", "@ctrlpcs"]
- name: smoking
  trait: smoker
  trait_type: binary
  reference: "no"
  covariates: [age, "@sex", "@cells", sentrix_row, "@ctrlpcs"]
  block: plate
```

| Field             | Required    | Description                                                                                                                             |
| ----------------- | ----------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| `name`            | yes         | Unique model name; lower-case letters, digits and `_`. Used in output file names                                                        |
| `trait`           | yes         | Column tested for association                                                                                                           |
| `trait_type`      | yes         | `continuous`, `binary` or `categorical`                                                                                                 |
| `transform`       | no          | Transformation of a continuous trait: `none` (default), `log`, `log2`, `log10`, `log1p`, `rint` (rank-based inverse normal) or `zscore` |
| `reference`       | conditional | Reference level, required for binary and categorical traits                                                                             |
| `contrasts`       | no          | Categorical traits only: a list of comparisons written `levelA-levelB`                                                                  |
| `covariates`      | no          | Columns or tokens to adjust for                                                                                                         |
| `block`           | no          | Column modelled as a random effect (limma `duplicateCorrelation`), for example `plate`                                                  |
| `subset`          | no          | Restrict the model to some samples: `{column: sex, equals: F}` or `{column: group, in: [a, b]}`                                         |
| `min_n`           | no          | Minimum number of complete cases (default 50, at least 10)                                                                              |
| `skip_dmr`        | no          | `true` to skip region analysis for this model                                                                                           |
| `skip_enrichment` | no          | `true` to skip enrichment for this model                                                                                                |
| `description`     | no          | Free text shown in the report                                                                                                           |

Columns can come from the samplesheet, the phenotype table or the covariates the pipeline derives from the methylation data:

| Name or token   | Meaning                                                      |
| --------------- | ------------------------------------------------------------ |
| `predicted_sex` | Sex predicted from the sex chromosomes                       |
| `dnam_smoking`  | DNA methylation smoking score (not available yet)            |
| `cc_<celltype>` | Estimated proportion of one cell type, for example `cc_Mono` |
| `ctrlpc<k>`     | Principal component `k` of the control-probe intensities     |
| `@cells`        | The configured cell-composition covariates                   |
| `@ctrlpcs`      | The selected control-probe principal components              |
| `@technical`    | `plate` and `sentrix_row`                                    |
| `@sex`          | `predicted_sex`                                              |
| `@smoking`      | `dnam_smoking` (not available yet)                           |

A model stops the run before any task starts when it names a column found in none of these sources, lists its trait among its covariates (directly or through a token), or uses its `block` column as a covariate as well. Without `--models` the pipeline runs preprocessing, QC, cell composition, biomarkers and the exploratory analysis, and no association models.

## Preprocessing

Select the method with `--preprocessing`. Every method writes the same SummarizedExperiment (see [docs/output.md](output.md)), so later steps do not depend on the choice.

| Method            | Steps                                                                                                                                                      | Masking                                              |
| ----------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------- |
| `minfi` (default) | minfi import, `preprocessNoob` or `--minfi_norm funnorm`, `quantile` or `swan`                                                                             | minfi detection p-value, bead count                  |
| `sesame`          | SeSAMe `prepSesame` per array (`--sesame_prep`, default `QCDPB`)                                                                                           | SeSAMe quality mask, pOOBAH detection, bead count    |
| `enmix`           | minfi import, ENmix out-of-band background with RELIC dye-bias correction, optional quantile normalisation (`--enmix_quantile`), RCP probe-type correction | ENmix detection p-value, bead count                  |
| `combined`        | the `enmix` steps, with the SeSAMe result computed on the same arrays                                                                                      | SeSAMe quality mask and pOOBAH detection, bead count |

SeSAMe processes each array on its own, so its results do not change when arrays are added to or removed from a run. The other methods read all arrays together. minfi detection p-values mask far fewer values than SeSAMe, which also removes probes with known design problems; the `combined` method keeps the ENmix correction and applies the SeSAMe masks.

| Parameter          | Default  | Meaning                                                                                                   |
| ------------------ | -------- | --------------------------------------------------------------------------------------------------------- |
| `--platform`       | `EPICv2` | Array platform. The platform of every array is detected from its IDATs and must match                     |
| `--sesame_prep`    | `QCDPB`  | `prepSesame` steps in order: quality mask, channel inference, non-linear dye bias, pOOBAH detection, noob |
| `--minfi_norm`     | `noob`   | Normalisation of the `minfi` method: `noob`, `funnorm`, `quantile` or `swan`                              |
| `--enmix_quantile` | `false`  | Add ENmix quantile normalisation to the `enmix` and `combined` methods                                    |
| `--detection_pval` | `0.05`   | Detection p-value threshold; values above it are set to `NA`                                              |
| `--min_beads`      | `3`      | Values measured by fewer beads are set to `NA`                                                            |

The masked beta values are used for association testing. The unmasked beta values are kept as a separate assay for steps that need complete data, such as cell-type deconvolution and epigenetic clocks. No array or probe is removed during preprocessing; ENmix quality flags are reported in the QC table (`enmix_low_quality`) instead of being applied.

### Offline runs

The `sesame` and `combined` methods need probe manifests, masks and platform signatures from Bioconductor's ExperimentHub (about 45 MB for EPIC v2). By default the first task of each run downloads them. The `minfi` and `enmix` methods use the manifest and annotation packages in their container and need no download. Probe filtering reads the probe manifest and mask from the network. To run without network access, save them once and pass them to later runs:

```bash
# once, with network access
nextflow run bhagesh-h/methylarray -profile docker --input samplesheet.csv --outdir results --save_reference

# later, offline
nextflow run bhagesh-h/methylarray -profile docker --input samplesheet.csv --outdir results2 --sesame_cache results/references/sesame_cache --probe_annotation results/references/probe_annotation.tsv.gz
```

### Remote IDAT files

IDAT paths may be HTTPS URLs, for example GEO supplementary files. Some servers, including NCBI GEO, refuse many parallel downloads. If staging fails with HTTP 503, limit the number of parallel transfers in a custom config passed with `-c`:

```groovy
threadPool.FileTransfer.maxThreads = 4
filePorter.maxRetries = 10
```

Nextflow prints an "Unrecognized config option" warning for these two settings, but applies them.

## Sample QC

After preprocessing, every array gets a QC status (`PASS`, `FLAG`, `FAIL` or `EXCLUDED`) with the reasons, in `qc/sample_qc.tsv`. No array is removed without a recorded reason. The checks:

- **Detection**: the fraction of probes with a detection p-value above `--detection_pval`. More than `--max_undetected_fraction` (default 10%) fails the array.
- **Sex**: the sex predicted by the preprocessing method and by ewastools from the normalised X and Y intensities, against `reported_sex`. If every prediction contradicts the reported sex, the array fails, since this usually means a sample swap. If the predictions disagree, or the X and Y intensities fall between the male and female clusters (possible sex chromosome aneuploidy or contamination), the array is flagged. The ewastools prediction needs at least one reported female and one reported male.
- **Control probes**: the 21 Illumina BeadArray Controls Reporter metrics (bisulfite conversion, extension, hybridisation, staining, specificity, target removal, non-polymorphic) against Illumina's thresholds. Arrays below a threshold are flagged, not failed: the thresholds are guidance for troubleshooting the assay, and arrays with good data can miss a single background ratio.
- **Contamination**: the ewastools SNP outlier score, the average log odds that the SNP probe beta values fit none of the three genotypes. Arrays above `--snp_outlier_threshold` (default -4, as recommended by the ewastools authors) fail.
- **Identity**: SNP genotypes of all arrays below the contamination threshold are compared. Arrays of one `replicate_group` with different genotypes, and arrays of different donors with the same genotypes (agreement above 0.9), are flagged.

| Parameter                   | Default                | Meaning                                                                |
| --------------------------- | ---------------------- | ---------------------------------------------------------------------- |
| `--max_undetected_fraction` | `0.1`                  | Largest fraction of undetected probes for a passing array              |
| `--snp_outlier_threshold`   | `-4`                   | Arrays with a higher SNP outlier score fail as possibly contaminated   |
| `--skip_ewastools`          | `false`                | Skip the control-probe, X/Y sex, contamination and identity checks     |
| `--ewastools_source`        | ewastools 1.8.2 source | ewastools source package (URL or local path), installed by the QC task |

ewastools is not available from conda or CRAN. The QC task installs it from its source package (a fixed GitHub commit) into the task directory, so it works with every profile. For runs without network access, download the file given by `--ewastools_source` once and pass its local path.

## Probe filtering

The arrays with QC status `PASS` or `FLAG` go on to probe filtering. Every removed probe is listed with its reasons in `probe_filter/` (see [docs/output.md](output.md#probe-filtering)).

| Parameter                    | Default      | Meaning                                                                                                     |
| ---------------------------- | ------------ | ----------------------------------------------------------------------------------------------------------- |
| `--probe_fail_fraction`      | `0.1`        | Remove probes with a detection p-value above `--detection_pval` in more than this fraction of arrays        |
| `--bead_fail_fraction`       | `0.05`       | Remove probes masked for another reason (bead count, quality mask) in more than this fraction of arrays     |
| `--snp_mask`                 | `common_1pt` | SNP mask: `common_1pt`, `common_5pt`, `afr`, `amr`, `eas`, `eur`, `sas` or `none`                           |
| `--crosshyb_filter`          | `wgbs`       | Cross-hybridising probes to remove: `wgbs` (WGBS evidence), `blat` (any predicted off-target hit) or `none` |
| `--keep_sex_chr`             | `false`      | Keep chrX and chrY probes                                                                                   |
| `--replicate_probe_strategy` | `mean`       | `mean`: one row per CpG, the mean of its replicate probes; `none`: one row per probe                        |
| `--probe_annotation`         |              | Saved `probe_annotation.tsv.gz`, used instead of building it                                                |

The probe annotation comes from two pinned sources: the Peters et al. 2024 extended EPIC v2 manifest (`--probe_manifest_source`, about 300 MB) and the Zhou lab InfiniumAnnotation mask (`--probe_mask_source`). Choose the SNP population that matches the study: the `common_*` masks use global allele frequencies. Values that failed detection or the bead count in only a few arrays stay in the data as `NA`.

Probes with manifest position mismatches (Peters et al. 2024 `MismatchPos`) and probes with poor titration (Zhou) are counted in the attrition table, not removed.

## Cell composition

With `--tissue blood` (default), leukocyte fractions are estimated with EpiDISH (robust partial correlations) from the filtered, CpG-level beta values without masking. The main reference is `cent12CT` (12 cell types, Luo et al. 2023); `centDHSbloodDMC` (7 types) is reported for comparison. Reference CpGs are matched through their EPIC v1 and 450K IDs as well as the CpG ID, because some CpGs have new names on EPIC v2. The run stops if less than `--min_reference_coverage` (default 90%) of a reference's CpGs are present. Use `--tissue other` for tissues without a blood reference.

## Epigenetic clocks

Clocks are computed with dnaMethyAge and DunedinPACE on the detection-masked beta values of the arrays that passed QC, averaged per CpG. The probe filter for association testing is not applied here, because each clock was trained on its own set of CpGs. Values that failed detection are set to the CpG's mean over the other arrays. CpGs that EPIC v2 does not measure take the reference means of dnaMethyAge, and DunedinPACE returns NA when less than `--dunedinpace_min_probes` of its CpGs are present. `clocks/clocks_coverage.tsv` reports the coverage of every clock.

| Parameter                  | Default                  | Meaning                                                   |
| -------------------------- | ------------------------ | --------------------------------------------------------- |
| `--clocks`                 | 16 clocks (see schema)   | Comma-separated dnaMethyAge clock names and `DunedinPACE` |
| `--dunedinpace_min_probes` | `0.8`                    | Minimum share of DunedinPACE CpGs in the data             |
| `--skip_clocks`            | `false`                  | Skip the clocks                                           |
| `--dnamethyage_source`     | dnaMethyAge 0.2.0 source | dnaMethyAge source package (URL or local path)            |
| `--dunedinpace_source`     | DunedinPACE 0.99.0       | DunedinPACE source package (URL or local path)            |

PCGrimAge needs the reported age and sex (`reported_age`, `reported_sex` in the samplesheet); arrays without them get NA. Age acceleration is computed when at least six arrays have a reported age. GrimAge, GrimAge2 and FitAge need licensed coefficients and are not available. BernabeuE2023c (cAge) is not in the default list: dnaMethyAge 0.2.0 fails for it when an array's estimate is below 20 years. A clock that fails gives NA, with the error in the notes of `clocks/clocks_coverage.tsv`. Intrinsic age acceleration needs at least five residual degrees of freedom after adjusting for age and six cell fractions, so at least 14 arrays with ages. Neither dnaMethyAge nor DunedinPACE is on conda or CRAN: the clock task installs them from their source packages at fixed GitHub commits. For runs without network access, download the two files once and pass their local paths.

## Exploratory analysis

The arrays that passed QC are summarised with two principal component analyses:

- Control-probe PCs (Lehne et al. 2015): the non-negative control probes of both channels, read by ewastools, describe technical variation such as bisulfite conversion, staining and hybridisation. `--n_ctrl_pcs` (default `auto`: the smallest number explaining at least 90% of their variance, at most 20) are kept as covariates.
- PCA of the `--n_variable_cpgs` (default 50,000) most variable CpGs of the filtered data (M-values). `--n_pcs` (default `auto`: the random matrix theory estimate of `isva::EstDimRMT`, at most 20) components are tested for association with plate, batch, Sentrix ID and row, detected fraction, control-probe PCs, predicted sex, age, cell lineages and every phenotype column. A strong association of a top component with plate or Sentrix row points to a batch effect that the models should include.

## Association models

Each model of `--models` is tested at every CpG of the filtered data with limma: a linear model on M-values with the trait and the covariates, moderated with `eBayes(trend = TRUE, robust = TRUE)`, and `duplicateCorrelation` when `block` is set. The same model on beta values gives the effect as a difference in beta (`delta_beta`). Test statistics are corrected for bias and inflation with bacon, and the genomic lambda is reported before and after.

- Covariate tokens: `@cells` adds five cell lineages (CD4 T, CD8 T, B, NK cells and monocytes; granulocytes are the reference), `@ctrlpcs` the kept control-probe PCs, `@technical` plate and Sentrix row, `@sex` the predicted sex. `@smoking` (DNA methylation smoking score) is not available yet, and a model that uses it stops the run at start-up.
- Arrays with a missing trait or covariate are left out of that model and listed with their reason. A model with fewer complete cases than its `min_n` is reported as not tested.
- Before testing, the trait is checked for association with plate and Sentrix row, and every covariate for collinearity (VIF). Confounding is flagged in `ewas_<model>_design_check.tsv`.
- Thresholds: FDR < 0.05, `--p_epigenome_wide` (default 9e-8, the EPIC threshold of Mansell et al. 2019) and `--p_suggestive` (default 1e-5).

## Differentially methylated regions

For every tested model, regions are called from the same statistics as the CpG-level analysis, by two methods:

- DMRcate (Peters et al. 2015): kernel smoothing of the moderated t-statistics along the genome (`--dmr_lambda` 1000 bp, `--dmr_c` 2). DMRcate's own EPIC v2 annotation step is not used, because it downloads the probe manifest and gene annotation during the run; the pipeline's probe annotation provides the positions and genes. Cross-hybridising probes, which DMRcate can remap to their off-target locus, are removed by the probe filter by default.
- dmrff (Suderman et al. 2018): regions of neighbouring CpGs (at most 500 bp apart) whose combined effect is significant after Bonferroni correction, accounting for the correlation between CpGs.

Region-level statistics are optimistic; the consensus of both methods is the more reliable result. `--skip_dmr` skips regions for every model, `skip_dmr: true` in the models file for one model. dmrff is not on conda or CRAN: the region task installs it from its source package at a fixed GitHub commit (`--dmrff_source`).

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run bhagesh-h/methylarray --input ./samplesheet.csv --outdir ./results  -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run bhagesh-h/methylarray -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull bhagesh-h/methylarray
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [bhagesh-h/methylarray releases page](https://github.com/bhagesh-h/methylarray/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
