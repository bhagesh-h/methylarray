# bhagesh-h/methylarray: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

- [Preprocessing](#preprocessing): one HDF5-backed SummarizedExperiment with the preprocessed arrays
- [Preprocessing QC](#preprocessing-qc): per-array QC statistics of the preprocessing method, with the sample metadata
- [Sample QC](#sample-qc): control metrics, sex check, SNP genotypes and one status per array
- [Reference resources](#reference-resources): the sesameData cache, when `--save_reference` is set
- [MultiQC](#multiqc): report with the sample QC table and software versions
- [Pipeline information](#pipeline-information): execution reports and parameters

### Preprocessing

- `preprocessing/<method>/<method>_se/`, where `<method>` is the value of `--preprocessing`
  - `assays.h5`: the assays, stored on disk.
  - `se.rds`: the SummarizedExperiment that refers to `assays.h5` by a relative path.

With `--preprocessing combined`, `preprocessing/sesame/` also holds the SeSAMe result whose masks were applied.

Load the directory in R, wherever it has been copied to:

```r
se <- HDF5Array::loadHDF5SummarizedExperiment("preprocessing/minfi/minfi_se")
```

Rows are probes (EPIC v2 probe IDs such as `cg00000029_TC21`; `rowData(se)$cpg_id` holds the CpG name without the suffix) and columns are arrays, sorted by sample name. Arrays marked `exclude` in the samplesheet are not included. Assays:

| Assay            | Content                                                                                                                                        |
| ---------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `beta`           | Beta values after preprocessing; values failing detection or the bead count (and the SeSAMe quality mask for `sesame` and `combined`) are `NA` |
| `beta_unmasked`  | The same beta values without masking, for methods that need complete data                                                                      |
| `detp`           | Detection p-values: minfi `detectionP` for `minfi`, ENmix `QCinfo` for `enmix`, SeSAMe pOOBAH for `sesame` and `combined`                      |
| `meth`, `unmeth` | Methylated and unmethylated intensities after preprocessing                                                                                    |

`colData(se)` holds the samplesheet metadata and the QC statistics described below. `metadata(se)` records the platform, the preprocessing method, its settings and the software versions.

### Preprocessing QC

- `qc/<method>/<method>_qc.tsv`: one row per array with the samplesheet metadata (`sample`, Sentrix ID, position, row and column, plate, well, batch, reported sex and age, replicate group) followed by:
  - `platform`, `prep`: detected platform and the `prepSesame` code used
  - `predicted_sex`: sex predicted by SeSAMe `inferSex`, or with the minfi `getSex` rule
  - `num_masked`, `frac_masked`: probes set to `NA` in `beta`
  - for `sesame`: the SeSAMe QC statistics from `sesameQC_calcStats` on the raw signal: detection (`frac_dt`, by probe type), intensity (`mean_intensity`, in-band and out-of-band means), channel and dye bias (`RGratio`, `RGdistort`, Infinium I channel switches) and the beta distribution (`mean_beta`, `frac_meth`, `frac_unmeth`, `frac_na`)
  - for `minfi`, `enmix` and `combined`: `normalisation`, `x_med` and `y_med` (median chrX and chrY copy number used for the sex call), `median_log2_meth` and `median_log2_unmeth` (minfi `getQC`), `frac_detected`, `num_low_beads`, `bisulfite_intensity` and `enmix_low_quality` (ENmix `QCinfo` flag, reported only)

### Sample QC

- `qc/sample_qc.tsv`: one row per samplesheet row, including excluded arrays:
  - `status`: `PASS`, `FLAG`, `FAIL` or `EXCLUDED`
  - `reasons`: why an array was flagged or failed, separated by `;`
  - `frac_detected`: fraction of probes with a detection p-value at or below `--detection_pval`
  - `reported_sex`, `predicted_sex_<method>` (from the preprocessing), `predicted_sex_ewastools` and `sex_check` (`match`, `mismatch`, `inconsistent` when the predictions disagree, `not_reported`, `not_predicted`)
  - `sex_x`, `sex_y`: normalised X and Y chromosome intensities
  - `snp_outlier`: ewastools SNP outlier score
  - `inferred_donor`: donor number from the SNP genotypes; technical replicates share a number
  - `replicate_group`, `failed_controls`
- `qc/sample_qc_retained.txt`: the arrays with status `PASS` or `FLAG`, one per line. Later steps use these arrays.
- `qc/ewastools/` (not with `--skip_ewastools`)
  - `ewastools_qc.tsv`: per array, the 21 Illumina BeadArray Controls Reporter metrics (`ctrl_*`), the metrics below Illumina's thresholds (`failed_controls`), X and Y intensities, the ewastools sex prediction, the SNP outlier score, the number of SNP probes, the inferred donor and whether the array is part of an identity conflict.
  - `ewastools_snp_conflicts.tsv`: pairs of arrays whose genotypes contradict the samplesheet: arrays from the same `replicate_group` with different genotypes, or arrays from different donors with the same genotypes. Each pair has its genotype agreement (0 to 1); pairs in the same `group` form one connected conflict.

How the status is set:

| Status     | Reason                                                                                                                                                             |
| ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `FAIL`     | more than `--max_undetected_fraction` of probes undetected; every sex prediction contradicts the reported sex; SNP outlier score above `--snp_outlier_threshold`   |
| `FLAG`     | control metrics below Illumina's thresholds; sex predictions disagree with each other; X and Y intensities between the male and female clusters; genotype conflict |
| `PASS`     | none of the above                                                                                                                                                  |
| `EXCLUDED` | `exclude` set in the samplesheet; the array was not processed                                                                                                      |

### Reference resources

- `references/` (only with `--save_reference`)
  - `sesame_cache/`: ExperimentHub cache with the sesameData resources for the platform. Pass it to `--sesame_cache` to run without network access.
  - `sesame_resources.tsv`: title, ExperimentHub ID, cached file and MD5 checksum of each resource.

### MultiQC

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The general statistics table shows, for every array, the QC status, the sex check, the SNP outlier score, the percentage of detected and masked probes, the mean intensity, the red/green ratio and the predicted sex. The Sample QC status table lists the reasons for every `FLAG` and `FAIL`, and the Sex check plot shows the X and Y intensities coloured by reported sex. The report also lists the software versions of every step and a methods description with references. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Parameters used by the pipeline run: `params.json`.

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
