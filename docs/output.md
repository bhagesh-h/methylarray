# bhagesh-h/methylarray: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

- [Preprocessing](#preprocessing): one HDF5-backed SummarizedExperiment with the preprocessed arrays
- [Preprocessing QC](#preprocessing-qc): per-array QC statistics of the preprocessing method, with the sample metadata
- [Sample QC](#sample-qc): control metrics, sex check, SNP genotypes and one status per array
- [Probe filtering](#probe-filtering): filtered data for analysis, one row per CpG and one row per probe, with the probe attrition
- [Cell composition](#cell-composition): blood cell-type fractions per array
- [Epigenetic clocks](#epigenetic-clocks): clock estimates, age acceleration and clock coverage per array
- [Exploratory analysis](#exploratory-analysis): control-probe PCs, PCA of the methylation data and its associations with technical and biological variables
- [Reference resources](#reference-resources): the sesameData cache and the probe annotation, when `--save_reference` is set
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

### Probe filtering

- `probe_filter/`
  - `<method>_filtered_se/`: the arrays with QC status `PASS` or `FLAG` and the probes that passed every filter, as an HDF5-backed SummarizedExperiment. With the default `--replicate_probe_strategy mean`, rows are CpGs (`cg########`, the mean of the replicate probes; `rowData` column `n_probes`). Assays: `beta` (failed values `NA`), `mvalue` (log2 of beta / (1 - beta), with beta clipped to [1e-4, 1 - 1e-4]) and `beta_unmasked`. `rowData` has the position, design, gene and CpG island context and the EPIC v1 and 450K IDs.
  - `<method>_filtered_uncollapsed_se/`: one row per probe, with the cross-hybridising probes kept so that DMR callers can remap them to their off-target locus. Assays `beta` and `beta_unmasked`.
  - `<method>_filtered_attrition.tsv`: probes removed at each step (`removed`), probes left (`remaining`), and probes carrying each reason in any order (`n_with_reason`). The `report_only_*` rows count flagged probes that were kept.
  - `<method>_filtered_exclusions.tsv.gz`: every removed probe with the step that removed it (`first_reason`), all reasons that apply, and its undetected and masked fractions.

The filters, in order: rs, nv and control probes; non-CG probes; probes without an hg38 position; chrX and chrY (kept with `--keep_sex_chr`); detection failure in more than `--probe_fail_fraction` of arrays; other masked values (bead count, quality mask) in more than `--bead_fail_fraction` of arrays; the SNP mask of `--snp_mask`; the Zhou general mask; cross-hybridisation (`--crosshyb_filter`). Manifest position mismatches (Peters et al. 2024 `MismatchPos`) and poor titration (Zhou) are counted but never removed.

### Cell composition

Only with `--tissue blood` (default).

- `cell_composition/cell_fractions.tsv`: one row per array that passed QC:
  - `lineage_B`, `lineage_CD4T`, `lineage_CD8T`, `lineage_NK`, `lineage_Mono`, `lineage_Gran`: the six lineages summed from the 12-type estimate (Gran: neutrophils, eosinophils, basophils)
  - `cent12_*`: EpiDISH RPC fractions for the 12 types of `cent12CT` (naive and memory B, CD4 T and CD8 T cells, regulatory T cells, NK cells, monocytes, neutrophils, eosinophils, basophils)
  - `dhs7_*`: fractions for the 7 types of `centDHSbloodDMC`, for comparison with earlier studies
  - `sum_cent12`, `sum_dhs7`: sum of the fractions of each reference
- `cell_composition/cell_fractions_reference_coverage.tsv`: for each reference, the number of reference CpGs, the missing ones by reason (`sex_chromosome`: left out on purpose, because X inactivation can mimic a cell-type signal; `not_on_array`: not measured by EPIC v2; `removed_by_filter`: removed by the probe filter), the CpGs found and how many of them only through their EPIC v1 or 450K ID, the coverage of the autosomal reference CpGs, the missing values set to the CpG median, and the status (`ok`, or `warn` below 98%; the run stops below `--min_reference_coverage`).

### Epigenetic clocks

Not with `--skip_clocks`.

- `clocks/clocks.tsv`: one row per array that passed QC, with the reported age (`age`) and one column per clock of `--clocks`. With at least six reported ages, `eaa_<clock>` holds the age acceleration (residual of the clock on age) of the clocks that estimate age, and, with at least 14 such arrays, `ieaa_<clock>` the intrinsic age acceleration (also adjusted for six of the seven `centDHSbloodDMC` cell fractions).
- `clocks/clocks_coverage.tsv`: for each clock, its CpGs, how many are in the data (`present`) and how many of those only under their EPIC v1 or 450K ID (`present_by_legacy_id`), the share present, the values that failed detection and were set to the CpG mean (`filled_values`), and notes.

### Exploratory analysis

- `exploration/pca_ctrl_pcs.tsv`: the control-probe principal components kept (`--n_ctrl_pcs`) per array, used as covariates by models with `@ctrlpcs`. Only the `sample` column without ewastools.
- `exploration/pca_ctrl_pca_variance.tsv`: share of the control-probe variance per component, and the cumulative share.
- `exploration/pca_scores.tsv`: the tested principal components (`--n_pcs`) of the most variable CpGs per array.
- `exploration/pca_variance.tsv`: share of variance per component, the cumulative share, and whether the component was tested.
- `exploration/pca_associations.tsv`: one row per variable and component: the variable group (`technical`, `biological`, `phenotype`), the test (linear regression for numeric variables, Kruskal-Wallis for categorical ones), the number of arrays, R-squared and the p-value. Technical variables: plate, batch, Sentrix ID and row, detected fraction, control-probe PCs. Biological variables: predicted sex, reported age and the cell lineages; every phenotype column is tested too.

The MultiQC report shows the associations as a heatmap of -log10 P, the variance per component, and PC1 against PC2 coloured by predicted sex.

### Reference resources

- `references/` (only with `--save_reference`)
  - `probe_annotation.tsv.gz`: one row per probe with the Peters et al. 2024 manifest columns and the Zhou mask flags used by the probe filter. Pass it to `--probe_annotation` to skip building it.
  - `probe_annotation_summary.tsv`, `probe_annotation_resources.tsv`: probes per flag, and the source files with their MD5 checksums.
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
