# bhagesh-h/methylarray: Output

## Introduction

This document describes the output produced by the pipeline. The MultiQC report and the HTML analysis report summarise the results at the end of the run.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

- [Preprocessing](#preprocessing): one HDF5-backed SummarizedExperiment with the preprocessed arrays
- [Preprocessing QC](#preprocessing-qc): per-array QC statistics of the preprocessing method, with the sample metadata
- [Sample QC](#sample-qc): control metrics, sex check, SNP genotypes and one status per array
- [Probe filtering](#probe-filtering): filtered data for analysis, one row per CpG and one row per probe, with the probe attrition
- [Cell composition](#cell-composition): blood cell-type fractions per array
- [Methylation QTLs](#methylation-qtls): CpGs with a cis-meQTL from a supplied meQTL table, with `--meqtl`
- [Epigenetic clocks](#epigenetic-clocks): clock estimates, age acceleration clock coverage and genotype-robust clocks per array
- [Exploratory analysis](#exploratory-analysis): control-probe PCs, PCA of the methylation data and its associations with technical and biological variables
- [Association models](#association-models): epigenome-wide association results per model of `--models`
- [Differentially methylated regions](#differentially-methylated-regions): DMRcate, dmrff and consensus regions per model
- [Gene set enrichment](#gene-set-enrichment): GO, KEGG and MSigDB enrichment, and KnowYourCG genomic context enrichment, of the significant CpGs per model
- [DNA methylation scores](#dna-methylation-scores): smoking, alcohol, BMI and IL-6 scores per array
- [Clocks, scores and cell fractions](#clocks-scores-and-cell-fractions): their agreement, their relation to the cell fractions, the batches and the phenotypes, and quartile groups of acceleration
- [GEO submission files](#geo-submission-files): processed matrix, signal intensities and a sample table, when `--export_geo` is set
- [Reference resources](#reference-resources): the sesameData cache and the probe annotation, when `--save_reference` is set
- [Analysis report](#analysis-report): one HTML report of the run
- [MultiQC](#multiqc): report with the QC and result summaries and software versions
- [Pipeline information](#pipeline-information): execution reports and parameters

### Preprocessing

- `preprocessing/<method>/<method>_se/`, where `<method>` is the value of `--preprocessing` (`minfi` by default)
  - `assays.h5`: the assays, stored on disk.
  - `se.rds`: the SummarizedExperiment that refers to `assays.h5` by a relative path.

With `--preprocessing combined`, `preprocessing/sesame/` also holds the SeSAMe result whose masks were applied, and `qc/sesame/` its QC table.

Load the directory in R, wherever it has been copied to:

```r
se <- HDF5Array::loadHDF5SummarizedExperiment("preprocessing/minfi/minfi_se")
```

Rows are probes (EPIC v2 probe IDs such as `cg00000029_TC21`; `rowData(se)$cpg_id` holds the CpG name without the suffix) and columns are arrays, sorted by sample name. Arrays marked `exclude` in the samplesheet are not included. Assays:

| Assay            | Content                                                                                                                                                                                                                       |
| ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `beta`           | Beta values after preprocessing; values failing detection or the bead count (and the SeSAMe quality mask for `sesame` and `combined`) are `NA`                                                                                |
| `beta_unmasked`  | The same beta values without masking, for methods that need complete data                                                                                                                                                     |
| `detp`           | Detection p-values: minfi `detectionP` for `minfi`, ENmix `QCinfo` for `enmix`, SeSAMe pOOBAH for `sesame` and `combined`                                                                                                     |
| `meth`, `unmeth` | Methylated and unmethylated intensities: after `prepSesame` for `sesame`, after Noob for every `--minfi_norm`, after ENmix background and RELIC correction (before quantile normalisation and RCP) for `enmix` and `combined` |

`colData(se)` holds the samplesheet metadata and the QC statistics described below. `metadata(se)` records the platform, the preprocessing method, its settings and the software versions.

### Preprocessing QC

- `qc/<method>/<method>_qc.tsv`: one row per array with the samplesheet metadata (`sample`, Sentrix ID, position, row and column, plate, well, batch, reported sex and age, replicate group) followed by:
  - `platform`: the array platform
  - for `sesame`: `prep`, the `prepSesame` code used; for `minfi`, `enmix` and `combined`: `method` and `normalisation`
  - `predicted_sex`: sex predicted by SeSAMe `inferSex`, or with the minfi `getSex` rule
  - `num_masked`, `frac_masked`: probes set to `NA` in `beta`
  - for `sesame`: the SeSAMe QC statistics from `sesameQC_calcStats` on the raw signal: detection (`frac_dt`, by probe type), intensity (`mean_intensity`, in-band and out-of-band means), channel and dye bias (`RGratio`, `RGdistort`, Infinium I channel switches) and the beta distribution (`mean_beta`, `frac_meth`, `frac_unmeth`, `frac_na`)
  - for `minfi`, `enmix` and `combined`: `x_med` and `y_med` (median chrX and chrY copy number used for the sex call), `median_log2_meth` and `median_log2_unmeth` (minfi `getQC`), `frac_detected`, `num_low_beads`, `bisulfite_intensity` and `enmix_low_quality` (ENmix `QCinfo` flag, reported only)

### Sample QC

- `qc/sample_qc.tsv`: one row per samplesheet row, including excluded arrays:
  - `status`: `PASS`, `FLAG`, `FAIL` or `EXCLUDED`
  - `reasons`: why an array was flagged or failed, separated by `;`
  - `frac_detected`: fraction of probes with a detection p-value at or below `--detection_pval`
  - `reported_sex`, `predicted_sex_<method>` (from the preprocessing), `predicted_sex_ewastools` and `sex_check` (`match`, `mismatch`, `inconsistent` when the predictions disagree, `not_reported`, `not_predicted`)
  - `sex_x`, `sex_y`: normalised X and Y chromosome intensities
  - `snp_outlier`: ewastools SNP outlier score
  - `relatives`: arrays called identical or first-degree relatives from the SNP genotypes (see `ewastools_relatedness.tsv.gz`); reported only, the status is set by the identity conflicts
  - `inferred_donor`: donor number from the SNP genotypes; technical replicates share a number
  - `replicate_group`, `failed_controls`
  - `ctrl_*` (with ewastools): the 21 Illumina BeadArray Controls Reporter metrics, also tested against the principal components
- `qc/sample_qc_retained.txt`: the arrays with status `PASS` or `FLAG`, one per line. Later steps use these arrays.
- `qc/ewastools/` (not with `--skip_ewastools`)
  - `ewastools_qc.tsv`: per array, the 21 Illumina BeadArray Controls Reporter metrics (`ctrl_*`), the metrics below Illumina's thresholds (`failed_controls`), X and Y intensities, the ewastools sex prediction, the SNP outlier score, the number of SNP probes, the inferred donor, whether the array is part of an identity conflict, and its related arrays (`relatives`, for example `B (identical); C (first-degree)`).
  - `ewastools_controls.tsv.gz`: raw intensities of every control probe in both channels: one row per probe and channel (`channel`, `group`, `name`), one column per array. The non-negative controls are the input of the control-probe PCs (see [Exploratory analysis](#exploratory-analysis)).
  - `ewastools_snp_summary.tsv`: per array, from the SNP probes: the sex used (reported, or predicted where none is reported), the SNP outlier score, the share of heterozygous calls at the autosomal SNPs and at the SNPs that behave as X-chromosome SNPs (men almost never heterozygous), the highest share of identical autosomal genotypes with any other array (near 1: the same person), and the number of autosomal and X-like SNPs. No genotypes are written.
  - `ewastools_snp_histograms.tsv`: counts of the autosomal SNP probe betas and of the pairwise genotype agreement, in bins of 0.02, for the report.
  - `ewastools_relatedness.tsv.gz`: every pair of arrays below the SNP outlier threshold, from the autosomal SNP genotypes: the SNPs called in both (`snps`), the mean and variance of identity by state (alleles shared per SNP, 0 to 2; the omicsPrint plot), the log10 likelihood ratios against unrelated for identical (`lod_identical`) and first-degree pairs (`lod_first_degree`, the better of parent-offspring and full siblings), the `relation` (`identical`, `first-degree`, `unrelated`, or `too few SNPs` below 20) and whether the samplesheet expects the same donor (`expected_same_donor`, same `replicate_group`). A pair is related when its likelihood ratio exceeds 60 times the number of pairs, which keeps the chance of any false call at or below 5%.
  - `ewastools_relatedness_summary.tsv`: the arrays and SNPs used, the number of pairs, the log10 threshold, the expected log10 likelihood ratio of a true first-degree pair (below the threshold, most relatives are missed: with the 60 or so SNP probes of the array this is the usual case, while identical pairs are always detected), and the counts of identical pairs, identical pairs from different expected donors, first-degree pairs and replicate pairs that are not identical.
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
  - `<method>_filtered_attrition.tsv`: probes removed at each step (`removed`), probes left (`remaining`), and probes carrying each reason in any order (`n_with_reason`). The `report_only_*` rows count flagged probes that were kept.
  - `<method>_filtered_exclusions.tsv.gz`: every removed probe with the step that removed it (`first_reason`), all reasons that apply, and its undetected and masked fractions.

The filters, in order: rs, nv and control probes; non-CG probes; probes without an hg38 position; chrX and chrY (kept with `--keep_sex_chr`); detection failure in more than `--probe_fail_fraction` of arrays; other masked values (bead count, quality mask) in more than `--bead_fail_fraction` of arrays; the SNP mask of `--snp_mask`; the Zhou general mask; cross-hybridisation (`--crosshyb_filter`). Manifest position mismatches (Peters et al. 2024 `MismatchPos`), poor titration (Zhou) and, when their lists are given, Illumina's flagged and mapping-inaccuracy probes are report-only: counted in the `report_only_*` rows of the attrition table, never removed.

### Cell composition

Only with `--tissue blood` (default).

- `cell_composition/cell_fractions.tsv`: one row per array that passed QC:
  - `lineage_B`, `lineage_CD4T`, `lineage_CD8T`, `lineage_NK`, `lineage_Mono`, `lineage_Gran`: the six lineages summed from the 12-type estimate (Gran: neutrophils, eosinophils, basophils)
  - `cent12_*`: EpiDISH RPC fractions for the 12 types of `cent12CT` (naive and memory B, CD4 T and CD8 T cells, regulatory T cells, NK cells, monocytes, neutrophils, eosinophils, basophils)
  - `dhs7_*`: fractions for the 7 types of `centDHSbloodDMC`, for comparison with earlier studies
  - `sum_cent12`, `sum_dhs7`: sum of the fractions of each reference
  - `idol_*` (not with `--skip_idol`): fractions of the six IDOL cell types (CD8T, CD4T, NK, Bcell, Mono, Neu), Houseman constrained projection
- `cell_composition/cell_fractions_reference_coverage.tsv`: for each reference, the number of reference CpGs, the missing ones by reason (`sex_chromosome`: left out on purpose, because X inactivation can mimic a cell-type signal; `not_on_array`: not measured by EPIC v2; `removed_by_filter`: removed by the probe filter), the CpGs found and how many of them only through their EPIC v1 or 450K ID, the coverage of the autosomal reference CpGs, the missing values set to the CpG median, and the status (`ok`, or `warn` below 98%; the run stops below `--min_reference_coverage`, except for IDOL, which is then skipped with the status `fail: skipped`).
- `cell_composition/cell_fractions_idol_agreement.tsv`: each cent12CT lineage against the matching IDOL cell type: Spearman rho and mean difference (lineage minus IDOL). Granulocytes include eosinophils and basophils, IDOL neutrophils do not.

### Methylation QTLs

- `meqtl/meqtl_cpgs.tsv` (with `--meqtl`): the EPIC v2 CpGs with a cis-meQTL in the supplied table: the name in the table, the number of SNPs, and the strongest SNP with its P and effect.
- `meqtl/meqtl_cpgs_summary.tsv`: rows read, cis rows, rows below `--meqtl_p`, CpGs in the table and CpGs matched on the array.

### Epigenetic clocks

Not with `--skip_clocks`.

- `clocks/clocks.tsv`: one row per array that passed QC, with the reported age (`age`) and one column per clock of `--clocks`. With at least six reported ages, `eaa_<clock>` holds the age acceleration (residual of the clock on age) of the clocks that estimate age, and, with at least 14 such arrays, `ieaa_<clock>` the intrinsic age acceleration (also adjusted for six of the seven `centDHSbloodDMC` cell fractions).
- `clocks/clocks_coverage.tsv`: for each clock, its CpGs, how many are in the data (`present`), how many have a common SNP (minor allele frequency 1% or more) at the CpG or extension base (`snp_cpgs`), how many are present only under their EPIC v1 or 450K ID (`present_by_legacy_id`), the share present, the values that failed detection and were set to the CpG mean (`filled_values`), and notes.
- `clocks/clocks_robust.tsv` (with `--meqtl`): the same columns as `clocks.tsv`, from the clocks computed with the CpGs that have a common SNP or a cis-meQTL held at their mean over the arrays (genotype-robust clocks). Compare their accelerations with those of `clocks.tsv`, not their ages. `clocks_coverage.tsv` then also gives the clock CpGs held per clock (`robust_removed`).
- `clocks/clocks_weights.tsv`: for each clock that is a weighted sum of CpGs (not the principal-component clocks or epiTOC2), every CpG with its gene, published weight, SD of its beta value in these arrays, standardised weight (weight times SD, in the clock's units per SD of methylation; sorted by its size) whether it has a common SNP, and, with `--meqtl`, whether it has a cis-meQTL (`meqtl`).

### Exploratory analysis

- `exploration/pca_ctrl_pcs.tsv`: the control-probe principal components kept (`--n_ctrl_pcs`) per array, used as covariates by models with `@ctrlpcs`. Only the `sample` column without ewastools.
- `exploration/ancestry_pcs.tsv` (not with `--skip_ancestry`): the genetic ancestry components kept (`--n_ancestry_pcs`) per array, used as covariates by models with `@ancestry`. Only the `sample` column when fewer than 10 arrays passed QC.
- `exploration/ancestry_pcs_variance.tsv`: share of variance per ancestry component (up to 10), the cumulative share, and whether the component is a covariate.
- `exploration/ancestry_pcs_summary.tsv`: the arrays, SNP CpGs and rs probes used, and the covariates the SNP CpGs were residualised for.
- `exploration/ancestry_pcs_associations.tsv`: each ancestry component against the DNA plate, BeadChip, Sentrix row, sex, age, the cell lineages and, when the phenotype table has one, an `ancestry`, `ethnicity` or `race` column (linear regression: R-squared and P). Technical variables should explain little.
- `exploration/pca_ctrl_pca_variance.tsv`: share of the control-probe variance per component, the cumulative share, and whether predicted sex was removed from the control probes before the PCA (`sex_removed`; when both sexes have at least two arrays).
- `exploration/pca_scores.tsv`: the tested principal components (`--n_pcs`) of the most variable CpGs per array.
- `exploration/pca_variance.tsv`: share of variance per component, the cumulative share, and whether the component was tested.
- `exploration/pca_outliers.tsv`: per array the robust Mahalanobis distance on PC1 and PC2 (each centred on its median and scaled by its MAD) and `pca_outlier` (`TRUE` above the chi-square 0.999 quantile with 2 degrees of freedom). Report only: flagged arrays stay in the analysis. Needs at least 5 arrays.
- `exploration/variance_partition.tsv.gz` (not with `--skip_variance_partition`): for each of the `--variance_partition_cpgs` most variable CpGs, the percentage of its variance explained by each variable and the residual.
- `exploration/variance_partition_summary.tsv`: per variable, its effect type (random or fixed), the median, mean and 90th percentile of the percentage over the CpGs, the CpGs above 10%, the arrays and CpGs used, and the model formula; a `note` says why it did not run (fewer than 10 arrays with complete variables).
- `exploration/pca_associations.tsv`: one row per variable and component: the variable group (`technical`, `control metrics`, `biological`, `phenotype`), the test (linear regression for numeric variables, Kruskal-Wallis for categorical ones), the number of arrays, R-squared, adjusted R-squared (corrected for the number of levels of a categorical variable) and the p-value. Technical variables: plate, batch, Sentrix ID and row, detected fraction, control-probe PCs. Control metrics: the Illumina BeadArray Controls Reporter metrics of ewastools (`control: <metric>`), as ChAMP's `champ.SVD` does with the internal controls. Biological variables: predicted sex, reported age and the cell lineages; every phenotype column is tested too.

The MultiQC report shows the associations as a heatmap of -log10 P, the variance per component, and PC1 against PC2 coloured by predicted sex.

### Association models

Only with `--models`.

- `ewas/ewas_summary_meqtl.tsv` (with `--meqtl`): per model and contrast, the tested CpGs and the CpGs at FDR < 0.05 with a cis-meQTL, their shares, the odds ratio and Fisher's exact P. Hits with a cis-meQTL are candidates for genetic confounding.
- `ewas/ewas_summary_comparison.tsv`: every pair of tested model contrasts: hits of each at FDR < 0.05, shared hits, the share of the first contrast's hits found by the second, the share of shared hits with the same direction, Fisher's odds ratio and P of the overlap against chance, and the Pearson correlation of the t statistics over all CpGs and over the hits of either.
- `ewas/ewas_summary_cells.tsv` (not with `--skip_cell_ewas`): one row per model and blood cell lineage: CpGs at FDR < 0.05 within the lineage, with higher and lower methylation, genomic lambda, the median fraction of the lineage, its `status` (`ok`, or `unreliable` with the reason: median fraction below 5% or lambda above 1.3) and the top CpG.
- `ewas/ewas_summary.tsv`: one row per model and contrast: status (`tested`, or why not), complete cases, CpGs tested, genomic lambda before and after bacon, bacon inflation and bias, CpGs at FDR < 0.05, at `--p_epigenome_wide` and at `--p_suggestive` (raw, and after bacon), the hits at FDR < 0.05 with higher and with lower methylation (`hits_higher`, `hits_lower`), the top CpG, and `design_flags`, the design checks of the model that were flagged (see `ewas_<model>_design_check.tsv`). With more than one tested model or contrast: `meff`, the effective number of independent tests across them (Galwey 2009, from the correlation of their t statistics over the CpGs), `p_study_wide` (`--p_epigenome_wide` / Meff) and `hits_study_wide`, the CpGs below it.
- `ewas/<model>/`
  - `ewas_<model>_stats.tsv.gz`: every CpG: `cpg_id`, position and gene context, `contrast`, `logfc_m` and `se_m` (coefficient and standard error on the M-value scale, per unit of the transformed trait or for the level against the reference), `delta_beta` (the same coefficient on beta values), moderated `t`, `p_value`, `fdr` (Benjamini-Hochberg), `p_bacon` and `fdr_bacon` (after bacon correction of bias and inflation), `mean_beta` and `n`.
  - `ewas_<model>_top.tsv`: the 1,000 CpGs with the smallest P per contrast.
  - `ewas_<model>_cells.tsv.gz`: the cell-type-specific association of every CpG (not with `--skip_cell_ewas`; header only when the model has a categorical trait with more than two levels or no cell fractions). The CellDMC design of EpiDISH (Zheng et al. 2018) models beta values on the six blood cell lineages, the trait times each lineage and the other covariates of the model, without an intercept; one row per CpG and lineage with `estimate` (change in beta within the lineage per unit of the trait, or between the groups), `se`, `t`, `p_value` (ordinary t-test, as CellDMC), `fdr` (Benjamini-Hochberg within the lineage) and `lineage_status`. Lineages with a median fraction below 5% or a genomic lambda above 1.3 are `unreliable`: with small fractions the interaction estimate is unstable and can give thousands of hits even for traits without a hit in the main model. Treat their hits as artefacts.
  - `ewas_<model>_design_check.tsv`: association of the trait with plate and Sentrix row (`confounded` below P = 0.01), the variance inflation factor of each design column (`collinear` above 10), and every numeric covariate against the trait (Spearman rho for a continuous trait, R-squared otherwise); methylation-derived covariates (cell fractions, control-probe and ancestry PCs, DNA methylation scores, PCs) are flagged `follows the trait` below P = 0.001, because adjusting for an axis that the trait drives can absorb its association.
  - `ewas_<model>_characteristics.tsv`: the arrays of the model described by trait group (binary or categorical trait) or overall (continuous trait): the number of arrays, mean (SD) of numeric covariates and n (%) of categorical ones, with a Kruskal-Wallis or Fisher's exact P between the groups. Aggregates only: counts from 1 to 4 are shown as `<5`.
  - `ewas_<model>_missing.tsv`: the arrays per combination of missing model variables (`missing`, separated by semicolons; empty for arrays with every value), after the subset and before complete cases are taken.
  - `ewas_<model>_resolved.json`: the model with its covariate tokens expanded, the number of complete cases, and every array left out with its reason.
  - `ewas_<model>_summary.json`: the per-model values of `ewas_summary.tsv`.

### Differentially methylated regions

Only with `--models`, and not with `--skip_dmr`. Every model gets a folder; for a model that was not tested or has `skip_dmr: true` the tables are empty and the summary says why.

- `dmr/<model>/`
  - `dmr_<model>_dmrcate.tsv`: DMRcate regions per contrast: position, number of CpGs, minimum smoothed FDR, Stouffer, harmonic mean FDR and Fisher P, maximum and mean delta-beta, and the genes of the CpGs in the region.
  - `dmr_<model>_dmrff.tsv`: dmrff regions per contrast at Bonferroni-adjusted P < 0.05 with at least `--dmr_min_cpgs` CpGs: position, number of CpGs, combined estimate (M-value scale), SE, z, P, adjusted P and genes.
  - `dmr_<model>_consensus.tsv` and `.bed`: regions found by both methods that overlap and change in the same direction. Regions found by one method only are in the method tables and have lower confidence.
  - `dmr_<model>_summary.json`: the number of CpGs, CpGs at FDR < 0.05 and regions per method and contrast, or why the analysis did not run.

### Gene set enrichment

Only with `--models`, and not with `--skip_enrichment`. Every model gets a folder; for a model that was not tested or has `skip_enrichment: true` the table is empty and the summary says why.

- `enrichment/<model>/enrichment_<model>.tsv`: one row per collection, contrast and gene set: `collection` (`GO`, `KEGG`, or the GMT file name such as `h.all.v2025.1.Hs`), `set` (GO or KEGG identifier, or gene set name), `term`, `ontology` (GO only: BP, CC, MF), `n` (genes in the set), `de` (genes with a significant CpG, weighted for CpGs in several genes), `p_value` and `fdr`, sorted by P.
- `enrichment/<model>/kycg_<model>.tsv` (not with `--skip_kycg`): KnowYourCG genomic context enrichment, one row per contrast, query (`all`, `higher`, `lower` methylation), knowledgebase (`database`) and feature (`set`): `log2_or` (log2 odds ratio; `Inf` or `-Inf` when every or no query probe is in the set), `p_value`, `fdr` (within the knowledgebase), `overlap` (query probes in the set), `query_probes`, `set_probes`.
- `enrichment/<model>/kycg_<model>_summary.json`: significant and tested CpGs per contrast and whether each query ran.
- `enrichment/<model>/enrichment_<model>_summary.json`: significant and tested CpGs per contrast, whether enrichment ran, and why a collection did not run (for example KEGG without network access).

### DNA methylation scores

Not with `--skip_dnam_scores`.

- `dnam_scores/dnam_scores.tsv`: one row per array that passed QC, one column `dnam_<score>` per score of `--dnam_scores`: the weighted sum of the beta values of the score's CpGs present on the array. Higher values mean more smoking, alcohol, BMI or IL-6 as predicted from methylation; the scale is that of the published weights, so compare scores within a study.
- `dnam_scores/dnam_scores_weights.tsv`: for each score, every CpG present with its weight, the SD of its beta value and the standardised weight.
- `dnam_scores/dnam_scores_coverage.tsv`: per score, its CpGs, how many are in the data (and how many only under their EPIC v1 or 450K ID), the share present, the values that failed detection and were set to the CpG mean, and the reference.

### Clocks, scores and cell fractions

Not with both `--skip_clocks` and `--skip_dnam_scores`. Only aggregates are written, never values of single arrays. Acceleration is a clock's residual on reported age; with cells held, it is also adjusted for the cell lineages other than granulocytes. Phenotype columns with more than two numeric values are continuous, columns with two values binary (reference `no`, `0`, `false` or the first in sort order); other columns and columns with fewer than 20 values are left out. P-values are adjusted with Benjamini-Hochberg (`q`).

- `biomarkers/biomarkers_agreement.tsv`: Spearman correlation of every pair of accelerations and scores.
- `biomarkers/biomarkers_cells_measures.tsv`: Spearman correlation of each EpiDISH cell lineage with each acceleration and score.
- `biomarkers/biomarkers_cells_phenotypes.tsv`: each cell lineage against age, sex and every phenotype, with age and sex held: the estimate in SD of the fraction per SD of a continuous trait, or with a binary one.
- `biomarkers/biomarkers_batch.tsv`: each acceleration and score against the DNA plate, the BeadChip and the Sentrix row (adjusted R-squared and F-test).
- `biomarkers/biomarkers_phenotypes.tsv`: each acceleration and score against every phenotype, plain, with cells held, with the ancestry components held (`ancestry held`) and, when some analysed arrays have a borderline SNP outlier score, without them (`low contamination`): Spearman correlation for a continuous phenotype, difference in SD between the groups of a binary one (t-test).
- `biomarkers/biomarkers_ancestry.tsv` (with the ancestry components): Spearman correlation of each acceleration and score with each genetic ancestry component; the phenotype associations with the components held are the `ancestry held` version of `biomarkers_phenotypes.tsv`.
- `biomarkers/biomarkers_robust.tsv` (with `--meqtl`): Spearman correlation of each clock acceleration with the one of the genotype-robust clock; the phenotype associations with the robust accelerations are the `genotype-robust` version of `biomarkers_phenotypes.tsv`.
- `biomarkers/biomarkers_missing.tsv` (with `--phenotypes`): the arrays per combination of missing phenotypes (`missing`, separated by semicolons; empty for arrays with every phenotype).
- `biomarkers/biomarkers_contamination.tsv` (with ewastools): Spearman correlation of the ewastools SNP outlier score with each acceleration and score on the analysed arrays, with the robust cutoff (median plus three MADs of the score) and the number of arrays above it. Contamination with another person's DNA pulls beta values towards the middle and can bias clocks; these arrays pass QC (their score is below `--snp_outlier_threshold`), so the `low contamination` version of the phenotype table shows whether they change the results.
- `biomarkers/biomarkers_quartiles.tsv`: for each acceleration, the bottom quarter (decelerated), middle half (normal) and top quarter (accelerated): the median of every characteristic (percentage for binary ones), its mean z per group, and a trend test across the ordered groups (Kendall's tau, equivalent to the Jonckheere-Terpstra test).

### GEO submission files

Only with `--export_geo`. For the arrays that passed sample QC, following the GEO rules for methylation arrays.

- `geo/geo_processed_matrix.tsv.gz`: `ID_REF`, then per array the normalised beta value of the preprocessing route (without detection masking) and `<array> Detection Pval`.
- `geo/geo_signal_intensities.tsv.gz`: `ID_REF`, then per array `<array> Unmethylated signal` and `<array> Methylated signal` (unnormalised, read from the IDAT files with minfi `preprocessRaw`) and `<array> Detection Pval`.
- `geo/geo_samples.tsv`: array name, red and green IDAT file names, Sentrix fields and the processing description, as a start for the GEO metadata sheet. No phenotype columns.

Columns are named after the samplesheet sample names. Rename them before submission if they identify participants, and submit the IDAT files themselves as raw data.

### Reference resources

- `references/` (only with `--save_reference`)
  - `probe_annotation.tsv.gz`: one row per probe with the Peters et al. 2024 manifest columns and the Zhou mask flags used by the probe filter. Pass it to `--probe_annotation` to skip building it.
  - `probe_annotation_summary.tsv`, `probe_annotation_resources.tsv`: probes per flag, and the source files with their MD5 checksums.
  - `sesame_cache/` (`sesame` and `combined` routes): ExperimentHub cache with the sesameData resources for the platform. Pass it to `--sesame_cache` to run without network access.
  - `kycg_cache/` and `kycg_resources.tsv` (with models, not with `--skip_kycg`): the KnowYourCG knowledgebases in an ExperimentHub cache, and their ExperimentHub IDs and MD5 checksums. Pass the folder to `--kycg_cache` for runs without network access.
  - `sesame_resources.tsv` (`sesame` and `combined` routes): title, ExperimentHub ID, cached file and MD5 checksum of each resource.

### Analysis report

Not with `--skip_report`.

- `report/methylarray_report.html`: a self-contained HTML report (Quarto) of the run. It shows aggregates and plots, not individual phenotype values; parts of steps that did not run are left out.
  - Summary: the main results, each with a link to its section, and the main settings of the run.
  - Samples: QC status with the flagged and failed arrays, sex check, SNP genotypes, relatedness and genetic ancestry.
  - Probes: the probes removed per filter step, and an UpSet plot of the reasons they share.
  - Cell composition: the lineage fractions and their comparison with the IDOL reference.
  - Clocks and scores: agreement with the reported age, age acceleration, methylation scores, the CpGs measured, the largest weights, the CpGs that clocks share and those with a SNP, a removed probe or an association, and with `--meqtl` the genetic influence.
  - Clocks, cells and phenotypes: agreement of the clocks and scores, cell fractions, batches, phenotypes, quartile groups of acceleration, missing phenotypes, ancestry and contamination.
  - Sources of variation: principal components with what they follow, and the variance partition.
  - Association models: overview, design checks and comparison of the models, then one section per model with its arrays, associations (Manhattan plots, volcano and QQ plot, the CpGs with the smallest P per direction), genomic context of the hits, cell types, regions, gene sets and KnowYourCG features.
  - Methods: the analysis steps with references, software versions and limitations.
- The contents on the left open for the part being read. The search box above them marks every match and opens the table or tab that holds it. Tables are closed until their name is clicked, or all at once with the buttons below the search box. Figures that belong together are tabs, are vector graphics (or 300 dpi images when they hold many points) and are numbered.
- `report/figures/`: every figure of the report as a 300 dpi PNG and a vector PDF, named by its number and description (for example `figure_34_<model>_cpgs_with_higher_and_lower_methylation.pdf`), for use in papers and slides.

### MultiQC

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The general statistics table shows, for every array, the QC status, the sex check, the SNP outlier score (with ewastools), and from the preprocessing route the percentage of detected and masked probes and the predicted sex, with the median log2 methylated intensity (minfi-based routes) or the mean intensity and red/green ratio (`sesame`). The Sample QC status table lists the reasons for every `FLAG` and `FAIL`, and the Sex check plot shows the X and Y intensities coloured by reported sex. Further sections follow the pipeline order: the SNP genotype table, the related and replicate pairs (SNP relatedness), the genetic ancestry components, probe attrition, cell composition and its agreement with IDOL, the meQTL table (with `--meqtl`), clock coverage and the clock values over the arrays (violin plot), score coverage and values, the agreement of the accelerations and scores and their relation to the cell fractions, the phenotypes, the ancestry components, the SNP outlier score and the technical batches, the PCA associations, variance and PC1 against PC2 (with the PCA outliers named), the variance partition, the association summary (with the hits that have a cis-meQTL and the flagged design checks), the cell-type-specific hits per model and lineage, the regions per DMR method, and one table of the top gene sets and KnowYourCG sets per model. The report also lists the software versions of every step and a methods description with references. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

- `pipeline_info/`
  - Nextflow reports: `execution_report_<timestamp>.html`, `execution_timeline_<timestamp>.html`, `execution_trace_<timestamp>.txt` and `pipeline_dag_<timestamp>.html`.
  - Pipeline reports: `pipeline_report.html` and `pipeline_report.txt`, only when `--email` or `--email_on_fail` is set.
  - Software versions of every task: `methylarray_software_mqc_versions.yml`.
  - Parameters of the run: `params_<timestamp>.json`.

The Nextflow reports help to troubleshoot a run and record the launch command, run times and resource usage of every task.
