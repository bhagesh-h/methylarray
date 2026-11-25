# bhagesh-h/methylarray: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### `Added`

- `-profile qi`: the sample QC and probe rules of Qi et al. 2025, with the EPIC v2 counterparts of their EPIC v1 probe lists.
- `--exclude_probes`: a list of probes or CpGs to remove, counted as its own step of the probe attrition.
- `--illumina_filter`: remove the probes of the Illumina lists (`--illumina_flagged_probes`, `--illumina_mapping_inaccuracy`) instead of only counting them.
- `--ewas_outlier_iqr`: set extreme M-values of a CpG to missing before the association tests (off by default).
- HTML report: UpSet plots of the reasons probes are removed for, the CpGs that clocks and scores share, the clock CpGs with a SNP, a cis-meQTL, a removed probe or an association, the missing phenotypes and model variables, and the hits that model contrasts share.
- `biomarkers_missing.tsv` and `ewas_<model>_missing.tsv`: the arrays per combination of missing phenotypes and of missing model variables.
- `--probe_filters`: apply only the named probe filter steps (or none); the others are counted as report-only.
- `--replicate_probe_strategy best`: keep one probe per CpG, with the others as their own step of the probe attrition.
- `--ewas_cells idol`: the IDOL cell fractions as the cell covariates of the association models.

### `Changed`

- HTML report: nine numbered parts with a contents list that opens for the part being read, a summary with links, one section per association model, tabs for figures that belong together, buttons to open or close all tables, and clocks under the names they are known by. Explanations are shorter, and the phenotype heatmaps have one row per phenotype.

### `Fixed`

- `--replicate_probe_strategy none`: cell fractions, gene sets, KnowYourCG, meQTL counts and regions now accept one row per probe.
- Association summary: a run in which no model has enough arrays no longer stops, and the model comparison is in a fixed order.
- HTML report: a phenotype table with many columns no longer stops the report. The quartile-group figure is split into parts, the phenotype heatmaps widen with the number of phenotypes, and no figure grows above the 50 inch limit of ggsave.

## v1.0.2 - [2026-11-24]

Fixes the bacon-corrected P values of the association models; no other result changes.

### `Fixed`

- EWAS: the bacon-corrected P value of a CpG with a very large positive t statistic (P below about 1e-16) was far too small, because the z value was computed through a probability that rounds to 1. It is now computed on the log scale, the same for both signs. Other bacon-corrected values can shift slightly, because the bacon sampler reacts to changes in the last decimals of its input.

## v1.0.1 - [2026-11-22]

Report readability release: no change to any result.

### `Added`

- Zenodo archive of each release, with DOI [10.5281/zenodo.24399752](https://doi.org/10.5281/zenodo.24399752) for all versions (1.0.2: [10.5281/zenodo.24516701](https://doi.org/10.5281/zenodo.24516701); 1.0.1: [10.5281/zenodo.24458958](https://doi.org/10.5281/zenodo.24458958); 1.0.0: [10.5281/zenodo.24399753](https://doi.org/10.5281/zenodo.24399753)); README badge and `CITATION.cff`.
- Zenodo DOI in the reports: the pipeline manifest carries the concept DOI, so the methods text of the MultiQC and the HTML report cites it, and the HTML report names it with the pipeline version.
- HTML report: a search box above the contents that highlights every match, steps through them and opens the tables that hold one.
- HTML report: a template of a gene and its CpG island (island, shore, shelf, open sea; TSS1500, TSS200, 5' and 3' UTR, first and other exons, introns, intergenic) before the first genomic context figure, to read those figures by.

### `Changed`

- HTML report: the contents are on the left and the report uses the rest of the page width.
- HTML report: every table is a numbered toggle with its name in bold, closed until clicked; this includes the software versions.
- HTML report: figures are embedded as vector graphics (SVG), or as the 300 dpi PNG when they hold more than 5,000 points, so they stay sharp when zoomed; they were 110 dpi images before.
- HTML report: the quartile-group figure has the trend stars in their own column inside each panel and the colour legend below the panels; the PCA association heatmap has vertical PC labels and grows with the number of variables; the methylation level figure of the hits is horizontal so its group labels do not overlap.
- Methods text (MultiQC and HTML report): one point per analysis step instead of one paragraph, the command line with one option per line and wrapped, and the references in a short APA style (first author and "et al." for three or more authors, DOI as a link).

## v1.0.0 - [2026-11-16]

Initial release of bhagesh-h/methylarray, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Samplesheet validation for IDAT input: matching red and green files, unique IDATs, Sentrix ID and position from the file names, manual exclusions with a reason.
- Optional phenotype table (`--phenotypes`) and model definitions (`--models`), checked against the samplesheet and each other before any task starts.
- SeSAMe preprocessing (`--preprocessing sesame`): sesameData cache (`SESAME_CACHE`, reusable offline with `--sesame_cache`), per-array `prepSesame` (`SESAME_PREPROCESS`, default `QCDPB`) and one HDF5-backed SummarizedExperiment with the sample QC statistics (`SESAME_BUILDSE`), plus a MultiQC general statistics table.
- minfi, ENmix and combined preprocessing (`--preprocessing minfi|enmix|combined`, #10): one `MINFI_PREPROCESS` module with Noob, Funnorm, quantile or SWAN normalisation (`--minfi_norm`), ENmix oob background with RELIC and RCP (optional `--enmix_quantile`), and the SeSAMe masks for the combined route. All routes write the same SummarizedExperiment and never drop arrays or probes.
- Sample QC (#12): ewastools control metrics against Illumina's thresholds, X/Y intensities and sex prediction, SNP genotypes with contamination scores and identity checks against `replicate_group` (`EWASTOOLS_QC`), and one status per array (`PASS`, `FLAG`, `FAIL`, `EXCLUDED`) with its reasons (`SAMPLEQC_SUMMARY`), in `qc/sample_qc.tsv`, `qc/sample_qc_retained.txt` and the MultiQC report. New parameters `--max_undetected_fraction`, `--snp_outlier_threshold`, `--skip_ewastools` and `--ewastools_source`.
- Probe annotation and filtering (#14): one EPIC v2 probe annotation table from the Peters et al. 2024 extended manifest and the Zhou lab InfiniumAnnotation mask (`PROBES_ANNOTATE`), and a probe filter on the arrays that passed sample QC (`PROBES_FILTER`): probe type, non-CG, no position, sex chromosomes, detection and bead-count failure fractions, SNP mask by population, Zhou general mask and cross-hybridisation, each counted in an attrition table and listed per probe. Manifest position mismatches and poor titration are reported, never removed. Replicate probes are collapsed to CpGs by mean; an uncollapsed matrix keeps the cross-hybridising probes for DMR remapping. M-values from clipped betas.
- Blood cell-type composition (#16): EpiDISH RPC with the 12-type `cent12CT` reference (six aggregated lineages) and the 7-type `centDHSbloodDMC` reference on the filtered CpG-level data (`EPIDISH_RPC`), reference CpGs matched through EPIC v1 and 450K IDs, coverage checked (`--min_reference_coverage`), MultiQC stacked bar plot. `--tissue other` skips it.
- Epigenetic clocks (#18): dnaMethyAge clocks (Horvath, Hannum, PhenoAge, skin and blood, Zhang, DNAmTL, mortality score, PC clocks, DNAmIC, epiTOC2) and DunedinPACE on the detection-masked CpG-level betas (`CLOCKS_COMPUTE`), clock CpGs matched through EPIC v1 and 450K IDs, coverage per clock, age acceleration and intrinsic age acceleration. Both packages are installed from pinned source packages (`--dnamethyage_source`, `--dunedinpace_source`).
- Exploratory analysis (#20): control-probe PCs from the ewastools control intensities (Lehne et al. 2015, `--n_ctrl_pcs`) and a PCA of the most variable CpGs (`--n_variable_cpgs`, `--n_pcs`, random matrix theory by default) with the association of each component with technical, biological and phenotype variables (`EXPLORATION_PCA`); MultiQC heatmap, scree plot and PC1/PC2 scatter. `EWASTOOLS_QC` also writes the control-probe intensities.
- Epigenome-wide association (#22): one limma model per entry of `--models` on the filtered CpG-level M-values (`EWAS_LIMMA`): covariate tokens (`@cells`, `@ctrlpcs`, `@technical`, `@sex`), subset, trait transforms, complete cases with reasons, design check (trait against plate and Sentrix row, VIF), `eBayes(trend, robust)`, `duplicateCorrelation` for blocks, delta-beta, genomic lambda and bacon, FDR, `--p_epigenome_wide` and `--p_suggestive` thresholds; one summary table across models (`EWAS_SUMMARY`).
- Differentially methylated regions (#24): DMRcate kernel smoothing on the EWAS statistics (no run-time annotation downloads) and dmrff (pinned source package) per tested model, and their consensus (overlap, same direction) with BED output (`DMR_CALL`); `--skip_dmr`, `--dmr_lambda`, `--dmr_c`, `--dmr_min_cpgs`, `--dmrff_source`.
- Gene set enrichment (#26): missMethyl `gometh` (GO; KEGG downloaded at run time, recorded as not run when unavailable) and `gsameth` on GMT collections (default MSigDB 2025.1 Hallmark and Reactome, staged files) for the CpGs at FDR < 0.05 of every tested model, against all tested CpGs, with the probe annotation passed as `anno` (`ENRICHMENT_MISSMETHYL`).
- DNA methylation scores (#28): smoking, alcohol and BMI (McCartney et al. 2018) and IL-6 (Stevenson et al. 2021) as weighted sums over the published CpG sets stored in `assets/dnam_scores` (CC BY 4.0, sources listed), on the detection-masked CpG-level betas with legacy-ID matching and coverage per score (`DNAM_SCORES`); `dnam_smoking` is now the `@smoking` model covariate. `--dnam_scores`, `--skip_dnam_scores`.
- HTML analysis report (#30): `REPORT_QUARTO` renders `report/methylarray_report.html` with Quarto from the result tables, software versions, methods text and parameters (`assets/report/methylarray_report.qmd`); `--skip_report`. MultiQC sections follow the pipeline order.
- GEO submission files (#32, `--export_geo`): processed matrix (normalised beta values with detection p-values), unnormalised methylated and unmethylated signal from the IDAT files, and a sample table with the IDAT file names, for the arrays that passed sample QC (`EXPORT_GEO`).
- Study-wide threshold across models (#32): `EWAS_SUMMARY` reports the effective number of independent tests over the tested models and contrasts (Meff, Galwey 2009, from the correlation of their t statistics), `p_study_wide` and its hits; shown in the report.
- PCA outlier flags (#32): `EXPLORATION_PCA` writes `pca_outliers.tsv`, a robust Mahalanobis distance on PC1 and PC2 with a flag above the chi-square 0.999 quantile. Report only; shown in MultiQC and the report.
- Clock, score and cell associations (`BIOMARKERS_ASSOCIATE`): agreement of the accelerations and scores, cell fractions against them and against age, sex and the phenotypes, accelerations against the DNA plate and BeadChip, every acceleration and score against every phenotype (plain and with cells held), and quartile groups of acceleration with trend tests; aggregates only (`biomarkers/`), with MultiQC heatmaps.
- Standardised weights of every weighted-sum clock and DNA methylation score (`clocks_weights.tsv`, `dnam_scores_weights.tsv`) and the number of clock CpGs with a common SNP.
- SNP genotype summaries from the ewastools SNP probes (`ewastools_snp_summary.tsv`, `ewastools_snp_histograms.tsv`): heterozygosity at autosomal and X-like SNPs, highest pairwise genotype agreement; MultiQC table.
- SNP relatedness (`EWASTOOLS_QC`): identity by state and likelihood ratios for identical and first-degree pairs of arrays from the SNP genotypes (`ewastools_relatedness.tsv.gz`, `ewastools_relatedness_summary.tsv`), the related arrays of each array in `sample_qc.tsv` (`relatives`), a MultiQC table of the related and replicate pairs and the identity-by-state plot in the report. Report only; the expected strength of evidence for a true first-degree pair is reported, since the SNP probes of the array detect identical pairs reliably but most relatives only with more SNPs.
- Contamination sensitivity (`BIOMARKERS_ASSOCIATE`): the ewastools SNP outlier score against every acceleration and score (`biomarkers_contamination.tsv`, MultiQC table), and the phenotype associations again without the analysed arrays above a robust cutoff (median plus three MADs), as the `low contamination` version in `biomarkers_phenotypes.tsv`; report plot of the score against the main clocks.
- Genetic ancestry components (`ANCESTRY_PCS`, after the EpiAnceR+ approach): PCA of the CpG probes with a common SNP, residualised for the control-probe PCs, cell fractions, sex and age, together with the rs-probe genotypes; `exploration/ancestry_pcs*.tsv`, the `@ancestry` model covariate token, `--n_ancestry_pcs` and `--skip_ancestry`, MultiQC scatter and scree plots. `BIOMARKERS_ASSOCIATE` relates every acceleration and score to the components (`biomarkers_ancestry.tsv`, MultiQC heatmap) and repeats the phenotype associations with them held (`ancestry held`).
- Genomic context enrichment (`KYCG_CACHE`, `ENRICHMENT_KYCG`): KnowYourCG on the CpGs at FDR < 0.05 of every tested model (all, higher and lower methylation) against all tested CpGs, in the EPIC v2 knowledgebases (chromatin states, TFBS, histone marks, CTCF, A/B compartments, imprinting DMRs, gene features) staged from ExperimentHub once per run; `enrichment/<model>/kycg_<model>.tsv`, MultiQC table and report dot plots; `--kycg_databases`, `--kycg_cache`, `--skip_kycg`.
- Methylation QTLs (`--meqtl`, `--meqtl_p`; `MEQTL_PREPARE`): a user-supplied meQTL table (for example GoDMC) gives the CpGs with a cis-meQTL, matched to EPIC v2 through the legacy IDs. The clock weights get a `meqtl` column, every clock is computed again with the CpGs that have a common SNP or a cis-meQTL held at their mean (`clocks_robust.tsv`, genotype-robust clocks) and the phenotype associations are repeated with them (`genotype-robust`, `biomarkers_robust.tsv`), and `EWAS_SUMMARY` reports the hits with a cis-meQTL per model (`ewas_summary_meqtl.tsv`). Report section on the genetic influence on the clocks.
- Variance partition (`VARIANCE_PARTITION`): variancePartition linear mixed models on the `--variance_partition_cpgs` most variable CpGs, with plate, BeadChip, Sentrix row, batch and sex as random effects and age and five cell lineages as fixed effects; `exploration/variance_partition*.tsv`, MultiQC bar graph and report violin plot; `--skip_variance_partition`.
- The Illumina control metrics of ewastools are in `sample_qc.tsv` (`ctrl_*`) and are tested against the principal components (group `control metrics` of `pca_associations.tsv`, the MultiQC heatmap and the report, where the heatmap is split by variable group).
- Illumina probe lists (`--illumina_flagged_probes`, `--illumina_mapping_inaccuracy`): user-supplied copies of Illumina's EPIC v2 flagged and mapping-inaccuracy probe files (not redistributable) are annotated in `probe_annotation.tsv.gz` and counted as report-only rows of the attrition table, never removed.
- IDOL cross-check of the cell fractions (`EPIDISH_RPC`, `--idol_reference`, `--skip_idol`): the six cell types of the FlowSorted.Blood.EPIC IDOL reference with Houseman's constrained projection (EpiDISH `CP`), as `idol_*` columns of `cell_fractions.tsv`, their agreement with the cent12CT lineages (`cell_fractions_idol_agreement.tsv`) and a report plot.
- PCA: the most variable CpGs are ranked by the variance of their beta values (the PCA itself stays on M-values), because ranking by M-value variance favours CpGs near 0 or 1; `pca_associations.tsv` adds the adjusted R-squared, which the report heatmap now shows, so that variables with many levels such as the BeadChip are not overstated.
- Report: SNP genotypes, clock and score weights, agreement and cell, batch, phenotype and quartile-group plots, R-squared heatmap of the principal components with their variance shares, and gene set enrichment dot plots; MultiQC tables of the top gene sets per model.
- Developer guide (`docs/developer.md`) and `dev/nfdev container`, which builds a module's images from its `environment.yml` with Wave.
- Every result reaches both reports. MultiQC adds the regions per DMR method, the meQTL table, the clock and score values (violin plots), the IDOL agreement, the batch tests of the accelerations and scores, and in the association summary the hits with a cis-meQTL and the flagged design checks. The HTML report adds the flagged design checks, the regions per DMR method, the age clocks against reported age, the distributions of the accelerations and scores, and PC1 against PC2.
- Cell-type-specific association for every model with a binary or continuous trait and blood cell fractions (`--skip_cell_ewas` to skip): the CellDMC design of EpiDISH, fitted for all CpGs at once with limma, in `ewas/<model>/ewas_<model>_cells.tsv.gz`, with hits per lineage in `ewas/ewas_summary_cells.tsv` and MultiQC.
- HTML report: a directional Manhattan plot (Uddin et al. 2022: -log10 P drawn upwards for higher and downwards for lower methylation, CpGs at FDR < 0.05 coloured) and tables of the top 20 CpGs with higher and with lower methylation, for every model and contrast and for every blood cell lineage of the cell-type-specific association.
- HTML report figures are publication-ready: one theme, labelled axes and legends, no titles on the images; each figure has a numbered caption with its description, and every figure is also saved as a 300 dpi PNG and a vector PDF in `report/figures/`.
- Interpretation of association results, after Uddin et al. (2022) and Kirmani et al. (2025): hits split by direction with a binomial test of the balance (`hits_higher`, `hits_lower`), volcano plots, the CpG island context and gene features of the hits by direction (Fisher's odds ratios), the methylation level of the hits, the characteristics of each model's arrays by trait group (`ewas_<model>_characteristics.tsv`, aggregates only), the overlap and t-statistic correlation of every pair of model contrasts (`ewas/ewas_summary_comparison.tsv` and MultiQC), and a design check that flags methylation-derived covariates that follow the trait.
- The task directories of a run and its cache of downloaded remote inputs are deleted once it has finished successfully (Nextflow `cleanup`); `--skip_cleanup` keeps them for `-resume`. The pipeline stops at the start when cleanup would break symlinked results (`--publish_dir_mode symlink` or `rellink`).
- Documentation website ([bhagesh-h.github.io/methylarray](https://bhagesh-h.github.io/methylarray/)), built with MkDocs Material from the README, `docs/`, `CITATIONS.md` and `CHANGELOG.md` (`mkdocs.yml`, `dev/docs/build_site.py`). Pull requests build it with broken-link checks; GitHub Pages publishes it only from `main`. Preview it locally with `dev/nfdev docs serve`.

### `Changed`

- Faster runs: `MINFI_PREPROCESS` processes its chunks on parallel worker processes (as many as the task's CPUs and memory allow); the per-module CPU and memory requests follow measured use, so the single-threaded steps (cell composition, PCA, probe annotation, association models, enrichment) request 2 CPUs and 12 GB and several models run side by side; with a `block`, the `duplicateCorrelation` consensus correlation is estimated on 50,000 random CpGs instead of every CpG.

- `MINFI_PREPROCESS` reads and corrects the arrays in chunks (hidden `--minfi_chunk_size`, default 40) and writes each chunk to HDF5 at once, which lowers its memory use. Results are unchanged: Noob, detection p-values, ENmix oob/RELIC/RCP and SWAN work array by array; Funnorm, quantile normalisation and the ENmix quantile step run on all arrays in a second pass, and the ENmix low-quality flag is computed over all arrays.
- minfi (Noob) is the default preprocessing route (`--preprocessing minfi`); SeSAMe stays available with `--preprocessing sesame`.
- CI tests run with Docker only, in up to 10 shards with a 30 minute limit per job; conda and Singularity jobs took longer and are checked locally with `dev/nfdev`. The pipeline download test runs only when started by hand. Stalled remote file reads fail after five minutes, and a shard whose image or file download broke off runs once more. CI downloads the public GEO test IDATs once per run (cached between runs) and the test jobs read local copies, because NCBI refuses part of the requests when many jobs fetch the same files.
- Task images are built from each module's `environment.yml`: by Wave in Seqera Containers (`community.wave.seqera.io`), with separate Docker and Singularity images, or, where Wave could not build them (the minfi and enrichment environments exceed its build time limit; the exploration, EWAS, DMR and report images were built locally while Wave was rate-limited), as Docker Hub images referenced by digest. `EXPORT_GEO` reuses the minfi image.

### `Fixed`

- Numeric Sentrix IDs, plates and batches in the samplesheet are accepted (#3).
- The test profiles stage the GEO test IDATs four at a time to avoid HTTP 503 responses (#7).
- MultiQC heatmaps with a single row or column (for example sex as the only phenotype, without a phenotype table or ages) are written with list categories, so MultiQC accepts them.
- The pipeline download test uses Apptainer 1.4.5, which converts task images whose OCI manifest holds a Docker-typed layer (Apptainer 1.3 refuses them).
- HTML report: wide tables (cell reference coverage, consensus DMRs) no longer run under the table of contents, which stays in the right margin; tables scroll sideways within the text column. Small P-values print in scientific notation instead of as 0.
- Control-probe PCs no longer carry sex. Some normalisation controls differ between women and men, which could give a PC that followed predicted sex. Models with both `@sex` and `@ctrlpcs` then had collinear covariates. Predicted sex is now removed from every control probe before the PCA; `pca_ctrl_pca_variance.tsv` records it in `sex_removed`.
- Cell-type-specific association: lineages with a median fraction below 5% or a genomic lambda above 1.3 are marked unreliable (`lineage_status`, `status`). With a small fraction, the interaction estimate of a lineage is unstable and can give thousands of false hits, even for traits without a hit in the main model. The HTML report lists them without plots or top tables. The limits are the `EWAS_LIMMA` options `--min_cell_fraction` and `--max_cell_lambda` (`ext.args`).

### `Dependencies`

- Nextflow 26.04.2 or later (#4).
- minfi 1.56.0, ENmix 1.46.2, IlluminaHumanMethylationEPICv2manifest 1.0.0 and IlluminaHumanMethylationEPICv2anno.20a1.hg38 1.0.0.
- ewastools 1.8.2 (GitHub commit d7d0a32, installed from source by the QC task), illuminaio 0.52.0 and data.table 1.18.6.1.
- SeSAMe 1.28.1, sesameData 1.28.0, ExperimentHub 3.0.0, SummarizedExperiment 1.40.0 and HDF5Array 1.38.0 on R 4.5.3, with the reference (netlib) BLAS.
- EpiDISH 2.26.0, isva 1.10, limma 3.68.5, bacon 1.38.0, DMRcate 3.6.0, missMethyl 1.44.0 with GO.db 3.22.0 and org.Hs.eg.db 3.22.0, Quarto 1.11.9 and MultiQC 1.35.
- knowYourCG 1.6.3 (AGPL-3.0) with ExperimentHub 3.0.0 and sesameData 1.28.0, image `bhagesh/methylarray-kycg`.
- variancePartition 1.40.2 with lme4, image `bhagesh/methylarray-variancepartition`.
- From pinned source packages: dnaMethyAge 0.2.0, DunedinPACE 0.99.0 and dmrff 1.1.2; MSigDB release 2025.1 GMT files.

### `Deprecated`
