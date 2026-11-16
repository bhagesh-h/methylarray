# bhagesh-h/methylarray: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
- Report: SNP genotypes, clock and score weights, agreement and cell, batch, phenotype and quartile-group plots, R-squared heatmap of the principal components with their variance shares, and gene set enrichment dot plots; MultiQC tables of the top gene sets per model.
- Developer guide (`docs/developer.md`) and `dev/nfdev container`, which builds a module's images from its `environment.yml` with Wave.

### `Changed`

- Faster runs: `MINFI_PREPROCESS` processes its chunks on parallel worker processes (as many as the task's CPUs and memory allow); the per-module CPU and memory requests follow measured use, so the single-threaded steps (cell composition, PCA, probe annotation, association models, enrichment) request 2 CPUs and 12 GB and several models run side by side; with a `block`, the `duplicateCorrelation` consensus correlation is estimated on 50,000 random CpGs instead of every CpG.

- `MINFI_PREPROCESS` reads and corrects the arrays in chunks (hidden `--minfi_chunk_size`, default 40) and writes each chunk to HDF5 at once, which lowers its memory use. Results are unchanged: Noob, detection p-values, ENmix oob/RELIC/RCP and SWAN work array by array; Funnorm, quantile normalisation and the ENmix quantile step run on all arrays in a second pass, and the ENmix low-quality flag is computed over all arrays.
- minfi (Noob) is the default preprocessing route (`--preprocessing minfi`); SeSAMe stays available with `--preprocessing sesame`.
- Task images are built from each module's `environment.yml`: by Wave in Seqera Containers (`community.wave.seqera.io`), with separate Docker and Singularity images, or, where Wave could not build them (the minfi and enrichment environments exceed its build time limit; the exploration, EWAS, DMR and report images were built locally while Wave was rate-limited), as Docker Hub images referenced by digest. `EXPORT_GEO` reuses the minfi image.

### `Fixed`

- Numeric Sentrix IDs, plates and batches in the samplesheet are accepted (#3).
- The test profiles stage the GEO test IDATs four at a time to avoid HTTP 503 responses (#7).
- The pipeline download test uses Apptainer 1.4.5, which converts task images whose OCI manifest holds a Docker-typed layer (Apptainer 1.3 refuses them).

### `Dependencies`

- Nextflow 26.04.2 or later (#4).
- minfi 1.56.0, ENmix 1.46.2, IlluminaHumanMethylationEPICv2manifest 1.0.0 and IlluminaHumanMethylationEPICv2anno.20a1.hg38 1.0.0.
- ewastools 1.8.2 (GitHub commit d7d0a32, installed from source by the QC task), illuminaio 0.52.0 and data.table 1.18.6.1.
- SeSAMe 1.28.1, sesameData 1.28.0, ExperimentHub 3.0.0, SummarizedExperiment 1.40.0 and HDF5Array 1.38.0 on R 4.5.3, with the reference (netlib) BLAS.
- EpiDISH 2.26.0, isva 1.10, limma 3.68.5, bacon 1.38.0, DMRcate 3.6.0, missMethyl 1.44.0 with GO.db 3.22.0 and org.Hs.eg.db 3.22.0, Quarto 1.11.9 and MultiQC 1.35.
- From pinned source packages: dnaMethyAge 0.2.0, DunedinPACE 0.99.0 and dmrff 1.1.2; MSigDB release 2025.1 GMT files.

### `Deprecated`
