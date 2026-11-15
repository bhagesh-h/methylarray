# bhagesh-h/methylarray: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

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
- DNA methylation scores (M11): smoking, alcohol and BMI (McCartney et al. 2018) and IL-6 (Stevenson et al. 2021) as weighted sums over the published CpG sets stored in `assets/dnam_scores` (CC BY 4.0, sources listed), on the detection-masked CpG-level betas with legacy-ID matching and coverage per score (`DNAM_SCORES`); `dnam_smoking` is now the `@smoking` model covariate. `--dnam_scores`, `--skip_dnam_scores`.
- Developer guide (`docs/developer.md`) and `dev/nfdev container`, which builds a module's images from its `environment.yml` with Wave.

### `Changed`

- minfi (Noob) is the default preprocessing route (`--preprocessing minfi`); SeSAMe stays available with `--preprocessing sesame`.
- Task images come from Seqera Containers (`community.wave.seqera.io`), built by Wave from each module's `environment.yml`, with separate Docker and Singularity images. The minfi module keeps its Docker Hub image, because installing the EPIC v2 annotation packages exceeds the Wave build time limit.

### `Fixed`

- Numeric Sentrix IDs, plates and batches in the samplesheet are accepted (#3).
- The test profiles stage the GEO test IDATs four at a time to avoid HTTP 503 responses (#7).

### `Dependencies`

- Nextflow 26.04.2 or later (#4).
- minfi 1.56.0, ENmix 1.46.2, IlluminaHumanMethylationEPICv2manifest 1.0.0 and IlluminaHumanMethylationEPICv2anno.20a1.hg38 1.0.0.
- ewastools 1.8.2 (GitHub commit d7d0a32, installed from source by the QC task), illuminaio 0.52.0 and data.table 1.18.6.1.
- SeSAMe 1.28.1, sesameData 1.28.0, ExperimentHub 3.0.0, SummarizedExperiment 1.40.0 and HDF5Array 1.38.0 on R 4.5.3, with the reference (netlib) BLAS.

### `Deprecated`
