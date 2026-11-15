#!/usr/bin/env Rscript
#
# Sample QC with ewastools on all arrays, following the package's recommended
# workflow:
#   - Illumina BeadArray Controls Reporter metrics against Illumina's thresholds
#     (raw intensities)
#   - normalised X and Y intensities and the predicted sex (dye-bias corrected
#     intensities before detection masking, so the Y probes are kept)
#   - SNP probe genotypes (detection-masked at p > 0.01, dye-bias corrected):
#     contamination score, inferred donors and conflicts with the expected donors
#     (replicate_group, otherwise the sample itself)
# Nothing is removed here: SAMPLEQC_SUMMARY applies the thresholds.
#

prefix <- "${prefix}"
samples_file <- "${samples}"

defaults <- list(platform = "EPICv2", snp_outlier_threshold = "-4")
parse_args <- function(x, defaults) {
    tokens <- strsplit(trimws(x), "[[:space:]]+")[[1]]
    tokens <- tokens[nzchar(tokens)]
    if (length(tokens) %% 2 != 0) stop("ext.args must be pairs of '--name value': ", x)
    opts <- defaults
    for (i in seq_len(length(tokens) / 2) * 2 - 1) {
        name <- sub("^--", "", tokens[i])
        if (!name %in% names(defaults)) stop("Unknown option in ext.args: ", tokens[i])
        opts[[name]] <- tokens[i + 1]
    }
    opts
}
opts <- parse_args("${task.ext.args ?: ''}", defaults)
snp_outlier_threshold <- as.numeric(opts[["snp_outlier_threshold"]])

# ewastools is not on conda or CRAN: install the staged source (pure R, a pinned
# GitHub commit) into the task directory
Sys.setenv(HOME = getwd())
library_dir <- file.path(getwd(), "rlib")
dir.create(library_dir, showWarnings = FALSE)
install.packages("${ewastools_source}", lib = library_dir, repos = NULL, type = "source", quiet = TRUE)
suppressPackageStartupMessages({
    library(ewastools, lib.loc = library_dir)
    library(data.table)
})

# Arrays, sorted by sample name for a stable row order
samples <- read.delim(samples_file, colClasses = "character", na.strings = c("", "NA"), check.names = FALSE)
samples <- samples[order(samples[["sample"]]), , drop = FALSE]
n <- nrow(samples)

raw <- read_idats(file.path("idat", samples[["idat_prefix"]]), quiet = TRUE)
if (raw[["platform"]] != opts[["platform"]]) {
    stop("The IDAT files are ", raw[["platform"]], " arrays, but the run is set to --platform ", opts[["platform"]])
}

# Control metrics: one column per metric, and the metrics below Illumina's threshold
metrics <- control_metrics(raw)
# control_metrics matches 'NP (C)' and 'NP (G)' as name prefixes, which on EPIC v2
# also selects the dark 'NP (G) 1' to 'NP (G) 5' controls and fails every array.
# Recompute both non-polymorphic metrics from the four classic controls, as the
# Illumina BeadArray Controls Reporter defines them.
controls <- raw[["controls"]]
np_cg <- controls[group == "NON-POLYMORPHIC" & name %in% c("NP (C)", "NP (G)"), index]
np_at <- controls[group == "NON-POLYMORPHIC" & name %in% c("NP (A)", "NP (T)"), index]
if (length(np_cg) == 2 && length(np_at) == 2) {
    green <- raw[["ctrlG"]]
    red <- raw[["ctrlR"]]
    metrics[["Non-polymorphic Green"]] <- structure(
        apply(green[np_cg, , drop = FALSE], 2, min) / apply(green[np_at, , drop = FALSE], 2, max), threshold = 5)
    metrics[["Non-polymorphic Red"]] <- structure(
        apply(red[np_at, , drop = FALSE], 2, min) / apply(red[np_cg, , drop = FALSE], 2, max), threshold = 5)
}
below <- do.call(cbind, lapply(metrics, function(metric) metric < attr(metric, "threshold")))
below[is.na(below)] <- FALSE
metric_names <- names(metrics)
failed_controls <- apply(below, 1, function(row) paste(metric_names[row], collapse = ";"))
metric_table <- as.data.frame(lapply(metrics, function(metric) signif(as.numeric(metric), 4)))
names(metric_table) <- paste0("ctrl_", gsub("[^a-z0-9]+", "_", tolower(metric_names)))

# Sex: normalised X and Y intensities; predict_sex needs reported males and females
xy <- check_sex(correct_dye_bias(raw))
male <- which(samples[["reported_sex"]] %in% "M")
female <- which(samples[["reported_sex"]] %in% "F")
predicted <- rep(NA_character_, n)
if (length(male) > 0 && length(female) > 0) {
    predicted <- c(m = "male", f = "female")[as.character(predict_sex(xy[["X"]], xy[["Y"]], male, female))]
}

# SNP genotypes
raw <- correct_dye_bias(mask(detectionP(raw), 0.01))
snp_rows <- raw[["manifest"]][probe_type == "rs", index]
snps <- beta_values(raw)[snp_rows, , drop = FALSE]
rm(raw)
invisible(gc())
genotypes <- call_genotypes(snps, learn = FALSE)
snp_outlier <- snp_outliers(genotypes)
inferred_donor <- enumerate_sample_donors(genotypes)

# Identity conflicts among the arrays below the contamination threshold
expected_donor <- ifelse(is.na(samples[["replicate_group"]]), samples[["sample"]], samples[["replicate_group"]])
keep <- which(snp_outlier <= snp_outlier_threshold)
conflicts <- NULL
if (length(keep) >= 2) {
    subset_genotypes <- genotypes
    subset_genotypes[["snps"]] <- genotypes[["snps"]][, keep, drop = FALSE]
    subset_genotypes[["outliers"]] <- genotypes[["outliers"]][, keep, drop = FALSE]
    subset_genotypes[["gamma"]] <- lapply(genotypes[["gamma"]], function(g) g[, keep, drop = FALSE])
    conflicts <- check_snp_agreement(subset_genotypes, expected_donor[keep], samples[["sample"]][keep])
}
conflict_table <- data.table(group = integer(), sample1 = character(), expected_donor1 = character(),
                             sample2 = character(), expected_donor2 = character(), agreement = numeric())
if (!is.null(conflicts)) {
    conflict_table <- rbindlist(conflicts, idcol = "group")
    setnames(conflict_table, c("donor1", "donor2"), c("expected_donor1", "expected_donor2"))
    setcolorder(conflict_table, names(conflict_table)[c(1, 3, 2, 5, 4, 6)])
    conflict_table[, agreement := round(agreement, 4)]
}
in_conflict <- samples[["sample"]] %in% c(conflict_table[["sample1"]], conflict_table[["sample2"]])

qc <- data.frame(
    sample = samples[["sample"]],
    platform = opts[["platform"]],
    control_failed = rowSums(below) > 0,
    failed_controls = failed_controls,
    sex_x = round(xy[["X"]], 4),
    sex_y = round(xy[["Y"]], 4),
    predicted_sex = unname(predicted),
    snp_outlier = round(snp_outlier, 3),
    num_snp_probes = colSums(!is.na(snps)),
    inferred_donor = inferred_donor,
    identity_conflict = in_conflict,
    metric_table,
    check.names = FALSE
)
write.table(qc, paste0(prefix, "_qc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, na = "")
write.table(conflict_table, paste0(prefix, "_snp_conflicts.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE)

# Software versions
writeLines(c('"${task.process}":',
    paste0("    ewastools: ", packageVersion("ewastools", lib.loc = library_dir)),
    paste0("    illuminaio: ", packageVersion("illuminaio"))), "versions.yml")
