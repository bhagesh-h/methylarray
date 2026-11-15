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

# Control-probe intensities of both channels, for the control-probe PCs (one row
# per probe and channel, one column per array)
control_info <- raw[["controls"]]
control_table <- data.table(
    channel = rep(c("green", "red"), each = nrow(control_info)),
    group = rep(control_info[["group"]], 2),
    name = rep(control_info[["name"]], 2),
    as.data.table(rbind(raw[["ctrlG"]], raw[["ctrlR"]]))
)
setnames(control_table, 4:ncol(control_table), samples[["sample"]])
fwrite(control_table, paste0(prefix, "_controls.tsv.gz"), sep = "\\t")

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
snp_ids <- raw[["manifest"]][probe_type == "rs", probe_id]
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

# SNP summaries (aggregates only, no genotypes are written): genotypes called
# from the SNP probe betas (below 0.2, 0.2 to 0.8, above 0.8 = 0, 1, 2 B
# alleles). SNPs at which men are almost never heterozygous (below 2%) and
# women often are (above 10%) behave as X-chromosome SNPs; sex is the reported
# sex, or the predicted sex where none is reported. Per array: share of
# heterozygous calls at the autosomal and the X-like SNPs, and the highest share
# of identical autosomal genotypes with any other array. Histograms of the SNP
# betas and of the pairwise agreement.
geno <- ifelse(snps < 0.2, 0L, ifelse(snps > 0.8, 2L, 1L))
rs <- sub("_.*", "", snp_ids)
sex <- ifelse(samples[["reported_sex"]] %in% c("M", "F"), samples[["reported_sex"]], c(male = "M", female = "F")[predicted])
x_like <- rep(FALSE, nrow(geno))
if (sum(sex %in% "M") >= 5 && sum(sex %in% "F") >= 5) {
    het_men <- rowMeans(geno[, sex %in% "M", drop = FALSE] == 1L, na.rm = TRUE)
    het_women <- rowMeans(geno[, sex %in% "F", drop = FALSE] == 1L, na.rm = TRUE)
    x_like <- !is.na(het_men) & !is.na(het_women) & het_men < 0.02 & het_women > 0.10
}
autosomal <- !x_like & !duplicated(rs)
g <- geno[autosomal, , drop = FALSE]
g[is.na(g)] <- -1L
agreement <- (crossprod(g == 0L) + crossprod(g == 1L) + crossprod(g == 2L)) / max(1, nrow(g))
diag(agreement) <- NA
max_agreement <- if (n > 1) apply(agreement, 2, max, na.rm = TRUE) else rep(NA_real_, n)
x_rows <- x_like & !duplicated(rs)
snp_summary <- data.frame(
    sample = samples[["sample"]],
    sex_used = sex,
    snp_outlier = round(snp_outlier, 3),
    het_share_autosomal = round(colMeans(geno[autosomal, , drop = FALSE] == 1L, na.rm = TRUE), 4),
    het_share_x = if (any(x_rows)) round(colMeans(geno[x_rows, , drop = FALSE] == 1L, na.rm = TRUE), 4) else NA_real_,
    max_agreement = round(max_agreement, 4),
    snps_autosomal = sum(autosomal),
    snps_x_like = sum(x_rows)
)
write.table(snp_summary, paste0(prefix, "_snp_summary.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, na = "")
histogram <- function(x, breaks, kind) {
    h <- hist(x[is.finite(x)], breaks = breaks, plot = FALSE)
    data.frame(histogram = kind, bin_start = head(breaks, -1), bin_end = tail(breaks, -1), count = h[["counts"]])
}
write.table(rbind(histogram(c(snps[autosomal, ]), seq(0, 1, by = 0.02), "snp_beta"),
                  histogram(agreement[upper.tri(agreement)], seq(0, 1, by = 0.02), "pair_agreement")),
            paste0(prefix, "_snp_histograms.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE)

# MultiQC: SNP summary per array
writeLines(c(
    "# id: 'snp_genotypes'",
    "# section_name: 'SNP genotypes'",
    "# description: 'From the SNP probes of the array (ewastools). Outlier: contamination score, higher (towards 0) is worse. Heterozygous: share of heterozygous calls at the autosomal SNPs and at the SNPs that behave as X-chromosome SNPs (near 0 in men). Max agreement: highest share of identical autosomal genotypes with any other array; values near 1 mean the same person.'",
    "# plot_type: 'table'"
), paste0(prefix, "_snp_mqc.tsv"))
write.table(data.frame(Sample = snp_summary[["sample"]], `SNP outlier` = snp_summary[["snp_outlier"]],
                       `Heterozygous (%)` = round(100 * snp_summary[["het_share_autosomal"]], 1),
                       `Heterozygous X (%)` = round(100 * snp_summary[["het_share_x"]], 1),
                       `Max agreement (%)` = round(100 * snp_summary[["max_agreement"]], 1), check.names = FALSE),
            paste0(prefix, "_snp_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE, na = "")
write.table(conflict_table, paste0(prefix, "_snp_conflicts.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE)

# Software versions
writeLines(c('"${task.process}":',
    paste0("    ewastools: ", packageVersion("ewastools", lib.loc = library_dir)),
    paste0("    illuminaio: ", packageVersion("illuminaio"))), "versions.yml")
