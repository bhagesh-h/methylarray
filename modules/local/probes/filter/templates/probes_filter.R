#!/usr/bin/env Rscript
#
# Probe filtering on the arrays that passed sample QC (PASS or FLAG).
# Exclusions, applied in this order and each counted in the attrition table:
#   1. probe type: rs (genotyping), nv (somatic variant) and control probes
#   2. non-CG (ch) probes
#   3. no hg38 position
#   4. chrX and chrY (unless --keep_sex_chr true)
#   5. detection: detection p-value above --detection_pval in more than
#      --probe_fail_fraction of the arrays
#   6. other masked values (bead count, quality mask) in more than
#      --bead_fail_fraction of the arrays
#   7. SNP mask of --snp_mask
#   8. Zhou general mask
#   9. cross-hybridisation (--crosshyb_filter wgbs: WGBS evidence; blat: any
#      BLAT prediction; none)
# Report-only flags (counted among the kept probes, never removed): manifest
# position mismatches (MismatchPos) and poor titration (Zhou).
# Remaining individual failures stay NA. Replicate probes are collapsed to the
# CpG (--replicate_probe_strategy mean, or none) and M-values are computed from
# betas clipped to [1e-4, 1 - 1e-4].
# The uncollapsed output skips step 9, so that DMR callers can remap
# cross-hybridising probes to their off-target locus.
#

prefix <- "${prefix}"
se_dir <- "${se}"
retained_file <- "${retained}"
annotation_file <- "${annotation}"

defaults <- list(detection_pval = "0.05", probe_fail_fraction = "0.1", bead_fail_fraction = "0.05",
                 crosshyb_filter = "wgbs", snp_mask = "common_1pt", keep_sex_chr = "false",
                 replicate_probe_strategy = "mean")
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
detection_pval <- as.numeric(opts[["detection_pval"]])
probe_fail_fraction <- as.numeric(opts[["probe_fail_fraction"]])
bead_fail_fraction <- as.numeric(opts[["bead_fail_fraction"]])
crosshyb_filter <- match.arg(opts[["crosshyb_filter"]], c("wgbs", "blat", "none"))
snp_mask <- match.arg(opts[["snp_mask"]], c("common_1pt", "common_5pt", "afr", "amr", "eas", "eur", "sas", "none"))
keep_sex_chr <- tolower(opts[["keep_sex_chr"]]) == "true"
replicate_strategy <- match.arg(opts[["replicate_probe_strategy"]], c("mean", "none"))

Sys.setenv(HOME = getwd())
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(HDF5Array)
    library(data.table)
})

se <- loadHDF5SummarizedExperiment(se_dir)
retained <- readLines(retained_file)
se <- se[, colnames(se) %in% retained]
if (ncol(se) == 0) stop("No array passed sample QC")
probes <- rownames(se)

# fread reads .gz files only with R.utils; base R decompresses instead
anno <- fread(text = readLines(gzfile(annotation_file)), na.strings = "", colClasses = list(character = c("chr", "suggested_offtarget", "epicv1_id", "hm450_id")))
anno <- anno[match(probes, anno[["probe_id"]])]
anno[, probe_id := probes]
anno[is.na(cpg_id), cpg_id := sub("_.*", "", probe_id)]
anno[is.na(probe_type), probe_type := substr(probe_id, 1, 2)]
for (column in names(anno)[vapply(anno, is.logical, logical(1))]) set(anno, which(is.na(anno[[column]])), column, FALSE)

beta <- as.matrix(assay(se, "beta"))
beta_unmasked <- as.matrix(assay(se, "beta_unmasked"))
detp <- as.matrix(assay(se, "detp"))
undetected <- rowMeans(is.na(detp) | detp > detection_pval)
other_masked <- rowMeans(is.na(beta) & !(is.na(detp) | detp > detection_pval))

# Every reason that applies to a probe, in filter order
snp_column <- if (snp_mask == "none") NULL else paste0("snp_", snp_mask)
crosshyb <- switch(crosshyb_filter, wgbs = anno[["ch_wgbs"]], blat = anno[["ch_blat"]], none = rep(FALSE, nrow(anno)))
reasons <- list(
    probe_type = !anno[["probe_type"]] %in% c("cg", "ch"),
    non_cg = anno[["probe_type"]] %in% "ch",
    no_position = is.na(anno[["pos"]]),
    sex_chromosome = if (keep_sex_chr) rep(FALSE, nrow(anno)) else anno[["chr"]] %in% c("chrX", "chrY"),
    detection = undetected > probe_fail_fraction,
    bead_or_quality_mask = other_masked > bead_fail_fraction,
    snp = if (is.null(snp_column)) rep(FALSE, nrow(anno)) else anno[[snp_column]],
    zhou_general_mask = anno[["zhou_general"]],
    cross_hybridisation = crosshyb
)
report_only <- list(mismatch_pos = anno[["mismatch_pos"]], zhou_uncorr_titration = anno[["zhou_uncorr_titration"]])

reason_matrix <- do.call(cbind, reasons)
excluded <- rowSums(reason_matrix) > 0
first_reason <- ifelse(excluded, colnames(reason_matrix)[max.col(reason_matrix, ties.method = "first")], NA_character_)

# Attrition: probes removed at each step, in order
remaining <- nrow(anno)
attrition <- data.table(step = "start", removed = 0L, remaining = remaining, n_with_reason = 0L)
for (step in names(reasons)) {
    removed <- sum(first_reason %in% step)
    remaining <- remaining - removed
    attrition <- rbind(attrition, data.table(step = step, removed = removed, remaining = remaining,
                                             n_with_reason = sum(reasons[[step]])))
}
kept <- !excluded
for (flag in names(report_only)) {
    attrition <- rbind(attrition, data.table(step = paste0("report_only_", flag), removed = 0L,
                                             remaining = remaining, n_with_reason = sum(report_only[[flag]] & kept)))
}
fwrite(attrition, paste0(prefix, "_attrition.tsv"), sep = "\\t")

exclusions <- data.table(
    probe_id = probes[excluded],
    first_reason = first_reason[excluded],
    reasons = apply(reason_matrix[excluded, , drop = FALSE], 1, function(row) paste(colnames(reason_matrix)[row], collapse = ";")),
    undetected_fraction = round(undetected[excluded], 4),
    other_masked_fraction = round(other_masked[excluded], 4)
)
fwrite(exclusions, paste0(prefix, "_exclusions.tsv.gz"), sep = "\\t")

col_data <- colData(se)
meta <- c(metadata(se), list(
    probe_filter = list(detection_pval = detection_pval, probe_fail_fraction = probe_fail_fraction,
                        bead_fail_fraction = bead_fail_fraction, crosshyb_filter = crosshyb_filter,
                        snp_mask = snp_mask, keep_sex_chr = keep_sex_chr,
                        replicate_probe_strategy = replicate_strategy)
))
row_columns <- c("probe_id", "cpg_id", "chr", "pos", "strand", "design_type", "replicate_group_size", "replicate_choice",
                 "ch_wgbs", "ch_blat", "suggested_offtarget", "mismatch_pos", "epicv1_id", "hm450_id", "gene", "gene_group", "cgi_relation")
row_columns <- intersect(row_columns, names(anno))

save_se <- function(assays, row_data, dir) {
    scratch <- paste0(dir, "_scratch.h5")
    on_disk <- lapply(names(assays), function(name) writeHDF5Array(assays[[name]], filepath = scratch, name = name))
    names(on_disk) <- names(assays)
    out <- SummarizedExperiment(assays = on_disk, colData = col_data, rowData = row_data, metadata = meta)
    saveHDF5SummarizedExperiment(out, dir = dir, replace = TRUE)
    unlink(scratch)
}

# Uncollapsed: every filter except cross-hybridisation, for DMR callers that remap those probes
keep_uncollapsed <- rowSums(reason_matrix[, setdiff(colnames(reason_matrix), "cross_hybridisation"), drop = FALSE]) == 0
save_se(
    list(beta = beta[keep_uncollapsed, , drop = FALSE], beta_unmasked = beta_unmasked[keep_uncollapsed, , drop = FALSE]),
    DataFrame(as.data.frame(anno[keep_uncollapsed, ..row_columns]), row.names = probes[keep_uncollapsed]),
    paste0(prefix, "_uncollapsed_se")
)

# Collapsed to the CpG: mean of the replicate probes
clip <- function(x) pmin(pmax(x, 1e-4), 1 - 1e-4)
group <- anno[["cpg_id"]][kept]
collapse_mean <- function(x) {
    present <- !is.na(x)
    sums <- rowsum(ifelse(present, x, 0), group, reorder = FALSE)
    counts <- rowsum(present * 1, group, reorder = FALSE)
    out <- sums / counts
    out[counts == 0] <- NA_real_
    out
}
if (replicate_strategy == "mean") {
    beta_out <- collapse_mean(beta[kept, , drop = FALSE])
    unmasked_out <- collapse_mean(beta_unmasked[kept, , drop = FALSE])
    row_anno <- anno[kept][, .SD[1], by = cpg_id][, n_probes := anno[kept][, .N, by = cpg_id][["N"]]]
    row_anno <- row_anno[match(rownames(beta_out), cpg_id)]
    row_ids <- rownames(beta_out)
} else {
    beta_out <- beta[kept, , drop = FALSE]
    unmasked_out <- beta_unmasked[kept, , drop = FALSE]
    row_anno <- anno[kept][, n_probes := 1L]
    row_ids <- probes[kept]
}
mvalue <- log2(clip(beta_out) / (1 - clip(beta_out)))
collapsed_columns <- c(setdiff(row_columns, c("probe_id", "replicate_choice")), "n_probes")
save_se(
    list(beta = beta_out, mvalue = mvalue, beta_unmasked = unmasked_out),
    DataFrame(as.data.frame(row_anno[, ..collapsed_columns]), row.names = row_ids),
    paste0(prefix, "_se")
)

# MultiQC: probe attrition as a bar plot
mqc <- attrition[step != "start" & !startsWith(step, "report_only_")]
writeLines(c(
    "# id: 'probe_attrition'",
    "# section_name: 'Probe filtering'",
    paste0("# description: 'Probes removed at each filter step, in order (", nrow(anno), " probes at the start, ",
           remaining, " kept; ", nrow(beta_out), " CpGs after collapsing replicate probes). Probes are counted at the first step that removes them.'"),
    "# plot_type: 'bargraph'",
    "# pconfig:",
    "#     id: 'probe_attrition_plot'",
    "#     title: 'Probe filtering: probes removed per step'",
    "#     ylab: 'Probes'"
), paste0(prefix, "_attrition_mqc.tsv"))
write.table(data.frame(step = mqc[["step"]], removed = mqc[["removed"]]), paste0(prefix, "_attrition_mqc.tsv"),
            sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE)

# Software versions
writeLines(c('"${task.process}":',
    paste0("    SummarizedExperiment: ", packageVersion("SummarizedExperiment")),
    paste0("    HDF5Array: ", packageVersion("HDF5Array")),
    paste0("    data.table: ", packageVersion("data.table"))), "versions.yml")
