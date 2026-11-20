#!/usr/bin/env Rscript
#
# Blood cell-type fractions with EpiDISH robust partial correlations (RPC) on
# the filtered, CpG-level, unmasked beta values:
#   - cent12CT.m: 12 leukocyte types (Luo et al. 2023), aggregated to six
#     lineages (B, CD4T, CD8T, NK, Mono, Gran)
#   - centDHSbloodDMC.m: 7 blood cell types (Teschendorff et al. 2017), for
#     comparison
#   - with an IDOL reference (--idol_reference): the 450 IDOL-optimised CpGs
#     of FlowSorted.Blood.EPIC (Salas et al. 2018) for six cell types, with
#     Houseman's constrained projection (EpiDISH method CP, inequality
#     constraint, as projectCellType_CP), as a cross-check of the cent12CT
#     lineages; a coverage below the minimum skips it with a note instead of
#     stopping the run
# Reference CpGs are matched through the EPIC v1 and 450K IDs of the probe
# annotation first, then by CpG ID, so that CpGs renamed on EPIC v2 are found.
# Coverage is computed on the autosomal reference CpGs: sex-chromosome CpGs are
# left out on purpose (X inactivation in women can mimic a cell-type signal).
# Every missing reference CpG is counted as on a sex chromosome, not on the
# array, or removed by the probe filter. Missing values are set to the CpG
# median over the arrays and counted.
#

prefix <- "${prefix}"
se_dir <- "${se}"
annotation_file <- "${annotation}"
idol_file <- "${idol_reference}"

defaults <- list(min_coverage = "0.9", warn_coverage = "0.98", maxit = "50")
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
min_coverage <- as.numeric(opts[["min_coverage"]])
warn_coverage <- as.numeric(opts[["warn_coverage"]])
maxit <- as.integer(opts[["maxit"]])

Sys.setenv(HOME = getwd())
suppressPackageStartupMessages({
    library(EpiDISH)
    library(SummarizedExperiment)
    library(HDF5Array)
    library(data.table)
})
set.seed(1)

se <- loadHDF5SummarizedExperiment(se_dir)
beta <- as.matrix(assay(se, "beta_unmasked"))
row_data <- rowData(se)
legacy <- function(column) if (column %in% names(row_data)) as.character(row_data[[column]]) else rep(NA_character_, nrow(se))
epicv1_id <- legacy("epicv1_id")
hm450_id <- legacy("hm450_id")

# Row of each reference CpG: EPIC v1 ID, then 450K ID, then CpG ID
match_reference <- function(ids) {
    by_epicv1 <- match(ids, epicv1_id)
    by_hm450 <- match(ids, hm450_id)
    by_cpg <- match(ids, rownames(se))
    row <- ifelse(!is.na(by_epicv1), by_epicv1, ifelse(!is.na(by_hm450), by_hm450, by_cpg))
    list(row = row, legacy_only = !is.na(row) & is.na(by_cpg))
}

# Probe annotation: which reference CpGs the array measures, and where
# (fread reads .gz files only with R.utils; base R decompresses instead)
annotation <- fread(text = readLines(gzfile(annotation_file)), na.strings = "",
                    select = c("cpg_id", "chr", "epicv1_id", "hm450_id"), colClasses = "character")
array_chr <- function(ids) {
    hit <- function(column) annotation[["chr"]][match(ids, annotation[[column]])]
    on_array <- !is.na(match(ids, annotation[["cpg_id"]])) | !is.na(match(ids, annotation[["epicv1_id"]])) |
        !is.na(match(ids, annotation[["hm450_id"]]))
    chr <- hit("epicv1_id")
    chr[is.na(chr)] <- hit("hm450_id")[is.na(chr)]
    chr[is.na(chr)] <- hit("cpg_id")[is.na(chr)]
    list(on_array = on_array, chr = chr)
}

data(cent12CT.m, package = "EpiDISH", envir = environment())
data(centDHSbloodDMC.m, package = "EpiDISH", envir = environment())
references <- list(cent12CT = cent12CT.m, centDHSbloodDMC = centDHSbloodDMC.m)
methods <- list(cent12CT = "RPC", centDHSbloodDMC = "RPC")
if (nzchar(idol_file)) {
    idol_env <- new.env()
    idol_name <- load(idol_file, envir = idol_env)[1]
    references[["IDOL"]] <- as.matrix(get(idol_name, envir = idol_env))
    methods[["IDOL"]] <- "CP"
}

coverage <- list()
fractions <- list()
for (name in names(references)) {
    ref <- references[[name]]
    matched <- match_reference(rownames(ref))
    found <- !is.na(matched[["row"]])
    placed <- array_chr(rownames(ref))
    sex_chr <- !found & placed[["chr"]] %in% c("chrX", "chrY")
    not_on_array <- !found & !sex_chr & !placed[["on_array"]]
    filtered <- !found & !sex_chr & placed[["on_array"]]
    ref_beta <- beta[matched[["row"]][found], , drop = FALSE]
    rownames(ref_beta) <- rownames(ref)[found]
    missing_values <- sum(is.na(ref_beta))
    if (missing_values > 0) {
        medians <- apply(ref_beta, 1, median, na.rm = TRUE)
        idx <- which(is.na(ref_beta), arr.ind = TRUE)
        ref_beta[idx] <- medians[idx[, 1]]
    }
    share <- sum(found) / (nrow(ref) - sum(sex_chr))
    coverage[[name]] <- data.table(
        reference = name, reference_cpgs = nrow(ref), sex_chromosome = sum(sex_chr),
        not_on_array = sum(not_on_array), removed_by_filter = sum(filtered), found = sum(found),
        found_by_legacy_id = sum(matched[["legacy_only"]]), coverage = round(share, 4),
        imputed_values = missing_values,
        status = if (share < min_coverage) "fail" else if (share < warn_coverage) "warn" else "ok"
    )
    if (share < min_coverage && name == "IDOL") {
        coverage[[name]][["status"]] <- "fail: skipped"
        next
    }
    if (share < min_coverage) {
        fwrite(rbindlist(coverage), paste0(prefix, "_reference_coverage.tsv"), sep = "\\t")
        stop(sprintf("Only %.1f%% of the autosomal %s reference CpGs are in the data (minimum %.0f%%)", 100 * share, name, 100 * min_coverage))
    }
    estimates <- if (methods[[name]] == "CP") {
        epidish(beta.m = ref_beta, ref.m = ref[found, , drop = FALSE], method = "CP", constraint = "inequality")[["estF"]]
    } else {
        epidish(beta.m = ref_beta, ref.m = ref[found, , drop = FALSE], method = "RPC", maxit = maxit)[["estF"]]
    }
    fractions[[name]] <- estimates[colnames(se), , drop = FALSE]
}

f12 <- fractions[["cent12CT"]]
lineages <- cbind(
    B = f12[, "Bnv"] + f12[, "Bmem"],
    CD4T = f12[, "CD4Tnv"] + f12[, "CD4Tmem"] + f12[, "Treg"],
    CD8T = f12[, "CD8Tnv"] + f12[, "CD8Tmem"],
    NK = f12[, "NK"],
    Mono = f12[, "Mono"],
    Gran = f12[, "Neu"] + f12[, "Eos"] + f12[, "Baso"]
)
out <- data.table(
    sample = colnames(se),
    setNames(as.data.table(round(lineages, 5)), paste0("lineage_", colnames(lineages))),
    setNames(as.data.table(round(f12, 5)), paste0("cent12_", colnames(f12))),
    setNames(as.data.table(round(fractions[["centDHSbloodDMC"]], 5)), paste0("dhs7_", colnames(fractions[["centDHSbloodDMC"]]))),
    sum_cent12 = round(rowSums(f12), 4),
    sum_dhs7 = round(rowSums(fractions[["centDHSbloodDMC"]]), 4)
)
# IDOL cross-check: fractions per array and their agreement with the lineages
idol_columns <- c(B = "Bcell", CD4T = "CD4T", CD8T = "CD8T", NK = "NK", Mono = "Mono", Gran = "Neu")
agreement <- data.table(lineage = character(), idol_type = character(), rho = numeric(), mean_difference = numeric(), n = integer())
if (!is.null(fractions[["IDOL"]])) {
    idol <- fractions[["IDOL"]]
    out <- cbind(out, setNames(as.data.table(round(idol, 5)), paste0("idol_", colnames(idol))))
    agreement <- rbindlist(lapply(names(idol_columns), function(l) {
        if (!idol_columns[[l]] %in% colnames(idol)) return(NULL)
        a <- lineages[, l]
        b <- idol[, idol_columns[[l]]]
        data.table(lineage = l, idol_type = idol_columns[[l]], rho = round(suppressWarnings(cor(a, b, method = "spearman")), 4),
                   mean_difference = round(mean(a - b), 4), n = length(a))
    }))
}
fwrite(agreement, paste0(prefix, "_idol_agreement.tsv"), sep = "\\t")
setorder(out, sample)
fwrite(out, paste0(prefix, ".tsv"), sep = "\\t")
fwrite(rbindlist(coverage), paste0(prefix, "_reference_coverage.tsv"), sep = "\\t")

# MultiQC: lineage fractions per array as a stacked bar plot
cov12 <- coverage[["cent12CT"]]
writeLines(c(
    "# id: 'cell_fractions'",
    "# section_name: 'Cell composition'",
    paste0("# description: 'Blood cell lineages from EpiDISH RPC with the 12-type reference cent12CT (",
           cov12[["found"]], " of ", cov12[["reference_cpgs"]] - cov12[["sex_chromosome"]], " autosomal reference CpGs in the data). Gran: neutrophils, eosinophils and basophils.'"),
    "# plot_type: 'bargraph'",
    "# pconfig:",
    "#     id: 'cell_fractions_plot'",
    "#     title: 'Cell composition: blood cell lineages'",
    "#     ylab: 'Fraction'",
    "#     cpswitch: false"
), paste0(prefix, "_mqc.tsv"))
mqc <- data.frame(Sample = out[["sample"]], out[, paste0("lineage_", colnames(lineages)), with = FALSE], check.names = FALSE)
names(mqc) <- sub("^lineage_", "", names(mqc))
write.table(mqc,
            paste0(prefix, "_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE)
if (nrow(agreement) > 0) {
    # MultiQC: agreement of the lineages with the IDOL cross-check
    writeLines(c(
        "# id: 'cell_fractions_idol'",
        "# section_name: 'Cell composition: IDOL cross-check'",
        "# description: 'Each cent12CT lineage against the matching IDOL cell type (FlowSorted.Blood.EPIC, constrained projection): Spearman rho over the arrays and the mean difference of the fractions. Gran includes eosinophils and basophils; IDOL neutrophils do not.'",
        "# plot_type: 'table'"
    ), paste0(prefix, "_idol_mqc.tsv"))
    write.table(data.frame(Lineage = agreement[["lineage"]], `IDOL type` = agreement[["idol_type"]], rho = agreement[["rho"]],
                           `Mean difference` = agreement[["mean_difference"]], Arrays = agreement[["n"]], check.names = FALSE),
                paste0(prefix, "_idol_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE)
}

# Software versions
writeLines(c('"${task.process}":',
    paste0("    EpiDISH: ", packageVersion("EpiDISH")),
    paste0("    SummarizedExperiment: ", packageVersion("SummarizedExperiment")),
    paste0("    HDF5Array: ", packageVersion("HDF5Array"))), "versions.yml")
