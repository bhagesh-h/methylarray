#!/usr/bin/env Rscript
#
# One table across the association models: N, genomic lambda, bacon inflation
# and bias, hits at FDR, the epigenome-wide and suggestive thresholds (raw and
# bacon-corrected), and the top CpG of each model and contrast.
# Study-wide threshold across the tested models and contrasts: the effective
# number of independent tests (Meff; Galwey 2009) from the eigenvalues of the
# correlation matrix of their t statistics over the CpGs, and
# p_study_wide = --p_epigenome_wide / Meff. Under the null the correlation of
# the statistics follows the correlation of the traits, so phenotype values are
# not needed.
# With a cis-meQTL list (--meqtl): per model and contrast, the CpGs at FDR <
# 0.05 that have a cis-meQTL, against the share among all tested CpGs (Fisher's
# exact test). Such hits are candidates for genetic confounding, especially
# without ancestry covariates.
#

prefix <- "${prefix}"
summary_files <- strsplit("${summaries}", " ")[[1]]
stats_files <- strsplit(trimws("${stats}"), "[[:space:]]+")[[1]]
stats_files <- stats_files[nzchar(stats_files)]
meqtl_file <- "${meqtl_cpgs}"

defaults <- list(p_epigenome_wide = "9e-8")
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
p_ewide <- as.numeric(opts[["p_epigenome_wide"]])

suppressPackageStartupMessages({
    library(jsonlite)
    library(data.table)
})
writeLines(c('"${task.process}":', paste0("    jsonlite: ", packageVersion("jsonlite"))), "versions.yml")

rows <- rbindlist(lapply(sort(summary_files), function(file) {
    s <- fromJSON(file, simplifyVector = FALSE)
    base <- list(model = s[["model"]], trait = s[["trait"]], n = s[["n"]], status = s[["status"]],
                 design_flags = if (is.null(s[["design_flags"]])) "" else s[["design_flags"]])
    if (is.null(s[["contrasts"]]) || length(s[["contrasts"]]) == 0) return(as.data.table(base))
    rbindlist(lapply(names(s[["contrasts"]]), function(contrast) {
        as.data.table(c(base, list(contrast = contrast), s[["contrasts"]][[contrast]]))
    }), fill = TRUE)
}), fill = TRUE)

# Cell-type-specific hits (CellDMC design), one row per model and lineage
cell_rows <- rbindlist(lapply(sort(summary_files), function(file) {
    s <- fromJSON(file, simplifyVector = FALSE)
    lineages <- s[["cells"]][["lineages"]]
    if (is.null(lineages) || length(lineages) == 0) return(NULL)
    rbindlist(lapply(names(lineages), function(cell) as.data.table(c(list(model = s[["model"]], cell = cell), lineages[[cell]]))), fill = TRUE)
}), fill = TRUE)
if (nrow(cell_rows) == 0) {
    cell_rows <- data.table(model = character(), cell = character(), hits_fdr = integer(), hits_higher = integer(), hits_lower = integer(),
                            lambda = numeric(), top_cpg = character(), top_p = numeric(), top_estimate = numeric())
}
fwrite(cell_rows, paste0(prefix, "_cells.tsv"), sep = "\\t", na = "")
if (nrow(cell_rows) > 0) {
    writeLines(c(
        "# id: 'ewas_cells'",
        "# section_name: 'Association models: cell types'",
        "# description: 'Cell-type-specific associations (CellDMC design: the trait times each blood cell lineage, with the other covariates of the model). One row per model and lineage: CpGs at FDR < 0.05 within the lineage, with higher and with lower methylation, genomic lambda, and the top CpG. Small lineages (B, NK) have little power.'",
        "# plot_type: 'table'"
    ), paste0(prefix, "_cells_mqc.tsv"))
    write.table(data.frame(Model = paste0(cell_rows[["model"]], " (", cell_rows[["cell"]], ")"), cell_rows[, !c("model", "cell")], check.names = FALSE),
                paste0(prefix, "_cells_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE, na = "")
}

# Meff and the study-wide threshold over the tested model contrasts
tests <- list()
for (file in stats_files) {
    model <- sub("_stats[.]tsv[.]gz\$", "", basename(file))
    s <- fread(text = readLines(gzfile(file)))
    if (!"fdr" %in% names(s)) s[, fdr := p.adjust(p_value, "BH"), by = contrast]
    s <- s[, c("cpg_id", "contrast", "t", "p_value", "fdr"), with = FALSE]
    if (nrow(s) > 0) tests[[model]] <- s
}
if (length(tests) > 0) {
    stats <- rbindlist(lapply(names(tests), function(prefix_model) {
        s <- tests[[prefix_model]]
        name <- rows[["model"]][paste0("ewas_", rows[["model"]]) == prefix_model | rows[["model"]] == prefix_model][1]
        s[, model := if (is.na(name)) prefix_model else name]
    }))
    stats[, key := paste(model, contrast, sep = " | ")]
    keys <- unique(stats[["key"]])
    meff <- 1
    if (length(keys) > 1) {
        t_wide <- dcast(stats, cpg_id ~ key, value.var = "t")
        correlation <- cor(as.matrix(t_wide[, keys, with = FALSE]), use = "pairwise.complete.obs")
        correlation[!is.finite(correlation)] <- 0
        diag(correlation) <- 1
        eigenvalues <- pmax(eigen(correlation, symmetric = TRUE, only.values = TRUE)[["values"]], 0)
        meff <- sum(sqrt(eigenvalues))^2 / sum(eigenvalues)
    }
    p_study <- p_ewide / meff
    hits <- stats[, .(hits_study_wide = sum(p_value < p_study, na.rm = TRUE)), by = .(model, contrast)]
    if (!"contrast" %in% names(rows)) rows[, contrast := NA_character_]
    rows <- merge(rows, hits, by = c("model", "contrast"), all.x = TRUE, sort = FALSE)
    rows[!is.na(hits_study_wide), `:=`(meff = round(meff, 3), p_study_wide = signif(p_study, 4))]
    setorder(rows, model)
}
fwrite(rows, paste0(prefix, ".tsv"), sep = "\\t", na = "")

# Hits with a cis-meQTL
meqtl <- data.table(model = character(), contrast = character(), tested = integer(), tested_meqtl = integer(), hits_fdr = integer(),
                    hits_meqtl = integer(), share_tested = numeric(), share_hits = numeric(), odds_ratio = numeric(), p_value = numeric())
if (nzchar(meqtl_file) && length(tests) > 0) {
    meqtl_cpgs <- fread(meqtl_file, colClasses = list(character = "cpg_id"))[["cpg_id"]]
    meqtl <- stats[, {
        has <- cpg_id %in% meqtl_cpgs
        hit <- !is.na(fdr) & fdr < 0.05
        test <- if (sum(hit) > 0 && sum(has) > 0) fisher.test(table(factor(hit, c(TRUE, FALSE)), factor(has, c(TRUE, FALSE)))) else NULL
        .(tested = .N, tested_meqtl = sum(has), hits_fdr = sum(hit), hits_meqtl = sum(hit & has),
          share_tested = round(mean(has), 4), share_hits = if (sum(hit) > 0) round(mean(has[hit]), 4) else NA_real_,
          odds_ratio = if (is.null(test)) NA_real_ else signif(unname(test[["estimate"]]), 4),
          p_value = if (is.null(test)) NA_real_ else signif(test[["p.value"]], 4))
    }, by = .(model, contrast)]
}
fwrite(meqtl, paste0(prefix, "_meqtl.tsv"), sep = "\\t", na = "")

# MultiQC: the summary, with the hits that have a cis-meQTL when a meQTL list was given
if (nrow(meqtl) > 0) rows <- merge(rows, meqtl[, .(model, contrast, hits_meqtl, share_hits_meqtl = share_hits)], by = c("model", "contrast"), all.x = TRUE, sort = FALSE)
columns <- intersect(c("model", "contrast", "n", "lambda", "bacon_inflation", "bacon_bias", "hits_fdr", "hits_epigenome_wide",
                       "hits_suggestive", "hits_fdr_bacon", "meff", "p_study_wide", "hits_study_wide", "hits_meqtl", "share_hits_meqtl",
                       "top_cpg", "top_p", "status", "design_flags"), names(rows))
writeLines(c(
    "# id: 'ewas_summary'",
    "# section_name: 'Association models'",
    "# description: 'One row per model and contrast: complete cases, genomic lambda, bacon inflation and bias, CpGs at FDR < 0.05, P < 9e-8 and P < 1e-5 (and at FDR < 0.05 after bacon), the effective number of tests across the models (Meff) with the study-wide threshold and its hits, the hits with a cis-meQTL (with --meqtl), the top CpG, and the design checks that were flagged (trait confounded with plate or Sentrix row at P < 0.01, or a variance inflation factor above 10).'",
    "# plot_type: 'table'"
), paste0(prefix, "_mqc.tsv"))
table <- rows[, ..columns]
table[, Model := if ("contrast" %in% names(table)) ifelse(is.na(contrast), model, paste0(model, " (", contrast, ")")) else model]
setcolorder(table, "Model")
table[, c("model") := NULL]
if ("contrast" %in% names(table)) table[, contrast := NULL]
write.table(table, paste0(prefix, "_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE, na = "")
