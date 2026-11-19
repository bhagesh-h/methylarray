#!/usr/bin/env Rscript
#
# DNA methylation scores: weighted sums of the beta values of published CpG sets
# (assets/dnam_scores; McCartney et al. 2018: smoking, alcohol, BMI; Stevenson
# et al. 2021: IL-6), on the arrays that passed sample QC.
# Input: the detection-masked beta values of the preprocessing route, averaged
# per CpG. The analysis probe filter is not applied, because each score was
# trained on its own CpG set. Score CpGs are matched through the EPIC v1 and
# 450K IDs of the probe annotation as well as the CpG ID. Values that failed
# detection are set to the CpG's mean over the other arrays; CpGs missing from
# the data are left out of the sum. Coverage is reported per score.
# Scores are named dnam_<score>; dnam_smoking is the @smoking model covariate.
#

prefix <- "${prefix}"
se_dir <- "${se}"
retained_file <- "${retained}"
annotation_file <- "${annotation}"
weight_files <- strsplit(trimws("${weights}"), "[[:space:]]+")[[1]]
sources_file <- "${sources}"

defaults <- list(scores = "smoking,alcohol,bmi,il6")
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
wanted <- strsplit(opts[["scores"]], ",")[[1]]

Sys.setenv(HOME = getwd())
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(HDF5Array)
    library(data.table)
})
writeLines(c('"${task.process}":',
    paste0("    SummarizedExperiment: ", packageVersion("SummarizedExperiment")),
    paste0("    data.table: ", packageVersion("data.table"))), "versions.yml")

sources <- fread(sources_file)
unknown <- setdiff(wanted, sources[["score"]])
if (length(unknown) > 0) stop("Unknown scores in --dnam_scores: ", paste(unknown, collapse = ", "), " (available: ", paste(sources[["score"]], collapse = ", "), ")")

# Detection-masked betas of the retained arrays, averaged per CpG
se <- loadHDF5SummarizedExperiment(se_dir)
se <- se[, colnames(se) %in% readLines(retained_file)]
if (ncol(se) == 0) stop("No array passed sample QC")
beta <- as.matrix(assay(se, "beta"))
cpg <- sub("_.*", "", rownames(beta))
present <- !is.na(beta)
beta <- rowsum(ifelse(present, beta, 0), cpg, reorder = FALSE) / rowsum(present * 1, cpg, reorder = FALSE)
beta[!is.finite(beta)] <- NA
rm(present, se)
invisible(gc())

# Old IDs of renamed CpGs (fread reads .gz files only with R.utils; base R decompresses instead)
annotation <- fread(text = readLines(gzfile(annotation_file)), na.strings = "",
                    select = c("cpg_id", "epicv1_id", "hm450_id"), colClasses = "character")
aliases <- unique(rbind(
    annotation[!is.na(epicv1_id) & epicv1_id != cpg_id, .(alias = epicv1_id, cpg_id)],
    annotation[!is.na(hm450_id) & hm450_id != cpg_id, .(alias = hm450_id, cpg_id)]
))
aliases <- aliases[cpg_id %in% rownames(beta) & !alias %in% rownames(beta)][!duplicated(alias)]
if (nrow(aliases) > 0) {
    extra <- beta[aliases[["cpg_id"]], , drop = FALSE]
    rownames(extra) <- aliases[["alias"]]
    beta <- rbind(beta, extra)
}

scores <- data.table(sample = colnames(beta))
coverage <- list()
weight_tables <- list()
for (name in wanted) {
    file <- weight_files[basename(weight_files) == sources[score == name][["file"]]]
    if (length(file) != 1) stop("Weight file for score ", name, " not staged")
    weights <- fread(file)
    found <- weights[cpg_id %in% rownames(beta)]
    values <- beta[found[["cpg_id"]], , drop = FALSE]
    filled <- sum(is.na(values))
    if (filled > 0) {
        means <- rowMeans(values, na.rm = TRUE)
        idx <- which(is.na(values), arr.ind = TRUE)
        values[idx] <- means[idx[, 1]]
    }
    usable <- rowSums(is.na(values)) == 0
    # standardised weight: the weight times the SD of the CpG's beta value in these arrays
    sds <- matrixStats::rowSds(values[usable, , drop = FALSE])
    weight_tables[[name]] <- data.table(score = name, cpg = found[["cpg_id"]][usable], weight = signif(found[["weight"]][usable], 6),
                                        sd_beta = signif(sds, 4), standardised = signif(found[["weight"]][usable] * sds, 4))[order(-abs(standardised))]
    scores[, (paste0("dnam_", name)) := round(colSums(values[usable, , drop = FALSE] * found[["weight"]][usable]), 6)]
    coverage[[name]] <- data.table(score = name, cpgs = nrow(weights), present = sum(usable),
                                   present_by_legacy_id = sum(found[["cpg_id"]][usable] %in% aliases[["alias"]]),
                                   share_present = round(sum(usable) / nrow(weights), 4), filled_values = filled,
                                   reference = sources[score == name][["reference"]], doi = sources[score == name][["doi"]])
}
setorder(scores, sample)
coverage <- rbindlist(coverage)
fwrite(rbindlist(weight_tables), paste0(prefix, "_weights.tsv"), sep = "\\t")
fwrite(scores, paste0(prefix, ".tsv"), sep = "\\t")
fwrite(coverage, paste0(prefix, "_coverage.tsv"), sep = "\\t")

# MultiQC: score coverage on the array
writeLines(c(
    "# id: 'dnam_scores_coverage'",
    "# section_name: 'DNA methylation scores'",
    "# description: 'CpGs of each score in the data. A score is the weighted sum of the beta values of its CpGs present on the array; values that failed detection take the CpG mean over the other arrays.'",
    "# plot_type: 'table'"
), paste0(prefix, "_coverage_mqc.tsv"))
write.table(data.frame(Score = coverage[["score"]], CpGs = coverage[["cpgs"]], Present = coverage[["present"]],
                       `Found by old ID` = coverage[["present_by_legacy_id"]], `Share present` = coverage[["share_present"]],
                       `Filled values` = coverage[["filled_values"]], check.names = FALSE),
            paste0(prefix, "_coverage_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE)

# MultiQC: the distribution of each score over the arrays
writeLines(c(
    "# id: 'dnam_score_values'",
    "# section_name: 'DNA methylation scores: values'",
    "# description: 'Each score over the analysed arrays. Scores are weighted sums on their own scale; compare arrays within a score, not across scores.'",
    "# plot_type: 'violin'"
), paste0(prefix, "_values_mqc.tsv"))
write.table(data.frame(Sample = scores[["sample"]], scores[, setdiff(names(scores), "sample"), with = FALSE], check.names = FALSE),
            paste0(prefix, "_values_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE, na = "")
