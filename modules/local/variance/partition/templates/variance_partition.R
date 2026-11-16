#!/usr/bin/env Rscript
#
# Variance partitioning (variancePartition; Hoffman and Schadt 2016) of the
# --variance_partition_cpgs most variable CpGs (M-values, complete in every
# array) of the filtered data: for each CpG, the share of its variance explained
# by each technical and biological variable, from one linear mixed model with
# the categorical variables as random effects (DNA plate, BeadChip, Sentrix row,
# batch, sex) and the continuous ones as fixed effects (age, five cell lineages
# with granulocytes as the reference), as the package recommends. A variable is
# used when it has at least two values; a categorical one also needs fewer
# levels than half the arrays. Arrays missing a used variable are left out, and
# a variable missing in more than 10% of the arrays is not used.
#

prefix <- "${prefix}"
se_dir <- "${se}"
samples_file <- "${samples}"
sample_qc_file <- "${sample_qc}"
cells_file <- "${cell_fractions}"

defaults <- list(n_cpgs = "20000", cores = "${task.cpus}")
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
n_cpgs <- as.integer(opts[["n_cpgs"]])
cores <- max(1L, as.integer(opts[["cores"]]))

Sys.setenv(HOME = getwd())
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(HDF5Array)
    library(variancePartition)
    library(BiocParallel)
    library(data.table)
})
writeLines(c('"${task.process}":',
    paste0("    variancePartition: ", packageVersion("variancePartition")),
    paste0("    lme4: ", packageVersion("lme4")),
    paste0("    HDF5Array: ", packageVersion("HDF5Array"))), "versions.yml")
set.seed(1)

read_table <- function(file) {
    if (!nzchar(file) || !file.exists(file) || file.info(file)[["size"]] == 0) return(NULL)
    fread(file, na.strings = c("", "NA"))
}

se <- loadHDF5SummarizedExperiment(se_dir)
arrays <- colnames(se)
mvalue <- as.matrix(assay(se, "mvalue"))
mvalue <- mvalue[rowSums(is.na(mvalue)) == 0, , drop = FALSE]
variances <- matrixStats::rowVars(mvalue)
top <- order(-variances)[seq_len(min(n_cpgs, nrow(mvalue)))]
mvalue <- mvalue[top, , drop = FALSE]

# Variables of every array
samples <- read_table(samples_file)
samples <- samples[match(arrays, sample)]
info <- data.frame(row.names = arrays)
categorical <- character()
continuous <- character()
add <- function(name, values, type) {
    values[values %in% ""] <- NA
    if (mean(is.na(values)) > 0.1) return(invisible(NULL))
    ok <- values[!is.na(values)]
    if (type == "categorical") {
        levels <- length(unique(ok))
        if (levels < 2 || levels >= length(arrays) / 2) return(invisible(NULL))
        info[[name]] <<- factor(values)
        categorical <<- c(categorical, name)
    } else {
        values <- suppressWarnings(as.numeric(values))
        if (sum(!is.na(values)) < 3 || var(values, na.rm = TRUE) == 0) return(invisible(NULL))
        info[[name]] <<- values
        continuous <<- c(continuous, name)
    }
}
if ("plate" %in% names(samples)) add("plate", as.character(samples[["plate"]]), "categorical")
if ("sentrix_id" %in% names(samples)) add("beadchip", as.character(samples[["sentrix_id"]]), "categorical")
if ("sentrix_position" %in% names(samples)) add("sentrix_row", substr(samples[["sentrix_position"]], 1, 3), "categorical")
if ("batch" %in% names(samples)) add("batch", as.character(samples[["batch"]]), "categorical")
sex <- if ("reported_sex" %in% names(samples)) as.character(samples[["reported_sex"]]) else rep(NA_character_, length(arrays))
sample_qc <- read_table(sample_qc_file)
predicted_column <- if (is.null(sample_qc)) NA else grep("^predicted_sex_", names(sample_qc), value = TRUE)[1]
if (!is.na(predicted_column)) {
    predicted <- c(male = "M", female = "F")[sample_qc[[predicted_column]][match(arrays, sample_qc[["sample"]])]]
    sex <- ifelse(sex %in% c("M", "F"), sex, predicted)
}
add("sex", sex, "categorical")
if ("reported_age" %in% names(samples)) add("age", samples[["reported_age"]], "continuous")
cells <- read_table(cells_file)
if (!is.null(cells)) {
    for (lineage in setdiff(grep("^lineage_", names(cells), value = TRUE), "lineage_Gran")) {
        add(sub("^lineage_", "cells_", lineage), cells[[lineage]][match(arrays, cells[["sample"]])], "continuous")
    }
}
keep <- complete.cases(info)
info <- info[keep, , drop = FALSE]
mvalue <- mvalue[, keep, drop = FALSE]
used <- c(continuous, categorical)
if (length(used) == 0 || nrow(info) < 10 || nrow(info) <= length(used) + 5) {
    # Too few arrays for the mixed models: empty tables with a note
    fwrite(data.table(cpg_id = character(), residual = numeric()), paste0(prefix, ".tsv.gz"), sep = "\\t")
    fwrite(data.table(variable = NA_character_, type = NA_character_, median = NA_real_, mean = NA_real_, q90 = NA_real_,
                      cpgs_over_10 = NA_integer_, arrays = nrow(info), cpgs = 0L, formula = NA_character_,
                      note = sprintf("not run: %d arrays with complete variables (%s); at least 10, and 6 more than the variables, are needed",
                                     nrow(info), if (length(used) == 0) "none" else paste(used, collapse = ", "))),
           paste0(prefix, "_summary.tsv"), sep = "\\t", na = "")
    quit(save = "no")
}
terms <- c(continuous, sprintf("(1|%s)", categorical))
form <- as.formula(paste("~", paste(terms, collapse = " + ")))
message("Formula: ", paste(deparse(form), collapse = ""), "; ", nrow(info), " arrays, ", nrow(mvalue), " CpGs")

param <- if (cores > 1) MulticoreParam(workers = cores, RNGseed = 1) else SerialParam()
varpart <- fitExtractVarPartModel(mvalue, form, info, BPPARAM = param, showWarnings = FALSE)
varpart <- as.data.frame(varpart)
table <- data.table(cpg_id = rownames(mvalue), round(100 * varpart, 3))
setnames(table, "Residuals", "residual", skip_absent = TRUE)
fwrite(table, paste0(prefix, ".tsv.gz"), sep = "\\t")

variables <- setdiff(names(table), "cpg_id")
summary <- data.table(variable = variables,
                      type = ifelse(variables %in% categorical, "random", ifelse(variables %in% continuous, "fixed", "residual")),
                      median = vapply(variables, function(v) median(table[[v]]), numeric(1)),
                      mean = vapply(variables, function(v) mean(table[[v]]), numeric(1)),
                      q90 = vapply(variables, function(v) unname(quantile(table[[v]], 0.9)), numeric(1)),
                      cpgs_over_10 = vapply(variables, function(v) sum(table[[v]] > 10), integer(1)))
summary[, c("median", "mean", "q90") := lapply(.SD, round, 3), .SDcols = c("median", "mean", "q90")]
summary[, `:=`(arrays = nrow(info), cpgs = nrow(table), formula = paste(deparse(form), collapse = ""))]
summary[, note := NA_character_]
fwrite(summary[order(-mean)], paste0(prefix, "_summary.tsv"), sep = "\\t", na = "")

# MultiQC: mean share of variance per variable
writeLines(c(
    "# id: 'variance_partition'",
    "# section_name: 'Variance partition'",
    paste0("# description: 'Mean percentage of the variance of the ", nrow(table), " most variable CpGs (M-values) explained by each variable, from one linear mixed model per CpG (variancePartition) on ",
           nrow(info), " arrays. Large technical shares (plate, BeadChip, row) call for these variables as model covariates.'"),
    "# plot_type: 'bargraph'",
    "# pconfig:",
    "#     id: 'variance_partition_plot'",
    "#     title: 'Variance partition: mean % of variance per variable'",
    "#     ylab: '% of variance'"
), paste0(prefix, "_mqc.tsv"))
fwrite(data.table(variable = "mean % of variance", t(setNames(summary[["mean"]], summary[["variable"]]))),
       paste0(prefix, "_mqc.tsv"), sep = "\\t", append = TRUE, col.names = TRUE)
