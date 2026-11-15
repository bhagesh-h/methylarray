#!/usr/bin/env Rscript
#
# Exploratory analysis of the arrays that passed QC:
#   - control-probe PCs (Lehne et al. 2015): PCA of the log2, standardised
#     intensities of the non-negative control probes in both channels; the
#     number kept is --n_ctrl_pcs, or with 'auto' the smallest number that
#     explains at least 90% of the variance (at most 20)
#   - PCA of the --n_variable_cpgs most variable CpGs (M-values; missing values
#     set to the CpG mean for the PCA only); --n_pcs, or with 'auto' the
#     random matrix theory estimate of isva::EstDimRMT (at most 20)
#   - association of each PC with technical and biological variables: linear
#     regression for numeric variables, Kruskal-Wallis for categorical ones,
#     with R-squared from the linear model
#   - outlier flags (report only, no array is removed): robust Mahalanobis
#     distance on PC1 and PC2, each centred on its median and scaled by its MAD
#     (the PCs are uncorrelated), flagged above the chi-square (2 df) 0.999
#     quantile
# Variables with a single value, and categorical variables with as many levels
# as arrays (for example the Sentrix ID in a small run), are left out.
#

prefix <- "${prefix}"
se_dir <- "${se}"
samples_file <- "${samples}"
sample_qc_file <- "${sample_qc}"
controls_file <- "${controls}"
cells_file <- "${cell_fractions}"
phenotypes_file <- "${phenotypes}"

defaults <- list(n_variable_cpgs = "50000", n_pcs = "auto", n_ctrl_pcs = "auto", max_pcs = "20")
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
n_variable <- as.integer(opts[["n_variable_cpgs"]])
max_pcs <- as.integer(opts[["max_pcs"]])

Sys.setenv(HOME = getwd())
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(HDF5Array)
    library(data.table)
    library(isva)
})
set.seed(1)

read_table <- function(file) {
    if (endsWith(file, ".gz")) fread(text = readLines(gzfile(file)), na.strings = c("", "NA"))
    else fread(file, na.strings = c("", "NA"))
}

se <- loadHDF5SummarizedExperiment(se_dir)
arrays <- colnames(se)
n <- length(arrays)
if (n < 3) stop("PCA needs at least three arrays, found ", n)
pc_limit <- min(max_pcs, n - 1)

# Control-probe PCs
ctrl_scores <- data.table(sample = arrays)
ctrl_variance <- data.table(pc = character(), variance_explained = numeric(), cumulative = numeric())
if (nzchar(controls_file)) {
    controls <- read_table(controls_file)
    ctrl <- as.matrix(controls[group != "NEGATIVE", ..arrays])
    ctrl <- ctrl[complete.cases(ctrl), , drop = FALSE]
    ctrl <- log2(pmax(ctrl, 1))
    ctrl <- ctrl[apply(ctrl, 1, var) > 0, , drop = FALSE]
    pca_ctrl <- prcomp(t(ctrl), center = TRUE, scale. = TRUE)
    share <- pca_ctrl[["sdev"]]^2 / sum(pca_ctrl[["sdev"]]^2)
    k <- if (opts[["n_ctrl_pcs"]] == "auto") which(cumsum(share) >= 0.9)[1] else as.integer(opts[["n_ctrl_pcs"]])
    k <- max(1L, min(k, pc_limit))
    ctrl_scores <- data.table(sample = arrays, round(pca_ctrl[["x"]][arrays, seq_len(k), drop = FALSE], 6))
    setnames(ctrl_scores, 2:ncol(ctrl_scores), paste0("ctrl_PC", seq_len(k)))
    shown <- seq_len(min(length(share), max_pcs))
    ctrl_variance <- data.table(pc = paste0("ctrl_PC", shown), variance_explained = round(share[shown], 5),
                                cumulative = round(cumsum(share)[shown], 5))
}
fwrite(ctrl_scores, paste0(prefix, "_ctrl_pcs.tsv"), sep = "\\t")
fwrite(ctrl_variance, paste0(prefix, "_ctrl_pca_variance.tsv"), sep = "\\t")

# PCA of the most variable CpGs
mvalue <- as.matrix(assay(se, "mvalue"))
variances <- apply(mvalue, 1, var, na.rm = TRUE)
top <- order(variances, decreasing = TRUE)[seq_len(min(n_variable, sum(is.finite(variances))))]
x <- mvalue[top, , drop = FALSE]
rm(mvalue)
missing <- which(is.na(x), arr.ind = TRUE)
if (nrow(missing) > 0) x[missing] <- rowMeans(x, na.rm = TRUE)[missing[, 1]]
centred <- x - rowMeans(x)
n_pcs <- if (opts[["n_pcs"]] == "auto") EstDimRMT(centred, plot = FALSE)[["dim"]] else as.integer(opts[["n_pcs"]])
n_pcs <- max(1L, min(n_pcs, pc_limit))
pca <- prcomp(t(centred), center = FALSE, scale. = FALSE)
share <- pca[["sdev"]]^2 / sum(pca[["sdev"]]^2)
scores <- pca[["x"]][arrays, seq_len(n_pcs), drop = FALSE]
fwrite(data.table(sample = arrays, round(scores, 6)), paste0(prefix, "_scores.tsv"), sep = "\\t")
shown <- seq_len(min(length(share), max_pcs))
variance <- data.table(pc = paste0("PC", shown), variance_explained = round(share[shown], 5),
                       cumulative = round(cumsum(share)[shown], 5), tested = shown <= n_pcs)
fwrite(variance, paste0(prefix, "_variance.tsv"), sep = "\\t")

# Outlier flags on PC1 and PC2 (report only)
outliers <- data.table(sample = arrays, robust_distance = NA_real_, pca_outlier = FALSE)
if (ncol(pca[["x"]]) >= 2 && length(arrays) >= 5) {
    pc12 <- pca[["x"]][arrays, 1:2, drop = FALSE]
    scaled <- sweep(sweep(pc12, 2, apply(pc12, 2, median)), 2, pmax(apply(pc12, 2, mad), .Machine[["double.eps"]]), "/")
    distance <- rowSums(scaled^2)
    outliers[, robust_distance := round(distance, 4)]
    outliers[, pca_outlier := distance > qchisq(0.999, df = 2)]
}
fwrite(outliers, paste0(prefix, "_outliers.tsv"), sep = "\\t")

# Variables to test, one row per array
samples <- read_table(samples_file)[match(arrays, sample)]
sample_qc <- read_table(sample_qc_file)[match(arrays, sample)]
variables <- list()
add <- function(name, values, group) variables[[name]] <<- list(values = values, group = group)
for (column in c("plate", "batch", "sentrix_id", "sentrix_row")) {
    if (column %in% names(samples)) add(column, as.character(samples[[column]]), "technical")
}
if ("frac_detected" %in% names(sample_qc)) add("frac_detected", as.numeric(sample_qc[["frac_detected"]]), "technical")
for (column in setdiff(names(ctrl_scores), "sample")) add(column, ctrl_scores[[column]], "technical")
sex_column <- grep("^predicted_sex_", names(sample_qc), value = TRUE)[1]
if (!is.na(sex_column)) add("predicted_sex", as.character(sample_qc[[sex_column]]), "biological")
if ("reported_age" %in% names(samples)) add("age", suppressWarnings(as.numeric(samples[["reported_age"]])), "biological")
if (nzchar(cells_file)) {
    cells <- read_table(cells_file)[match(arrays, sample)]
    for (column in grep("^lineage_", names(cells), value = TRUE)) add(sub("^lineage_", "", column), cells[[column]], "biological")
}
if (nzchar(phenotypes_file)) {
    phenotypes <- read_table(phenotypes_file)
    phenotypes <- phenotypes[match(arrays, phenotypes[["sample"]])]
    for (column in setdiff(names(phenotypes), "sample")) {
        values <- phenotypes[[column]]
        add(column, if (is.numeric(values)) values else as.character(values), "phenotype")
    }
}

testable <- function(values) {
    present <- !is.na(values)
    if (sum(present) < 3) return(FALSE)
    if (is.numeric(values)) return(length(unique(values[present])) >= 3)
    levels <- length(unique(values[present]))
    levels >= 2 && levels < sum(present)
}
associations <- rbindlist(lapply(names(variables), function(name) {
    values <- variables[[name]][["values"]]
    if (!testable(values)) return(NULL)
    rbindlist(lapply(seq_len(n_pcs), function(i) {
        y <- scores[, i]
        keep <- !is.na(values)
        if (is.numeric(values)) {
            fit <- summary(lm(y[keep] ~ values[keep]))
            p <- pf(fit[["fstatistic"]][1], fit[["fstatistic"]][2], fit[["fstatistic"]][3], lower.tail = FALSE)
            test <- "linear regression"
        } else {
            group <- factor(values[keep])
            fit <- summary(lm(y[keep] ~ group))
            p <- kruskal.test(y[keep], group)[["p.value"]]
            test <- "Kruskal-Wallis"
        }
        data.table(variable = name, group = variables[[name]][["group"]], pc = paste0("PC", i), test = test,
                   n = sum(keep), r_squared = round(fit[["r.squared"]], 4), p_value = signif(unname(p), 4))
    }))
}))
fwrite(associations, paste0(prefix, "_associations.tsv"), sep = "\\t")

# MultiQC: association heatmap (-log10 P), scree plot, PC1 and PC2
pcs <- paste0("PC", seq_len(n_pcs))
heatmap <- list()
for (name in unique(associations[["variable"]])) {
    rows <- associations[variable == name]
    heatmap[[name]] <- as.list(setNames(round(-log10(rows[["p_value"]]), 2), rows[["pc"]]))[pcs]
}
if (length(heatmap) > 0) {
    jsonlite::write_json(list(
        id = "pca_associations",
        section_name = "PCA associations",
        description = paste0("Association of the top ", n_pcs, " principal components of the ", nrow(x),
                             " most variable CpGs with technical and biological variables: -log10 P of a linear regression (numeric) or Kruskal-Wallis test (categorical). R-squared values are in pca_associations.tsv."),
        plot_type = "heatmap",
        pconfig = list(id = "pca_associations_plot", title = "PCA: association with variables (-log10 P)",
                       xlab = "Principal component", ylab = "Variable", min = 0),
        xcats = pcs,
        ycats = names(heatmap),
        data = heatmap
    ), paste0(prefix, "_associations_mqc.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
}
writeLines(c(
    "# id: 'pca_scree'",
    "# section_name: 'PCA variance'",
    paste0("# description: 'Share of the variance of the ", nrow(x), " most variable CpGs explained by each principal component. The first ", n_pcs, " are tested for associations.'"),
    "# plot_type: 'bargraph'",
    "# pconfig:",
    "#     id: 'pca_scree_plot'",
    "#     title: 'PCA: variance explained'",
    "#     ylab: 'Share of variance'",
    "#     cpswitch: false"
), paste0(prefix, "_scree_mqc.tsv"))
write.table(data.frame(PC = variance[["pc"]], variance = variance[["variance_explained"]]),
            paste0(prefix, "_scree_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE)
if (ncol(pca[["x"]]) >= 2) {
    sex <- if (is.na(sex_column)) rep(NA_character_, n) else as.character(sample_qc[[sex_column]])
    colours <- c(female = "#d95f02", male = "#1b9e77")
    points <- list()
    for (i in seq_len(n)) {
        colour <- if (is.na(sex[i]) || !sex[i] %in% names(colours)) "#999999" else colours[[sex[i]]]
        points[[arrays[i]]] <- list(x = round(pca[["x"]][arrays[i], 1], 3), y = round(pca[["x"]][arrays[i], 2], 3), color = colour)
    }
    jsonlite::write_json(list(
        id = "pca_scatter",
        section_name = "PCA",
        description = paste0("First two principal components of the most variable CpGs (PC1 ", round(100 * share[1], 1),
                             "%, PC2 ", round(100 * share[2], 1), "% of the variance), coloured by predicted sex (orange female, green male). Flagged as outliers (robust Mahalanobis distance, report only): ",
                             if (any(outliers[["pca_outlier"]])) paste(outliers[["sample"]][outliers[["pca_outlier"]]], collapse = ", ") else "none", "."),
        plot_type = "scatter",
        pconfig = list(id = "pca_scatter_plot", title = "PCA: PC1 and PC2", xlab = "PC1", ylab = "PC2"),
        data = points
    ), paste0(prefix, "_scatter_mqc.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
}

# Software versions
writeLines(c('"${task.process}":',
    paste0("    isva: ", packageVersion("isva")),
    paste0("    SummarizedExperiment: ", packageVersion("SummarizedExperiment")),
    paste0("    HDF5Array: ", packageVersion("HDF5Array"))), "versions.yml")
