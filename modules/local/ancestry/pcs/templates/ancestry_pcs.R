#!/usr/bin/env Rscript
#
# Genetic ancestry principal components from the array itself, after the
# EpiAnceR+ approach (Höffler et al. 2025), on the arrays that passed sample QC:
#   - CpG probes with a common SNP (minor allele frequency 5% or more) within
#     5 bp of the 3' end, at the target CpG or the extension base (Zhou lab
#     mask, the snp_common_5pt column of the probe annotation), on the autosomes;
#     their beta values carry the genotype of the SNP
#   - probes called in fewer than 90% of the arrays are left out; other missing
#     values are set to the probe mean
#   - the beta values are residualised for the control-probe PCs, the PCs of the
#     blood cell lineage fractions (at most five), sex (reported, otherwise
#     predicted) and age (when reported for at least 90% of the arrays; missing
#     ages set to the mean), so that technical and biological variation does not
#     drive the components
#   - the genotypes of the rs probes (k-means on the beta values with centres
#     0.2, 0.5 and 0.8, coded 0, 0.5 and 1) are added unchanged
#   - PCA of the combined matrix (centred, not scaled)
# The first --n_ancestry_pcs components are the ancestry covariates (@ancestry).
# Each component is tested against the technical and biological variables
# (linear regression, R-squared) to show what it captures; when the phenotype
# table has an ancestry, ethnicity or race column, against that column too.
# The code is written from the published description; the EpiAnceR scripts are
# not used.
#

prefix <- "${prefix}"
se_dir <- "${se}"
annotation_file <- "${annotation}"
samples_file <- "${samples}"
sample_qc_file <- "${sample_qc}"
ctrl_file <- "${ctrl_pcs}"
cells_file <- "${cell_fractions}"
phenotypes_file <- "${phenotypes}"

defaults <- list(n_ancestry_pcs = "4", max_pcs = "10", min_call_rate = "0.9")
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
n_keep <- as.integer(opts[["n_ancestry_pcs"]])
max_pcs <- as.integer(opts[["max_pcs"]])
min_call_rate <- as.numeric(opts[["min_call_rate"]])

Sys.setenv(HOME = getwd())
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(HDF5Array)
    library(data.table)
    library(jsonlite)
})
writeLines(c('"${task.process}":',
    paste0("    R: ", getRversion()),
    paste0("    SummarizedExperiment: ", packageVersion("SummarizedExperiment")),
    paste0("    HDF5Array: ", packageVersion("HDF5Array"))), "versions.yml")
set.seed(1)

read_table <- function(file) {
    if (!nzchar(file) || !file.exists(file) || file.info(file)[["size"]] == 0) return(NULL)
    if (endsWith(file, ".gz")) fread(text = readLines(gzfile(file)), na.strings = c("", "NA"))
    else fread(file, na.strings = c("", "NA"), sep = if (endsWith(file, ".csv")) "," else "\\t")
}

# Arrays that passed sample QC
sample_qc <- read_table(sample_qc_file)
se <- loadHDF5SummarizedExperiment(se_dir)
arrays <- sort(intersect(colnames(se), sample_qc[status %in% c("PASS", "FLAG"), sample]))
n <- length(arrays)
if (n < 10) {
    # Too few arrays for population structure: empty tables, so models with
    # @ancestry stop with a clear message
    fwrite(data.table(sample = arrays), paste0(prefix, ".tsv"), sep = "\t")
    fwrite(data.table(pc = character(), variance_explained = numeric(), cumulative = numeric(), covariate = logical()), paste0(prefix, "_variance.tsv"), sep = "\t")
    fwrite(data.table(metric = c("arrays", "note"), value = c(n, "fewer than 10 arrays passed QC; no ancestry PCs computed")), paste0(prefix, "_summary.tsv"), sep = "\t")
    fwrite(data.table(variable = character(), pc = character(), n = integer(), r_squared = numeric(), p_value = numeric()), paste0(prefix, "_associations.tsv"), sep = "\t")
    quit(save = "no")
}

# SNP CpGs and rs probes
anno <- read_table(annotation_file)
snp_probes <- anno[probe_type == "cg" & snp_common_5pt == TRUE & !is.na(chr) & !chr %in% c("chrX", "chrY", "chrM"), probe_id]
snp_probes <- intersect(rownames(se), snp_probes)
rs_probes <- grep("^rs", rownames(se), value = TRUE)
beta <- as.matrix(assay(se[snp_probes, arrays], "beta"))
call_rate <- rowMeans(!is.na(beta))
beta <- beta[call_rate >= min_call_rate & apply(beta, 1, var, na.rm = TRUE) > 0, , drop = FALSE]
if (nrow(beta) == 0) stop("No CpG probe with a common SNP passed the call-rate filter")
missing <- which(is.na(beta), arr.ind = TRUE)
if (nrow(missing) > 0) beta[missing] <- rowMeans(beta, na.rm = TRUE)[missing[, 1]]

rs_beta <- as.matrix(assay(se[rs_probes, arrays], "beta_unmasked"))
call_genotypes <- function(x) {
    ok <- is.finite(x)
    g <- rep(NA_real_, length(x))
    fit <- tryCatch(kmeans(x[ok], centers = c(0.2, 0.5, 0.8)), error = function(e) NULL)
    g[ok] <- if (is.null(fit)) ifelse(x[ok] < 0.2, 0, ifelse(x[ok] > 0.8, 2, 1)) else fit[["cluster"]] - 1
    g / 2
}
genotypes <- t(apply(rs_beta, 1, call_genotypes))
genotypes <- genotypes[rowMeans(is.na(genotypes)) <= 1 - min_call_rate, , drop = FALSE]
genotypes <- genotypes[apply(genotypes, 1, var, na.rm = TRUE) > 0, , drop = FALSE]
missing <- which(is.na(genotypes), arr.ind = TRUE)
if (nrow(missing) > 0) genotypes[missing] <- rowMeans(genotypes, na.rm = TRUE)[missing[, 1]]

# Covariates for the residualisation
samples <- read_table(samples_file)
samples <- samples[match(arrays, sample)]
covariates <- data.table(sample = arrays)
used <- character()
ctrl <- read_table(ctrl_file)
ctrl_columns <- if (is.null(ctrl)) character() else grep("^ctrl_PC", names(ctrl), value = TRUE)
if (length(ctrl_columns) > 0) {
    covariates <- cbind(covariates, ctrl[match(arrays, sample), ..ctrl_columns])
    used <- c(used, sprintf("%d control-probe PCs", length(ctrl_columns)))
}
cells <- read_table(cells_file)
lineages <- if (is.null(cells)) character() else grep("^lineage_", names(cells), value = TRUE)
if (length(lineages) >= 2) {
    fractions <- as.matrix(cells[match(arrays, sample), ..lineages])
    fractions <- fractions[, apply(fractions, 2, var, na.rm = TRUE) > 0, drop = FALSE]
    if (ncol(fractions) >= 2 && all(complete.cases(fractions))) {
        k <- min(5, ncol(fractions) - 1)
        cell_pcs <- prcomp(fractions)[["x"]][, seq_len(k), drop = FALSE]
        colnames(cell_pcs) <- paste0("cell_PC", seq_len(k))
        covariates <- cbind(covariates, cell_pcs)
        used <- c(used, sprintf("%d PCs of the cell lineage fractions", k))
    }
}
sex <- samples[["reported_sex"]]
predicted_column <- grep("^predicted_sex_", names(sample_qc), value = TRUE)[1]
if (!is.na(predicted_column)) {
    predicted <- c(male = "M", female = "F")[sample_qc[[predicted_column]][match(arrays, sample_qc[["sample"]])]]
    sex <- ifelse(sex %in% c("M", "F"), sex, predicted)
}
if (all(sex %in% c("M", "F")) && length(unique(sex)) == 2) {
    covariates[["men"]] <- as.numeric(sex == "M")
    used <- c(used, "sex")
}
age <- suppressWarnings(as.numeric(samples[["reported_age"]]))
if (mean(is.finite(age)) >= 0.9 && var(age, na.rm = TRUE) > 0) {
    age[!is.finite(age)] <- mean(age, na.rm = TRUE)
    covariates[["age"]] <- age
    used <- c(used, "age")
}
design <- cbind(1, as.matrix(covariates[, -1, with = FALSE]))
design <- design[, colSums(!is.finite(design)) == 0, drop = FALSE]
if (ncol(design) >= n - 2) {
    design <- design[, 1, drop = FALSE]
    used <- "none (more covariates than arrays allow)"
}
residuals <- qr.resid(qr(design), t(beta))

# PCA of the residualised SNP CpGs and the rs genotypes
combined <- cbind(residuals, t(genotypes))
pca <- prcomp(combined, center = TRUE, scale. = FALSE)
k_out <- min(max_pcs, n - 1, ncol(combined))
share <- pca[["sdev"]]^2 / sum(pca[["sdev"]]^2)
n_keep <- min(n_keep, k_out)
scores <- pca[["x"]][, seq_len(k_out), drop = FALSE]
colnames(scores) <- paste0("ancestry_pc", seq_len(k_out))
fwrite(data.table(sample = arrays, round(scores[, seq_len(n_keep), drop = FALSE], 6)), paste0(prefix, ".tsv"), sep = "\\t")
fwrite(data.table(pc = paste0("ancestry_pc", seq_len(k_out)), variance_explained = round(share[seq_len(k_out)], 5),
                  cumulative = round(cumsum(share)[seq_len(k_out)], 5), covariate = seq_len(k_out) <= n_keep),
       paste0(prefix, "_variance.tsv"), sep = "\\t")
fwrite(data.table(metric = c("arrays", "snp_cpgs", "rs_snps", "residualised_for"),
                  value = c(n, nrow(beta), nrow(genotypes), paste(if (length(used) == 0) "none" else used, collapse = "; "))),
       paste0(prefix, "_summary.tsv"), sep = "\\t")

# What each component captures: technical and biological variables, and the
# ancestry column of the phenotype table when there is one
variables <- list()
if ("plate" %in% names(samples)) variables[["DNA plate"]] <- samples[["plate"]]
if ("sentrix_id" %in% names(samples)) variables[["BeadChip"]] <- samples[["sentrix_id"]]
if ("sentrix_position" %in% names(samples)) variables[["Sentrix row"]] <- substr(samples[["sentrix_position"]], 1, 3)
if (any(sex %in% c("M", "F"))) variables[["Sex"]] <- sex
if (any(is.finite(suppressWarnings(as.numeric(samples[["reported_age"]]))))) variables[["Age"]] <- suppressWarnings(as.numeric(samples[["reported_age"]]))
for (column in lineages) variables[[sub("^lineage_", "", column)]] <- cells[[column]][match(arrays, cells[["sample"]])]
phenotypes <- read_table(phenotypes_file)
ancestry_column <- if (is.null(phenotypes)) NA else grep("^(ancestry|ethnicity|race)\$", names(phenotypes), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(ancestry_column)) variables[[ancestry_column]] <- as.character(phenotypes[[ancestry_column]][match(arrays, phenotypes[["sample"]])])
associations <- rbindlist(lapply(names(variables), function(v) rbindlist(lapply(seq_len(k_out), function(i) {
    x <- variables[[v]]
    ok <- !is.na(x)
    if (is.character(x)) x <- factor(x)
    levels_ok <- if (is.factor(x)) nlevels(droplevels(x[ok])) >= 2 && nlevels(droplevels(x[ok])) < sum(ok) - 1 else var(x[ok]) > 0
    if (sum(ok) < 5 || !isTRUE(levels_ok)) return(NULL)
    fit <- summary(lm(scores[ok, i] ~ x[ok]))
    f <- fit[["fstatistic"]]
    data.table(variable = v, pc = colnames(scores)[i], n = sum(ok), r_squared = round(fit[["r.squared"]], 4),
               p_value = signif(unname(pf(f[1], f[2], f[3], lower.tail = FALSE)), 4))
}))))
if (nrow(associations) == 0) associations <- data.table(variable = character(), pc = character(), n = integer(), r_squared = numeric(), p_value = numeric())
fwrite(associations, paste0(prefix, "_associations.tsv"), sep = "\\t")

# MultiQC: PC1 against PC2 (coloured by the ancestry column when there is one,
# otherwise by DNA plate), and the variance of each component
group <- if (!is.na(ancestry_column)) variables[[ancestry_column]] else if ("DNA plate" %in% names(variables)) as.character(variables[["DNA plate"]]) else rep("all", n)
group[is.na(group)] <- "unknown"
palette <- c("#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3", "#937860", "#da8bc3", "#8c8c8c", "#ccb974", "#64b5cd")
colours <- setNames(rep(palette, length.out = length(unique(group))), sort(unique(group)))
if (k_out >= 2) {
    points <- lapply(seq_len(n), function(i) list(x = round(scores[i, 1], 4), y = round(scores[i, 2], 4), color = colours[[group[i]]]))
    names(points) <- arrays
    write_json(list(id = "ancestry_scatter", section_name = "Ancestry PCs",
                    description = paste0("Genetic ancestry components from the CpGs with a common SNP and the rs probes (EpiAnceR+ approach), residualised for ",
                                         paste(used, collapse = ", "), ". Coloured by ", if (!is.na(ancestry_column)) ancestry_column else "DNA plate",
                                         ": ", paste(sprintf("%s %s", names(colours), colours), collapse = ", "),
                                         ". The first ", n_keep, " components are the @ancestry covariates."),
                    plot_type = "scatter",
                    pconfig = list(id = "ancestry_scatter_plot", title = "Ancestry PCs: PC1 against PC2", xlab = "Ancestry PC1", ylab = "Ancestry PC2"),
                    data = points),
               paste0(prefix, "_scatter_mqc.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
}
writeLines(c("# id: 'ancestry_scree'", "# section_name: 'Ancestry PCs: variance'",
             "# description: 'Share of the variance of the residualised SNP CpGs and rs genotypes explained by each ancestry component.'",
             "# plot_type: 'bargraph'", "# pconfig:", "#     id: 'ancestry_scree_plot'", "#     title: 'Ancestry PCs: variance explained (%)'"),
           paste0(prefix, "_scree_mqc.tsv"))
fwrite(data.table(pc = paste0("PC", seq_len(k_out)), variance = round(100 * share[seq_len(k_out)], 2)),
       paste0(prefix, "_scree_mqc.tsv"), sep = "\\t", append = TRUE, col.names = TRUE)
