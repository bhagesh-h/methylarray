#!/usr/bin/env Rscript
#
# One epigenome-wide association model on the filtered, CpG-level data:
#   1. Resolve the model: covariate tokens, subset, trait transform, complete
#      cases (min_n); every dropped array is recorded with its reason
#   2. Design check: trait against plate and Sentrix row, VIF of the covariates
#   3. limma on M-values: lmFit (+ duplicateCorrelation with a block), then
#      eBayes(trend = TRUE, robust = TRUE); delta-beta from the same design on
#      beta values
#   4. Inflation: genomic lambda, bacon bias and inflation with corrected P/FDR
#   5. Thresholds: FDR < 0.05, P < 9e-8 (EPIC) and P < 1e-5
# Token expansion: @cells = cc_CD4T, cc_CD8T, cc_B, cc_NK, cc_Mono (granulocytes
# as reference); @ctrlpcs = the kept control-probe PCs (ctrlpc1..k); @technical =
# plate, sentrix_row; @sex = predicted_sex; @smoking = dnam_smoking (McCartney et al. 2018).
#

prefix <- "${prefix}"
model_file <- "${model}"
se_dir <- "${se}"
samples_file <- "${samples}"
sample_qc_file <- "${sample_qc}"
phenotypes_file <- "${phenotypes}"
cells_file <- "${cell_fractions}"
ctrl_file <- "${ctrl_pcs}"
scores_file <- "${dnam_scores}"

defaults <- list(p_epigenome_wide = "9e-8", p_suggestive = "1e-5", fdr = "0.05", n_top = "1000")
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
p_suggestive <- as.numeric(opts[["p_suggestive"]])
fdr_level <- as.numeric(opts[["fdr"]])
n_top <- as.integer(opts[["n_top"]])

Sys.setenv(HOME = getwd())
suppressPackageStartupMessages({
    library(limma)
    library(bacon)
    library(SummarizedExperiment)
    library(HDF5Array)
    library(data.table)
    library(jsonlite)
})
set.seed(1)

# Software versions
writeLines(c('"${task.process}":',
    paste0("    limma: ", packageVersion("limma")),
    paste0("    bacon: ", packageVersion("bacon")),
    paste0("    SummarizedExperiment: ", packageVersion("SummarizedExperiment"))), "versions.yml")

model <- fromJSON(model_file, simplifyVector = TRUE)
read_table <- function(file) fread(file, na.strings = c("", "NA"))
se <- loadHDF5SummarizedExperiment(se_dir)
arrays <- colnames(se)

# One row per array with every column a model may use
covariates <- data.table(sample = arrays)
join <- function(table) {
    table <- table[match(arrays, table[["sample"]])]
    new <- setdiff(names(table), names(covariates))
    if (length(new) > 0) covariates[, (new) := table[, ..new]]
}
join(read_table(samples_file))
sample_qc <- read_table(sample_qc_file)
sex_column <- grep("^predicted_sex_", names(sample_qc), value = TRUE)[1]
if (!is.na(sex_column)) join(data.table(sample = sample_qc[["sample"]], predicted_sex = sample_qc[[sex_column]]))
if (nzchar(cells_file)) {
    cells <- read_table(cells_file)
    lineages <- grep("^lineage_", names(cells), value = TRUE)
    join(setnames(cells[, c("sample", lineages), with = FALSE], lineages, sub("^lineage_", "cc_", lineages)))
}
if (nzchar(ctrl_file)) {
    ctrl <- read_table(ctrl_file)
    pcs <- grep("^ctrl_PC", names(ctrl), value = TRUE)
    if (length(pcs) > 0) join(setnames(ctrl[, c("sample", pcs), with = FALSE], pcs, sub("^ctrl_PC", "ctrlpc", pcs)))
}
if (nzchar(scores_file)) join(read_table(scores_file))
if (nzchar(phenotypes_file)) join(read_table(phenotypes_file))

# Covariate tokens
expand <- function(covariate) {
    switch(covariate,
        "@cells" = {
            cols <- c("cc_CD4T", "cc_CD8T", "cc_B", "cc_NK", "cc_Mono")
            if (!all(cols %in% names(covariates))) stop("Model '", model[["name"]], "' uses @cells, but no cell fractions are available (--tissue other)")
            cols
        },
        "@ctrlpcs" = {
            cols <- grep("^ctrlpc[0-9]+\$", names(covariates), value = TRUE)
            if (length(cols) == 0) stop("Model '", model[["name"]], "' uses @ctrlpcs, but no control-probe PCs are available (--skip_ewastools)")
            cols
        },
        "@technical" = c("plate", "sentrix_row"),
        "@sex" = "predicted_sex",
        "@smoking" = {
            if (!"dnam_smoking" %in% names(covariates)) stop("Model '", model[["name"]], "' uses @smoking, but no DNA methylation smoking score is available (--skip_dnam_scores, or smoking not in --dnam_scores)")
            "dnam_smoking"
        },
        covariate)
}
covariate_names <- unique(unlist(lapply(as.character(unlist(model[["covariates"]])), expand)))
trait <- model[["trait"]]
block <- if (is.null(model[["block"]])) NULL else model[["block"]]
needed <- unique(c(trait, covariate_names, block))
absent <- setdiff(needed, names(covariates))
if (length(absent) > 0) stop("Model '", model[["name"]], "' uses columns without data: ", paste(absent, collapse = ", "))

# Subset, then complete cases; record why arrays are dropped
dropped <- data.table(sample = character(), reason = character())
keep <- rep(TRUE, length(arrays))
if (!is.null(model[["subset"]])) {
    column <- model[["subset"]][["column"]]
    values <- as.character(covariates[[column]])
    allowed <- as.character(unlist(if (!is.null(model[["subset"]][["equals"]])) model[["subset"]][["equals"]] else model[["subset"]][["in"]]))
    outside <- keep & !(values %in% allowed)
    dropped <- rbind(dropped, data.table(sample = arrays[outside], reason = paste0("subset: ", column, " not in ", paste(allowed, collapse = ", "))))
    keep <- keep & !outside
}
for (column in needed) {
    missing <- keep & is.na(covariates[[column]])
    if (any(missing)) dropped <- rbind(dropped, data.table(sample = arrays[missing], reason = paste0("missing ", column)))
    keep <- keep & !missing
}
n <- sum(keep)
data <- as.data.frame(covariates[keep])
rownames(data) <- data[["sample"]]

report <- list(model = model[["name"]], trait = trait, trait_type = model[["trait_type"]], n = n,
                n_dropped = sum(!keep), covariates = covariate_names, block = block)
write_outputs <- function(stats, top, design_check, status) {
    report[["status"]] <<- status
    fwrite(stats, paste0(prefix, "_stats.tsv.gz"), sep = "\\t", na = "")
    fwrite(top, paste0(prefix, "_top.tsv"), sep = "\\t", na = "")
    fwrite(design_check, paste0(prefix, "_design_check.tsv"), sep = "\\t", na = "")
    write_json(list(model = model, covariates = covariate_names, n = n, dropped = dropped),
               paste0(prefix, "_resolved.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
    write_json(report, paste0(prefix, "_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
}
if (n < model[["min_n"]]) {
    message("Model '", model[["name"]], "': ", n, " complete cases, fewer than min_n = ", model[["min_n"]], "; not tested")
    write_outputs(data.table(), data.table(), data.table(), sprintf("not tested: %d complete cases, fewer than min_n %d", n, model[["min_n"]]))
    quit(save = "no", status = 0)
}

# Trait coding and transform
if (model[["trait_type"]] == "continuous") {
    y <- as.numeric(data[[trait]])
    transform <- model[["transform"]]
    y <- switch(transform,
        none = y, log = log(y), log2 = log2(y), log10 = log10(y), log1p = log1p(y),
        zscore = as.numeric(scale(y)),
        rint = qnorm((rank(y, ties.method = "average") - 0.5) / length(y)))
    if (any(!is.finite(y))) stop("Model '", model[["name"]], "': the ", transform, " transform gives non-finite values")
    data[[trait]] <- y
} else {
    levels_present <- unique(as.character(data[[trait]]))
    if (!as.character(model[["reference"]]) %in% levels_present) {
        stop("Model '", model[["name"]], "': reference level '", model[["reference"]], "' is not among the trait values (", paste(sort(levels_present), collapse = ", "), ")")
    }
    data[[trait]] <- relevel(factor(as.character(data[[trait]])), ref = as.character(model[["reference"]]))
    if (nlevels(data[[trait]]) < 2) stop("Model '", model[["name"]], "': the trait has a single level after complete cases")
    if (model[["trait_type"]] == "binary" && nlevels(data[[trait]]) != 2) stop("Model '", model[["name"]], "': a binary trait needs exactly two levels, found ", nlevels(data[[trait]]))
}
used_covariates <- character()
for (column in covariate_names) {
    values <- data[[column]]
    if (!is.numeric(values)) values <- factor(as.character(values))
    if (length(unique(values)) < 2) {
        message("Model '", model[["name"]], "': covariate ", column, " has a single value and is left out")
        next
    }
    data[[column]] <- values
    used_covariates <- c(used_covariates, column)
}
report[["covariates_used"]] <- used_covariates
formula <- as.formula(paste("~", paste(c(paste0("`", trait, "`"), paste0("`", used_covariates, "`")), collapse = " + ")))
design <- model.matrix(formula, data = data)
if (qr(design)[["rank"]] < ncol(design)) stop("Model '", model[["name"]], "': the design matrix is not of full rank (confounded covariates)")

# Design check: trait against plate and Sentrix row; VIF of the numeric design columns
design_check <- rbindlist(lapply(intersect(c("plate", "sentrix_row"), names(data)), function(batch) {
    groups <- factor(as.character(data[[batch]]))
    if (nlevels(groups) < 2) return(NULL)
    p <- if (model[["trait_type"]] == "continuous") kruskal.test(data[[trait]], groups)[["p.value"]]
         else suppressWarnings(fisher.test(table(data[[trait]], groups), simulate.p.value = TRUE, B = 10000)[["p.value"]])
    data.table(check = paste0("trait ~ ", batch), value = nlevels(groups), p_value = signif(p, 4),
               flag = if (p < 0.01) "confounded" else "")
}))
term_columns <- setdiff(colnames(design), "(Intercept)")
if (length(term_columns) > 1) {
    vif <- sapply(term_columns, function(column) {
        others <- setdiff(term_columns, column)
        r2 <- summary(lm(design[, column] ~ design[, others]))[["r.squared"]]
        1 / (1 - r2)
    })
    design_check <- rbind(design_check, data.table(check = paste0("VIF ", sub("^`|`\$", "", term_columns)), value = round(vif, 3),
                                                   p_value = NA_real_, flag = ifelse(vif > 10, "collinear", "")), fill = TRUE)
}

# limma on M-values, with the coefficients or contrasts of the trait
mvalue <- as.matrix(assay(se, "mvalue")[, rownames(data), drop = FALSE])
beta <- as.matrix(assay(se, "beta")[, rownames(data), drop = FALSE])
correlation <- NULL
if (!is.null(block)) {
    correlation <- duplicateCorrelation(mvalue, design, block = data[[block]])[["consensus.correlation"]]
    report[["block_correlation"]] <- correlation
}
fit_m <- lmFit(mvalue, design, block = if (is.null(block)) NULL else data[[block]], correlation = correlation)
fit_b <- lmFit(beta, design, block = if (is.null(block)) NULL else data[[block]], correlation = correlation)
trait_columns <- grep(paste0("^`?", trait, "`?"), colnames(design), value = TRUE)
if (model[["trait_type"]] == "categorical" && !is.null(model[["contrasts"]])) {
    level_of <- function(level) {
        column <- paste0("`", trait, "`", level)
        if (column %in% colnames(design)) column else if (paste0(trait, level) %in% colnames(design)) paste0(trait, level) else NA
    }
    contrasts <- lapply(as.character(unlist(model[["contrasts"]])), function(contrast) {
        parts <- strsplit(contrast, "-", fixed = TRUE)[[1]]
        vector <- setNames(rep(0, ncol(design)), colnames(design))
        for (i in 1:2) {
            if (parts[i] == as.character(model[["reference"]])) next
            column <- level_of(parts[i])
            if (is.na(column)) stop("Model '", model[["name"]], "': contrast level ", parts[i], " is not a level of ", trait)
            vector[column] <- if (i == 1) 1 else -1
        }
        vector
    })
    contrast_matrix <- do.call(cbind, contrasts)
    colnames(contrast_matrix) <- as.character(unlist(model[["contrasts"]]))
} else {
    contrast_matrix <- diag(ncol(design))[, match(trait_columns, colnames(design)), drop = FALSE]
    rownames(contrast_matrix) <- colnames(design)
    colnames(contrast_matrix) <- if (model[["trait_type"]] == "continuous") trait
                                 else paste0(sub(paste0("^`?", trait, "`?"), "", trait_columns), "-", model[["reference"]])
}
fit_m <- eBayes(contrasts.fit(fit_m, contrast_matrix), trend = TRUE, robust = TRUE)
fit_b <- contrasts.fit(fit_b, contrast_matrix)

genomic_lambda <- function(p) median(qchisq(p, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)
row_info <- as.data.frame(rowData(se))
annotation <- data.table(cpg_id = rownames(se),
                         chr = row_info[["chr"]], pos = row_info[["pos"]], gene = row_info[["gene"]],
                         gene_group = row_info[["gene_group"]], cgi_relation = row_info[["cgi_relation"]])
n_per_cpg <- rowSums(!is.na(mvalue))
mean_beta <- rowMeans(beta, na.rm = TRUE)
stats <- list()
tops <- list()
report[["contrasts"]] <- list()
for (k in seq_len(ncol(contrast_matrix))) {
    contrast <- colnames(contrast_matrix)[k]
    t_stat <- fit_m[["t"]][, k]
    p <- fit_m[["p.value"]][, k]
    df <- fit_m[["df.total"]]
    z <- qnorm(pt(t_stat, df))
    z[!is.finite(z)] <- sign(t_stat[!is.finite(z)]) * 37
    bc <- bacon(teststatistics = z, na.exclude = TRUE)
    p_bacon <- as.numeric(pval(bc))
    table <- data.table(
        annotation,
        contrast = contrast,
        logfc_m = signif(fit_m[["coefficients"]][, k], 5),
        se_m = signif(fit_m[["stdev.unscaled"]][, k] * sqrt(fit_m[["s2.post"]]), 5),
        delta_beta = signif(fit_b[["coefficients"]][, k], 5),
        t = signif(t_stat, 5),
        p_value = signif(p, 5),
        fdr = signif(p.adjust(p, "BH"), 5),
        p_bacon = signif(p_bacon, 5),
        fdr_bacon = signif(p.adjust(p_bacon, "BH"), 5),
        mean_beta = round(mean_beta, 4),
        n = n_per_cpg
    )
    setorder(table, p_value)
    stats[[contrast]] <- table
    tops[[contrast]] <- head(table, n_top)
    report[["contrasts"]][[contrast]] <- list(
        n_cpgs = sum(!is.na(p)),
        lambda = round(genomic_lambda(p), 4),
        lambda_bacon = round(genomic_lambda(p_bacon), 4),
        bacon_inflation = round(as.numeric(inflation(bc)), 4),
        bacon_bias = round(as.numeric(bias(bc)), 4),
        hits_fdr = sum(table[["fdr"]] < fdr_level, na.rm = TRUE),
        hits_epigenome_wide = sum(table[["p_value"]] < p_ewide, na.rm = TRUE),
        hits_suggestive = sum(table[["p_value"]] < p_suggestive, na.rm = TRUE),
        hits_fdr_bacon = sum(table[["fdr_bacon"]] < fdr_level, na.rm = TRUE),
        hits_epigenome_wide_bacon = sum(table[["p_bacon"]] < p_ewide, na.rm = TRUE),
        top_cpg = table[["cpg_id"]][1],
        top_p = table[["p_value"]][1],
        top_delta_beta = table[["delta_beta"]][1]
    )
}
write_outputs(rbindlist(stats), rbindlist(tops), design_check, "tested")
