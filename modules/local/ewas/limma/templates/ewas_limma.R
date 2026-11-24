#!/usr/bin/env Rscript
#
# One epigenome-wide association model on the filtered, CpG-level data:
#   1. Resolve the model: covariate tokens, subset, trait transform, complete
#      cases (min_n); every dropped array is recorded with its reason
#   2. Design check: trait against plate and Sentrix row, VIF of the covariates
#   3. limma on M-values: lmFit (+ duplicateCorrelation with a block, its
#      consensus correlation estimated on --block_cpgs random CpGs), then
#      eBayes(trend = TRUE, robust = TRUE); delta-beta from the same design on
#      beta values
#   4. Inflation: genomic lambda, bacon bias and inflation with corrected P/FDR
#   5. Thresholds: FDR < 0.05, P < 9e-8 (EPIC) and P < 1e-5
#   6. Cell-type-specific associations (--cell_ewas true, binary or continuous
#      trait, blood cell fractions available): the CellDMC design of EpiDISH
#      (Zheng et al. 2018), beta ~ lineages + trait:lineages + the other
#      covariates without an intercept, fitted for all CpGs at once with limma
#      lmFit (the same estimates as CellDMC's per-CpG lm) and an ordinary t-test
#      of each interaction; BH FDR within each lineage. A lineage is marked
#      unreliable when its median fraction is below --min_cell_fraction (5%)
#      or its genomic lambda is above --max_cell_lambda (1.3): with a small
#      fraction, the interaction estimate is unstable and can give thousands of
#      false hits, even for traits without a hit in the main model
# Token expansion: @cells = cc_CD4T, cc_CD8T, cc_B, cc_NK, cc_Mono (granulocytes
# as reference); @ctrlpcs = the kept control-probe PCs (ctrlpc1..k); @technical =
# plate, sentrix_row; @sex = predicted_sex; @smoking = dnam_smoking (McCartney et al. 2018);
# @ancestry = the genetic ancestry PCs (ancestry_pc1..k, --n_ancestry_pcs).
#

prefix <- "${prefix}"
model_file <- "${model}"
se_dir <- "${se}"
samples_file <- "${samples}"
sample_qc_file <- "${sample_qc}"
phenotypes_file <- "${phenotypes}"
cells_file <- "${cell_fractions}"
ctrl_file <- "${ctrl_pcs}"
ancestry_file <- "${ancestry_pcs}"
scores_file <- "${dnam_scores}"

defaults <- list(p_epigenome_wide = "9e-8", p_suggestive = "1e-5", fdr = "0.05", n_top = "1000", block_cpgs = "50000", cell_ewas = "true", min_cell_fraction = "0.05", max_cell_lambda = "1.3")
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
    if (length(new) > 0) covariates[, (new) := table[, new, with = FALSE]]
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
if (nzchar(ancestry_file)) {
    ancestry <- read_table(ancestry_file)
    join(ancestry[, c("sample", grep("^ancestry_pc[0-9]+\$", names(ancestry), value = TRUE)), with = FALSE])
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
        "@ancestry" = {
            cols <- grep("^ancestry_pc[0-9]+\$", names(covariates), value = TRUE)
            if (length(cols) == 0) stop("Model '", model[["name"]], "' uses @ancestry, but no ancestry PCs are available (--skip_ancestry)")
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
characteristics <- data.table(variable = character(), level = character(), p_value = numeric())
empty_cells <- data.table(cpg_id = character(), chr = character(), pos = integer(), gene = character(), contrast = character(), cell = character(),
                          estimate = numeric(), se = numeric(), t = numeric(), p_value = numeric(), fdr = numeric(), lineage_status = character())
# fwrite leaves a header-only .gz without its gzip trailer: write those through a connection
write_gz <- function(table, file) {
    if (nrow(table) > 0) return(fwrite(table, file, sep = "\\t", na = ""))
    connection <- gzfile(file, "w")
    if (ncol(table) > 0) writeLines(paste(names(table), collapse = "\\t"), connection)
    close(connection)
}
write_outputs <- function(stats, top, design_check, status, cells = empty_cells) {
    report[["status"]] <<- status
    write_gz(stats, paste0(prefix, "_stats.tsv.gz"))
    write_gz(cells, paste0(prefix, "_cells.tsv.gz"))
    fwrite(top, paste0(prefix, "_top.tsv"), sep = "\\t", na = "")
    fwrite(design_check, paste0(prefix, "_design_check.tsv"), sep = "\\t", na = "")
    fwrite(characteristics, paste0(prefix, "_characteristics.tsv"), sep = "\\t", na = "")
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
# Covariates that follow the trait (Uddin et al. 2022: adjusting for methylation
# axes that the trait itself drives can mask its association). Every numeric
# covariate is tested against the trait; methylation-derived axes (cell
# fractions, control-probe and ancestry PCs, DNA methylation scores and other
# PC-like columns) are flagged at P < 0.001.
derived <- function(column) grepl("^(cc_|ctrlpc[0-9]|ancestry_pc[0-9]|dnam_|PC[0-9]|pc[0-9])", column)
trait_checks <- rbindlist(lapply(used_covariates, function(column) {
    values <- data[[column]]
    if (!is.numeric(values)) return(NULL)
    test <- if (model[["trait_type"]] == "continuous") suppressWarnings(cor.test(values, data[[trait]], method = "spearman"))
            else kruskal.test(values, data[[trait]])
    p_cov <- test[["p.value"]]
    value <- if (model[["trait_type"]] == "continuous") unname(test[["estimate"]])
             else summary(lm(values ~ data[[trait]]))[["r.squared"]]
    data.table(check = paste0(column, " ~ trait"), value = round(value, 4), p_value = signif(p_cov, 4),
               flag = if (p_cov < 0.001 && derived(column)) "follows the trait" else "")
}))
if (nrow(trait_checks) > 0) design_check <- rbind(design_check, trait_checks, fill = TRUE)
flagged <- design_check[!is.na(flag) & nzchar(flag)]

# Baseline characteristics of the model's arrays (Table 1 of Uddin et al. 2022
# and Kirmani et al. 2025): by trait group for a binary or categorical trait,
# overall for a continuous one. Aggregates only; any count from 1 to 4 is shown
# as "<5", as are means of fewer than 5 arrays.
describe_columns <- c(if (model[["trait_type"]] == "continuous") trait, used_covariates[!grepl("^(ctrlpc|ancestry_pc|PC|pc)[0-9]", used_covariates)])
group <- if (model[["trait_type"]] == "continuous") rep("All arrays", n) else as.character(data[[trait]])
group_levels <- if (model[["trait_type"]] == "continuous") "All arrays" else levels(data[[trait]])
count_text <- function(k, total) if (k == 0) "0" else if (k < 5) "<5" else sprintf("%d (%.1f%%)", k, 100 * k / total)
characteristics <- rbindlist(c(
    list(as.data.table(c(list(variable = "Arrays", level = ""), setNames(lapply(group_levels, function(g) {
        k <- sum(group == g)
        if (k < 5) "<5" else as.character(k)
    }), group_levels), list(p_value = NA_real_)))),
    lapply(describe_columns, function(column) {
        values <- data[[column]]
        compare <- model[["trait_type"]] != "continuous" && column != trait
        if (is.numeric(values)) {
            cells <- setNames(lapply(group_levels, function(g) {
                v <- values[group == g]
                if (length(v) < 5) "<5" else sprintf("%.3g (%.3g)", mean(v), sd(v))
            }), group_levels)
            p_cell <- if (compare) kruskal.test(values, factor(group))[["p.value"]] else NA_real_
            list(as.data.table(c(list(variable = column, level = "mean (SD)"), cells, list(p_value = signif(p_cell, 4)))))
        } else {
            values <- as.character(values)
            p_cell <- if (compare && length(unique(values)) > 1) suppressWarnings(fisher.test(table(values, group), simulate.p.value = TRUE, B = 2000)[["p.value"]]) else NA_real_
            lapply(sort(unique(values)), function(level) as.data.table(c(list(variable = column, level = level),
                setNames(lapply(group_levels, function(g) count_text(sum(values == level & group == g), sum(group == g))), group_levels),
                list(p_value = if (level == sort(unique(values))[1]) signif(p_cell, 4) else NA_real_))))
        }
    }) |> unlist(recursive = FALSE)
), fill = TRUE)
report[["design_flags"]] <- if (nrow(flagged) > 0) paste(paste0(flagged[["check"]], " (", flagged[["flag"]], ")"), collapse = "; ") else ""

# limma on M-values, with the coefficients or contrasts of the trait
mvalue <- as.matrix(assay(se, "mvalue")[, rownames(data), drop = FALSE])
beta <- as.matrix(assay(se, "beta")[, rownames(data), drop = FALSE])
correlation <- NULL
if (!is.null(block)) {
    # The consensus correlation is one trimmed mean over CpGs, so a fixed random
    # subset gives the same estimate at a fraction of the cost (0 = every CpG)
    block_cpgs <- as.integer(opts[["block_cpgs"]])
    rows <- if (block_cpgs > 0 && nrow(mvalue) > block_cpgs) sort(sample.int(nrow(mvalue), block_cpgs)) else seq_len(nrow(mvalue))
    correlation <- duplicateCorrelation(mvalue[rows, , drop = FALSE], design, block = data[[block]])[["consensus.correlation"]]
    report[["block_correlation"]] <- correlation
    report[["block_correlation_cpgs"]] <- length(rows)
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
    # log-scale lower tail of both signs: pt() of a large positive t rounds to 1 and its z to infinity
    z <- -sign(t_stat) * qnorm(pt(-abs(t_stat), df, log.p = TRUE), log.p = TRUE)
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
        hits_higher = sum(table[["fdr"]] < fdr_level & table[["delta_beta"]] > 0, na.rm = TRUE),
        hits_lower = sum(table[["fdr"]] < fdr_level & table[["delta_beta"]] < 0, na.rm = TRUE),
        hits_epigenome_wide = sum(table[["p_value"]] < p_ewide, na.rm = TRUE),
        hits_suggestive = sum(table[["p_value"]] < p_suggestive, na.rm = TRUE),
        hits_fdr_bacon = sum(table[["fdr_bacon"]] < fdr_level, na.rm = TRUE),
        hits_epigenome_wide_bacon = sum(table[["p_bacon"]] < p_ewide, na.rm = TRUE),
        top_cpg = table[["cpg_id"]][1],
        top_p = table[["p_value"]][1],
        top_delta_beta = table[["delta_beta"]][1]
    )
}
# Cell-type-specific associations (CellDMC design)
cells_table <- empty_cells
lineage_columns <- grep("^cc_", names(data), value = TRUE)
cell_note <- (if (opts[["cell_ewas"]] != "true") "skipped: --skip_cell_ewas"
              else if (length(lineage_columns) < 2) "not run: no cell fractions"
              else if (anyNA(data[, lineage_columns])) "not run: cell fractions missing for some arrays"
              else if (ncol(contrast_matrix) != 1 || !model[["trait_type"]] %in% c("continuous", "binary")) "not run: needs a binary or continuous trait"
              else "")
if (!nzchar(cell_note)) {
    fractions <- as.matrix(data[, lineage_columns])
    colnames(fractions) <- sub("^cc_", "", lineage_columns)
    pheno <- if (model[["trait_type"]] == "continuous") as.numeric(data[[trait]]) else as.numeric(data[[trait]] != as.character(model[["reference"]]))
    cell_design <- model.matrix(~ fractions + pheno:fractions)[, -1, drop = FALSE]
    others <- setdiff(used_covariates, lineage_columns)
    if (length(others) > 0) {
        cell_design <- cbind(cell_design, model.matrix(as.formula(paste("~", paste0("`", others, "`", collapse = " + "))), data = data)[, -1, drop = FALSE])
    }
    interaction_columns <- grep(":pheno\$", colnames(cell_design), value = TRUE)
    if (qr(cell_design)[["rank"]] < ncol(cell_design)) {
        cell_note <- "not run: the cell-type design is not of full rank (a lineage is constant or follows another)"
    } else if (n - ncol(cell_design) < 10) {
        cell_note <- sprintf("not run: %d arrays for %d design columns", n, ncol(cell_design))
    } else {
        fit_cells <- lmFit(beta, cell_design, block = if (is.null(block)) NULL else data[[block]], correlation = correlation)
        cells_table <- rbindlist(lapply(interaction_columns, function(column) {
            estimate <- fit_cells[["coefficients"]][, column]
            se_cell <- fit_cells[["stdev.unscaled"]][, column] * fit_cells[["sigma"]]
            t_cell <- estimate / se_cell
            p_cell <- 2 * pt(-abs(t_cell), fit_cells[["df.residual"]])
            table <- data.table(annotation[, .(cpg_id, chr, pos, gene)], contrast = colnames(contrast_matrix)[1],
                                cell = sub("^fractions", "", sub(":pheno\$", "", column)),
                                estimate = signif(estimate, 5), se = signif(se_cell, 5), t = signif(t_cell, 5),
                                p_value = signif(p_cell, 5), fdr = signif(p.adjust(p_cell, "BH"), 5))
            setorder(table, p_value)
            table
        }))
        # lineages too small, or too inflated, to trust
        median_fraction <- apply(fractions, 2, median)
        lineage_status <- vapply(names(median_fraction), function(l) {
            p <- cells_table[["p_value"]][cells_table[["cell"]] == l]
            lambda <- genomic_lambda(p)
            reasons <- c(if (median_fraction[[l]] < as.numeric(opts[["min_cell_fraction"]])) sprintf("median fraction %.1f%% below %.0f%%", 100 * median_fraction[[l]], 100 * as.numeric(opts[["min_cell_fraction"]])),
                         if (is.finite(lambda) && lambda > as.numeric(opts[["max_cell_lambda"]])) sprintf("lambda %.2f above %s", lambda, opts[["max_cell_lambda"]]))
            if (length(reasons)) paste0("unreliable: ", paste(reasons, collapse = "; ")) else "ok"
        }, character(1))
        cells_table[, lineage_status := lineage_status[cell]]
        cell_note <- "run"
    }
}
report[["cells"]] <- list(status = cell_note)
if (nrow(cells_table) > 0) {
    report[["cells"]][["lineages"]] <- lapply(split(cells_table, by = "cell"), function(x) list(
        hits_fdr = sum(x[["fdr"]] < fdr_level, na.rm = TRUE),
        hits_higher = sum(x[["fdr"]] < fdr_level & x[["estimate"]] > 0, na.rm = TRUE),
        hits_lower = sum(x[["fdr"]] < fdr_level & x[["estimate"]] < 0, na.rm = TRUE),
        lambda = round(genomic_lambda(x[["p_value"]]), 4),
        median_fraction = round(median_fraction[[x[["cell"]][1]]], 4),
        status = x[["lineage_status"]][1],
        top_cpg = x[["cpg_id"]][1], top_p = x[["p_value"]][1], top_estimate = x[["estimate"]][1]))
}
write_outputs(rbindlist(stats), rbindlist(tops), design_check, "tested", cells_table)
