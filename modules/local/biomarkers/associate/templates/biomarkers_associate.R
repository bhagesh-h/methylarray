#!/usr/bin/env Rscript
#
# Epigenetic clocks, DNA methylation scores and blood cell fractions, related to
# each other, to the technical batches and to the phenotypes, on the arrays that
# passed sample QC. Only aggregates are written (correlations, test results,
# group medians), never values of single arrays.
#   - acceleration: a clock's residual on reported age (the EAA of the clock
#     table where it exists); with cells held: also adjusted for the six blood
#     lineages but granulocytes (the IEAA of the clock table where it exists)
#   - agreement: Spearman correlation of every pair of accelerations and scores
#   - cells: Spearman correlation of each cell lineage with each acceleration
#     and score; each lineage against age, sex and every phenotype, with age and
#     sex held (estimate in SD of the fraction)
#   - batch: each acceleration against DNA plate, BeadChip and Sentrix row
#     (adjusted R-squared and F-test)
#   - phenotypes: each acceleration and score against every phenotype, plain and
#     with cells held: difference in SD of acceleration between the two groups of
#     a binary phenotype (t-test), Spearman correlation for a continuous one
#   - contamination: Spearman correlation of the ewastools SNP outlier score with
#     each acceleration and score, and the phenotype associations again without
#     the arrays above the median plus three MADs of the score (borderline
#     contamination below the FAIL threshold), as the 'low contamination' version
#   - ancestry: Spearman correlation of each acceleration and score with each
#     genetic ancestry PC, and the phenotype associations again with the
#     measures residualised on the ancestry PCs ('ancestry held'). Clocks are
#     less accurate in admixed individuals (Cruz-Gonzalez et al. 2026), so a
#     difference in acceleration between ancestries can be calibration error
#   - genotype-robust clocks (with --meqtl): each acceleration against the one
#     from the same clock without the CpGs with a common SNP or cis-meQTL
#     (Spearman), and the phenotype associations with the robust accelerations
#     ('genotype-robust')
#   - quartile groups: bottom quarter (decelerated), middle half (normal) and top
#     quarter (accelerated) of each acceleration, with the median (or percentage)
#     of every characteristic per group and a trend test across the ordered
#     groups (Kendall's tau, equivalent to the Jonckheere-Terpstra test)
# P-values are adjusted with Benjamini-Hochberg within each table and version.
# Phenotypes: numeric columns with more than two values are continuous, columns
# with two values are binary (the reference is "no", "0", "false" or the first
# value in sort order); columns with more levels, or fewer than --min_n values,
# are left out.
#

prefix <- "${prefix}"
samples_file <- "${samples}"
sample_qc_file <- "${sample_qc}"
clocks_file <- "${clocks}"
cells_file <- "${cell_fractions}"
scores_file <- "${dnam_scores}"
phenotypes_file <- "${phenotypes}"
ancestry_file <- "${ancestry_pcs}"
robust_file <- "${clocks_robust}"

defaults <- list(min_n = "20")
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
min_n <- as.integer(opts[["min_n"]])

suppressPackageStartupMessages(library(jsonlite))
writeLines(c('"${task.process}":', paste0("    R: ", getRversion()), paste0("    jsonlite: ", packageVersion("jsonlite"))), "versions.yml")

read_table <- function(file) {
    if (!nzchar(file) || !file.exists(file) || file.info(file)[["size"]] == 0) return(NULL)
    sep <- if (grepl("[.]csv\$", file)) "," else "\\t"
    read.delim(file, sep = sep, check.names = FALSE, na.strings = c("", "NA"), stringsAsFactors = FALSE)
}

# Arrays that passed sample QC, with the samplesheet fields
sample_qc <- read_table(sample_qc_file)
retained <- sample_qc[["sample"]][sample_qc[["status"]] %in% c("PASS", "FLAG")]
data <- read_table(samples_file)
data <- data[data[["sample"]] %in% retained, , drop = FALSE]
data <- data[order(data[["sample"]]), , drop = FALSE]
merge_in <- function(table, columns = NULL) {
    if (is.null(table)) return(invisible(NULL))
    table <- table[match(data[["sample"]], table[["sample"]]), , drop = FALSE]
    keep <- setdiff(if (is.null(columns)) names(table) else intersect(columns, names(table)), names(data))
    for (column in keep) data[[column]] <<- table[[column]]
}
clocks <- read_table(clocks_file)
cells <- read_table(cells_file)
scores <- read_table(scores_file)
phenotypes <- read_table(phenotypes_file)
merge_in(clocks)
merge_in(cells, grep("^lineage_", names(cells), value = TRUE))
merge_in(scores)
phenotype_columns <- if (is.null(phenotypes)) character() else setdiff(names(phenotypes), c("sample", names(data)))
merge_in(phenotypes)

age <- suppressWarnings(as.numeric(data[["reported_age"]]))
men <- if ("reported_sex" %in% names(data)) ifelse(data[["reported_sex"]] %in% "M", 1, ifelse(data[["reported_sex"]] %in% "F", 0, NA)) else rep(NA_real_, nrow(data))
lineages <- grep("^lineage_", names(data), value = TRUE)
held_cells <- setdiff(lineages, "lineage_Gran")

# Ageing measures: clocks (as acceleration when ages are known) and scores
clock_columns <- if (is.null(clocks)) character() else setdiff(names(clocks), c("sample", "age", grep("^i?eaa_", names(clocks), value = TRUE)))
score_columns <- grep("^dnam_", names(data), value = TRUE)
residual <- function(y, design) {
    d <- data.frame(y = y, design)
    ok <- complete.cases(d)
    out <- rep(NA_real_, length(y))
    if (sum(ok) - ncol(d) < 5) return(out)
    out[ok] <- residuals(lm(y ~ ., data = d[ok, , drop = FALSE]))
    out
}
measures <- list()
measures_cells <- list()
for (clock in clock_columns) {
    y <- suppressWarnings(as.numeric(data[[clock]]))
    if (sum(is.finite(y)) < min_n) next
    has_age <- sum(is.finite(age)) >= min_n
    name <- if (has_age) paste(clock, "acceleration") else clock
    measures[[name]] <- if (paste0("eaa_", clock) %in% names(data)) as.numeric(data[[paste0("eaa_", clock)]])
                        else if (has_age) residual(y, data.frame(age = age)) else y
    if (length(held_cells) > 0) {
        measures_cells[[name]] <- if (paste0("ieaa_", clock) %in% names(data)) as.numeric(data[[paste0("ieaa_", clock)]])
                                  else residual(y, cbind(data.frame(age = if (has_age) age else 0), data[, held_cells, drop = FALSE]))
    }
}
for (score in score_columns) {
    y <- as.numeric(data[[score]])
    if (sum(is.finite(y)) < min_n) next
    measures[[score]] <- y
    if (length(held_cells) > 0) measures_cells[[score]] <- residual(y, data[, held_cells, drop = FALSE])
}

# Phenotypes: continuous or binary (coded 0 and 1)
reference_level <- function(levels) {
    known <- levels[tolower(levels) %in% c("no", "0", "false", "f", "n", "never", "control", "none")]
    if (length(known) > 0) known[1] else sort(levels)[1]
}
traits <- list()
if (sum(is.finite(age)) >= min_n) traits[["Age"]] <- list(values = age, type = "continuous")
if (sum(!is.na(men)) >= min_n && length(unique(na.omit(men))) == 2) traits[["Men"]] <- list(values = men, type = "binary")
for (column in phenotype_columns) {
    x <- data[[column]]
    values <- unique(na.omit(x))
    if (sum(!is.na(x)) < min_n || length(values) < 2) next
    numeric <- suppressWarnings(as.numeric(as.character(x)))
    if (length(values) == 2) {
        reference <- reference_level(as.character(values))
        other <- setdiff(as.character(values), reference)
        traits[[paste0(column, " (", other, ")")]] <- list(values = ifelse(is.na(x), NA, as.numeric(as.character(x) == other)), type = "binary")
    } else if (all(is.na(numeric) == is.na(x))) {
        traits[[column]] <- list(values = numeric, type = "continuous")
    }
}
bh <- function(p) p.adjust(p, "BH")
spearman <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < min_n) return(c(rho = NA, p = NA, n = sum(ok)))
    test <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    c(rho = unname(test[["estimate"]]), p = test[["p.value"]], n = sum(ok))
}
write_tsv <- function(table, file) write.table(table, file, sep = "\\t", quote = FALSE, row.names = FALSE, na = "")
round_table <- function(table) {
    for (column in names(table)) if (is.numeric(table[[column]])) table[[column]] <- signif(table[[column]], 4)
    table
}
empty <- function(columns) as.data.frame(setNames(replicate(length(columns), character(), simplify = FALSE), columns))

# Agreement between the ageing measures
agreement <- empty(c("measure1", "measure2", "rho", "p", "n"))
if (length(measures) > 1) {
    pairs <- t(combn(names(measures), 2))
    agreement <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
        r <- spearman(measures[[pairs[i, 1]]], measures[[pairs[i, 2]]])
        data.frame(measure1 = pairs[i, 1], measure2 = pairs[i, 2], rho = r[["rho"]], p = r[["p"]], n = r[["n"]])
    }))
}
write_tsv(round_table(agreement), paste0(prefix, "_agreement.tsv"))

# Cell lineages against the ageing measures
cells_measures <- empty(c("cell", "measure", "rho", "p", "n", "q"))
if (length(lineages) > 0 && length(measures) > 0) {
    cells_measures <- do.call(rbind, lapply(lineages, function(cell) do.call(rbind, lapply(names(measures), function(m) {
        r <- spearman(as.numeric(data[[cell]]), measures[[m]])
        data.frame(cell = sub("^lineage_", "", cell), measure = m, rho = r[["rho"]], p = r[["p"]], n = r[["n"]])
    }))))
    cells_measures[["q"]] <- bh(cells_measures[["p"]])
}
write_tsv(round_table(cells_measures), paste0(prefix, "_cells_measures.tsv"))

# Cell lineages against age, sex and the phenotypes, with age and sex held
cells_traits <- empty(c("cell", "trait", "type", "estimate", "p", "n", "q"))
if (length(lineages) > 0 && length(traits) > 0) {
    cells_traits <- do.call(rbind, lapply(lineages, function(cell) do.call(rbind, lapply(names(traits), function(t) {
        y <- as.numeric(data[[cell]])
        y <- (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)
        x <- traits[[t]][["values"]]
        if (traits[[t]][["type"]] == "continuous") x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
        covariates <- data.frame(age = age, men = men)
        if (t == "Age") covariates[["age"]] <- NULL
        if (t == "Men") covariates[["men"]] <- NULL
        covariates <- covariates[, vapply(covariates, function(v) sum(!is.na(v)) >= min_n && length(unique(na.omit(v))) > 1, logical(1)), drop = FALSE]
        d <- data.frame(y = y, x = x, covariates)
        d <- d[complete.cases(d), , drop = FALSE]
        if (nrow(d) < min_n || length(unique(d[["x"]])) < 2) return(NULL)
        fit <- summary(lm(y ~ ., data = d))[["coefficients"]]
        data.frame(cell = sub("^lineage_", "", cell), trait = t, type = traits[[t]][["type"]], estimate = fit["x", 1], p = fit["x", 4], n = nrow(d))
    }))))
    if (!is.null(cells_traits) && nrow(cells_traits) > 0) cells_traits[["q"]] <- bh(cells_traits[["p"]]) else cells_traits <- empty(c("cell", "trait", "type", "estimate", "p", "n", "q"))
}
write_tsv(round_table(cells_traits), paste0(prefix, "_cells_phenotypes.tsv"))

# Accelerations against the technical batches
f_test <- function(y, g) {
    ok <- is.finite(y) & !is.na(g)
    g <- factor(g[ok])
    if (sum(ok) < min_n || nlevels(g) < 2 || nlevels(g) >= sum(ok) - 1) return(NULL)
    fit <- summary(lm(y[ok] ~ g))
    f <- fit[["fstatistic"]]
    c(adj_r2 = fit[["adj.r.squared"]], p = unname(pf(f[1], f[2], f[3], lower.tail = FALSE)), levels = nlevels(g), n = sum(ok))
}
batches <- list()
if ("plate" %in% names(data)) batches[["DNA plate"]] <- data[["plate"]]
if ("sentrix_id" %in% names(data)) batches[["BeadChip"]] <- data[["sentrix_id"]]
if ("sentrix_position" %in% names(data)) batches[["Sentrix row"]] <- substr(data[["sentrix_position"]], 1, 3)
batch <- empty(c("measure", "batch", "adj_r2", "p", "levels", "n"))
if (length(batches) > 0 && length(measures) > 0) {
    batch <- do.call(rbind, lapply(names(measures), function(m) do.call(rbind, lapply(names(batches), function(b) {
        r <- f_test(measures[[m]], batches[[b]])
        if (is.null(r)) return(NULL)
        data.frame(measure = m, batch = b, adj_r2 = r[["adj_r2"]], p = r[["p"]], levels = r[["levels"]], n = r[["n"]])
    }))))
    if (is.null(batch)) batch <- empty(c("measure", "batch", "adj_r2", "p", "levels", "n"))
}
write_tsv(round_table(batch), paste0(prefix, "_batch.tsv"))
if (nrow(batch) > 0) {
    # MultiQC: one row per measure, the F-test P of each batch variable
    wide <- reshape(batch[, c("measure", "batch", "p")], idvar = "measure", timevar = "batch", direction = "wide")
    names(wide) <- sub("^p[.]", "P ", names(wide))
    names(wide)[1] <- "Measure"
    writeLines(c(
        "# id: 'biomarker_batch'",
        "# section_name: 'Clocks and scores: technical batches'",
        "# description: 'F-test P of each acceleration and score against the DNA plate, BeadChip and Sentrix row. A small P means the batch explains part of the measure; then include it in the models.'",
        "# plot_type: 'table'"
    ), paste0(prefix, "_batch_mqc.tsv"))
    write.table(round_table(wide), paste0(prefix, "_batch_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE, na = "")
}

# Contamination: DNA from another person pulls beta values towards the middle
# and can bias clocks. The SNP outlier score against every measure, and the
# arrays above a robust cutoff among the analysed arrays
snp_outlier <- if ("snp_outlier" %in% names(sample_qc)) {
    suppressWarnings(as.numeric(sample_qc[["snp_outlier"]][match(data[["sample"]], sample_qc[["sample"]])]))
} else rep(NA_real_, nrow(data))
contamination <- empty(c("measure", "rho", "p", "n", "q", "cutoff", "arrays_above_cutoff"))
above_cutoff <- rep(FALSE, nrow(data))
if (sum(is.finite(snp_outlier)) >= min_n && length(measures) > 0) {
    spread <- mad(snp_outlier, na.rm = TRUE)
    cutoff <- median(snp_outlier, na.rm = TRUE) + 3 * spread
    above_cutoff <- spread > 0 & is.finite(snp_outlier) & snp_outlier > cutoff
    contamination <- do.call(rbind, lapply(names(measures), function(m) {
        r <- spearman(snp_outlier, measures[[m]])
        data.frame(measure = m, rho = r[["rho"]], p = r[["p"]], n = r[["n"]])
    }))
    contamination[["q"]] <- bh(contamination[["p"]])
    contamination[["cutoff"]] <- cutoff
    contamination[["arrays_above_cutoff"]] <- sum(above_cutoff)
}
write_tsv(round_table(contamination), paste0(prefix, "_contamination.tsv"))

# Ancestry: the genetic ancestry PCs against every measure
ancestry <- read_table(ancestry_file)
ancestry_pcs <- if (is.null(ancestry)) character() else grep("^ancestry_pc[0-9]+\$", names(ancestry), value = TRUE)
ancestry_table <- empty(c("measure", "pc", "rho", "p", "n", "q"))
measures_ancestry <- list()
if (length(ancestry_pcs) > 0 && length(measures) > 0) {
    pcs <- as.data.frame(ancestry[match(data[["sample"]], ancestry[["sample"]]), ancestry_pcs, drop = FALSE])
    ancestry_table <- do.call(rbind, lapply(names(measures), function(m) do.call(rbind, lapply(ancestry_pcs, function(pc) {
        r <- spearman(measures[[m]], pcs[[pc]])
        data.frame(measure = m, pc = pc, rho = r[["rho"]], p = r[["p"]], n = r[["n"]])
    }))))
    ancestry_table[["q"]] <- bh(ancestry_table[["p"]])
    measures_ancestry <- lapply(measures, function(a) residual(a, pcs))
}
write_tsv(round_table(ancestry_table), paste0(prefix, "_ancestry.tsv"))

# Genotype-robust clocks: the same measure from the clocks without the CpGs
# with a common SNP or a cis-meQTL
robust_clocks <- read_table(robust_file)
measures_robust <- list()
robust_table <- empty(c("measure", "rho", "n"))
if (!is.null(robust_clocks)) {
    robust_clocks <- robust_clocks[match(data[["sample"]], robust_clocks[["sample"]]), , drop = FALSE]
    has_age <- sum(is.finite(age)) >= min_n
    for (clock in intersect(clock_columns, names(robust_clocks))) {
        name <- if (has_age) paste(clock, "acceleration") else clock
        if (!name %in% names(measures)) next
        y <- suppressWarnings(as.numeric(robust_clocks[[clock]]))
        measures_robust[[name]] <- if (paste0("eaa_", clock) %in% names(robust_clocks)) as.numeric(robust_clocks[[paste0("eaa_", clock)]])
                                   else if (has_age) residual(y, data.frame(age = age)) else y
    }
    if (length(measures_robust) > 0) {
        robust_table <- do.call(rbind, lapply(names(measures_robust), function(m) {
            r <- spearman(measures[[m]], measures_robust[[m]])
            data.frame(measure = m, rho = r[["rho"]], n = r[["n"]])
        }))
    }
}
write_tsv(round_table(robust_table), paste0(prefix, "_robust.tsv"))

# Ageing measures against the phenotypes, plain, with cells held, with the
# ancestry PCs held and without the arrays above the contamination cutoff
association <- function(a, trait) {
    a <- a / sd(a, na.rm = TRUE)
    x <- trait[["values"]]
    if (trait[["type"]] == "binary") {
        ok <- is.finite(a) & !is.na(x)
        if (sum(ok) < min_n || min(table(x[ok])) < 3) return(NULL)
        test <- t.test(a[ok & x == 1], a[ok & x == 0], var.equal = TRUE)
        return(c(estimate = unname(diff(rev(test[["estimate"]]))), p = test[["p.value"]], n = sum(ok)))
    }
    r <- spearman(a, x)
    if (is.na(r[["p"]])) return(NULL)
    c(estimate = r[["rho"]], p = r[["p"]], n = r[["n"]])
}
versions <- list(plain = measures, `cells held` = measures_cells)
if (length(measures_ancestry) > 0) versions[["ancestry held"]] <- measures_ancestry
if (length(measures_robust) > 0) versions[["genotype-robust"]] <- measures_robust
if (any(above_cutoff)) versions[["low contamination"]] <- lapply(measures, function(a) replace(a, above_cutoff, NA))
measure_traits <- empty(c("measure", "trait", "type", "version", "estimate", "p", "n", "q"))
rows <- list()
for (version in names(versions)) for (m in names(versions[[version]])) for (t in names(traits)) {
    if (t == "Age" && grepl("acceleration\$", m)) next
    r <- association(versions[[version]][[m]], traits[[t]])
    if (is.null(r)) next
    rows[[length(rows) + 1]] <- data.frame(measure = m, trait = t, type = traits[[t]][["type"]], version = version,
                                           estimate = r[["estimate"]], p = r[["p"]], n = r[["n"]])
}
if (length(rows) > 0) {
    measure_traits <- do.call(rbind, rows)
    measure_traits[["q"]] <- ave(measure_traits[["p"]], measure_traits[["version"]], measure_traits[["trait"]], FUN = bh)
}
write_tsv(round_table(measure_traits), paste0(prefix, "_phenotypes.tsv"))

# Quartile groups of acceleration
group_names <- c("decelerated", "normal", "accelerated")
quartiles <- empty(c("measure", "trait", "type", "group", "value", "mean_z", "p_trend", "q"))
rows <- list()
for (m in grep("acceleration\$", names(measures), value = TRUE)) {
    a <- measures[[m]]
    if (sum(is.finite(a)) < 4 * min_n / 2) next
    cut_low <- quantile(a, 0.25, na.rm = TRUE)
    cut_high <- quantile(a, 0.75, na.rm = TRUE)
    group <- ifelse(is.na(a), NA, ifelse(a <= cut_low, 1L, ifelse(a >= cut_high, 3L, 2L)))
    for (t in names(traits)) {
        x <- traits[[t]][["values"]]
        ok <- !is.na(group) & is.finite(x)
        if (sum(ok) < min_n || length(unique(x[ok])) < 2) next
        test <- suppressWarnings(cor.test(group[ok], x[ok], method = "kendall", exact = FALSE))
        z <- (x - mean(x[ok])) / sd(x[ok])
        for (k in 1:3) {
            in_group <- ok & group == k
            value <- if (traits[[t]][["type"]] == "binary") 100 * mean(x[in_group]) else median(x[in_group])
            rows[[length(rows) + 1]] <- data.frame(measure = m, trait = t, type = traits[[t]][["type"]], group = group_names[k],
                                                   value = value, mean_z = mean(z[in_group]), p_trend = test[["p.value"]])
        }
    }
}
if (length(rows) > 0) {
    quartiles <- do.call(rbind, rows)
    first <- !duplicated(paste(quartiles[["measure"]], quartiles[["trait"]]))
    q <- ave(ifelse(first, quartiles[["p_trend"]], NA), quartiles[["measure"]], FUN = function(p) p.adjust(p, "BH"))
    quartiles[["q"]] <- ave(q, quartiles[["measure"]], quartiles[["trait"]], FUN = function(v) v[!is.na(v)][1])
}
write_tsv(round_table(quartiles), paste0(prefix, "_quartiles.tsv"))

# MultiQC heatmaps
heatmap_json <- function(id, section, description, title, rows, columns, values, file, min = -1, max = 1) {
    write_json(list(id = id, section_name = section, description = description, plot_type = "heatmap",
                    pconfig = list(id = paste0(id, "_plot"), title = title, min = min, max = max, square = FALSE, decimalPlaces = 2),
                    xcats = I(columns), ycats = I(rows), data = values),
               file, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
}
if (nrow(agreement) > 0) {
    names_m <- names(measures)
    lookup <- function(a, b) {
        if (a == b) return(1)
        i <- which((agreement[["measure1"]] == a & agreement[["measure2"]] == b) | (agreement[["measure1"]] == b & agreement[["measure2"]] == a))
        if (length(i) == 0) NA else round(agreement[["rho"]][i[1]], 2)
    }
    values <- lapply(names_m, function(a) lapply(names_m, function(b) lookup(a, b)))
    heatmap_json("biomarker_agreement", "Agreement of the ageing measures",
                 "Spearman correlation between the accelerations of the clocks (residual on age) and the DNA methylation scores, on the arrays that passed QC.",
                 "Agreement of accelerations and scores (Spearman)", names_m, names_m, values, paste0(prefix, "_agreement_mqc.json"))
}
if (nrow(cells_measures) > 0) {
    cell_names <- unique(cells_measures[["cell"]])
    values <- lapply(cell_names, function(cell) lapply(names(measures), function(m) {
        v <- cells_measures[["rho"]][cells_measures[["cell"]] == cell & cells_measures[["measure"]] == m]
        if (length(v) == 0 || is.na(v)) NA else round(v, 2)
    }))
    heatmap_json("biomarker_cells", "Cell fractions and ageing measures",
                 "Spearman correlation between the EpiDISH cell lineage fractions and the accelerations and scores.",
                 "Cell fractions against accelerations and scores (Spearman)", cell_names, names(measures), values, paste0(prefix, "_cells_mqc.json"))
}
plain <- measure_traits[measure_traits[["version"]] == "plain", , drop = FALSE]
if (nrow(plain) > 0) {
    trait_names <- unique(plain[["trait"]])
    measure_names <- unique(plain[["measure"]])
    values <- lapply(measure_names, function(m) lapply(trait_names, function(t) {
        i <- which(plain[["measure"]] == m & plain[["trait"]] == t)
        if (length(i) == 0) NA else round(plain[["estimate"]][i[1]], 2)
    }))
    heatmap_json("biomarker_phenotypes", "Ageing measures and phenotypes",
                 "Each acceleration and score against every phenotype: Spearman correlation for continuous phenotypes, difference in SD of acceleration between the groups for binary ones. Tested values with q are in biomarkers_phenotypes.tsv.",
                 "Accelerations and scores against the phenotypes", measure_names, trait_names, values, paste0(prefix, "_phenotypes_mqc.json"))
}
if (nrow(contamination) > 0) {
    write_json(list(id = "biomarker_contamination", section_name = "Contamination and ageing measures",
                    description = paste0("Spearman correlation of the ewastools SNP outlier score (higher: more likely contaminated) with each acceleration and score, on the arrays that passed QC. ",
                                         contamination[["arrays_above_cutoff"]][1], " arrays lie above the robust cutoff of ", signif(contamination[["cutoff"]][1], 3),
                                         " (median plus three MADs); the phenotype associations without them are the 'low contamination' version in biomarkers_phenotypes.tsv."),
                    plot_type = "table", pconfig = list(id = "biomarker_contamination_table", title = "Contamination against accelerations and scores"),
                    data = setNames(lapply(seq_len(nrow(contamination)), function(i) list(rho = round(contamination[["rho"]][i], 3), q = signif(contamination[["q"]][i], 3),
                                                                                         n = contamination[["n"]][i])), contamination[["measure"]])),
               paste0(prefix, "_contamination_mqc.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
}
if (nrow(ancestry_table) > 0) {
    measure_names <- unique(ancestry_table[["measure"]])
    values <- lapply(measure_names, function(m) lapply(ancestry_pcs, function(pc) {
        v <- ancestry_table[["rho"]][ancestry_table[["measure"]] == m & ancestry_table[["pc"]] == pc]
        if (length(v) == 0 || is.na(v)) NA else round(v, 2)
    }))
    heatmap_json("biomarker_ancestry", "Ancestry and ageing measures",
                 "Spearman correlation of each acceleration and score with the genetic ancestry PCs. Clocks are less accurate in admixed individuals, so differences between ancestries can be calibration error rather than biology. The phenotype associations with the ancestry PCs held are the 'ancestry held' version in biomarkers_phenotypes.tsv.",
                 "Accelerations and scores against the ancestry PCs (Spearman)", measure_names, ancestry_pcs, values, paste0(prefix, "_ancestry_mqc.json"))
}
