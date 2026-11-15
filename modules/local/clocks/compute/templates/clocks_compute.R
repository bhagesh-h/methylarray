#!/usr/bin/env Rscript
#
# Epigenetic clocks on the arrays that passed sample QC:
#   - dnaMethyAge (Wang et al. 2024) clocks, selected with --clocks
#   - DunedinPACE (Belsky et al. 2022) with PACEProjector
# Input: the detection-masked beta values of the preprocessing route, collapsed
# to CpGs by mean. The analysis probe filter is not applied, because each clock
# was trained on its own CpG set.
# Clock CpGs are matched through the EPIC v1 and 450K IDs of the probe annotation
# as well as the CpG ID: CpGs renamed on EPIC v2 are added under their old name.
# Missing values are set to the CpG's mean over the other arrays before the
# clocks are computed (dnaMethyAge's own filling fails when a single CpG is
# missing); CpGs missing in every array are left to dnaMethyAge, which uses its
# reference means. Each clock's coverage and filled values are reported.
# With reported ages: age acceleration (EAA, residual of clock on age) and, with
# cell fractions, intrinsic EAA (also adjusted for six of the seven EpiDISH
# centDHSbloodDMC fractions).
#

prefix <- "${prefix}"
se_dir <- "${se}"
retained_file <- "${retained}"
annotation_file <- "${annotation}"
samples_file <- "${samples}"
cells_file <- "${cell_fractions}"
dnamethyage_source <- "${dnamethyage_source}"
dunedinpace_source <- "${dunedinpace_source}"

defaults <- list(clocks = "HorvathS2013,HannumG2013,LevineM2018,HorvathS2018,DunedinPACE",
                 dunedinpace_min_probes = "0.8", cores = "${task.cpus}")
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
clocks <- strsplit(opts[["clocks"]], ",")[[1]]
dunedinpace_min <- as.numeric(opts[["dunedinpace_min_probes"]])
# dnaMethyAge forks one worker per core, each with a copy of the data: cap at 3
cores <- max(1L, min(3L, as.integer(opts[["cores"]])))

# dnaMethyAge and DunedinPACE are not on conda or CRAN: install the staged
# sources (pure R, pinned GitHub commits) into the task directory
Sys.setenv(HOME = getwd())
options(repos = character(0))
library_dir <- file.path(getwd(), "rlib")
dir.create(library_dir, showWarnings = FALSE)
.libPaths(c(library_dir, .libPaths()))
install.packages(c(dnamethyage_source, dunedinpace_source), lib = library_dir, repos = NULL, type = "source", quiet = TRUE)
suppressPackageStartupMessages({
    library(dnaMethyAge, lib.loc = library_dir)
    library(DunedinPACE, lib.loc = library_dir)
    library(SummarizedExperiment)
    library(HDF5Array)
    library(data.table)
})
set.seed(1)
unknown <- setdiff(clocks, c(suppressMessages(availableClock()), "DunedinPACE"))
if (length(unknown) > 0) stop("Unknown clocks in --clocks: ", paste(unknown, collapse = ", "))

# Detection-masked betas of the retained arrays, collapsed to CpGs by mean
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
annotation_header <- names(fread(text = readLines(gzfile(annotation_file), n = 2), nrows = 0))
annotation <- fread(text = readLines(gzfile(annotation_file)), na.strings = "",
                    select = intersect(c("cpg_id", "epicv1_id", "hm450_id", "gene", "snp_common_1pt"), annotation_header),
                    colClasses = "character")
for (column in setdiff(c("gene", "snp_common_1pt"), names(annotation))) annotation[, (column) := NA_character_]
# Per CpG: its first gene and whether any of its probes has a common SNP
# (minor allele frequency 1% or more) at the CpG or the extension base
cpg_info <- annotation[, .(gene = sub(";.*", "", na.omit(gene)[1]),
                           snp = if (all(is.na(snp_common_1pt))) NA else any(toupper(snp_common_1pt) %in% "TRUE")),
                       by = cpg_id]
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

# CpGs of each clock, from the packages' coefficient tables
coefficient_cpgs <- function(clock) {
    if (clock == "DunedinPACE") return(DunedinPACE:::mPACE_Models[["model_probes"]][[1]])
    env <- new.env()
    if (grepl("^PC[A-Z]", clock)) {
        data(list = "PC-clocks", package = "dnaMethyAge", envir = env)
    } else {
        data(list = clock, package = "dnaMethyAge", envir = env)
    }
    coefs <- get("coefs", env)
    ids <- if (is.data.frame(coefs) && "Probe" %in% names(coefs)) as.character(coefs[["Probe"]])
           else if (is.matrix(coefs) && !is.null(rownames(coefs))) rownames(coefs)
           else as.character(as.data.frame(coefs)[[1]])
    setdiff(unique(sub("_2\$", "", ids)), c("Intercept", "(Intercept)"))
}
to_cpg <- function(ids) ifelse(ids %in% aliases[["alias"]], aliases[["cpg_id"]][match(ids, aliases[["alias"]])], ids)
coverage <- rbindlist(lapply(clocks, function(clock) {
    ids <- coefficient_cpgs(clock)
    on_data <- ids %in% rownames(beta)
    data.table(clock = clock, cpgs = length(ids), present = sum(on_data),
               snp_cpgs = sum(cpg_info[["snp"]][match(to_cpg(ids), cpg_info[["cpg_id"]])] %in% TRUE),
               present_by_legacy_id = sum(ids %in% aliases[["alias"]]),
               share_present = round(mean(on_data), 4),
               filled_values = sum(is.na(beta[intersect(ids, rownames(beta)), , drop = FALSE])))
}))

# Fill missing values with the CpG mean over the other arrays; drop CpGs missing everywhere
beta <- beta[rowSums(!is.na(beta)) > 0, , drop = FALSE]
missing <- which(is.na(beta), arr.ind = TRUE)
if (nrow(missing) > 0) beta[missing] <- rowMeans(beta, na.rm = TRUE)[missing[, 1]]

# Standardised weights: each CpG's coefficient times the SD of its beta value in
# these arrays, in the clock's own units per SD (the clocks that are weighted
# sums of CpGs; principal-component clocks and epiTOC2 have no per-CpG weights)
clock_weights <- function(clock) {
    if (clock == "DunedinPACE") {
        w <- DunedinPACE:::mPACE_Models[["model_weights"]][[1]]
        return(data.table(cpg = names(w), weight = unname(w)))
    }
    if (grepl("^PC[A-Z]", clock)) return(NULL)
    env <- new.env()
    data(list = clock, package = "dnaMethyAge", envir = env)
    coefs <- get("coefs", env)
    if (!is.data.frame(coefs) || !all(c("Probe", "Coefficient") %in% names(coefs))) return(NULL)
    data.table(cpg = sub("_2\$", "", as.character(coefs[["Probe"]])), weight = as.numeric(coefs[["Coefficient"]]))[!cpg %in% c("Intercept", "(Intercept)")]
}
weights <- rbindlist(lapply(clocks, function(clock) {
    w <- clock_weights(clock)
    if (is.null(w) || nrow(w) == 0) return(NULL)
    w <- w[weight != 0]
    sds <- matrixStats::rowSds(beta[intersect(w[["cpg"]], rownames(beta)), , drop = FALSE])
    w[, sd_beta := unname(sds[match(cpg, names(sds))])]
    w[, standardised := weight * sd_beta]
    info <- cpg_info[match(to_cpg(w[["cpg"]]), cpg_id)]
    w[, `:=`(clock = clock, gene = info[["gene"]], snp = info[["snp"]])]
    setcolorder(w, c("clock", "cpg", "gene", "weight", "sd_beta", "standardised", "snp"))
    w[order(-abs(standardised))]
}), fill = TRUE)
if (nrow(weights) == 0) weights <- data.table(clock = character(), cpg = character(), gene = character(), weight = numeric(),
                                              sd_beta = numeric(), standardised = numeric(), snp = logical())
weights[, `:=`(weight = signif(weight, 6), sd_beta = signif(sd_beta, 4), standardised = signif(standardised, 4))]
fwrite(weights, paste0(prefix, "_weights.tsv"), sep = "\\t", na = "")

# Reported age and sex, and the cell fractions
samples <- read.delim(samples_file, colClasses = "character", na.strings = c("", "NA"), check.names = FALSE)
samples <- samples[match(colnames(beta), samples[["sample"]]), , drop = FALSE]
age <- suppressWarnings(as.numeric(samples[["reported_age"]]))
sex <- c(F = "Female", M = "Male")[samples[["reported_sex"]]]
age_info <- data.frame(Sample = colnames(beta), Age = age, Sex = unname(sex))
complete_age_sex <- !is.na(age_info[["Age"]]) & !is.na(age_info[["Sex"]])

result <- data.table(sample = colnames(beta), age = age)
notes <- setNames(rep("", length(clocks)), clocks)
# A clock that fails gives NA with the error in the coverage notes, so that one
# clock does not stop the others (dnaMethyAge 0.2.0, for example, fails for
# BernabeuE2023c when an array's linear estimate is below 20 years)
run_clock <- function(clock, ...) {
    tryCatch(suppressMessages(suppressWarnings(methyAge(..., clock = clock, do_plot = FALSE, use_cores = cores))),
             error = function(e) { notes[[clock]] <<- paste("failed:", trimws(conditionMessage(e))); NULL })
}
for (clock in clocks) {
    values <- rep(NA_real_, ncol(beta))
    if (clock == "DunedinPACE") {
        pace <- tryCatch(PACEProjector(beta, proportionOfProbesRequired = dunedinpace_min)[["DunedinPACE"]],
                         error = function(e) { notes[[clock]] <<- conditionMessage(e); NULL })
        if (!is.null(pace)) values <- unname(pace[colnames(beta)])
    } else if (clock == "PCGrimAge") {
        if (any(complete_age_sex)) {
            fit <- run_clock(clock, beta[, complete_age_sex, drop = FALSE], age_info = age_info[complete_age_sex, ])
            if (!is.null(fit)) values <- fit[["mAge"]][match(colnames(beta), fit[["Sample"]])]
        }
        if (!all(complete_age_sex) && !nzchar(notes[[clock]])) notes[[clock]] <- "needs reported age and sex; NA for arrays without them"
    } else {
        fit <- run_clock(clock, beta)
        if (!is.null(fit)) values <- fit[["mAge"]][match(colnames(beta), fit[["Sample"]])]
    }
    result[, (clock) := round(values, 5)]
    message(clock, ": ", if (nzchar(notes[[clock]])) notes[[clock]] else "done")
}
coverage[, note := ifelse(nzchar(notes[clock]), notes[clock], NA_character_)]

# Age acceleration for the clocks that estimate age, when ages are given
age_clocks <- intersect(clocks, c("HorvathS2013", "HannumG2013", "LevineM2018", "HorvathS2018", "ZhangQ2019",
                                  "BernabeuE2023c", "PCGrimAge", "PCPhenoAge", "PCHorvathS2013", "PCHannumG2013", "PCHorvathS2018"))
cells <- NULL
if (nzchar(cells_file)) {
    cells <- fread(cells_file)
    cells <- cells[match(colnames(beta), sample)]
}
if (sum(!is.na(age)) >= 6) {
    for (clock in age_clocks) {
        y <- result[[clock]]
        if (sum(!is.na(y) & !is.na(age)) < 6) next
        result[, (paste0("eaa_", clock)) := round(residuals(lm(y ~ age, na.action = na.exclude)), 5)]
        if (!is.null(cells)) {
            # six of the seven fractions: they sum to about one, so one is left out
            cell_columns <- setdiff(grep("^dhs7_", names(cells), value = TRUE), "dhs7_Neutro")
            design <- data.frame(y = y, age = age, cells[, ..cell_columns])
            # at least five residual degrees of freedom, so that the residuals are not forced to zero
            if (sum(complete.cases(design)) - (ncol(design) + 1) >= 5) {
                result[, (paste0("ieaa_", clock)) := round(residuals(lm(y ~ ., data = design, na.action = na.exclude)), 5)]
            }
        }
    }
}
setorder(result, sample)
fwrite(result, paste0(prefix, ".tsv"), sep = "\\t", na = "")
fwrite(coverage, paste0(prefix, "_coverage.tsv"), sep = "\\t", na = "")

# MultiQC: clock coverage on the array
writeLines(c(
    "# id: 'clock_coverage'",
    "# section_name: 'Epigenetic clocks'",
    "# description: 'CpGs of each clock in the data (EPIC v2 does not measure every CpG of clocks trained on 450K or EPIC v1). Missing CpGs take the reference means of dnaMethyAge; values that failed detection take the CpG mean over the other arrays. SNP CpGs: clock CpGs with a common SNP at the CpG or extension base.'",
    "# plot_type: 'table'"
), paste0(prefix, "_coverage_mqc.tsv"))
write.table(data.frame(Clock = coverage[["clock"]], CpGs = coverage[["cpgs"]], Present = coverage[["present"]],
                       `Found by old ID` = coverage[["present_by_legacy_id"]], `Share present` = coverage[["share_present"]], `SNP CpGs` = coverage[["snp_cpgs"]],
                       `Filled values` = coverage[["filled_values"]], check.names = FALSE),
            paste0(prefix, "_coverage_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE)

# Software versions
writeLines(c('"${task.process}":',
    paste0("    dnaMethyAge: ", packageVersion("dnaMethyAge", lib.loc = library_dir)),
    paste0("    DunedinPACE: ", packageVersion("DunedinPACE", lib.loc = library_dir)),
    paste0("    preprocessCore: ", packageVersion("preprocessCore")),
    paste0("    impute: ", packageVersion("impute")),
    paste0("    RPMM: ", packageVersion("RPMM"))), "versions.yml")
