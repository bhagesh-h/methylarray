#!/usr/bin/env Rscript
#
# minfi-based preprocessing of all arrays into one HDF5-backed
# SummarizedExperiment (probes x samples), for three routes:
#   minfi:    Noob (default), Funnorm, quantile or SWAN; minfi detection p-values
#   enmix:    ENmix oob background + RELIC dye-bias correction (+ optional
#             quantile normalisation) + RCP probe-type correction
#   combined: the ENmix route with the SeSAMe quality and pOOBAH masks
# Values that fail detection or the bead count are set to NA in 'beta';
# 'beta_unmasked' keeps them. No array or probe is removed.
#

prefix <- "${prefix}"
samples_file <- "${samples}"
sesame_se <- "${sesame_se}"

defaults <- list(method = "minfi", minfi_norm = "noob", enmix_quantile = "false",
                 min_beads = "3", detection_pval = "0.05", platform = "EPICv2")
parse_args <- function(x, defaults) {
    tokens <- strsplit(trimws(x), "[[:space:]]+")[[1]]
    tokens <- tokens[nzchar(tokens)]
    if (length(tokens) %% 2 != 0) stop("ext.args must be pairs of '--name value': ", x)
    opts <- defaults
    for (i in seq(1, length(tokens), by = 2)) {
        name <- sub("^--", "", tokens[i])
        if (!name %in% names(defaults)) stop("Unknown option in ext.args: ", tokens[i])
        opts[[name]] <- tokens[i + 1]
    }
    opts
}
opts <- parse_args("${task.ext.args ?: ''}", defaults)
method <- match.arg(opts[["method"]], c("minfi", "enmix", "combined"))
minfi_norm <- match.arg(opts[["minfi_norm"]], c("noob", "funnorm", "quantile", "swan"))
enmix_quantile <- tolower(opts[["enmix_quantile"]]) == "true"
min_beads <- as.integer(opts[["min_beads"]])
detection_pval <- as.numeric(opts[["detection_pval"]])
arrays <- c(EPICv2 = "IlluminaHumanMethylationEPICv2")
if (!opts[["platform"]] %in% names(arrays)) stop("Unsupported platform: ", opts[["platform"]])
if (method == "combined" && !nzchar(sesame_se)) stop("The combined route needs the SeSAMe SummarizedExperiment")

Sys.setenv(HOME = getwd())
suppressPackageStartupMessages({
    library(minfi)
    library(ENmix)
    library(SummarizedExperiment)
    library(HDF5Array)
})
set.seed(1)

# Arrays, sorted by sample name for a stable column order
samples <- read.delim(samples_file, colClasses = "character", na.strings = c("", "NA"), check.names = FALSE)
samples <- samples[order(samples[["sample"]]), , drop = FALSE]
rownames(samples) <- samples[["sample"]]
targets <- data.frame(Basename = file.path("idat", samples[["idat_prefix"]]), Sample_Name = samples[["sample"]])

rg <- read.metharray.exp(targets = targets, extended = TRUE)
sampleNames(rg) <- samples[["sample"]]
if (annotation(rg)[["array"]] != arrays[[opts[["platform"]]]]) {
    stop("The IDAT files are ", annotation(rg)[["array"]], " arrays, but the run is set to --platform ", opts[["platform"]])
}

# Detection p-values and bead counts per probe (ENmix QCinfo; its sample and
# CpG flags are reported, not used to drop data). distplot = TRUE because
# QCinfo 1.46.2 fails with outlier = FALSE and distplot = FALSE ('mdat' not found);
# the plots stay in the task directory.
qci <- QCinfo(rg, detPthre = detection_pval, nbthre = min_beads, outlier = FALSE, distplot = TRUE)

noob <- preprocessNoob(rg)
probes <- rownames(noob)

# Put any probes x samples matrix on the common probe order, NA where absent
align <- function(m) {
    out <- matrix(NA_real_, nrow = length(probes), ncol = nrow(samples), dimnames = list(probes, samples[["sample"]]))
    common <- intersect(probes, rownames(m))
    out[common, ] <- as.matrix(m[common, samples[["sample"]], drop = FALSE])
    out
}

nbead <- align(qci[["nbead"]])
if (method == "minfi") {
    beta_unmasked <- switch(minfi_norm,
        noob = getBeta(noob),
        funnorm = getBeta(preprocessFunnorm(rg)),
        quantile = getBeta(preprocessQuantile(rg)),
        swan = getBeta(preprocessSWAN(rg, mSet = noob))
    )
    beta_unmasked <- align(beta_unmasked)
    meth <- align(getMeth(noob))
    unmeth <- align(getUnmeth(noob))
    detp <- align(detectionP(rg))
    normalisation <- minfi_norm
} else {
    mdat <- preprocessENmix(rg, bgParaEst = "oob", dyeCorr = "RELIC", nCores = 1)
    if (enmix_quantile) mdat <- norm.quantile(mdat, method = "quantile1")
    beta_unmasked <- align(rcp(mdat))
    meth <- align(getMeth(mdat))
    unmeth <- align(getUnmeth(mdat))
    detp <- align(qci[["detP"]])
    normalisation <- paste0("ENmix oob + RELIC", if (enmix_quantile) " + quantile" else "", " + RCP")
}

# Masks: detection and bead count, or for the combined route the SeSAMe masks
# (quality mask, pOOBAH, bead count) plus the bead count of this import
if (method == "combined") {
    sesame <- loadHDF5SummarizedExperiment(sesame_se)
    missing <- setdiff(samples[["sample"]], colnames(sesame))
    if (length(missing)) stop("Arrays missing from the SeSAMe result: ", paste(missing, collapse = ", "))
    sesame_beta <- align(assay(sesame, "beta")[intersect(probes, rownames(sesame)), samples[["sample"]], drop = FALSE])
    detp <- align(assay(sesame, "detp")[intersect(probes, rownames(sesame)), samples[["sample"]], drop = FALSE])
    mask <- is.na(sesame_beta) | is.na(nbead) | nbead < min_beads
    rm(sesame_beta)
} else {
    mask <- is.na(detp) | detp > detection_pval | is.na(nbead) | nbead < min_beads
}
beta <- beta_unmasked
beta[mask] <- NA

# Per-array QC
raw_qc <- getQC(preprocessRaw(rg))
# Sex with the getSex rule (median chrY minus median chrX copy number below -2
# means female). getSex itself also runs a k-means consistency check that only
# warns, and fails when there are very few arrays, so the rule is applied here.
genome <- mapToGenome(noob)
copy_number <- getCN(genome)
chromosome <- as.character(seqnames(genome))
sex <- data.frame(
    xMed = matrixStats::colMedians(copy_number[chromosome == "chrX", , drop = FALSE], na.rm = TRUE),
    yMed = matrixStats::colMedians(copy_number[chromosome == "chrY", , drop = FALSE], na.rm = TRUE),
    row.names = colnames(copy_number)
)
sex[["predictedSex"]] <- ifelse(sex[["yMed"]] - sex[["xMed"]] < -2, "F", "M")
rm(genome, copy_number)
qc <- data.frame(
    sample = samples[["sample"]],
    platform = opts[["platform"]],
    method = method,
    normalisation = normalisation,
    predicted_sex = ifelse(sex[samples[["sample"]], "predictedSex"] == "F", "female", "male"),
    x_med = sex[samples[["sample"]], "xMed"],
    y_med = sex[samples[["sample"]], "yMed"],
    median_log2_meth = raw_qc[samples[["sample"]], "mMed"],
    median_log2_unmeth = raw_qc[samples[["sample"]], "uMed"],
    frac_detected = colMeans(detp <= detection_pval, na.rm = TRUE),
    num_low_beads = colSums(is.na(nbead) | nbead < min_beads),
    num_masked = colSums(mask),
    frac_masked = colMeans(mask),
    bisulfite_intensity = qci[["bisul"]][samples[["sample"]]],
    enmix_low_quality = samples[["sample"]] %in% qci[["badsample"]],
    check.names = FALSE
)

# One assay at a time into a scratch HDF5 file
scratch_h5 <- "scratch_assays.h5"
assays <- list()
for (name in c("beta", "beta_unmasked", "detp", "meth", "unmeth")) {
    assays[[name]] <- writeHDF5Array(get(name), filepath = scratch_h5, name = name)
}
rm(beta, beta_unmasked, detp, meth, unmeth, mask, nbead, rg, noob)
invisible(gc())

col_data <- cbind(samples[, setdiff(colnames(samples), "idat_prefix"), drop = FALSE], qc[, setdiff(colnames(qc), "sample"), drop = FALSE])
se <- SummarizedExperiment(
    assays = assays,
    colData = DataFrame(col_data, check.names = FALSE),
    rowData = DataFrame(probe_id = probes, cpg_id = sub("_.*", "", probes), row.names = probes),
    metadata = list(
        platform = opts[["platform"]],
        preprocessing = method,
        normalisation = normalisation,
        detection = if (method == "combined") "SeSAMe masks (quality mask, pOOBAH)" else if (method == "minfi") "minfi detectionP" else "ENmix QCinfo",
        minfi_version = as.character(packageVersion("minfi")),
        enmix_version = as.character(packageVersion("ENmix"))
    )
)
saveHDF5SummarizedExperiment(se, dir = paste0(prefix, "_se"), replace = TRUE)
unlink(scratch_h5)

write.table(col_data, paste0(prefix, "_qc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE)

# MultiQC general statistics: one row per array
mqc <- data.frame(
    Sample = qc[["sample"]],
    frac_detected = round(100 * qc[["frac_detected"]], 2),
    frac_masked = round(100 * qc[["frac_masked"]], 2),
    median_log2_meth = round(qc[["median_log2_meth"]], 2),
    predicted_sex = qc[["predicted_sex"]],
    check.names = FALSE
)
names(mqc)[-1] <- paste0(method, "_", names(mqc)[-1])
header <- c(
    paste0("# id: '", method, "_qc'"),
    paste0("# section_name: '", method, " sample QC'"),
    "# plot_type: 'generalstats'",
    "# pconfig:",
    paste0("#     ", method, "_frac_detected:"),
    "#         title: 'Detected (%)'",
    "#         description: 'Probes with a detection p-value below the threshold'",
    "#         min: 0",
    "#         max: 100",
    "#         suffix: '%'",
    paste0("#     ", method, "_frac_masked:"),
    "#         title: 'Masked (%)'",
    "#         description: 'Values set to NA by the detection and bead-count masks'",
    "#         min: 0",
    "#         max: 100",
    "#         suffix: '%'",
    paste0("#     ", method, "_median_log2_meth:"),
    "#         title: 'Median log2 M'",
    "#         description: 'Median log2 methylated intensity (minfi getQC)'",
    paste0("#     ", method, "_predicted_sex:"),
    "#         title: 'Predicted sex'",
    "#         description: 'Sex predicted with the minfi getSex rule'"
)
writeLines(header, paste0(prefix, "_qc_mqc.tsv"))
write.table(mqc, paste0(prefix, "_qc_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE)

# Software versions
writeLines(c('"${task.process}":',
    paste0("    minfi: ", packageVersion("minfi")),
    paste0("    ENmix: ", packageVersion("ENmix")),
    paste0("    SummarizedExperiment: ", packageVersion("SummarizedExperiment")),
    paste0("    HDF5Array: ", packageVersion("HDF5Array"))), "versions.yml")
