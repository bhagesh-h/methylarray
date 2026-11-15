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
# Memory and time: the arrays are read and corrected in chunks of --chunk_size
# arrays, on parallel worker processes when the task has several CPUs
# (at least 2; Noob, detection p-values, ENmix oob/RELIC/RCP and SWAN work array by array),
# and each chunk is written to HDF5 at once. Funnorm, quantile normalisation and
# the ENmix quantile step need all arrays together and run in a second pass.
# The ENmix low-quality flag (QCinfo rule: more than 5% of values failing
# detection or bead count, or bisulfite intensity below mean - 3 SD) is computed
# over all arrays.
#

prefix <- "${prefix}"
samples_file <- "${samples}"
sesame_se <- "${sesame_se}"

defaults <- list(method = "minfi", minfi_norm = "noob", enmix_quantile = "false",
                 min_beads = "3", detection_pval = "0.05", platform = "EPICv2", chunk_size = "40")
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
chunk_size <- max(2L, as.integer(opts[["chunk_size"]]))
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

# Arrays, sorted by sample name for a stable column order, in chunks of about
# equal size and at least --chunk_size arrays (QCinfo needs two or more)
samples <- read.delim(samples_file, colClasses = "character", na.strings = c("", "NA"), check.names = FALSE)
samples <- samples[order(samples[["sample"]]), , drop = FALSE]
rownames(samples) <- samples[["sample"]]
task_dir <- getwd()
targets <- data.frame(Basename = file.path(task_dir, "idat", samples[["idat_prefix"]]), Sample_Name = samples[["sample"]])
n_chunks <- max(1, floor(nrow(samples) / chunk_size))
chunks <- split(seq_len(nrow(samples)), (seq_len(nrow(samples)) - 1) %/% ceiling(nrow(samples) / n_chunks))
between_arrays <- (method == "minfi" && minfi_norm %in% c("funnorm", "quantile")) || (method != "minfi" && enmix_quantile)

read_arrays <- function(idx, extended) {
    rg <- read.metharray.exp(targets = targets[idx, , drop = FALSE], extended = extended)
    sampleNames(rg) <- samples[["sample"]][idx]
    if (annotation(rg)[["array"]] != arrays[[opts[["platform"]]]]) {
        stop("The IDAT files are ", annotation(rg)[["array"]], " arrays, but the run is set to --platform ", opts[["platform"]])
    }
    rg
}

# Put any probes x arrays matrix on the common probe order, NA where absent
probes <- NULL
align <- function(m, columns) {
    out <- matrix(NA_real_, nrow = length(probes), ncol = length(columns), dimnames = list(probes, columns))
    common <- intersect(probes, rownames(m))
    out[common, ] <- as.matrix(m[common, columns, drop = FALSE])
    out
}

# One chunk: read, correct and write its matrices to its own scratch HDF5 file,
# in its own directory (QCinfo writes its plots to the working directory).
# Returns the per-array QC values and the names of the datasets written.
process_chunk <- function(k) {
    idx <- chunks[[k]]
    columns <- samples[["sample"]][idx]
    chunk_dir <- file.path(task_dir, sprintf("chunk_%03d", k))
    dir.create(chunk_dir, showWarnings = FALSE)
    owd <- setwd(chunk_dir)
    on.exit(setwd(owd))
    h5 <- file.path(chunk_dir, "assays.h5")
    written <- character()
    store <- function(name, m) {
        writeHDF5Array(m, filepath = h5, name = name)
        written <<- c(written, name)
    }
    message(sprintf("Chunk %d of %d: %d arrays", k, length(chunks), length(idx)))
    rg <- read_arrays(idx, extended = TRUE)

    # Detection p-values and bead counts per probe (ENmix QCinfo; no array or
    # CpG is flagged here). distplot = TRUE because QCinfo 1.46.2 fails with
    # outlier = FALSE and distplot = FALSE ('mdat' not found); the plots stay
    # in the chunk directory.
    qci <- QCinfo(rg, detPthre = detection_pval, nbthre = min_beads, samplethre = 1, CpGthre = 1, bisulthre = 0,
                  outlier = FALSE, distplot = TRUE)
    noob <- preprocessNoob(rg)
    if (is.null(probes)) {
        probes <<- rownames(noob)
        array_annotation <<- annotation(noob)
    }
    store("nbead", align(qci[["nbead"]], columns))
    if (method == "minfi") {
        store("detp", align(detectionP(rg), columns))
        store("meth", align(getMeth(noob), columns))
        store("unmeth", align(getUnmeth(noob), columns))
        if (minfi_norm == "noob") store("beta_unmasked", align(getBeta(noob), columns))
        if (minfi_norm == "swan") store("beta_unmasked", align(getBeta(preprocessSWAN(rg, mSet = noob)), columns))
    } else {
        store("detp", align(qci[["detP"]], columns))
        mdat <- preprocessENmix(rg, bgParaEst = "oob", dyeCorr = "RELIC", nCores = 1)
        store("meth", align(getMeth(mdat), columns))
        store("unmeth", align(getUnmeth(mdat), columns))
        if (!enmix_quantile) store("beta_unmasked", align(rcp(mdat), columns))
        rm(mdat)
    }

    # Per-array QC. Sex with the getSex rule (median chrY minus median chrX
    # copy number below -2 means female). getSex itself also runs a k-means
    # consistency check that only warns, and fails when there are very few
    # arrays, so the rule is applied here.
    raw_qc <- getQC(preprocessRaw(rg))
    genome <- mapToGenome(noob)
    copy_number <- getCN(genome)
    chromosome <- as.character(seqnames(genome))
    bad_values <- qci[["nbead"]] < min_beads | qci[["detP"]] > detection_pval
    qc <- data.frame(
        sample = columns,
        x_med = matrixStats::colMedians(copy_number[chromosome == "chrX", , drop = FALSE], na.rm = TRUE),
        y_med = matrixStats::colMedians(copy_number[chromosome == "chrY", , drop = FALSE], na.rm = TRUE),
        median_log2_meth = raw_qc[columns, "mMed"],
        median_log2_unmeth = raw_qc[columns, "uMed"],
        bisulfite_intensity = qci[["bisul"]][columns],
        enmix_bad_fraction = colSums(bad_values[, columns, drop = FALSE], na.rm = TRUE) / nrow(bad_values)
    )
    rm(rg, qci, noob, raw_qc, genome, copy_number, chromosome, bad_values)
    invisible(gc())
    list(qc = qc, h5 = h5, datasets = written)
}

# The first chunk runs here and sets the probe order; the others run in
# parallel worker processes, as many as the CPUs and the memory allow
# (about 0.25 GB per array of a chunk, plus 2 GB)
results <- list(process_chunk(1))
if (length(chunks) > 1) {
    memory_gb <- ${task.memory ? task.memory.toGiga() : 8}
    per_chunk_gb <- 0.25 * max(lengths(chunks)) + 2
    workers <- max(1L, min(${task.cpus}, length(chunks) - 1, floor(memory_gb / per_chunk_gb)))
    message("Remaining chunks on ", workers, " worker(s)")
    rest <- parallel::mclapply(seq_along(chunks)[-1], process_chunk, mc.cores = workers, mc.preschedule = FALSE)
    failed <- vapply(rest, function(r) inherits(r, "try-error") || is.null(r), logical(1))
    if (any(failed)) stop("Chunk processing failed: ", paste(unique(vapply(rest[failed], function(r) if (is.null(r)) "worker ended without a result (out of memory?)" else as.character(r), character(1))), collapse = "; "))
    results <- c(results, rest)
}
qc_parts <- lapply(results, `[[`, "qc")

# Chunks bound by column at the end; arrays computed across all chunks go to a
# scratch file of their own
parts <- list()
for (k in seq_along(results)) {
    for (name in results[[k]][["datasets"]]) parts[[name]][[k]] <- HDF5Array(results[[k]][["h5"]], name)
}
scratch_h5 <- file.path(task_dir, "scratch_assays.h5")
store <- function(name, k, m) {
    parts[[name]][[k]] <<- writeHDF5Array(m, filepath = scratch_h5, name = paste0(name, "_", k))
}
bound <- function(name) do.call(cbind, parts[[name]])


# Normalisations across all arrays
if (between_arrays) {
    message("Normalisation across all arrays")
    if (method == "minfi") {
        rg <- read_arrays(seq_len(nrow(samples)), extended = FALSE)
        beta_all <- switch(minfi_norm, funnorm = getBeta(preprocessFunnorm(rg)), quantile = getBeta(preprocessQuantile(rg)))
        rm(rg)
    } else {
        mset <- MethylSet(Meth = as.matrix(bound("meth")), Unmeth = as.matrix(bound("unmeth")), annotation = array_annotation)
        beta_all <- rcp(norm.quantile(mset, method = "quantile1"))
        rm(mset)
    }
    parts[["beta_unmasked"]] <- list(writeHDF5Array(align(beta_all, samples[["sample"]]), filepath = scratch_h5, name = "beta_unmasked_all"))
    rm(beta_all)
    invisible(gc())
}
normalisation <- if (method == "minfi") minfi_norm else paste0("ENmix oob + RELIC", if (enmix_quantile) " + quantile" else "", " + RCP")

# Masks: detection and bead count, or for the combined route the SeSAMe masks
# (quality mask, pOOBAH, bead count) plus the bead count of this import
if (method == "combined") {
    sesame <- loadHDF5SummarizedExperiment(sesame_se)
    missing <- setdiff(samples[["sample"]], colnames(sesame))
    if (length(missing)) stop("Arrays missing from the SeSAMe result: ", paste(missing, collapse = ", "))
    sesame_probes <- intersect(probes, rownames(sesame))
}
beta_unmasked_all <- bound("beta_unmasked")
detp_all <- bound("detp")
nbead_all <- bound("nbead")
mask_parts <- list()
for (k in seq_along(chunks)) {
    columns <- samples[["sample"]][chunks[[k]]]
    beta <- as.matrix(beta_unmasked_all[, columns, drop = FALSE])
    nbead <- as.matrix(nbead_all[, columns, drop = FALSE])
    low_beads <- is.na(nbead) | nbead < min_beads
    if (method == "combined") {
        detp <- align(as.matrix(assay(sesame, "detp")[sesame_probes, columns, drop = FALSE]), columns)
        store("detp_sesame", k, detp)
        mask <- is.na(align(as.matrix(assay(sesame, "beta")[sesame_probes, columns, drop = FALSE]), columns)) | low_beads
    } else {
        detp <- as.matrix(detp_all[, columns, drop = FALSE])
        mask <- is.na(detp) | detp > detection_pval | low_beads
    }
    beta[mask] <- NA
    store("beta", k, beta)
    mask_parts[[k]] <- data.frame(
        sample = columns,
        frac_detected = colMeans(detp <= detection_pval, na.rm = TRUE),
        num_low_beads = colSums(low_beads),
        num_masked = colSums(mask),
        frac_masked = colMeans(mask)
    )
    rm(beta, nbead, low_beads, detp, mask)
    invisible(gc())
}
if (method == "combined") parts[["detp"]] <- parts[["detp_sesame"]]

qc_values <- merge(do.call(rbind, qc_parts), do.call(rbind, mask_parts), by = "sample")
qc_values <- qc_values[match(samples[["sample"]], qc_values[["sample"]]), ]
bisul <- qc_values[["bisulfite_intensity"]]
bisul_threshold <- if (length(bisul) > 1) mean(bisul, na.rm = TRUE) - 3 * sd(bisul, na.rm = TRUE) else -Inf
qc <- data.frame(
    sample = samples[["sample"]],
    platform = opts[["platform"]],
    method = method,
    normalisation = normalisation,
    predicted_sex = ifelse(qc_values[["y_med"]] - qc_values[["x_med"]] < -2, "female", "male"),
    x_med = qc_values[["x_med"]],
    y_med = qc_values[["y_med"]],
    median_log2_meth = qc_values[["median_log2_meth"]],
    median_log2_unmeth = qc_values[["median_log2_unmeth"]],
    frac_detected = qc_values[["frac_detected"]],
    num_low_beads = qc_values[["num_low_beads"]],
    num_masked = qc_values[["num_masked"]],
    frac_masked = qc_values[["frac_masked"]],
    bisulfite_intensity = bisul,
    enmix_low_quality = qc_values[["enmix_bad_fraction"]] > 0.05 | (!is.na(bisul) & bisul < bisul_threshold),
    check.names = FALSE
)

# The chunks bound by column, written block by block into the final directory
assays <- lapply(c(beta = "beta", beta_unmasked = "beta_unmasked", detp = "detp", meth = "meth", unmeth = "unmeth"), bound)
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
unlink(c(scratch_h5, file.path(task_dir, sprintf("chunk_%03d", seq_along(chunks)))), recursive = TRUE)


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
