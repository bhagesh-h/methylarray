#!/usr/bin/env Rscript
#
# Gather the per-sample SeSAMe results into one HDF5-backed
# SummarizedExperiment (probes x samples) with the sample metadata as colData.
#

prefix <- "${prefix}"
samples_file <- "${samples}"

suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(HDF5Array)
})

# Sample metadata, one row per array, sorted by sample name for a stable column order
samples <- read.delim(samples_file, colClasses = "character", na.strings = c("", "NA"), check.names = FALSE)
samples <- samples[order(samples[["sample"]]), , drop = FALSE]
rownames(samples) <- samples[["sample"]]

rds_files <- list.files("rds", pattern = "[.]sesame[.]rds\$", full.names = TRUE)
results <- lapply(rds_files, readRDS)
names(results) <- vapply(results, function(r) r[["sample"]], character(1))

missing <- setdiff(samples[["sample"]], names(results))
extra <- setdiff(names(results), samples[["sample"]])
if (length(missing) || length(extra)) {
    stop("Sample sheet and SeSAMe results differ. Missing results: ", paste(missing, collapse = ", "),
         "; results without metadata: ", paste(extra, collapse = ", "))
}
results <- results[samples[["sample"]]]

platforms <- unique(vapply(results, function(r) r[["platform"]], character(1)))
if (length(platforms) != 1) {
    stop("All arrays must come from one platform; found: ", paste(platforms, collapse = ", "))
}
probes <- results[[1]][["probes"]]
same_probes <- vapply(results, function(r) identical(r[["probes"]], probes), logical(1))
if (!all(same_probes)) {
    stop("Probe sets differ between samples: ", paste(names(results)[!same_probes], collapse = ", "))
}

# One assay at a time into a scratch HDF5 file, so that only one
# probes x samples matrix is in memory
assay_names <- c("beta", "beta_unmasked", "detp", "meth", "unmeth")
scratch_h5 <- "scratch_assays.h5"
assays <- list()
for (name in assay_names) {
    values <- vapply(results, function(r) r[[name]], numeric(length(probes)))
    dimnames(values) <- list(probes, names(results))
    assays[[name]] <- writeHDF5Array(values, filepath = scratch_h5, name = name)
    rm(values)
    invisible(gc())
}
rm(results)
invisible(gc())

# Sample QC statistics joined to the sample metadata
qc_files <- list.files("qc", pattern = "[.]qc[.]tsv\$", full.names = TRUE)
qc <- do.call(rbind, lapply(qc_files, read.delim, check.names = FALSE))
rownames(qc) <- qc[["sample"]]
qc <- qc[samples[["sample"]], , drop = FALSE]
col_data <- cbind(samples, qc[, setdiff(colnames(qc), "sample"), drop = FALSE])

row_data <- DataFrame(probe_id = probes, cpg_id = sub("_.*", "", probes), row.names = probes)
se <- SummarizedExperiment(
    assays = assays,
    colData = DataFrame(col_data, check.names = FALSE),
    rowData = row_data,
    metadata = list(
        platform = platforms,
        preprocessing = "sesame",
        prep = unique(qc[["prep"]]),
        sesame_version = as.character(packageVersion("sesame"))
    )
)
# Self-contained directory (assays.h5 + se.rds, relative paths);
# load it with HDF5Array::loadHDF5SummarizedExperiment()
saveHDF5SummarizedExperiment(se, dir = paste0(prefix, "_se"), replace = TRUE)
unlink(scratch_h5)

write.table(col_data, paste0(prefix, "_qc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE)

# MultiQC general statistics: one row per sample
mqc <- data.frame(
    Sample = qc[["sample"]],
    sesame_frac_detected = round(100 * qc[["frac_dt"]], 2),
    sesame_frac_masked = round(100 * qc[["frac_masked"]], 2),
    sesame_mean_intensity = round(qc[["mean_intensity"]], 0),
    sesame_rg_ratio = round(qc[["RGratio"]], 3),
    sesame_predicted_sex = qc[["predicted_sex"]],
    check.names = FALSE
)
header <- c(
    "# id: 'sesame_qc'",
    "# section_name: 'SeSAMe sample QC'",
    "# plot_type: 'generalstats'",
    "# pconfig:",
    "#     sesame_frac_detected:",
    "#         title: 'Detected (%)'",
    "#         description: 'Probes with a pOOBAH detection p-value below the threshold'",
    "#         min: 0",
    "#         max: 100",
    "#         suffix: '%'",
    "#     sesame_frac_masked:",
    "#         title: 'Masked (%)'",
    "#         description: 'Probes set to NA by the SeSAMe quality and detection masks'",
    "#         min: 0",
    "#         max: 100",
    "#         suffix: '%'",
    "#     sesame_mean_intensity:",
    "#         title: 'Mean intensity'",
    "#         description: 'Mean in-band signal intensity'",
    "#         format: '{:,.0f}'",
    "#     sesame_rg_ratio:",
    "#         title: 'R/G ratio'",
    "#         description: 'Ratio of red to green in-band median intensity'",
    "#     sesame_predicted_sex:",
    "#         title: 'Predicted sex'",
    "#         description: 'Sex predicted by SeSAMe inferSex'"
)
writeLines(header, paste0(prefix, "_qc_mqc.tsv"))
write.table(mqc, paste0(prefix, "_qc_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE)

# Software versions
writeLines(c('"${task.process}":', paste0("    SummarizedExperiment: ", packageVersion("SummarizedExperiment")), paste0("    HDF5Array: ", packageVersion("HDF5Array"))), "versions.yml")
