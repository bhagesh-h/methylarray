#!/usr/bin/env Rscript
#
# GEO submission files for the arrays that passed sample QC, following the GEO
# rules for methylation arrays:
#   - processed matrix: ID_REF, then per array the normalised beta value of the
#     preprocessing route (no detection masking) and its detection p-value
#   - signal intensities: ID_REF, then per array the unnormalised methylated and
#     unmethylated signal (minfi preprocessRaw on the IDAT files) and the
#     detection p-value
#   - samples: array name, IDAT file names and Sentrix fields, as a start for
#     the GEO metadata sheet (no phenotype columns)
# Columns are named after the samplesheet sample names; rename them before
# submission if those names identify participants. Rows are written in blocks,
# so only one block of each matrix is in memory.
#

prefix <- "${prefix}"
samples_file <- "${samples}"
se_dir <- "${se}"
retained_file <- "${retained}"

defaults <- list(block_size = "100000")
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
block_size <- as.integer(opts[["block_size"]])

Sys.setenv(HOME = getwd())
suppressPackageStartupMessages({
    library(minfi)
    library(SummarizedExperiment)
    library(HDF5Array)
    library(data.table)
})
writeLines(c('"${task.process}":',
    paste0("    minfi: ", packageVersion("minfi")),
    paste0("    data.table: ", packageVersion("data.table"))), "versions.yml")

# Arrays that passed sample QC, sorted by sample name
samples <- read.delim(samples_file, colClasses = "character", na.strings = c("", "NA"), check.names = FALSE)
samples <- samples[samples[["sample"]] %in% readLines(retained_file), , drop = FALSE]
samples <- samples[order(samples[["sample"]]), , drop = FALSE]
if (nrow(samples) == 0) stop("No array passed sample QC")
arrays <- samples[["sample"]]

se <- loadHDF5SummarizedExperiment(se_dir)
missing <- setdiff(arrays, colnames(se))
if (length(missing) > 0) stop("Arrays missing from the preprocessing result: ", paste(missing, collapse = ", "))
se <- se[, arrays]
probes <- rownames(se)

# Unnormalised intensities from the IDAT files
rg <- read.metharray.exp(targets = data.frame(Basename = file.path("idat", samples[["idat_prefix"]])))
sampleNames(rg) <- arrays
raw <- preprocessRaw(rg)
rm(rg)
invisible(gc())
common <- intersect(probes, rownames(raw))
if (length(common) < 0.9 * length(probes)) {
    stop("Only ", length(common), " of ", length(probes), " probes of the preprocessing result are in the raw data")
}
meth <- getMeth(raw)[common, arrays, drop = FALSE]
unmeth <- getUnmeth(raw)[common, arrays, drop = FALSE]
rm(raw)
invisible(gc())

interleave <- function(blocks, suffixes) {
    out <- do.call(cbind, blocks)
    order <- as.vector(t(matrix(seq_len(ncol(out)), ncol = length(blocks))))
    out <- out[, order, drop = FALSE]
    colnames(out) <- paste0(rep(arrays, each = length(blocks)), suffixes)
    out
}
write_block <- function(id, values, file, first) {
    fwrite(data.table(ID_REF = id, values), file, sep = "\\t", na = "", append = !first, col.names = first, compress = "gzip")
}

processed_file <- paste0(prefix, "_processed_matrix.tsv.gz")
intensity_file <- paste0(prefix, "_signal_intensities.tsv.gz")
unlink(c(processed_file, intensity_file))
starts <- seq(1, length(probes), by = block_size)
for (s in starts) {
    rows <- probes[s:min(s + block_size - 1, length(probes))]
    beta <- round(as.matrix(assay(se, "beta_unmasked")[rows, , drop = FALSE]), 6)
    detp <- signif(as.matrix(assay(se, "detp")[rows, , drop = FALSE]), 6)
    write_block(rows, interleave(list(beta, detp), c("", " Detection Pval")), processed_file, s == 1)
    in_raw <- rows[rows %in% common]
    write_block(in_raw, interleave(list(unmeth[in_raw, , drop = FALSE], meth[in_raw, , drop = FALSE], detp[in_raw, , drop = FALSE]),
                                   c(" Unmethylated signal", " Methylated signal", " Detection Pval")), intensity_file, s == 1)
}

# Start of the GEO metadata sheet
meta <- data.frame(sample = arrays, idat_red = paste0(samples[["idat_prefix"]], "_Red.idat"),
                   idat_grn = paste0(samples[["idat_prefix"]], "_Grn.idat"), check.names = FALSE)
for (column in intersect(c("sentrix_id", "sentrix_position"), colnames(samples))) meta[[column]] <- samples[[column]]
meta[["processing"]] <- paste0("preprocessing route ", metadata(se)[["preprocessing"]],
                               if (!is.null(metadata(se)[["normalisation"]])) paste0(" (", metadata(se)[["normalisation"]], ")") else "",
                               "; beta values without detection masking; detection: ",
                               if (!is.null(metadata(se)[["detection"]])) metadata(se)[["detection"]] else "pOOBAH")
fwrite(meta, paste0(prefix, "_samples.tsv"), sep = "\\t")
cat(length(probes), "probes,", length(arrays), "arrays\n")
