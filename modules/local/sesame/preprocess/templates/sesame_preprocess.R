#!/usr/bin/env Rscript
#
# SeSAMe preprocessing of one array: read the IDAT pair, run prepSesame and
# write beta values, detection p-values, intensities, QC statistics and
# control-probe intensities.
#

sample <- "${meta.id}"
prefix <- "${task.ext.prefix ?: meta.id}"
idat_red <- "${idat_red}"
idat_grn <- "${idat_grn}"
cache_dir <- normalizePath("${sesame_cache}")

# Options from ext.args, written as '--name value'
defaults <- list(prep = "QCDPB", min_beads = "3", detection_pval = "0.05", platform = "EPICv2")
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
min_beads <- as.integer(opts[["min_beads"]])
detection_pval <- as.numeric(opts[["detection_pval"]])
prep <- opts[["prep"]]

# Offline: every resource comes from the cache prepared by SESAME_CACHE
Sys.setenv(
    EXPERIMENT_HUB_CACHE = cache_dir,
    HOME = getwd(),
    XDG_CACHE_HOME = file.path(getwd(), ".cache")
)
suppressPackageStartupMessages({
    library(ExperimentHub)
    setExperimentHubOption("CACHE", cache_dir)
    setExperimentHubOption("LOCAL", TRUE)
    library(sesame)
})
set.seed(1)

# readIDATpair expects <prefix>_Red.idat(.gz) and <prefix>_Grn.idat(.gz)
idat_prefix <- sub("_Red[.]idat([.]gz)?\$", "", idat_red)
if (sub("_Grn[.]idat([.]gz)?\$", "", idat_grn) != idat_prefix) {
    stop("The red and green IDAT files of ", sample, " do not share a prefix")
}

sdf <- readIDATpair(idat_prefix, min_beads = min_beads)
platform <- sdfPlatform(sdf)
if (platform != opts[["platform"]]) {
    stop("Sample ", sample, " is a ", platform, " array, but the run is set to --platform ", opts[["platform"]])
}

# Preprocessing (default QCDPB: quality mask, channel inference, non-linear dye
# bias correction, pOOBAH detection masking, noob background correction)
prep_args <- list(P = list(pval.threshold = detection_pval))
processed <- prepSesame(sdf, prep = prep, prep_args = prep_args)

# Detection p-values from the signal that the P step sees: after the steps before it
detp_input <- if (grepl("P", prep, fixed = TRUE)) prepSesame(sdf, prep = sub("P.*", "", prep), prep_args = prep_args) else sdf
detp <- pOOBAH(detp_input, return.pval = TRUE)

beta <- getBetas(processed, mask = TRUE)
beta_unmasked <- getBetas(processed, mask = FALSE)
mu <- signalMU(processed, mask = FALSE)
probes <- names(beta)
rownames(mu) <- mu[["Probe_ID"]]

saveRDS(
    list(
        sample = sample,
        platform = platform,
        probes = probes,
        beta = unname(beta),
        beta_unmasked = unname(beta_unmasked[probes]),
        detp = unname(detp[probes]),
        meth = unname(mu[probes, "M"]),
        unmeth = unname(mu[probes, "U"])
    ),
    paste0(prefix, ".sesame.rds")
)

# Per-sample QC: SeSAMe QC statistics on the raw signal, plus the masking result
stats <- as.data.frame(sesameQC_calcStats(sdf))
qc <- data.frame(
    sample = sample,
    platform = platform,
    prep = prep,
    predicted_sex = tolower(inferSex(beta, platform = platform)),
    num_masked = sum(is.na(beta)),
    frac_masked = mean(is.na(beta)),
    stats,
    check.names = FALSE
)
write.table(qc, paste0(prefix, ".qc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE)

# Control-probe intensities (for control metrics and control-probe PCA)
controls <- controls(sdf)
controls <- data.frame(sample = sample, controls[, c("Probe_ID", "UG", "UR")])
gz <- gzfile(paste0(prefix, ".controls.tsv.gz"), "w")
write.table(controls, gz, sep = "\\t", quote = FALSE, row.names = FALSE)
close(gz)

# Software versions
writeLines(c('"${task.process}":', paste0("    sesame: ", packageVersion("sesame")), paste0("    sesameData: ", packageVersion("sesameData"))), "versions.yml")
