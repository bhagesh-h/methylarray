#!/usr/bin/env Rscript
#
# Download the sesameData resources that SeSAMe needs for one array platform
# into a local ExperimentHub cache, so that every later task runs offline.
#

platform <- "${platform}"

# Resources used by readIDATpair (platform detection, probe addresses),
# prepSesame (Q: quality mask) and inferSex (lift-over to HM450 probe IDs)
resources <- list(
    EPICv2 = c("idatSignature", "probeIDSignature", "EPICv2.address", "KYCG.EPICv2.Mask.20230314", "HM450.address")
)
if (!platform %in% names(resources)) {
    stop("SeSAMe resources are defined for these platforms only: ", paste(names(resources), collapse = ", "))
}
titles <- resources[[platform]]

# Keep every cache inside the task directory: the home directory may not be writable
cache_dir <- file.path(getwd(), "sesame_cache")
dir.create(cache_dir, showWarnings = FALSE)
Sys.setenv(
    EXPERIMENT_HUB_CACHE = cache_dir,
    HOME = getwd(),
    XDG_CACHE_HOME = file.path(getwd(), ".cache")
)

suppressPackageStartupMessages({
    library(ExperimentHub)
    library(sesameData)
})

sesameDataCache(data_titles = titles)

# Record what was downloaded: ExperimentHub ID, cached file and its MD5
eh <- ExperimentHub(cache = cache_dir, localHub = TRUE)
hub <- mcols(eh)
rows <- lapply(titles, function(title) {
    eh_id <- rownames(hub)[hub[["title"]] == title]
    if (length(eh_id) != 1) {
        stop("Expected one ExperimentHub record for ", title, ", found ", length(eh_id))
    }
    path <- unname(cache(eh[eh_id]))
    data.frame(title = title, eh_id = eh_id, file = basename(path), md5 = unname(tools::md5sum(path)))
})
write.table(do.call(rbind, rows), "sesame_resources.tsv", sep = "\\t", quote = FALSE, row.names = FALSE)

# Software versions
writeLines(c('"${task.process}":', paste0("    sesame: ", packageVersion("sesame")), paste0("    sesameData: ", packageVersion("sesameData")), paste0("    ExperimentHub: ", packageVersion("ExperimentHub"))), "versions.yml")
