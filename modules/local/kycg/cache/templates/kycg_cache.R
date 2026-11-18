#!/usr/bin/env Rscript
#
# Download the KnowYourCG knowledgebases of one array platform (ExperimentHub,
# sesameData) into a local ExperimentHub cache, so that the enrichment tasks
# run offline. --databases names the groups without platform and date, for
# example ChromHMM,TFBS (see knowYourCG::listDBGroups()).
#

platform <- "${platform}"
defaults <- list(databases = "ChromHMM,TFBS,HM,CTCFbind,ABCompartment,ImprintingDMR,MetagenePC")
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
wanted <- trimws(strsplit(opts[["databases"]], ",")[[1]])

# Keep every cache inside the task directory: the home directory may not be writable
cache_dir <- file.path(getwd(), "kycg_cache")
dir.create(cache_dir, showWarnings = FALSE)
Sys.setenv(EXPERIMENT_HUB_CACHE = cache_dir, HOME = getwd(), XDG_CACHE_HOME = file.path(getwd(), ".cache"))
suppressPackageStartupMessages({
    library(ExperimentHub)
    library(knowYourCG)
})
writeLines(c('"${task.process}":', paste0("    knowYourCG: ", packageVersion("knowYourCG")),
             paste0("    ExperimentHub: ", packageVersion("ExperimentHub"))), "versions.yml")

# Group titles of the platform, e.g. KYCG.EPICv2.ChromHMM.20220303
groups <- listDBGroups(paste0("^KYCG[.]", platform, "[.]"))[["Title"]]
short <- sub("[.][0-9]{8}\$", "", sub(paste0("^KYCG[.]", platform, "[.]"), "", groups))
unknown <- setdiff(wanted, short)
if (length(unknown) > 0) {
    stop("Unknown KnowYourCG databases for ", platform, ": ", paste(unknown, collapse = ", "),
         ". Available: ", paste(sort(short), collapse = ", "))
}
titles <- groups[match(wanted, short)]
kycgDataCache(data_titles = titles)

# Record what was downloaded: group, ExperimentHub ID and the MD5 of the cached
# file (not its name: ExperimentHub prefixes it with a random ID)
eh <- ExperimentHub(cache = cache_dir, localHub = TRUE)
hub <- mcols(eh)
rows <- lapply(titles, function(title) {
    eh_id <- rownames(hub)[hub[["title"]] == title]
    if (length(eh_id) != 1) stop("Expected one ExperimentHub record for ", title, ", found ", length(eh_id))
    path <- unname(cache(eh[eh_id]))
    data.frame(database = short[match(title, groups)], title = title, eh_id = eh_id, md5 = unname(tools::md5sum(path)))
})
write.table(do.call(rbind, rows), "kycg_resources.tsv", sep = "\\t", quote = FALSE, row.names = FALSE)
