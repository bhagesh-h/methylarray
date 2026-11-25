#!/usr/bin/env Rscript
#
# Genomic context enrichment of the CpGs significant in one association model,
# with KnowYourCG (Zhou et al.): chromatin states (ChromHMM), transcription
# factor binding sites, histone marks, CTCF binding, A/B compartments,
# imprinting DMRs and gene features (metagene), from the staged knowledgebases.
# Per contrast, three queries: all CpGs at FDR < --fdr, those with higher and
# those with lower methylation (sign of delta_beta), each against every tested
# CpG (the universe); queries with fewer than --min_hits CpGs are not run.
# The CpG IDs are mapped to the EPIC v2 probe IDs of the knowledgebases through
# the probe annotation (replicate probes of a CpG all count).
# Test: Fisher's exact test, log2 odds ratio as the estimate, FDR within each
# knowledgebase group (knowYourCG::testEnrichment).
#

prefix <- "${prefix}"
stats_file <- "${stats}"
resolved_file <- "${resolved}"
cache_dir <- "${cache}"
annotation_file <- "${annotation}"

defaults <- list(fdr = "0.05", min_hits = "10", platform = "EPICv2")
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
fdr_level <- as.numeric(opts[["fdr"]])
min_hits <- as.integer(opts[["min_hits"]])

# The staged cache is read only: copy it, ExperimentHub writes to its index
local_cache <- file.path(getwd(), "kycg_cache_local")
dir.create(local_cache, showWarnings = FALSE)
file.copy(list.files(cache_dir, full.names = TRUE), local_cache, recursive = TRUE)
Sys.setenv(EXPERIMENT_HUB_CACHE = local_cache, HOME = getwd(), XDG_CACHE_HOME = file.path(getwd(), ".cache"))
suppressPackageStartupMessages({
    library(ExperimentHub)
    library(knowYourCG)
    library(data.table)
    library(jsonlite)
})
setExperimentHubOption("LOCAL", TRUE)
writeLines(c('"${task.process}":', paste0("    knowYourCG: ", packageVersion("knowYourCG")),
             paste0("    ExperimentHub: ", packageVersion("ExperimentHub"))), "versions.yml")

resolved <- fromJSON(resolved_file, simplifyVector = TRUE)
model <- resolved[["model"]]
report <- list(model = model[["name"]], contrasts = list())
empty <- data.table(contrast = character(), query = character(), database = character(), set = character(), log2_or = numeric(),
                    p_value = numeric(), fdr = numeric(), overlap = integer(), query_probes = integer(), set_probes = integer())
write_outputs <- function(table, status) {
    report[["status"]] <<- status
    fwrite(table, paste0(prefix, ".tsv"), sep = "\\t", na = "")
    write_json(report, paste0(prefix, "_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
    # MultiQC: the 10 sets with the smallest P per contrast and query, enriched sets only
    top <- table[log2_or > 0][order(p_value)][, head(.SD, 10), by = .(contrast, query)]
    if (nrow(top) > 0) {
        writeLines(c(
            paste0("# id: 'kycg_", model[["name"]], "'"),
            paste0("# section_name: 'Genomic context: ", model[["name"]], "'"),
            "# description: 'KnowYourCG enrichment of the CpGs at FDR < 0.05 (all, higher and lower methylation) against all tested CpGs: the 10 enriched sets with the smallest P per contrast and query. Log2 OR: log2 odds ratio; Overlap: significant probes in the set.'",
            "# plot_type: 'table'"
        ), paste0(prefix, "_mqc.tsv"))
        write.table(data.frame(Set = paste(top[["contrast"]], top[["query"]], top[["database"]], top[["set"]], sep = " | "),
                               Contrast = top[["contrast"]], Query = top[["query"]], Database = top[["database"]], Set = top[["set"]],
                               `Log2 OR` = round(top[["log2_or"]], 2), Overlap = top[["overlap"]], P = signif(top[["p_value"]], 3),
                               FDR = signif(top[["fdr"]], 3), check.names = FALSE),
                    paste0(prefix, "_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE)
    }
}

stats <- fread(text = readLines(gzfile(stats_file)), na.strings = c("", "NA"))
if (isTRUE(model[["skip_enrichment"]]) || nrow(stats) == 0) {
    write_outputs(empty, if (isTRUE(model[["skip_enrichment"]])) "skipped: skip_enrichment" else "not run: the model was not tested")
    quit(save = "no", status = 0)
}

# Knowledgebases in the staged cache
eh <- ExperimentHub(localHub = TRUE)
titles <- grep(paste0("^KYCG[.]", opts[["platform"]], "[.]"), mcols(eh)[["title"]], value = TRUE)
titles <- setdiff(titles, grep("[.]Mask[.]", titles, value = TRUE))
if (length(titles) == 0) stop("No KnowYourCG knowledgebase for ", opts[["platform"]], " in the staged cache")
databases <- getDBs(titles, silent = TRUE)
group_name <- function(g) sub("[.][0-9]{8}\$", "", sub("^KYCG[.][^.]+[.]", "", g))

# CpG to probe IDs (probe IDs, with --replicate_probe_strategy none, stay as they are)
anno <- fread(text = readLines(gzfile(annotation_file)), select = c("probe_id", "cpg_id"))
probes_of <- function(cpgs) anno[["probe_id"]][anno[["cpg_id"]] %in% cpgs | anno[["probe_id"]] %in% cpgs]

tables <- list()
for (k in unique(stats[["contrast"]])) {
    s <- stats[stats[["contrast"]] == k & !is.na(stats[["p_value"]])]
    universe <- probes_of(unique(s[["cpg_id"]]))
    sig <- s[s[["fdr"]] < fdr_level]
    queries <- list(all = unique(sig[["cpg_id"]]))
    if ("delta_beta" %in% names(sig)) {
        queries[["higher"]] <- unique(sig[["cpg_id"]][sig[["delta_beta"]] > 0])
        queries[["lower"]] <- unique(sig[["cpg_id"]][sig[["delta_beta"]] < 0])
    }
    status <- list()
    for (q in names(queries)) {
        if (length(queries[[q]]) < min_hits) {
            status[[q]] <- sprintf("not run: %d significant CpGs, fewer than %d", length(queries[[q]]), min_hits)
            next
        }
        result <- testEnrichment(probes_of(queries[[q]]), databases = databases, universe = universe,
                                 platform = opts[["platform"]], silent = TRUE)
        result <- as.data.table(result)
        # testEnrichment reports an infinite odds ratio (no overlap, or every query probe in the set) as +-1024
        result[abs(estimate) >= 1000, estimate := sign(estimate) * Inf]
        tables[[paste(k, q)]] <- data.table(contrast = k, query = q, database = group_name(result[["group"]]), set = result[["dbname"]],
                                            log2_or = signif(result[["estimate"]], 4), p_value = signif(result[["p.value"]], 4),
                                            fdr = signif(result[["FDR"]], 4), overlap = result[["overlap"]], query_probes = result[["nQ"]],
                                            set_probes = result[["nD"]])
        status[[q]] <- sprintf("run: %d significant CpGs", length(queries[[q]]))
    }
    report[["contrasts"]][[k]] <- list(significant_cpgs = length(queries[["all"]]), tested_cpgs = length(unique(s[["cpg_id"]])), queries = status)
}
table <- if (length(tables) > 0) rbindlist(tables) else empty
setorder(table, contrast, query, p_value)
write_outputs(table, if (length(tables) > 0) "run" else "not run: too few significant CpGs")
