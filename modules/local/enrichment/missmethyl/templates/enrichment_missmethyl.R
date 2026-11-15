#!/usr/bin/env Rscript
#
# Gene set enrichment of the CpGs significant in one association model, with
# missMethyl (Phipson et al. 2016; Maksimovic et al. 2021), which corrects for
# the number of CpGs per gene:
#   - GO (gometh; org.Hs.eg.db and GO.db, offline)
#   - KEGG (gometh; the pathways are downloaded at run time: if that fails, the
#     reason is recorded and the other collections still run)
#   - every staged MSigDB GMT file (Entrez IDs; gsameth), e.g. Hallmark, Reactome
# Significant CpGs: FDR < --fdr per contrast; background: every tested CpG.
# The gene names and gene regions of the probe annotation are passed as 'anno',
# so the CpG IDs match directly and the EPIC v2 annotation package is not read.
#

prefix <- "${prefix}"
stats_file <- "${stats}"
resolved_file <- "${resolved}"
gmt_files <- strsplit(trimws("${gmts}"), "[[:space:]]+")[[1]]
gmt_files <- gmt_files[nzchar(gmt_files)]

defaults <- list(fdr = "0.05", min_hits = "10", collections = "GO,KEGG")
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
collections <- intersect(strsplit(opts[["collections"]], ",")[[1]], c("GO", "KEGG"))

Sys.setenv(HOME = getwd())
suppressPackageStartupMessages({
    library(missMethyl)
    library(data.table)
    library(jsonlite)
})
set.seed(1)
writeLines(c('"${task.process}":',
    paste0("    missMethyl: ", packageVersion("missMethyl")),
    paste0("    org.Hs.eg.db: ", packageVersion("org.Hs.eg.db")),
    paste0("    GO.db: ", packageVersion("GO.db"))), "versions.yml")

resolved <- fromJSON(resolved_file, simplifyVector = TRUE)
model <- resolved[["model"]]
report <- list(model = model[["name"]], contrasts = list())
empty <- data.table(collection = character(), contrast = character(), set = character(), term = character(),
                    ontology = character(), n = integer(), de = integer(), p_value = numeric(), fdr = numeric())
write_outputs <- function(table, status) {
    report[["status"]] <<- status
    fwrite(table, paste0(prefix, ".tsv"), sep = "\\t", na = "")
    write_json(report, paste0(prefix, "_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
}

stats <- fread(text = readLines(gzfile(stats_file)), na.strings = c("", "NA"))
if (isTRUE(model[["skip_enrichment"]]) || nrow(stats) == 0) {
    write_outputs(empty, if (isTRUE(model[["skip_enrichment"]])) "skipped: skip_enrichment" else "not run: the model was not tested")
    quit(save = "no", status = 0)
}

# Annotation of every tested CpG in the form missMethyl expects
genes <- unique(stats[, .(cpg_id, gene, gene_group)], by = "cpg_id")
anno <- data.frame(Name = genes[["cpg_id"]],
                   UCSC_RefGene_Name = ifelse(is.na(genes[["gene"]]), "", genes[["gene"]]),
                   UCSC_RefGene_Group = ifelse(is.na(genes[["gene_group"]]), "", genes[["gene_group"]]),
                   row.names = genes[["cpg_id"]])

# MSigDB collections from GMT files (set name, description, Entrez IDs)
read_gmt <- function(file) {
    lines <- strsplit(readLines(file), "\\t")
    setNames(lapply(lines, function(x) x[-(1:2)]), vapply(lines, `[`, character(1), 1))
}
gmts <- setNames(lapply(gmt_files, read_gmt), sub("\\\\.(entrez\\\\.)?gmt\$", "", basename(gmt_files)))

tables <- list()
for (k in unique(stats[["contrast"]])) {
    s <- stats[stats[["contrast"]] == k & !is.na(stats[["p_value"]])]
    all_cpg <- unique(s[["cpg_id"]])
    sig_cpg <- unique(s[["cpg_id"]][s[["fdr"]] < fdr_level])
    notes <- character()
    if (length(sig_cpg) < min_hits) {
        report[["contrasts"]][[k]] <- list(significant_cpgs = length(sig_cpg), status = sprintf("not run: %d significant CpGs, fewer than %d", length(sig_cpg), min_hits))
        next
    }
    for (collection in collections) {
        result <- tryCatch(suppressMessages(suppressWarnings(
            gometh(sig.cpg = sig_cpg, all.cpg = all_cpg, collection = collection, array.type = "EPIC_V2", anno = anno))),
            error = function(e) { notes <<- c(notes, paste0(collection, " not run: ", trimws(conditionMessage(e)))); NULL })
        if (is.null(result)) next
        tables[[paste(k, collection)]] <- data.table(collection = collection, contrast = k, set = rownames(result),
            term = if ("TERM" %in% names(result)) result[["TERM"]] else if ("Description" %in% names(result)) result[["Description"]] else NA_character_,
            ontology = if ("ONTOLOGY" %in% names(result)) result[["ONTOLOGY"]] else NA_character_,
            n = result[["N"]], de = result[["DE"]], p_value = signif(result[["P.DE"]], 4), fdr = signif(result[["FDR"]], 4))
    }
    for (name in names(gmts)) {
        result <- tryCatch(suppressMessages(suppressWarnings(
            gsameth(sig.cpg = sig_cpg, all.cpg = all_cpg, collection = gmts[[name]], array.type = "EPIC_V2", anno = anno))),
            error = function(e) { notes <<- c(notes, paste0(name, " not run: ", trimws(conditionMessage(e)))); NULL })
        if (is.null(result)) next
        tables[[paste(k, name)]] <- data.table(collection = name, contrast = k, set = rownames(result), term = rownames(result),
            ontology = NA_character_, n = result[["N"]], de = result[["DE"]], p_value = signif(result[["P.DE"]], 4), fdr = signif(result[["FDR"]], 4))
    }
    report[["contrasts"]][[k]] <- list(significant_cpgs = length(sig_cpg), tested_cpgs = length(all_cpg), status = "run",
                                       notes = if (length(notes) > 0) notes else NULL)
}
table <- if (length(tables) > 0) rbindlist(tables, fill = TRUE) else empty
setorder(table, collection, contrast, p_value)
write_outputs(table, "run")
