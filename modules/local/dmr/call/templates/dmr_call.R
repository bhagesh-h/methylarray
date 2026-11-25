#!/usr/bin/env Rscript
#
# Differentially methylated regions for one association model, from the same
# statistics as the CpG-level analysis (EWAS_LIMMA):
#   - DMRcate (Peters et al. 2015): the kernel smoothing of dmrcate() on a
#     CpGannotated object built from the moderated t, delta-beta, P and FDR at
#     each CpG's hg38 position. cpg.annotate() and extractRanges() are not used:
#     for EPIC v2 they download the Peters manifest and gene annotation from
#     AnnotationHub and ExperimentHub during the task. Genes come from the probe
#     annotation of the CpGs in each region.
#   - dmrff (Suderman et al. 2018): regions from the M-value estimates, standard
#     errors and P-values, with the CpG correlation from the M-values of the
#     arrays in the model; Bonferroni-adjusted P across candidate regions.
#   - Consensus: regions of both methods that overlap and agree in direction.
#

prefix <- "${prefix}"
stats_file <- "${stats}"
resolved_file <- "${resolved}"
se_dir <- "${se}"
dmrff_source <- "${dmrff_source}"

defaults <- list(lambda = "1000", C = "2", min_cpgs = "2", fdr = "0.05", maxgap = "500", p_cutoff = "0.05")
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
lambda <- as.numeric(opts[["lambda"]])
C <- as.numeric(opts[["C"]])
min_cpgs <- as.integer(opts[["min_cpgs"]])
fdr_level <- as.numeric(opts[["fdr"]])
maxgap <- as.numeric(opts[["maxgap"]])
p_cutoff <- as.numeric(opts[["p_cutoff"]])

# dmrff is not on conda or CRAN: install the staged source (pure R, a pinned
# GitHub commit) into the task directory; never install from the network
Sys.setenv(HOME = getwd())
options(repos = character(0), mc.cores = as.integer("${task.cpus}"))
library_dir <- file.path(getwd(), "rlib")
dir.create(library_dir, showWarnings = FALSE)
.libPaths(c(library_dir, .libPaths()))
install.packages(dmrff_source, lib = library_dir, repos = NULL, type = "source", quiet = TRUE)
suppressPackageStartupMessages({
    library(DMRcate)
    library(GenomicRanges)
    library(SummarizedExperiment)
    library(HDF5Array)
    library(data.table)
    library(jsonlite)
    library(dmrff, lib.loc = library_dir)
})
set.seed(1)
writeLines(c('"${task.process}":',
    paste0("    DMRcate: ", packageVersion("DMRcate")),
    paste0("    dmrff: ", packageVersion("dmrff", lib.loc = library_dir))), "versions.yml")

resolved <- fromJSON(resolved_file, simplifyVector = TRUE)
model <- resolved[["model"]]
report <- list(model = model[["name"]], contrasts = list())
empty_dmrcate <- data.table(contrast = character(), chr = character(), start = integer(), end = integer(), n_cpgs = integer(),
                            min_smoothed_fdr = numeric(), stouffer = numeric(), hmfdr = numeric(), fisher = numeric(),
                            maxdiff = numeric(), meandiff = numeric(), genes = character())
empty_dmrff <- data.table(contrast = character(), chr = character(), start = integer(), end = integer(), n_cpgs = integer(),
                          estimate = numeric(), se = numeric(), z = numeric(), p_value = numeric(), p_adjusted = numeric(),
                          genes = character())
write_outputs <- function(dmrcate_table, dmrff_table, consensus, status) {
    report[["status"]] <<- status
    fwrite(dmrcate_table, paste0(prefix, "_dmrcate.tsv"), sep = "\\t", na = "")
    fwrite(dmrff_table, paste0(prefix, "_dmrff.tsv"), sep = "\\t", na = "")
    fwrite(consensus, paste0(prefix, "_consensus.tsv"), sep = "\\t", na = "")
    bed <- if (nrow(consensus) > 0) consensus[, .(chr, start = start - 1L, end, name = paste0(contrast, ":", genes), score = 0, strand = ".")]
           else data.table()
    fwrite(bed, paste0(prefix, "_consensus.bed"), sep = "\\t", col.names = FALSE)
    write_json(report, paste0(prefix, "_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
    # MultiQC: one row per model and contrast; the files of all models share the
    # section id, so MultiQC joins them into one table
    rows <- if (length(report[["contrasts"]]) == 0) {
        data.table(Model = report[["model"]], cpgs = NA_integer_, cpgs_fdr = NA_integer_, dmrcate = NA_integer_, dmrff = NA_integer_,
                   consensus = NA_integer_, status = status)
    } else {
        rbindlist(lapply(names(report[["contrasts"]]), function(k) {
            x <- report[["contrasts"]][[k]]
            data.table(Model = paste0(report[["model"]], " (", k, ")"), cpgs = x[["cpgs"]], cpgs_fdr = x[["cpgs_fdr"]], dmrcate = x[["dmrcate"]],
                       dmrff = x[["dmrff"]], consensus = x[["consensus"]], status = if (is.null(x[["note"]])) status else x[["note"]])
        }))
    }
    writeLines(c(
        "# id: 'dmr_summary'",
        "# section_name: 'Differentially methylated regions'",
        "# description: 'One row per model and contrast: CpGs tested, CpGs at FDR < 0.05, regions found by DMRcate and by dmrff, and the consensus regions found by both (same direction, overlapping by at least 1 bp).'",
        "# plot_type: 'table'"
    ), paste0(prefix, "_mqc.tsv"))
    write.table(rows, paste0(prefix, "_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE, na = "")
}

stats <- fread(text = readLines(gzfile(stats_file)), na.strings = c("", "NA"))
if (isTRUE(model[["skip_dmr"]]) || nrow(stats) == 0) {
    write_outputs(empty_dmrcate, empty_dmrff, data.table(), if (isTRUE(model[["skip_dmr"]])) "skipped: skip_dmr" else "not run: the model was not tested")
    quit(save = "no", status = 0)
}

# M-values of the arrays in the model, for the dmrff correlation
se <- loadHDF5SummarizedExperiment(se_dir)
dropped <- if (length(resolved[["dropped"]]) == 0) character() else as.character(resolved[["dropped"]][["sample"]])
arrays <- setdiff(colnames(se), dropped)
mvalue <- as.matrix(assay(se, "mvalue")[, arrays, drop = FALSE])

genes_in <- function(table, chr, start, end) {
    inside <- table[["chr"]] == chr & table[["pos"]] >= start & table[["pos"]] <= end
    genes <- unique(unlist(strsplit(na.omit(table[["gene"]][inside]), ";")))
    paste(sort(genes[nzchar(genes)]), collapse = ";")
}

dmrcate_tables <- list()
dmrff_tables <- list()
consensus_tables <- list()
for (k in unique(stats[["contrast"]])) {
    s <- stats[stats[["contrast"]] == k & !is.na(stats[["p_value"]]) & !is.na(stats[["chr"]]) & !is.na(stats[["pos"]])]
    setorder(s, chr, pos)
    # One row per probe (--replicate_probe_strategy none): the first probe of each position, so that
    # replicate probes of one CpG do not make a region
    if (any(grepl("_[TB][CO][0-9]+\$", s[["cpg_id"]]))) s <- unique(s, by = c("chr", "pos"))

    # DMRcate kernel on the model statistics
    ranges <- GRanges(s[["chr"]], IRanges(s[["pos"]], s[["pos"]]), stat = s[["t"]], diff = s[["delta_beta"]],
                      rawpval = s[["p_value"]], ind.fdr = s[["fdr"]], is.sig = s[["fdr"]] < fdr_level)
    names(ranges) <- s[["cpg_id"]]
    ranges <- sort(ranges)
    note <- ""
    dmrcate_table <- empty_dmrcate
    if (sum(mcols(ranges)[["is.sig"]]) == 0) {
        note <- "DMRcate: no CpG at FDR < threshold"
    } else {
        found <- tryCatch(suppressMessages(dmrcate(new("CpGannotated", ranges = ranges, betas = matrix()),
                                                   lambda = lambda, C = C, min.cpgs = min_cpgs)),
                          error = function(e) { note <<- paste("DMRcate:", conditionMessage(e)); NULL })
        if (!is.null(found) && length(found@coord) > 0) {
            coords <- tstrsplit(found@coord, "[:-]")
            dmrcate_table <- data.table(contrast = k, chr = coords[[1]], start = as.integer(coords[[2]]), end = as.integer(coords[[3]]),
                                        n_cpgs = found@no.cpgs, min_smoothed_fdr = signif(found@min_smoothed_fdr, 4),
                                        stouffer = signif(found@Stouffer, 4), hmfdr = signif(found@HMFDR, 4), fisher = signif(found@Fisher, 4),
                                        maxdiff = signif(found@maxdiff, 4), meandiff = signif(found@meandiff, 4))
            dmrcate_table[, genes := mapply(function(c, a, b) genes_in(s, c, a, b), chr, start, end)]
        }
    }

    # dmrff on the M-value estimates
    m <- mvalue[s[["cpg_id"]], , drop = FALSE]
    candidates <- suppressMessages(dmrff(estimate = s[["logfc_m"]], se = s[["se_m"]], p.value = s[["p_value"]],
                                         methylation = m, chr = s[["chr"]], pos = s[["pos"]],
                                         maxgap = maxgap, p.cutoff = p_cutoff, verbose = FALSE))
    dmrff_table <- empty_dmrff
    if (!is.null(candidates) && nrow(candidates) > 0) {
        hits <- as.data.table(candidates)[p.adjust < fdr_level & n >= min_cpgs]
        if (nrow(hits) > 0) {
            dmrff_table <- hits[, .(contrast = k, chr = as.character(chr), start = as.integer(start), end = as.integer(end),
                                    n_cpgs = as.integer(n), estimate = signif(estimate, 4), se = signif(se, 4), z = signif(z, 4),
                                    p_value = signif(p.value, 4), p_adjusted = signif(p.adjust, 4))]
            dmrff_table[, genes := mapply(function(c, a, b) genes_in(s, c, a, b), chr, start, end)]
        }
    }

    # Consensus: overlapping regions with the same direction
    consensus <- data.table()
    if (nrow(dmrcate_table) > 0 && nrow(dmrff_table) > 0) {
        a <- GRanges(dmrcate_table[["chr"]], IRanges(dmrcate_table[["start"]], dmrcate_table[["end"]]))
        b <- GRanges(dmrff_table[["chr"]], IRanges(dmrff_table[["start"]], dmrff_table[["end"]]))
        hit <- findOverlaps(a, b)
        same <- sign(dmrcate_table[["meandiff"]][queryHits(hit)]) == sign(dmrff_table[["estimate"]][subjectHits(hit)])
        hit <- hit[same]
        if (length(hit) > 0) {
            consensus <- data.table(contrast = k,
                                    chr = dmrcate_table[["chr"]][queryHits(hit)],
                                    start = pmin(dmrcate_table[["start"]][queryHits(hit)], dmrff_table[["start"]][subjectHits(hit)]),
                                    end = pmax(dmrcate_table[["end"]][queryHits(hit)], dmrff_table[["end"]][subjectHits(hit)]),
                                    dmrcate_min_smoothed_fdr = dmrcate_table[["min_smoothed_fdr"]][queryHits(hit)],
                                    dmrcate_meandiff = dmrcate_table[["meandiff"]][queryHits(hit)],
                                    dmrff_p_adjusted = dmrff_table[["p_adjusted"]][subjectHits(hit)],
                                    dmrff_estimate = dmrff_table[["estimate"]][subjectHits(hit)])
            consensus[, genes := mapply(function(c, a, b) genes_in(s, c, a, b), chr, start, end)]
        }
    }

    dmrcate_tables[[k]] <- dmrcate_table
    dmrff_tables[[k]] <- dmrff_table
    consensus_tables[[k]] <- consensus
    report[["contrasts"]][[k]] <- list(cpgs = nrow(s), cpgs_fdr = sum(s[["fdr"]] < fdr_level, na.rm = TRUE),
                                              dmrcate = nrow(dmrcate_table), dmrff = nrow(dmrff_table),
                                              consensus = nrow(consensus), note = if (nzchar(note)) note else NULL)
}
write_outputs(rbindlist(dmrcate_tables, fill = TRUE), rbindlist(dmrff_tables, fill = TRUE),
              rbindlist(consensus_tables, fill = TRUE), "run")
