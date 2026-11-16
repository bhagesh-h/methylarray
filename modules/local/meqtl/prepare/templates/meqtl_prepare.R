#!/usr/bin/env Rscript
#
# CpGs with a cis-meQTL, from a user-supplied meQTL table (--meqtl), for example
# the GoDMC meta-analysis (Min et al. 2021). Only a 'cpg' column is required;
# when present, 'cistrans' keeps the cis associations (TRUE) and 'pval' keeps
# P < --meqtl_p (default 1e-8, the GoDMC cis threshold). The CpG names (450K or
# EPIC v1) are matched to the EPIC v2 CpGs directly and through the EPIC v1 and
# 450K IDs of the probe annotation. Output: one row per EPIC v2 CpG with its
# number of meQTL SNPs and the strongest one.
#

meqtl_file <- "${meqtl}"
annotation_file <- "${annotation}"
prefix <- "${prefix}"

defaults <- list(meqtl_p = "1e-8")
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
p_max <- as.numeric(opts[["meqtl_p"]])

suppressPackageStartupMessages(library(data.table))
writeLines(c('"${task.process}":', paste0("    R: ", getRversion()), paste0("    data.table: ", packageVersion("data.table"))), "versions.yml")

# Read only the needed columns; compressed tables are streamed through gzip
reader <- function(...) if (grepl("[.]gz\$", meqtl_file)) fread(cmd = paste("gzip -dc", shQuote(meqtl_file)), ...) else fread(meqtl_file, ...)
header <- names(reader(nrows = 0))
if (!"cpg" %in% header) stop("The meQTL table needs a 'cpg' column; found: ", paste(header, collapse = ", "))
columns <- intersect(c("cpg", "snp", "pval", "cistrans", "beta_a1"), header)
meqtl <- reader(select = columns, colClasses = list(character = intersect(c("cpg", "snp", "cistrans"), columns)))
rows_read <- nrow(meqtl)
if ("cistrans" %in% columns) meqtl <- meqtl[toupper(cistrans) %in% "TRUE"]
rows_cis <- nrow(meqtl)
if ("pval" %in% columns) meqtl <- meqtl[!is.na(pval) & pval < p_max]
rows_kept <- nrow(meqtl)
for (column in setdiff(c("snp", "pval", "beta_a1"), columns)) meqtl[, (column) := NA]

# Per CpG name: the number of SNPs and the strongest association
per_name <- meqtl[order(pval), .(snps = .N, top_snp = snp[1], top_p = pval[1], top_beta = beta_a1[1]), by = cpg]

# Match the names to the EPIC v2 CpGs
anno <- fread(cmd = paste("gzip -dc", shQuote(annotation_file)), select = c("cpg_id", "epicv1_id", "hm450_id"),
              colClasses = "character", na.strings = "")
anno <- unique(anno)
names_map <- unique(rbind(
    data.table(name = anno[["cpg_id"]], cpg_id = anno[["cpg_id"]]),
    anno[!is.na(epicv1_id), .(name = epicv1_id, cpg_id)],
    anno[!is.na(hm450_id), .(name = hm450_id, cpg_id)]
))
matched <- merge(per_name, names_map, by.x = "cpg", by.y = "name")
matched <- matched[order(top_p)][!duplicated(cpg_id)]
result <- matched[, .(cpg_id, meqtl_name = cpg, snps, top_snp, top_p = signif(top_p, 4), top_beta = signif(top_beta, 4))][order(cpg_id)]
fwrite(result, paste0(prefix, ".tsv"), sep = "\\t", na = "")
fwrite(data.table(metric = c("rows_read", "rows_cis", "rows_kept", "p_threshold", "cpgs_in_table", "cpgs_on_array"),
                  value = c(rows_read, rows_cis, rows_kept, if ("pval" %in% columns) p_max else NA, nrow(per_name), nrow(result))),
       paste0(prefix, "_summary.tsv"), sep = "\\t", na = "")
