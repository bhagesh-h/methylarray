#!/usr/bin/env Rscript
#
# One table across the association models: N, genomic lambda, bacon inflation
# and bias, hits at FDR, the epigenome-wide and suggestive thresholds (raw and
# bacon-corrected), and the top CpG of each model and contrast.
#

prefix <- "${prefix}"
summary_files <- strsplit("${summaries}", " ")[[1]]

suppressPackageStartupMessages({
    library(jsonlite)
    library(data.table)
})
writeLines(c('"${task.process}":', paste0("    jsonlite: ", packageVersion("jsonlite"))), "versions.yml")

rows <- rbindlist(lapply(sort(summary_files), function(file) {
    s <- fromJSON(file, simplifyVector = FALSE)
    base <- list(model = s[["model"]], trait = s[["trait"]], n = s[["n"]], status = s[["status"]])
    if (is.null(s[["contrasts"]]) || length(s[["contrasts"]]) == 0) return(as.data.table(base))
    rbindlist(lapply(names(s[["contrasts"]]), function(contrast) {
        as.data.table(c(base, list(contrast = contrast), s[["contrasts"]][[contrast]]))
    }), fill = TRUE)
}), fill = TRUE)
fwrite(rows, paste0(prefix, ".tsv"), sep = "\\t", na = "")

columns <- intersect(c("model", "contrast", "n", "lambda", "bacon_inflation", "bacon_bias", "hits_fdr", "hits_epigenome_wide",
                       "hits_suggestive", "hits_fdr_bacon", "top_cpg", "top_p", "status"), names(rows))
writeLines(c(
    "# id: 'ewas_summary'",
    "# section_name: 'Association models'",
    "# description: 'One row per model and contrast: complete cases, genomic lambda, bacon inflation and bias, CpGs at FDR < 0.05, P < 9e-8 and P < 1e-5 (and at FDR < 0.05 after bacon), and the top CpG.'",
    "# plot_type: 'table'"
), paste0(prefix, "_mqc.tsv"))
table <- rows[, ..columns]
table[, Model := if ("contrast" %in% names(table)) ifelse(is.na(contrast), model, paste0(model, " (", contrast, ")")) else model]
setcolorder(table, "Model")
table[, c("model") := NULL]
if ("contrast" %in% names(table)) table[, contrast := NULL]
write.table(table, paste0(prefix, "_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE, na = "")
