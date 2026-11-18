# Builds the test fixture in this folder (test_se/): the rows of a filtered,
# CpG-level SummarizedExperiment (minfi route, -profile test, six public GEO
# arrays) that match a reference CpG of cent12CT.m, centDHSbloodDMC.m or the
# IDOL table of FlowSorted.Blood.EPIC (IDOLOptimizedCpGs.compTable) by CpG,
# EPIC v1 or 450K ID.
#   Rscript make_fixture.R <filtered_se> <reference ID list> <this folder>
# annotation.tsv.gz holds the rows of probe_annotation.tsv.gz whose cpg_id,
# epicv1_id or hm450_id is a reference CpG (columns probe_id, cpg_id, chr, pos,
# epicv1_id, hm450_id).
suppressPackageStartupMessages({ library(SummarizedExperiment); library(HDF5Array) })
args <- commandArgs(TRUE)
se <- loadHDF5SummarizedExperiment(args[1])
ids <- readLines(args[2])
rd <- rowData(se)
keep <- rownames(se) %in% ids | rd[["epicv1_id"]] %in% ids | rd[["hm450_id"]] %in% ids
small <- se[keep, ]
assays(small) <- lapply(assays(small)[c("beta_unmasked")], as.matrix)
rowData(small) <- rowData(small)[, c("cpg_id", "chr", "pos", "epicv1_id", "hm450_id")]
colData(small) <- colData(small)[, "sample", drop = FALSE]
metadata(small) <- list()
saveHDF5SummarizedExperiment(small, dir = file.path(args[3], "test_se"), replace = TRUE)
cat(nrow(small), "rows x", ncol(small), "arrays\n")
