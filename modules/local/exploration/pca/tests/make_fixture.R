# Builds test_se/ in this folder: 5,000 CpGs (every 170th row) of the filtered,
# CpG-level SummarizedExperiment of the minfi route of -profile test (six public
# GEO arrays), with the mvalue assay only. samples.tsv, sample_qc.tsv,
# controls.tsv.gz and cell_fractions.tsv are the matching tables of that run.
#   Rscript make_fixture.R <filtered_se> <this folder>
suppressPackageStartupMessages({ library(SummarizedExperiment); library(HDF5Array) })
args <- commandArgs(TRUE)
se <- loadHDF5SummarizedExperiment(args[1])
small <- se[seq(1, nrow(se), by = 170)[1:5000], ]
assays(small) <- list(mvalue = as.matrix(assay(small, "mvalue")))
rowData(small) <- rowData(small)[, c("cpg_id", "chr", "pos")]
colData(small) <- colData(small)[, "sample", drop = FALSE]
metadata(small) <- list()
saveHDF5SummarizedExperiment(small, dir = file.path(args[2], "test_se"), replace = TRUE)
cat(nrow(small), "CpGs x", ncol(small), "arrays\n")
