# Builds the test fixture in this folder from the minfi route of -profile test
# (six public GEO arrays): the probe rows that the clocks HorvathS2013,
# HannumG2013, LevineM2018 and DunedinPACE use, with Horvath's 21k reference CpGs
# and the DunedinPACE background CpGs,
# and the matching rows of the probe annotation.
#   Rscript make_fixture.R <minfi_se> <probe_annotation.tsv.gz> <dnaMethyAge.tar.gz> <DunedinPACE.tar.gz> <this folder>
suppressPackageStartupMessages({ library(SummarizedExperiment); library(HDF5Array) })
args <- commandArgs(TRUE)
lib <- file.path(tempdir(), "rlib"); dir.create(lib)
install.packages(args[3:4], lib = lib, repos = NULL, type = "source", quiet = TRUE)
library(dnaMethyAge, lib.loc = lib); library(DunedinPACE, lib.loc = lib)
cpgs_of <- function(clock) { env <- new.env(); data(list = clock, package = "dnaMethyAge", envir = env); as.character(get("coefs", env)[["Probe"]]) }
# Horvath's preprocessing normalises against its 21k reference CpGs
env <- new.env(); data(list = "27k_reference", package = "dnaMethyAge", envir = env)
models <- DunedinPACE:::mPACE_Models
needed <- unique(c(unlist(lapply(c("HorvathS2013", "HannumG2013", "LevineM2018"), cpgs_of)),
                   as.character(get("probeAnnotation21kdatMethUsed", env)[["Name"]]),
                   unlist(models[["model_probes"]]), unlist(models[["gold_standard_probes"]])))
anno <- read.delim(gzfile(args[2]), colClasses = "character", na.strings = "")
anno <- anno[anno[["cpg_id"]] %in% needed | anno[["epicv1_id"]] %in% needed | anno[["hm450_id"]] %in% needed,
             c("probe_id", "cpg_id", "epicv1_id", "hm450_id")]
se <- loadHDF5SummarizedExperiment(args[1])
small <- se[rownames(se) %in% anno[["probe_id"]], ]
assays(small) <- list(beta = as.matrix(assay(small, "beta")))
colData(small) <- colData(small)[, "sample", drop = FALSE]
rowData(small) <- rowData(small)[, c("probe_id", "cpg_id")]
metadata(small) <- list()
saveHDF5SummarizedExperiment(small, dir = file.path(args[5], "test_se"), replace = TRUE)
con <- gzfile(file.path(args[5], "annotation.tsv.gz"), "w")
write.table(anno, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
close(con)
cat(nrow(small), "probe rows,", nrow(anno), "annotation rows\n")
