# Builds the simulated test fixture in this folder: a probe-level
# SummarizedExperiment of 6 arrays with the CpGs of the four DNA methylation
# scores (assets/dnam_scores), where five smoking CpGs are missing, one IL-6 CpG
# is measured only under a new EPIC v2 name (its old ID is in annotation.tsv.gz),
# and some values failed detection; S6 failed sample QC. Random values.
#   Rscript make_fixture.R <assets/dnam_scores> <this folder>
suppressPackageStartupMessages({ library(SummarizedExperiment); library(HDF5Array); library(data.table) })
args <- commandArgs(TRUE)
set.seed(5)
weights <- rbindlist(lapply(list.files(args[1], pattern = "^(mccartney|stevenson).*tsv$", full.names = TRUE), fread))
cpgs <- unique(weights[["cpg_id"]])
smoking <- fread(file.path(args[1], "mccartney2018_smoking.tsv"))[["cpg_id"]]
il6 <- fread(file.path(args[1], "stevenson2021_il6.tsv"))[["cpg_id"]]
cpgs <- setdiff(cpgs, smoking[1:5])
renamed <- setdiff(il6, smoking)[1]
ids <- ifelse(cpgs == renamed, "cg99999999", cpgs)
probes <- paste0(ids, "_TC11")
samples <- paste0("S", 1:6)
beta <- matrix(runif(length(probes) * 6, 0.05, 0.95), length(probes), 6, dimnames = list(probes, samples))
beta[sample(length(beta), 30)] <- NA
se <- SummarizedExperiment(assays = list(beta = beta), colData = DataFrame(sample = samples, row.names = samples),
                           rowData = DataFrame(probe_id = probes, row.names = probes))
saveHDF5SummarizedExperiment(se, dir = file.path(args[2], "test_se"), replace = TRUE)
writeLines(samples[1:5], file.path(args[2], "retained.txt"))
anno <- data.table(probe_id = probes, cpg_id = ids, epicv1_id = ifelse(ids == "cg99999999", renamed, ids), hm450_id = NA)
con <- gzfile(file.path(args[2], "annotation.tsv.gz"), "w")
write.table(anno, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
close(con)
cat(length(probes), "probes; renamed", renamed, "\n")
