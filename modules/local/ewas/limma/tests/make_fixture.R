# Builds the simulated test fixture in this folder: 40 arrays x 2,000 CpGs with a
# planted association of the first 20 CpGs with 'score' (continuous) and of CpGs
# 21-40 with 'smoker' (binary), plus the covariate tables a model can use. The
# values are random; they do not describe any study.
#   Rscript make_fixture.R <this folder>
suppressPackageStartupMessages({ library(SummarizedExperiment); library(HDF5Array) })
out <- commandArgs(TRUE)[1]
set.seed(42)
n <- 40; p <- 2000
samples <- sprintf("S%02d", seq_len(n))
age <- round(runif(n, 30, 70), 1)
sex <- rep(c("female", "male"), length.out = n)
plate <- rep(c("P1", "P2", "P3", "P4"), each = 10)
score <- rnorm(n)
smoker <- rep(c("no", "no", "yes", "yes"), times = 10)
cells <- matrix(runif(n * 6), n, 6, dimnames = list(samples, c("B", "CD4T", "CD8T", "NK", "Mono", "Gran")))
cells <- cells / rowSums(cells)
mvalue <- matrix(rnorm(n * p, mean = rep(rnorm(p, 0, 2), n), sd = 0.4), p, n,
                 dimnames = list(sprintf("cg%08d", seq_len(p)), samples))
mvalue[1:20, ] <- mvalue[1:20, ] + outer(rep(0.5, 20), score)
mvalue[21:40, ] <- mvalue[21:40, ] + outer(rep(0.6, 20), smoker == "yes")
mvalue[sample(length(mvalue), 200)] <- NA
beta <- 2^mvalue / (1 + 2^mvalue)
se <- SummarizedExperiment(
    assays = list(beta = beta, mvalue = mvalue),
    colData = DataFrame(sample = samples, row.names = samples),
    rowData = DataFrame(cpg_id = rownames(mvalue), chr = paste0("chr", rep(1:20, each = 100)), pos = rep(seq(1000, 100000, length.out = 100), 20),
                        gene = NA_character_, gene_group = NA_character_, cgi_relation = NA_character_, row.names = rownames(mvalue))
)
saveHDF5SummarizedExperiment(se, dir = file.path(out, "test_se"), replace = TRUE)
write.table(data.frame(sample = samples, plate = plate, sentrix_row = rep(paste0("R0", 1:8), 5), reported_age = age),
            file.path(out, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = samples, status = "PASS", predicted_sex_minfi = sex),
            file.path(out, "sample_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pheno <- data.frame(sample = samples, score = round(score, 4), age = age, smoker = smoker)
pheno$score[c(3, 17)] <- NA
write.table(pheno, file.path(out, "phenotypes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(data.frame(sample = samples, setNames(as.data.frame(round(cells, 5)), paste0("lineage_", colnames(cells)))),
            file.path(out, "cell_fractions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = samples, ctrl_PC1 = round(rnorm(n), 4), ctrl_PC2 = round(rnorm(n), 4)),
            file.path(out, "ctrl_pcs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("fixture written\n")
