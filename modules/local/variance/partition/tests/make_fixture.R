# Builds the simulated fixtures in this folder (run in the module image):
#   Rscript make_fixture.R
# 48 arrays on four DNA plates and six BeadChips, half of them men, with ages
# and five cell lineages. 600 CpGs (M-values): a plate effect in the first 200,
# a sex effect in the next 200 and an age effect in the last 200, each about
# 50% of the variance. Random values.
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(HDF5Array)
})
set.seed(5)
n <- 48
arrays <- sprintf("S%02d", seq_len(n))
plate <- sprintf("P%d", (seq_len(n) - 1) %% 4 + 1)
chip <- sprintf("2000%02d", (seq_len(n) - 1) %/% 8)
row <- sprintf("R0%dC01", (seq_len(n) - 1) %% 8 + 1)
sex <- rep(c("F", "M"), length.out = n)
age <- round(runif(n, 50, 80), 1)
plate_effect <- c(P1 = -1, P2 = -1, P3 = 1, P4 = 1)[plate]  # balanced against sex
noise <- matrix(rnorm(600 * n), 600)
m <- noise
m[1:200, ] <- noise[1:200, ] + rep(plate_effect, each = 200)
m[201:400, ] <- noise[201:400, ] + rep(ifelse(sex == "M", 1, -1), each = 200)
m[401:600, ] <- noise[401:600, ] + rep(as.numeric(scale(age)), each = 200)
dimnames(m) <- list(sprintf("cg%08d", 1:600), arrays)
unlink("test_se", recursive = TRUE)
saveHDF5SummarizedExperiment(SummarizedExperiment(assays = list(mvalue = m)), "test_se")
write.table(data.frame(sample = arrays, plate = plate, sentrix_id = chip, sentrix_position = row, reported_sex = sex, reported_age = age),
            "samples.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = arrays, status = "PASS", predicted_sex_minfi = ifelse(sex == "M", "male", "female")),
            "sample_qc.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
fractions <- matrix(runif(n * 6), n)
fractions <- round(fractions / rowSums(fractions), 5)
colnames(fractions) <- paste0("lineage_", c("B", "CD4T", "CD8T", "NK", "Mono", "Gran"))
write.table(data.frame(sample = arrays, fractions), "cell_fractions.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
