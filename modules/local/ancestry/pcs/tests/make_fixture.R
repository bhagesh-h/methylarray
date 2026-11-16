# Builds the simulated fixtures in this folder (run in the module image):
#   Rscript make_fixture.R
# 60 arrays from two populations (A and B) with different allele frequencies, on
# three DNA plates. 300 autosomal CpG probes with a common SNP and 40 rs probes
# carry the genotypes; the SNP CpGs also carry a strong plate effect, which the
# control-probe PC of ctrl_pcs.tsv tracks. 5 SNP CpGs on chrX and 100 CpGs
# without a SNP must be left out; one SNP CpG fails in half of the arrays.
# Two arrays failed sample QC. Random values.
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(HDF5Array)
})
set.seed(7)
n <- 60
arrays <- sprintf("S%02d", seq_len(n))
population <- rep(c("A", "B"), length.out = n)
plate <- sprintf("P%d", (seq_len(n) - 1) %% 3 + 1)
plate_index <- (seq_len(n) - 1) %% 3

simulate <- function(m, slope, base) {
    p_a <- runif(m, 0.1, 0.9)
    p_b <- pmin(pmax(p_a + rnorm(m, 0, 0.3), 0.05), 0.95)
    g <- sapply(seq_len(n), function(j) rbinom(m, 2, if (population[j] == "A") p_a else p_b))
    pmin(pmax(base + slope * g + matrix(rnorm(m * n, 0, 0.03), m), 0.01), 0.99)
}
snp <- simulate(300, 0.3, 0.15)
snp <- snp + outer(rnorm(300, 0, 1), plate_index) * 0.08
snp <- pmin(pmax(snp, 0.01), 0.99)
snp_x <- simulate(5, 0.3, 0.15)
plain <- matrix(runif(100 * n, 0.05, 0.95), 100)
rs <- simulate(40, 0.45, 0.05)
beta <- rbind(snp, snp_x, plain, rs)
ids <- c(sprintf("cg%08d_BC21", 1:300), sprintf("cg%08d_BC21", 301:305), sprintf("cg%08d_TC21", 1:100 + 1000), sprintf("rs%d_BC11", 1:40))
dimnames(beta) <- list(ids, arrays)
masked <- beta
masked[sample(length(masked), 0.01 * length(masked))] <- NA
masked[1, seq(1, n, 2)] <- NA
unlink("test_se", recursive = TRUE)
saveHDF5SummarizedExperiment(SummarizedExperiment(assays = list(beta = masked, beta_unmasked = beta)), "test_se")

anno <- data.frame(probe_id = ids, probe_type = substr(ids, 1, 2),
                   chr = c(rep("chr1", 300), rep("chrX", 5), rep("chr2", 100), rep(NA, 40)),
                   snp_common_5pt = c(rep(TRUE, 305), rep(FALSE, 100), rep(FALSE, 40)))
gz <- gzfile("annotation.tsv.gz", "w")
write.table(anno, gz, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
close(gz)
write.table(data.frame(sample = arrays, plate = plate, sentrix_id = sprintf("2000%02d", (seq_len(n) - 1) %/% 8),
                       sentrix_position = sprintf("R0%dC01", (seq_len(n) - 1) %% 8 + 1),
                       reported_sex = rep(c("F", "F", "M"), length.out = n), reported_age = round(runif(n, 50, 80), 1)),
            "samples.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = arrays, status = c("FAIL", "FAIL", rep("PASS", n - 2)),
                       predicted_sex_minfi = ifelse(rep(c("F", "F", "M"), length.out = n) == "M", "male", "female")),
            "sample_qc.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = arrays, ctrl_PC1 = round(plate_index + rnorm(n, 0, 0.05), 4), ctrl_PC2 = round(rnorm(n), 4)),
            "ctrl_pcs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
fractions <- matrix(runif(n * 6), n)
fractions <- round(fractions / rowSums(fractions), 5)
colnames(fractions) <- paste0("lineage_", c("B", "CD4T", "CD8T", "NK", "Mono", "Gran"))
write.table(data.frame(sample = arrays, fractions), "cell_fractions.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = arrays, ancestry = population, bmi = round(rnorm(n, 28, 4), 1)),
            "phenotypes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
