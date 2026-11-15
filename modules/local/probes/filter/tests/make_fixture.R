# Builds the test fixture in this folder (test_se/, retained.txt, annotation.tsv.gz):
#   Rscript make_fixture.R <this folder>
# Each probe is built to hit one filter rule; S4 failed sample QC and is not retained.
suppressPackageStartupMessages({ library(SummarizedExperiment); library(HDF5Array) })
out <- commandArgs(TRUE)[1]
probes <- c("cg00000001_TC11", "cg00000002_TC11", "cg00000002_TC12", "rs00000003_TC11", "ch.1.1000F_BC11",
            "cg00000005_TC11", "cg00000006_TC11", "cg00000007_TC11", "cg00000008_TC11", "cg00000009_TC11",
            "cg00000010_TC11", "cg00000011_TC11", "cg00000012_TC11", "cg00000013_TC11", "cg00000014_TC11",
            "cg00000015_TC11")
samples <- c("S1", "S2", "S3", "S4")
beta <- matrix(rep(c(0.1, 0.2, 0.3, 0.4), each = length(probes)) + seq(0, 0.3, length.out = length(probes)),
               nrow = length(probes), dimnames = list(probes, samples))
beta["cg00000002_TC12", ] <- beta["cg00000002_TC11", ] + 0.1
detp <- matrix(0.001, nrow = length(probes), ncol = 4, dimnames = list(probes, samples))
detp["cg00000007_TC11", c("S1", "S2")] <- 0.2      # detection failure in 2 of 3 retained arrays
detp["cg00000015_TC11", c("S1", "S2")] <- 0.2      # chrX and detection failure
detp["cg00000001_TC11", "S4"] <- 0.9               # S4 failed sample QC: ignored
beta_unmasked <- beta
beta[detp > 0.05] <- NA
beta["cg00000008_TC11", "S3"] <- NA                # bead-count mask in 1 of 3 retained arrays
se <- SummarizedExperiment(
    assays = list(beta = beta, beta_unmasked = beta_unmasked, detp = detp),
    colData = DataFrame(sample = samples, row.names = samples),
    rowData = DataFrame(probe_id = probes, cpg_id = sub("_.*", "", probes), row.names = probes)
)
saveHDF5SummarizedExperiment(se, dir = file.path(out, "test_se"), replace = TRUE)
writeLines(c("S1", "S2", "S3"), file.path(out, "retained.txt"))
anno <- data.frame(
    probe_id = probes, cpg_id = sub("_.*", "", probes),
    probe_type = substr(probes, 1, 2),
    chr = c("chr1", "chr2", "chr2", "chr3", "chr1", NA, "chrX", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chrX"),
    pos = c(100L, 200L, 200L, 300L, 1000L, NA, 400L, 500L, 600L, 700L, 800L, 900L, 1000L, 1100L, 1200L, 1300L),
    replicate_group_size = c(1L, 2L, 2L, rep(1L, 13)),
    ch_wgbs = probes == "cg00000011_TC11",
    ch_blat = probes %in% c("cg00000011_TC11", "cg00000012_TC11"),
    mismatch_pos = probes == "cg00000013_TC11",
    zhou_general = probes == "cg00000010_TC11",
    zhou_uncorr_titration = FALSE,
    snp_common_1pt = probes == "cg00000009_TC11",
    snp_eur = probes %in% c("cg00000009_TC11", "cg00000014_TC11")
)
con <- gzfile(file.path(out, "annotation.tsv.gz"), "w")
write.table(anno, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
close(con)
