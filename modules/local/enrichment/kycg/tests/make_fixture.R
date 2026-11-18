# Builds the simulated fixtures in this folder (run in the module image, with
# the ChromHMM knowledgebase in the ExperimentHub cache of EXPERIMENT_HUB_CACHE):
#   Rscript make_fixture.R
# 3,000 tested CpGs: 300 in the active TSS state (TssA), the rest from other
# states. 100 CpGs are significant for 'score': 80 TssA CpGs with lower
# methylation and 20 others with higher methylation. Random values.
suppressPackageStartupMessages({
    library(knowYourCG)
    library(data.table)
})
set.seed(3)
dbs <- getDBs("KYCG.EPICv2.ChromHMM.20220303", silent = TRUE)
tss <- unique(dbs[["TssA"]])
other <- setdiff(unique(unlist(dbs[setdiff(names(dbs), c("TssA", "TssFlnk", "TssFlnkU", "TssFlnkD"))])), tss)
probes <- c(sample(tss, 300), sample(other, 2700))
cpgs <- sub("_.*", "", probes)
keep <- !duplicated(cpgs)
probes <- probes[keep]
cpgs <- cpgs[keep]
in_tss <- probes %in% tss
sig <- c(sample(which(in_tss), 80), sample(which(!in_tss), 20))
p <- runif(length(cpgs), 0.01, 1)
p[sig] <- 10^-runif(length(sig), 8, 12)
delta <- rnorm(length(cpgs), 0, 0.005)
delta[sig] <- ifelse(in_tss[sig], -0.05, 0.05)
stats <- data.table(cpg_id = cpgs, contrast = "score", delta_beta = round(delta, 4), p_value = signif(p, 4), fdr = signif(p.adjust(p, "BH"), 4))
gz <- gzfile("stats.tsv.gz", "w")
write.table(stats, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)
gz <- gzfile("annotation.tsv.gz", "w")
write.table(data.table(probe_id = probes, cpg_id = cpgs), gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)
