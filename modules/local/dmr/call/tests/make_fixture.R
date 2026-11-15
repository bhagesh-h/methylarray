# Builds the simulated test fixture in this folder: 40 arrays x 2,000 CpGs at
# 100 bp spacing on four chromosomes, with a planted region (CpGs 101-115)
# associated with 'score'; stats.tsv.gz holds the limma results in the columns
# of EWAS_LIMMA, resolved.json the model with two arrays left out. The values are
# random and describe no study.
#   Rscript make_fixture.R <this folder>
suppressPackageStartupMessages({ library(SummarizedExperiment); library(HDF5Array); library(limma); library(data.table) })
out <- commandArgs(TRUE)[1]
set.seed(7)
n <- 40; p <- 2000
samples <- sprintf("S%02d", seq_len(n))
score <- rnorm(n)
ids <- sprintf("cg%08d", seq_len(p))
chr <- paste0("chr", rep(1:4, each = 500))
pos <- rep(seq(10000, by = 100, length.out = 500), 4)
base <- rnorm(p, 0, 2)
# neighbouring CpGs share a common component, as in real data
shared <- matrix(rnorm(n * p / 5, sd = 0.3), p / 5, n)[rep(seq_len(p / 5), each = 5), ]
mvalue <- base + shared + matrix(rnorm(n * p, sd = 0.3), p, n)
dimnames(mvalue) <- list(ids, samples)
mvalue[101:115, ] <- mvalue[101:115, ] + outer(rep(0.6, 15), score)
beta <- 2^mvalue / (1 + 2^mvalue)
se <- SummarizedExperiment(assays = list(beta = beta, mvalue = mvalue),
                           colData = DataFrame(sample = samples, row.names = samples),
                           rowData = DataFrame(cpg_id = ids, chr = chr, pos = pos, row.names = ids))
saveHDF5SummarizedExperiment(se, dir = file.path(out, "test_se"), replace = TRUE)
keep <- samples[-c(3, 17)]
design <- model.matrix(~ score[-c(3, 17)])
fit <- eBayes(lmFit(mvalue[, keep], design), trend = TRUE, robust = TRUE)
fit_b <- lmFit(beta[, keep], design)
stats <- data.table(cpg_id = ids, chr = chr, pos = pos, gene = ifelse(seq_len(p) %in% 95:120, "GENEA", NA), gene_group = NA, cgi_relation = NA,
                    contrast = "score", logfc_m = fit$coefficients[, 2], se_m = fit$stdev.unscaled[, 2] * sqrt(fit$s2.post),
                    delta_beta = fit_b$coefficients[, 2], t = fit$t[, 2], p_value = fit$p.value[, 2],
                    fdr = p.adjust(fit$p.value[, 2], "BH"), p_bacon = fit$p.value[, 2], fdr_bacon = p.adjust(fit$p.value[, 2], "BH"),
                    mean_beta = rowMeans(beta[, keep]), n = length(keep))
fwrite(stats[order(p_value)], file.path(out, "stats.tsv.gz"), sep = "\t")
writeLines('{"model": {"name": "score", "trait": "score", "trait_type": "continuous", "skip_dmr": false}, "covariates": [], "n": 38, "dropped": [{"sample": "S03", "reason": "missing score"}, {"sample": "S17", "reason": "missing score"}]}',
           file.path(out, "resolved.json"))
writeLines('{"model": {"name": "skipped", "trait": "score", "trait_type": "continuous", "skip_dmr": true}, "covariates": [], "n": 38, "dropped": []}',
           file.path(out, "resolved_skip.json"))
cat(sum(stats$fdr < 0.05), "CpGs at FDR < 0.05\n")
