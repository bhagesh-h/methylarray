# Builds the simulated test fixture in this folder: EWAS statistics of 3,000 CpGs
# in 600 real human genes (5 CpGs each, gene symbols from org.Hs.eg.db), where
# every CpG of 30 planted genes and 20 random CpGs are significant, and a GMT
# file with the planted set and two control sets (Entrez IDs). The statistics
# are random and describe no study.
#   Rscript make_fixture.R <this folder>
suppressPackageStartupMessages({ library(org.Hs.eg.db); library(data.table) })
out <- commandArgs(TRUE)[1]
set.seed(11)
eg <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID"),
                                             columns = c("ENTREZID", "SYMBOL", "GENETYPE"), keytype = "ENTREZID"))
eg <- eg[eg$GENETYPE == "protein-coding" & !is.na(eg$SYMBOL), ]
genes <- eg[sample(nrow(eg), 600), ]
cpg_gene <- rep(seq_len(600), each = 5)
n <- length(cpg_gene)
p <- runif(n)
planted <- cpg_gene %in% 1:30
p[planted] <- runif(sum(planted), 0, 1e-6)
p[sample(which(!planted), 20)] <- runif(20, 0, 1e-6)
stats <- data.table(cpg_id = sprintf("cg%08d", seq_len(n)), chr = "chr1", pos = seq_len(n) * 100L,
                    gene = genes$SYMBOL[cpg_gene], gene_group = "TSS200", cgi_relation = NA, contrast = "score",
                    logfc_m = 0, se_m = 1, delta_beta = 0, t = qnorm(p / 2, lower.tail = FALSE), p_value = p,
                    fdr = p.adjust(p, "BH"), p_bacon = p, fdr_bacon = p.adjust(p, "BH"), mean_beta = 0.5, n = 40L)
fwrite(stats, file.path(out, "stats.tsv.gz"), sep = "\t")
writeLines(c(paste(c("PLANTED", "planted genes", genes$ENTREZID[1:30]), collapse = "\t"),
             paste(c("CONTROL_A", "control genes", genes$ENTREZID[101:130]), collapse = "\t"),
             paste(c("CONTROL_B", "control genes", genes$ENTREZID[201:230]), collapse = "\t")),
           file.path(out, "sets.entrez.gmt"))
writeLines('{"model": {"name": "score", "trait": "score", "trait_type": "continuous", "skip_enrichment": false}, "covariates": [], "n": 40, "dropped": []}',
           file.path(out, "resolved.json"))
writeLines('{"model": {"name": "skipped", "trait": "score", "trait_type": "continuous", "skip_enrichment": true}, "covariates": [], "n": 40, "dropped": []}',
           file.path(out, "resolved_skip.json"))
cat(sum(stats$fdr < 0.05), "CpGs at FDR < 0.05\n")
