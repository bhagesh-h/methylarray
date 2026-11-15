#!/usr/bin/env Rscript
#
# One probe annotation table for EPIC v2, the single source of the probe flags
# used by the probe filter and later steps:
#   - Peters et al. 2024 extended manifest (Bioconductor EPICv2manifest,
#     AnnotationHub AH116484, staged as an .Rda file): position, design,
#     replicate probes and the best replicate, cross-hybridisation (BLAT
#     prediction and WGBS evidence), manifest mismatches, EPIC v1 and 450K
#     IDs, gene and CpG island context
#   - Zhou lab InfiniumAnnotation EPIC v2 mask (staged file): general mask,
#     mapping and non-unique probes, SNP masks by population and allele
#     frequency, poor titration
# The array's rs, nv and control probes are in neither source; the probe
# filter gives them their type from the ID prefix.
#

prefix <- "${prefix}"
manifest_file <- "${manifest}"
mask_file <- "${probe_mask}"

defaults <- list(platform = "EPICv2")
parse_args <- function(x, defaults) {
    tokens <- strsplit(trimws(x), "[[:space:]]+")[[1]]
    tokens <- tokens[nzchar(tokens)]
    if (length(tokens) %% 2 != 0) stop("ext.args must be pairs of '--name value': ", x)
    opts <- defaults
    for (i in seq_len(length(tokens) / 2) * 2 - 1) {
        name <- sub("^--", "", tokens[i])
        if (!name %in% names(defaults)) stop("Unknown option in ext.args: ", tokens[i])
        opts[[name]] <- tokens[i + 1]
    }
    opts
}
opts <- parse_args("${task.ext.args ?: ''}", defaults)
if (opts[["platform"]] != "EPICv2") stop("Probe annotation is available for EPICv2 only, not ", opts[["platform"]])

Sys.setenv(HOME = getwd())
suppressPackageStartupMessages(library(data.table))

blank_to_na <- function(x) {
    x <- as.character(x)
    x[!is.na(x) & !nzchar(x)] <- NA_character_
    x
}

# Peters et al. 2024 manifest: one data.frame named EPICv2manifest
loaded <- new.env()
if (!"EPICv2manifest" %in% load(manifest_file, envir = loaded)) stop(manifest_file, " does not contain the EPICv2manifest table")
manifest <- as.data.table(loaded[["EPICv2manifest"]])
rm(loaded)

anno <- manifest[, .(
    probe_id = IlmnID,
    cpg_id = Name,
    probe_type = Probe_Type,
    chr = blank_to_na(CHR),
    pos = as.integer(MAPINFO),
    strand = blank_to_na(Strand_FR),
    design_type = blank_to_na(Infinium_Design_Type),
    color_channel = blank_to_na(Color_Channel),
    replicate_group_size = .N,
    replicate_choice = blank_to_na(Rep_results_by_LOCATION),
    ch_blat = CH_BLAT %in% "Y",
    ch_wgbs = CH_WGBS_evidence %in% "Y",
    suggested_offtarget = blank_to_na(Suggested_offtarget),
    mismatch_pos = MismatchPos %in% "Y",
    epicv1_id = blank_to_na(EPICv1_Loci),
    hm450_id = blank_to_na(Methyl450_Loci),
    gene = blank_to_na(UCSC_RefGene_Name),
    gene_group = blank_to_na(UCSC_RefGene_Group),
    cgi_relation = blank_to_na(Relation_to_UCSC_CpG_Island)
), by = Name][, Name := NULL]
anno[chr %in% "chr0", c("chr", "pos") := list(NA_character_, NA_integer_)]

# Zhou lab mask: one row per masked probe, the mask names in maskUniq
# fread reads .gz files only with R.utils; base R decompresses instead
mask <- fread(text = readLines(gzfile(mask_file)), colClasses = "character", na.strings = c("", "NA"))
mask_sets <- strsplit(ifelse(is.na(mask[["maskUniq"]]), "", mask[["maskUniq"]]), ",")
has_mask <- function(names) vapply(mask_sets, function(set) any(set %in% names), logical(1))
snp_names <- function(set) {
    if (startsWith(set, "common")) {
        paste0(c("M_SNP", "M_1baseSwitchSNP", "M_2extBase_SNP"), set)
    } else {
        paste0(c("M_SNP_", "M_1baseSwitchSNP_", "M_2extBase_SNP_"), set, "_1pt")
    }
}
flags <- data.table(
    probe_id = mask[["Probe_ID"]],
    zhou_general = toupper(mask[["M_general"]]) %in% "TRUE",
    zhou_mapping = has_mask("M_mapping"),
    zhou_nonunique = has_mask("M_nonuniq"),
    zhou_uncorr_titration = has_mask("M_uncorr_titration")
)
snp_sets <- c(common_1pt = "common_1pt", common_5pt = "common_5pt", AFR = "AFR", AMR = "AMR", EAS = "EAS", EUR = "EUR", SAS = "SAS")
for (set in names(snp_sets)) {
    flags[, (paste0("snp_", tolower(set))) := has_mask(snp_names(snp_sets[[set]]))]
}

# Probes in the mask but not in the manifest (none in the current sources)
extra <- setdiff(flags[["probe_id"]], anno[["probe_id"]])
if (length(extra) > 0) {
    anno <- rbind(anno, data.table(probe_id = extra, cpg_id = sub("_.*", "", extra),
                                   probe_type = substr(extra, 1, 2)), fill = TRUE)
}
anno <- merge(anno, flags, by = "probe_id", all.x = TRUE, sort = FALSE)
flag_columns <- setdiff(names(flags), "probe_id")
for (column in flag_columns) set(anno, which(is.na(anno[[column]])), column, FALSE)
anno[is.na(replicate_group_size), replicate_group_size := 1L]
for (column in c("ch_blat", "ch_wgbs", "mismatch_pos")) set(anno, which(is.na(anno[[column]])), column, FALSE)
setorder(anno, probe_id)

fwrite(anno, paste0(prefix, ".tsv.gz"), sep = "\\t", na = "")

# Number of probes per flag
logical_columns <- names(anno)[vapply(anno, is.logical, logical(1))]
summary <- rbind(
    data.table(flag = "probes", n = nrow(anno)),
    data.table(flag = paste0("type_", sort(unique(anno[["probe_type"]]))),
               n = as.integer(table(anno[["probe_type"]])[sort(unique(anno[["probe_type"]]))])),
    data.table(flag = "no_position", n = sum(is.na(anno[["pos"]]))),
    data.table(flag = "replicate_probes", n = sum(anno[["replicate_group_size"]] > 1)),
    data.table(flag = logical_columns, n = vapply(logical_columns, function(column) sum(anno[[column]]), integer(1)))
)
fwrite(summary, paste0(prefix, "_summary.tsv"), sep = "\\t")

# Provenance of the two sources
resources <- data.table(
    resource = c("Peters et al. 2024 EPIC v2 manifest", "Zhou lab EPIC v2 mask"),
    file = c(basename(manifest_file), basename(mask_file)),
    md5 = unname(tools::md5sum(c(manifest_file, mask_file)))
)
fwrite(resources, paste0(prefix, "_resources.tsv"), sep = "\\t")

# Software versions
writeLines(c('"${task.process}":',
    paste0("    data.table: ", packageVersion("data.table"))), "versions.yml")
