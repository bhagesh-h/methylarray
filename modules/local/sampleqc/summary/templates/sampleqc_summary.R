#!/usr/bin/env Rscript
#
# Sample QC summary: joins the samplesheet, the preprocessing QC and the
# ewastools QC (when it ran) and gives every array a status with its reasons:
#   FAIL      detection failure, sex mismatch or SNP contamination
#   FLAG      control metrics below Illumina's thresholds, sex predictors that
#             disagree or fall between the clusters, genotype conflicts
#   PASS      none of the above
#   EXCLUDED  'exclude' set in the samplesheet (not preprocessed)
# No array is dropped here; samples_retained.txt lists the PASS and FLAG arrays.
#

prefix <- "${prefix}"
samples_file <- "${samples}"
preprocessing_qc_file <- "${preprocessing_qc}"
ewastools_qc_file <- "${ewastools_qc}"
conflicts_file <- "${snp_conflicts}"

defaults <- list(method = "sesame", max_undetected_fraction = "0.1", snp_outlier_threshold = "-4")
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
method <- opts[["method"]]
max_undetected <- as.numeric(opts[["max_undetected_fraction"]])
snp_outlier_threshold <- as.numeric(opts[["snp_outlier_threshold"]])

read_tsv <- function(file) {
    read.delim(file, colClasses = "character", na.strings = c("", "NA"), check.names = FALSE)
}

samples <- read_tsv(samples_file)
samples <- samples[order(samples[["sample"]]), , drop = FALSE]
excluded <- tolower(samples[["exclude"]]) %in% "true"

pre <- read_tsv(preprocessing_qc_file)
detected_column <- intersect(c("frac_detected", "frac_dt"), names(pre))[1]
if (is.na(detected_column)) stop("The preprocessing QC table has no detected-fraction column")
pre <- pre[match(samples[["sample"]], pre[["sample"]]), , drop = FALSE]
frac_detected <- as.numeric(pre[[detected_column]])
pre_sex <- pre[["predicted_sex"]]

has_ewastools <- nzchar(ewastools_qc_file)
ew <- NULL
conflicts <- NULL
if (has_ewastools) {
    ew <- read_tsv(ewastools_qc_file)
    ew <- ew[match(samples[["sample"]], ew[["sample"]]), , drop = FALSE]
    conflicts <- read_tsv(conflicts_file)
}
ew_value <- function(column) if (is.null(ew) || !column %in% names(ew)) rep(NA_character_, nrow(samples)) else ew[[column]]
ew_sex <- ew_value("predicted_sex")
# predict_sex runs only with reported males and females; then NA means the X and
# Y intensities fall outside both clusters
ew_sex_ran <- any(!is.na(ew_sex))
snp_outlier <- as.numeric(ew_value("snp_outlier"))
reported_sex <- c(M = "male", F = "female")[samples[["reported_sex"]]]

statuses <- character(nrow(samples))
reasons <- character(nrow(samples))
sex_check <- character(nrow(samples))
for (i in seq_len(nrow(samples))) {
    fail <- character()
    flag <- character()

    # Sex: every available prediction against the reported sex
    predictions <- c(pre_sex[i], ew_sex[i])
    predictions <- predictions[!is.na(predictions)]
    sex_check[i] <- if (is.na(reported_sex[i])) {
        "not_reported"
    } else if (length(predictions) == 0) {
        "not_predicted"
    } else if (all(predictions == reported_sex[i])) {
        "match"
    } else if (all(predictions != reported_sex[i])) {
        "mismatch"
    } else {
        "inconsistent"
    }

    if (excluded[i]) {
        statuses[i] <- "EXCLUDED"
        reasons[i] <- "excluded in the samplesheet"
        sex_check[i] <- NA
        next
    }

    undetected <- 1 - frac_detected[i]
    if (!is.na(undetected) && undetected > max_undetected) {
        fail <- c(fail, sprintf("%.1f%% of probes undetected (maximum %.1f%%)", 100 * undetected, 100 * max_undetected))
    }
    if (sex_check[i] == "mismatch") {
        fail <- c(fail, sprintf("predicted sex %s, reported %s", predictions[1], reported_sex[i]))
    }
    if (sex_check[i] == "inconsistent") {
        flag <- c(flag, sprintf("sex predictions disagree (%s %s, ewastools %s)", method, pre_sex[i], ew_sex[i]))
    }
    if (ew_sex_ran && is.na(ew_sex[i])) {
        flag <- c(flag, "X and Y intensities outside the male and female clusters")
    }
    if (!is.na(snp_outlier[i]) && snp_outlier[i] > snp_outlier_threshold) {
        fail <- c(fail, sprintf("SNP outlier score %.2f above %g (possible contamination)", snp_outlier[i], snp_outlier_threshold))
    }
    if (has_ewastools && ew_value("control_failed")[i] %in% "TRUE") {
        flag <- c(flag, paste0("control metrics below Illumina thresholds: ", gsub(";", ", ", ew_value("failed_controls")[i])))
    }
    if (!is.null(conflicts) && nrow(conflicts) > 0) {
        mine <- conflicts[conflicts[["sample1"]] == samples[["sample"]][i] | conflicts[["sample2"]] == samples[["sample"]][i], , drop = FALSE]
        for (j in seq_len(nrow(mine))) {
            other <- setdiff(c(mine[["sample1"]][j], mine[["sample2"]][j]), samples[["sample"]][i])
            same_donor <- mine[["expected_donor1"]][j] == mine[["expected_donor2"]][j]
            flag <- c(flag, sprintf(
                if (same_donor) "genotypes differ from %s in the same replicate group (agreement %s)"
                else "genotypes match %s from a different donor (agreement %s)",
                other, mine[["agreement"]][j]))
        }
    }

    statuses[i] <- if (length(fail) > 0) "FAIL" else if (length(flag) > 0) "FLAG" else "PASS"
    reasons[i] <- paste(c(fail, flag), collapse = "; ")
}

qc <- data.frame(
    sample = samples[["sample"]],
    status = statuses,
    reasons = reasons,
    frac_detected = round(frac_detected, 4),
    reported_sex = unname(reported_sex),
    predicted_sex = pre_sex,
    predicted_sex_ewastools = ew_sex,
    sex_check = sex_check,
    sex_x = ew_value("sex_x"),
    sex_y = ew_value("sex_y"),
    snp_outlier = snp_outlier,
    inferred_donor = ew_value("inferred_donor"),
    relatives = ew_value("relatives"),
    replicate_group = samples[["replicate_group"]],
    failed_controls = ew_value("failed_controls"),
    check.names = FALSE
)
names(qc)[names(qc) == "predicted_sex"] <- paste0("predicted_sex_", method)
# The Illumina control metrics of ewastools, for the PCA associations
if (!is.null(ew)) for (column in grep("^ctrl_", names(ew), value = TRUE)) qc[[column]] <- ew[[column]]
write.table(qc, paste0(prefix, ".tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, na = "")
writeLines(qc[["sample"]][qc[["status"]] %in% c("PASS", "FLAG")], paste0(prefix, "_retained.txt"))

# MultiQC general statistics
mqc <- data.frame(Sample = qc[["sample"]], qc_status = qc[["status"]], sex_check = qc[["sex_check"]])
pconfig <- c(
    "# id: 'sample_qc'",
    "# section_name: 'Sample QC'",
    "# plot_type: 'generalstats'",
    "# pconfig:",
    "#     qc_status:",
    "#         title: 'QC status'",
    "#         description: 'PASS, FLAG, FAIL or EXCLUDED (see the Sample QC status table)'",
    "#     sex_check:",
    "#         title: 'Sex check'",
    "#         description: 'Predicted against reported sex'"
)
if (has_ewastools) {
    mqc[["snp_outlier"]] <- qc[["snp_outlier"]]
    pconfig <- c(pconfig,
        "#     snp_outlier:",
        "#         title: 'SNP outlier'",
        paste0("#         description: 'ewastools SNP outlier score; above ", snp_outlier_threshold, " suggests contamination'"),
        "#         format: '{:,.2f}'")
}
writeLines(pconfig, paste0(prefix, "_mqc.tsv"))
write.table(mqc, paste0(prefix, "_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE, na = "")

# MultiQC table with the status and reasons
writeLines(c(
    "# id: 'sample_qc_status'",
    "# section_name: 'Sample QC status'",
    "# description: 'Status of every array with the reasons for FAIL and FLAG. Arrays marked FAIL or EXCLUDED do not go on to the analysis.'",
    "# plot_type: 'table'"
), paste0(prefix, "_status_mqc.tsv"))
write.table(data.frame(Sample = qc[["sample"]], Status = qc[["status"]], Reasons = qc[["reasons"]]),
    paste0(prefix, "_status_mqc.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE, append = TRUE)

# MultiQC scatter plot of the ewastools X and Y intensities
if (has_ewastools && any(!is.na(qc[["sex_x"]]))) {
    colours <- c(female = "#d95f02", male = "#1b9e77")
    points <- list()
    for (i in which(!is.na(qc[["sex_x"]]))) {
        colour <- if (is.na(qc[["reported_sex"]][i])) "#999999" else colours[[qc[["reported_sex"]][i]]]
        points[[qc[["sample"]][i]]] <- list(x = as.numeric(qc[["sex_x"]][i]), y = as.numeric(qc[["sex_y"]][i]), color = colour)
    }
    jsonlite::write_json(list(
        id = "sample_qc_sex",
        section_name = "Sex check",
        description = "Normalised X and Y chromosome intensities from ewastools, coloured by reported sex (orange female, green male, grey not reported). Females have high X and low Y intensities, males the reverse.",
        plot_type = "scatter",
        pconfig = list(id = "sample_qc_sex_plot", title = "Sample QC: X and Y intensities",
                       xlab = "Normalised X intensity", ylab = "Normalised Y intensity"),
        data = points
    ), paste0(prefix, "_sex_mqc.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
}

# Software versions
writeLines(c('"${task.process}":', paste0("    R: ", getRversion())), "versions.yml")
