//
// Subworkflow with functionality specific to the bhagesh-h/methylarray pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    phenotypes        //  string: Path to the phenotype table (optional)
    models            //  string: Path to the YAML model definitions (optional)
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Work directory cleanup: published results must not point into the work
    // directory that is deleted at the end of a successful run
    //
    if (!params.skip_cleanup && params.publish_dir_mode in ['symlink', 'rellink']) {
        error("--publish_dir_mode ${params.publish_dir_mode} links the results to the work directory, which is deleted after a successful run. Use --publish_dir_mode copy (the default) or link, or keep the work directory with --skip_cleanup.")
    }

    //
    // Create channel from input file provided through params.input:
    // one element per array, [ meta, idat_red, idat_grn ]
    //
    def samplesheet_rows = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
    validateUniqueIdats(samplesheet_rows)
    def validated_rows = samplesheet_rows.collect { meta, idat_red, idat_grn -> validateIdatPair(meta, idat_red, idat_grn) }
    ch_samplesheet = channel.fromList(validated_rows)

    //
    // Phenotype table and model definitions: check that they agree with the
    // samplesheet and with each other before any task starts
    //
    def phenotype_columns = phenotypes ? validatePhenotypes(phenotypes, validated_rows) : []
    def model_list = models ? validateModels(models, validated_rows, phenotype_columns) : []
    ch_phenotypes = channel.value(phenotypes ? file(phenotypes, checkIfExists: true) : [])
    ch_models = channel.fromList(model_list)

    emit:
    samplesheet = ch_samplesheet // channel: [ meta, idat_red, idat_grn ]
    phenotypes  = ch_phenotypes  // value channel: phenotype table, or [] when not given
    models      = ch_models      // channel: one map per model, defaults filled in
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.isBound() ? multiqc_reports.getVal() : [],
            )
        }

        completionSummary(monochrome_logs)
        cleanupWorkDirectory()

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// After a successful run: Nextflow removes the task directories of the run
// (cleanup = true); remove also the run's cache of downloaded remote inputs
// (work/stage-<session>), which Nextflow keeps
//
def cleanupWorkDirectory() {
    if (!workflow.success || params.skip_cleanup) return
    def stage = workflow.workDir.resolve("stage-${workflow.sessionId}")
    if (stage.exists()) stage.deleteDir()
    log.info "The task directories and downloaded inputs of this run in ${workflow.workDir} are being removed (Nextflow cleanup). Use --skip_cleanup to keep them for -resume."
}

//
// Validate one samplesheet row and add the Sentrix ID, position and row to meta
//
def validateIdatPair(meta_in, idat_red, idat_grn) {
    // nf-schema stores missing optional values as []: use null instead, and keep
    // numeric-looking identifiers (read as integers) as strings
    def cleaned = meta_in.collectEntries { key, value -> [ key, (value == [] || value == '') ? null : value ] }
    def meta = cleaned + [
        sentrix_id: cleaned.sentrix_id?.toString(),
        plate: cleaned.plate?.toString(),
        batch: cleaned.batch?.toString(),
        replicate_group: cleaned.replicate_group?.toString()
    ]

    // Both files must belong to the same array: identical prefix before _Red/_Grn
    def red_prefix = idat_red.name.replaceFirst(/_Red\.idat(\.gz)?$/, '')
    def grn_prefix = idat_grn.name.replaceFirst(/_Grn\.idat(\.gz)?$/, '')
    if (red_prefix != grn_prefix) {
        error("Samplesheet: the IDAT files of sample '${meta.id}' do not belong to the same array (${idat_red.name}, ${idat_grn.name})")
    }

    // Sentrix ID and position from the standard file name <...>_<sentrix_id>_<position>
    def tokens = red_prefix.tokenize('_')
    def file_sentrix_id = (tokens.size() >= 2 && tokens[-2] ==~ /[0-9]{10,14}/) ? tokens[-2] : null
    def file_position = (tokens.size() >= 1 && tokens[-1] ==~ /R[0-9]{2}C[0-9]{2}/) ? tokens[-1] : null
    if (meta.sentrix_id && file_sentrix_id && meta.sentrix_id != file_sentrix_id) {
        error("Samplesheet: sample '${meta.id}' has sentrix_id ${meta.sentrix_id} but its IDAT files name ${file_sentrix_id}")
    }
    if (meta.sentrix_position && file_position && meta.sentrix_position != file_position) {
        error("Samplesheet: sample '${meta.id}' has sentrix_position ${meta.sentrix_position} but its IDAT files name ${file_position}")
    }
    def sentrix_id = meta.sentrix_id ?: file_sentrix_id
    def sentrix_position = meta.sentrix_position ?: file_position
    def sentrix_row = sentrix_position ? sentrix_position.substring(0, 3) : null
    def sentrix_col = sentrix_position ? sentrix_position.substring(3) : null

    def new_meta = meta + [
        sentrix_id: sentrix_id,
        sentrix_position: sentrix_position,
        sentrix_row: sentrix_row,
        sentrix_col: sentrix_col,
        exclude: meta.exclude ?: false
    ]
    return [ new_meta, idat_red, idat_grn ]
}

//
// Each IDAT file may appear in only one samplesheet row
//
def validateUniqueIdats(samplesheet) {
    def files = samplesheet.collectMany { row -> [ row[1].toUriString(), row[2].toUriString() ] }
    def duplicated = files.countBy { f -> f }.findAll { _f, n -> n > 1 }.keySet()
    if (duplicated) {
        error("Samplesheet: these IDAT files appear in more than one row: ${duplicated.join(', ')}")
    }
}

//
// Covariates that the pipeline derives from the data. Phenotype columns may not
// reuse these names, and models may refer to them without a phenotype table.
//
def isDerivedCovariate(name) {
    return name in ['predicted_sex', 'dnam_smoking'] || name ==~ /cc_[A-Za-z0-9_]+/ || name ==~ /ctrlpc[0-9]+/ || name ==~ /ancestry_pc[0-9]+/
}

//
// Columns behind a model covariate entry: a token with a fixed expansion, or the column itself
//
def tokenExpansion(covariate) {
    def tokens = [ '@technical': ['plate', 'sentrix_row'], '@sex': ['predicted_sex'], '@smoking': ['dnam_smoking'] ]
    return tokens.containsKey(covariate) ? tokens[covariate] : [ covariate ]
}

//
// Samplesheet columns that models may use: every technical column with at least
// one value, plus the Sentrix row and column derived from the file names
//
def samplesheetColumns(samplesheet_rows) {
    def columns = samplesheet_rows.collectMany { row -> row[0].findAll { _key, value -> value != null }.keySet() } as Set
    return columns - ['id', 'exclude', 'exclude_reason']
}

//
// Phenotype table: one row per sample, column names that cannot be confused with
// samplesheet or derived columns. Only column names and counts are logged, never values.
// Returns the phenotype column names (without 'sample').
//
def validatePhenotypes(phenotypes, samplesheet_rows) {
    def path = file(phenotypes, checkIfExists: true)
    def separator = path.name.toLowerCase().endsWith('.tsv') ? '\t' : ','
    def rows = path.splitCsv(header: true, sep: separator)
    if (!rows) {
        error("Phenotype table: ${path.name} has no data rows")
    }
    def columns = (rows[0].keySet() as List) - ['sample']

    def reserved = samplesheet_rows[0][0].keySet() + ['sample']
    def clashes = columns.findAll { column -> column in reserved || isDerivedCovariate(column) }
    if (clashes) {
        error("Phenotype table: these columns reuse samplesheet or derived covariate names; rename them: ${clashes.join(', ')}")
    }

    def phenotype_ids = rows.collect { row -> row.sample }
    def duplicated = phenotype_ids.countBy { id -> id }.findAll { _id, n -> n > 1 }.keySet()
    if (duplicated) {
        error("Phenotype table: these samples appear in more than one row: ${duplicated.join(', ')}")
    }

    def sheet_ids = samplesheet_rows.collect { row -> row[0].id } as Set
    def shared = sheet_ids.intersect(phenotype_ids as Set)
    if (!shared) {
        error("Phenotype table: none of its samples are in the samplesheet")
    }
    def unknown = phenotype_ids.findAll { id -> !(id in sheet_ids) }
    if (unknown) {
        def shown = unknown.take(10).join(', ') + (unknown.size() > 10 ? ', ...' : '')
        log.warn("Phenotype table: ${unknown.size()} sample(s) are not in the samplesheet and are ignored: ${shown}")
    }
    def without = sheet_ids.size() - shared.size()
    if (without) {
        log.info("Phenotype table: ${without} of ${sheet_ids.size()} samplesheet samples have no phenotype row. They are processed, but left out of models that use phenotype columns.")
    }
    log.info("Phenotype table: ${shared.size()} samples, ${columns.size()} columns")
    return columns
}

//
// Model definitions: every column a model names must exist in the samplesheet,
// the phenotype table or the derived covariates. Returns the models with defaults filled in.
//
def validateModels(models, samplesheet_rows, phenotype_columns) {
    def path = file(models, checkIfExists: true)
    def entries = new org.yaml.snakeyaml.Yaml().load(path.text)
    if (!(entries instanceof List) || !entries) {
        error("Models: ${path.name} must be a non-empty YAML list of models")
    }

    def available = samplesheetColumns(samplesheet_rows) + phenotype_columns
    def defaults = [ transform: 'none', covariates: [], min_n: 50, skip_dmr: false, skip_enrichment: false ]
    def problems = []
    def resolved = entries.collect { entry ->
        def model = defaults + entry
        def named = [ model.trait ] + model.covariates.findAll { covariate -> !covariate.startsWith('@') } + [ model.block, model.subset?.column ]
        def missing = named.findAll { column -> column != null && !(column in available) && !isDerivedCovariate(column) }.unique()
        if (missing) {
            problems << "model '${model.name}' uses columns that are not in the samplesheet, the phenotype table or the derived covariates: ${missing.join(', ')}"
        }
        // Tokens with a fixed expansion; @cells, @ctrlpcs and @ancestry depend on the data
        def fixed = model.covariates.collectMany { covariate -> tokenExpansion(covariate) }
        if (model.trait in fixed) {
            problems << "model '${model.name}' lists its trait '${model.trait}' as a covariate"
        }
        if (model.block && model.block in fixed) {
            problems << "model '${model.name}' uses '${model.block}' both as block and as covariate"
        }
        return model
    }
    def duplicated = resolved.collect { model -> model.name }.countBy { name -> name }.findAll { _name, n -> n > 1 }.keySet()
    if (duplicated) {
        problems << "these model names are used more than once: ${duplicated.join(', ')}"
    }
    if (problems) {
        error("Models:\n  - " + problems.join('\n  - '))
    }
    log.info("Models: ${resolved.size()} defined (${resolved.collect { model -> model.name }.join(', ')})")
    return resolved
}

//
// Sample metadata written to the SummarizedExperiment colData, in this order
//
def sampleColumns() {
    return [
        'id', 'sentrix_id', 'sentrix_position', 'sentrix_row', 'sentrix_col', 'plate', 'well',
        'batch', 'reported_sex', 'reported_age', 'replicate_group'
    ]
}

//
// One TSV file with a row per array ('sample' plus sampleColumns); with
// with_idat_prefix, a last column holds the IDAT file prefix used by minfi
//
def sampleTable(ch_arrays, with_idat_prefix) {
    def columns = sampleColumns()
    def header = (['sample'] + columns.drop(1) + (with_idat_prefix ? ['idat_prefix'] : [])).join('\t')
    return ch_arrays
        .map { meta, idat_red, _idat_grn ->
            def values = columns.collect { column -> meta[column] ?: '' }
            if (with_idat_prefix) {
                values << idat_red.name.replaceFirst(/_Red\.idat(\.gz)?$/, '')
            }
            values.join('\t')
        }
        .collectFile(name: 'samples.tsv', seed: header, newLine: true, sort: true)
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def norm_text = [
        noob: "Noob background and dye-bias correction (Triche et al. 2013)",
        funnorm: "functional normalisation (Fortin et al. 2014)",
        quantile: "stratified quantile normalisation (Touleimat and Tost 2012)",
        swan: "SWAN normalisation (Maksimovic et al. 2012) of the Noob-corrected intensities (Triche et al. 2013)"
    ]
    def enmix_text = "ENmix out-of-band background correction (Xu et al. 2016) with RELIC dye-bias correction (Xu et al. 2017)${params.enmix_quantile ? ', quantile normalisation' : ''} and RCP probe-type correction (Niu et al. 2016)"
    def masks = "values with a detection p-value above ${params.detection_pval} or fewer than ${params.min_beads} beads were set to missing"
    def preprocessing = [
        sesame: "Arrays were preprocessed per sample with SeSAMe (Zhou et al. 2018) using prepSesame code ${params.sesame_prep}, with pOOBAH detection masking at p > ${params.detection_pval}, noob background correction (Triche et al. 2013) and probes with fewer than ${params.min_beads} beads masked.",
        minfi: "Arrays were read with minfi (Aryee et al. 2014; Fortin et al. 2017) and preprocessed with ${norm_text[params.minfi_norm]}; ${masks}.",
        enmix: "Arrays were read with minfi (Aryee et al. 2014) and preprocessed with ${enmix_text}; ${masks}.",
        combined: "Arrays were read with minfi (Aryee et al. 2014) and preprocessed with ${enmix_text}. Values masked by SeSAMe (Zhou et al. 2018; quality mask and pOOBAH detection at p > ${params.detection_pval}) or measured by fewer than ${params.min_beads} beads were set to missing."
    ]
    def sample_qc = params.skip_ewastools
        ? "Arrays with more than ${String.format("%.1f", (params.max_undetected_fraction as double) * 100)}% undetected probes or a predicted sex that contradicted the reported sex failed quality control."
        : "Control-probe metrics, normalised X and Y intensities and SNP probe genotypes were computed with ewastools (Heiss and Just 2018). Arrays with more than ${String.format("%.1f", (params.max_undetected_fraction as double) * 100)}% undetected probes, a predicted sex that contradicted the reported sex or a SNP outlier score above ${params.snp_outlier_threshold} failed quality control; arrays with control metrics below Illumina's thresholds or genotype conflicts with the expected donors were flagged."
    def crosshyb = [ wgbs: "cross-hybridising probes with WGBS evidence of off-target binding (Peters et al. 2024)", blat: "probes with a predicted off-target BLAT hit (Peters et al. 2024)", none: null ][params.crosshyb_filter]
    def snp = params.snp_mask == 'none' ? null : "probes with a SNP (${params.snp_mask.startsWith('common') ? 'global minor allele frequency of at least ' + params.snp_mask.replace('common_', '').replace('pt', '%') : 'allele frequency of at least 1% in the ' + params.snp_mask.toUpperCase() + ' population'}) in the probe body, at the CpG or at the extension base"
    def removed = [ "non-CG and genotyping probes", "probes without an hg38 position", params.keep_sex_chr ? null : "probes on the sex chromosomes",
        "probes that failed detection in more than ${String.format('%.0f', (params.probe_fail_fraction as double) * 100)}% or other masks in more than ${String.format('%.0f', (params.bead_fail_fraction as double) * 100)}% of arrays",
        snp, "probes in the general mask of Zhou et al. (2017)", crosshyb ].findAll { item -> item }
    def probe_filter = "Probes were annotated with the extended EPIC v2 manifest (Peters et al. 2024) and the InfiniumAnnotation masks (Zhou et al. 2017). The following were removed: ${removed.join('; ')}.${params.replicate_probe_strategy == 'mean' ? ' Replicate probes were averaged per CpG.' : ''} Association tests use M-values."
    def cells = params.tissue == 'blood' ? "Leukocyte fractions were estimated with EpiDISH robust partial correlations (Teschendorff et al. 2017) using the 12-cell-type reference of Luo et al. (2023)." : ''
    def clocks = params.skip_clocks ? '' : "Epigenetic clocks were computed with dnaMethyAge (Wang et al. 2024)${params.clocks.tokenize(',').contains('DunedinPACE') ? ' and DunedinPACE (Belsky et al. 2022)' : ''} on the detection-masked beta values averaged per CpG; values that failed detection were set to the CpG mean over the other arrays."
    def exploration = "Technical variation was summarised with principal components of the control-probe intensities (Lehne et al. 2015), and the principal components of the most variable CpGs were tested for association with technical and biological variables${params.n_pcs == 'auto' ? ', their number estimated by random matrix theory (isva, Teschendorff et al. 2011)' : ''}."
    def ewas = params.models ? "Each association model was fitted at every CpG with limma (Ritchie et al. 2015) on M-values with eBayes moderation (trend and robust), effect sizes were taken from the same model on beta values, and test statistics were corrected for bias and inflation with bacon (van Iterson et al. 2017). CpGs were reported at FDR < 0.05 and at P < ${params.p_epigenome_wide} (Mansell et al. 2019); across several models, this threshold was divided by the effective number of independent tests, estimated from the correlation of the model statistics (Galwey 2009)." : ''
    def cell_ewas = params.models && params.tissue == 'blood' && !params.skip_cell_ewas ? "Cell-type-specific associations were tested with the CellDMC design (Zheng et al. 2018): beta values on the six leukocyte lineages, the trait times each lineage and the other covariates, fitted with limma, with an ordinary t-test of each interaction and FDR within each lineage." : ''
    def dmr = params.models && !params.skip_dmr ? "Differentially methylated regions were called with DMRcate (Peters et al. 2015; lambda ${params.dmr_lambda}, C ${params.dmr_c}) and dmrff (Suderman et al. 2018) from the same model statistics, and regions found by both methods in the same direction were reported as consensus regions." : ''
    def enrichment = params.models && !params.skip_enrichment ? "Gene set enrichment of the CpGs at FDR < 0.05 against all tested CpGs was tested with missMethyl (Maksimovic et al. 2021), which accounts for the number of CpGs per gene, for Gene Ontology and KEGG terms and the MSigDB gene set collections." : ''
    def scores = params.skip_dnam_scores ? '' : "DNA methylation scores for smoking, alcohol and BMI (McCartney et al. 2018) and IL-6 (Stevenson et al. 2021) were computed as weighted sums of the published CpG sets."
    // One point per step, so that the methods read as a list rather than one long paragraph
    def steps = [
        'Preprocessing': preprocessing[params.preprocessing] ?: '', 'Sample quality control': sample_qc, 'Probe filtering': probe_filter,
        'Cell composition': cells, 'Epigenetic clocks': clocks, 'DNA methylation scores': scores, 'Exploration': exploration,
        'Association models': ewas, 'Cell-type-specific associations': cell_ewas, 'Regions': dmr, 'Enrichment': enrichment,
        'Reporting': "Results were summarised with MultiQC (Ewels et al. 2016) and a Quarto report."
    ].findAll { _name, text -> text }
    return "<ul>" + steps.collect { name, text -> "<li><strong>${name}.</strong> ${text}</li>" }.join('') + "</ul>"
}

//
// One reference as an HTML list item for the MultiQC methods description
//
def referenceHtml(authors, year, title, journal, doi) {
    // APA style, shortened: one or two authors in full, the first author and "et al." for more
    def names = authors.split(/(?<=\.),\s+(?:&\s+)?/)
    def shown = names.size() <= 2 ? authors : "${names[0]}, et al."
    def journal_parts = journal.split(', ', 2)
    def source = journal_parts.size() == 2 ? "<em>${journal_parts[0]}</em>, ${journal_parts[1]}" : "<em>${journal}</em>"
    return "<li>${shown} (${year}). ${title}. ${source}. <a href=\"https://doi.org/${doi}\">https://doi.org/${doi}</a></li>"
}

def toolBibliographyText() {
    def refs = [
        sesame: referenceHtml('Zhou, W., Triche, T. J., Laird, P. W., & Shen, H.', 2018, 'SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions', 'Nucleic Acids Research, 46(20), e123', '10.1093/nar/gky691'),
        noob: referenceHtml('Triche, T. J., Weisenberger, D. J., Van Den Berg, D., Laird, P. W., & Siegmund, K. D.', 2013, 'Low-level processing of Illumina Infinium DNA Methylation BeadArrays', 'Nucleic Acids Research, 41(7), e90', '10.1093/nar/gkt090'),
        minfi: referenceHtml('Aryee, M. J., Jaffe, A. E., Corrada-Bravo, H., et al.', 2014, 'Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays', 'Bioinformatics, 30(10), 1363-1369', '10.1093/bioinformatics/btu049'),
        minfi_epic: referenceHtml('Fortin, J. P., Triche, T. J., & Hansen, K. D.', 2017, 'Preprocessing, normalization and integration of the Illumina HumanMethylationEPIC array with minfi', 'Bioinformatics, 33(4), 558-560', '10.1093/bioinformatics/btw691'),
        funnorm: referenceHtml('Fortin, J. P., Labbe, A., Lemire, M., et al.', 2014, 'Functional normalization of 450k methylation array data improves replication in large cancer studies', 'Genome Biology, 15, 503', '10.1186/s13059-014-0503-2'),
        quantile: referenceHtml('Touleimat, N., & Tost, J.', 2012, 'Complete pipeline for Infinium Human Methylation 450K BeadChip data processing using subset quantile normalization for accurate DNA methylation estimation', 'Epigenomics, 4(3), 325-341', '10.2217/epi.12.21'),
        swan: referenceHtml('Maksimovic, J., Gordon, L., & Oshlack, A.', 2012, 'SWAN: Subset-quantile within array normalization for Illumina Infinium HumanMethylation450 BeadChips', 'Genome Biology, 13(6), R44', '10.1186/gb-2012-13-6-r44'),
        enmix: referenceHtml('Xu, Z., Niu, L., Li, L., & Taylor, J. A.', 2016, 'ENmix: a novel background correction method for Illumina HumanMethylation450 BeadChip', 'Nucleic Acids Research, 44(3), e20', '10.1093/nar/gkv907'),
        relic: referenceHtml('Xu, Z., Langie, S. A., De Boever, P., Taylor, J. A., & Niu, L.', 2017, 'RELIC: a novel dye-bias correction method for Illumina Methylation BeadChip', 'BMC Genomics, 18, 4', '10.1186/s12864-016-3426-3'),
        rcp: referenceHtml('Niu, L., Xu, Z., & Taylor, J. A.', 2016, 'RCP: a novel probe design bias correction method for Illumina Methylation BeadChip', 'Bioinformatics, 32(17), 2659-2663', '10.1093/bioinformatics/btw285'),
        ewastools: referenceHtml('Heiss, J. A., & Just, A. C.', 2018, 'Identifying mislabeled and contaminated DNA methylation microarray data: an extended quality control toolset with examples from GEO', 'Clinical Epigenetics, 10, 73', '10.1186/s13148-018-0504-1'),
        peters2024: referenceHtml('Peters, T. J., Meyer, B., Ryan, L., et al.', 2024, 'Characterisation and reproducibility of the HumanMethylationEPIC v2.0 BeadChip for DNA methylation profiling', 'BMC Genomics, 25, 251', '10.1186/s12864-024-10027-5'),
        zhou2017: referenceHtml('Zhou, W., Laird, P. W., & Shen, H.', 2017, 'Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes', 'Nucleic Acids Research, 45(4), e22', '10.1093/nar/gkw967'),
        epidish: referenceHtml('Teschendorff, A. E., Breeze, C. E., Zheng, S. C., & Beck, S.', 2017, 'A comparison of reference-based algorithms for correcting cell-type heterogeneity in Epigenome-Wide Association Studies', 'BMC Bioinformatics, 18, 105', '10.1186/s12859-017-1511-5'),
        luo2023: referenceHtml('Luo, Q., Dwaraka, V. B., Chen, Q., et al.', 2023, 'A meta-analysis of immune-cell fractions at high resolution reveals novel associations with common phenotypes and health outcomes', 'Genome Medicine, 15, 59', '10.1186/s13073-023-01211-5'),
        dnamethyage: referenceHtml('Wang, Y., Grant, O. A., Zhai, X., McDonald-Maier, K. D., & Schalkwyk, L. C.', 2024, 'Insights into ageing rates comparison across tissues from recalibrating cerebellum DNA methylation clock', 'GeroScience, 46, 39-56', '10.1007/s11357-023-00871-w'),
        dunedinpace: referenceHtml('Belsky, D. W., Caspi, A., Corcoran, D. L., et al.', 2022, 'DunedinPACE, a DNA methylation biomarker of the pace of aging', 'eLife, 11, e73420', '10.7554/eLife.73420'),
        lehne2015: referenceHtml('Lehne, B., Drong, A. W., Loh, M., et al.', 2015, 'A coherent approach for analysis of the Illumina HumanMethylation450 BeadChip improves data quality and performance in epigenome-wide association studies', 'Genome Biology, 16, 37', '10.1186/s13059-015-0600-x'),
        rmt: referenceHtml('Teschendorff, A. E., Zhuang, J., & Widschwendter, M.', 2011, 'Independent surrogate variable analysis to deconvolve confounding factors in large-scale microarray profiling studies', 'Bioinformatics, 27(11), 1496-1505', '10.1093/bioinformatics/btr171'),
        celldmc: referenceHtml('Zheng, S. C., Breeze, C. E., Beck, S., & Teschendorff, A. E.', 2018, 'Identification of differentially methylated cell types in epigenome-wide association studies', 'Nature Methods, 15(12), 1059-1066', '10.1038/s41592-018-0213-x'),
        limma: referenceHtml('Ritchie, M. E., Phipson, B., Wu, D., et al.', 2015, 'limma powers differential expression analyses for RNA-sequencing and microarray studies', 'Nucleic Acids Research, 43(7), e47', '10.1093/nar/gkv007'),
        bacon: referenceHtml('van Iterson, M., van Zwet, E. W., the BIOS Consortium, & Heijmans, B. T.', 2017, 'Controlling bias and inflation in epigenome- and transcriptome-wide association studies using the empirical null distribution', 'Genome Biology, 18, 19', '10.1186/s13059-016-1131-9'),
        mansell2019: referenceHtml('Mansell, G., Gorrie-Stone, T. J., Bao, Y., et al.', 2019, 'Guidance for DNA methylation studies: statistical insights from the Illumina EPIC array', 'BMC Genomics, 20, 366', '10.1186/s12864-019-5761-7'),
        galwey2009: referenceHtml('Galwey, N. W.', 2009, 'A new measure of the effective number of tests, a practical tool for comparing families of non-independent significance tests', 'Genetic Epidemiology, 33(7), 559-568', '10.1002/gepi.20408'),
        dmrcate: referenceHtml('Peters, T. J., Buckley, M. J., Statham, A. L., et al.', 2015, 'De novo identification of differentially methylated regions in the human genome', 'Epigenetics & Chromatin, 8, 6', '10.1186/1756-8935-8-6'),
        dmrff: referenceHtml('Suderman, M., Staley, J. R., French, R., Arathimos, R., Simpkin, A., & Tilling, K.', 2018, 'dmrff: identifying differentially methylated regions efficiently with power and control', 'bioRxiv', '10.1101/508556'),
        missmethyl: referenceHtml('Maksimovic, J., Oshlack, A., & Phipson, B.', 2021, 'Gene set enrichment analysis for genome-wide DNA methylation data', 'Genome Biology, 22, 173', '10.1186/s13059-021-02388-x'),
        mccartney2018: referenceHtml('McCartney, D. L., Hillary, R. F., Stevenson, A. J., et al.', 2018, 'Epigenetic prediction of complex traits and death', 'Genome Biology, 19, 136', '10.1186/s13059-018-1514-1'),
        stevenson2021: referenceHtml('Stevenson, A. J., Gadd, D. A., Hillary, R. F., et al.', 2021, 'Creating and validating a DNA methylation-based proxy for interleukin-6', 'The Journals of Gerontology: Series A, 76(12), 2284-2292', '10.1093/gerona/glab046'),
        multiqc: referenceHtml('Ewels, P., Magnusson, M., Lundin, S., & Käller, M.', 2016, 'MultiQC: summarize analysis results for multiple tools and samples in a single report', 'Bioinformatics, 32(19), 3047-3048', '10.1093/bioinformatics/btw354')
    ]
    def minfi_norm_refs = [ noob: ['noob'], funnorm: ['funnorm'], quantile: ['quantile'], swan: ['swan', 'noob'] ]
    def keys = [
        sesame: ['sesame', 'noob'],
        minfi: ['minfi', 'minfi_epic'] + minfi_norm_refs[params.minfi_norm],
        enmix: ['minfi', 'enmix', 'relic', 'rcp'],
        combined: ['minfi', 'enmix', 'relic', 'rcp', 'sesame']
    ][params.preprocessing] ?: []
    def qc_keys = params.skip_ewastools ? [] : ['ewastools']
    def cell_keys = params.tissue == 'blood' ? ['epidish', 'luo2023'] : []
    def clock_keys = params.skip_clocks ? [] : ['dnamethyage'] + (params.clocks.tokenize(',').contains('DunedinPACE') ? ['dunedinpace'] : [])
    def exploration_keys = (params.skip_dnam_scores ? [] : ['mccartney2018', 'stevenson2021']) + ['lehne2015'] + (params.n_pcs == 'auto' ? ['rmt'] : [])
    def ewas_keys = params.models ? ['limma', 'bacon', 'mansell2019', 'galwey2009'] + (params.tissue == 'blood' && !params.skip_cell_ewas ? ['celldmc'] : []) + (params.skip_dmr ? [] : ['dmrcate', 'dmrff']) + (params.skip_enrichment ? [] : ['missmethyl']) : []
    return (keys + qc_keys + ['peters2024', 'zhou2017'] + cell_keys + clock_keys + exploration_keys + ewas_keys + ['multiqc']).unique().collect { key -> refs[key] }.join(' ')
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"
    // The command with one option per line, so that it wraps in the reports
    meta["command_text"] = workflow.commandLine.replaceAll(/ (--?[A-Za-z])/, ' \\\\\n    $1')

    // Tool references
    meta["tool_citations"] = toolCitationText()
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
