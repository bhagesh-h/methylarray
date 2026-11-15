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
    return name in ['predicted_sex', 'dnam_smoking'] || name ==~ /cc_[A-Za-z0-9_]+/ || name ==~ /ctrlpc[0-9]+/
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
        // Tokens with a fixed expansion; @cells and @ctrlpcs depend on the data
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
    return [ preprocessing[params.preprocessing] ?: '', "Results were summarised with MultiQC (Ewels et al. 2016)." ].join(' ').trim()
}

//
// One reference as an HTML list item for the MultiQC methods description
//
def referenceHtml(authors, year, title, journal, doi) {
    return "<li>${authors} (${year}). ${title}. ${journal}. doi: <a href=\"https://doi.org/${doi}\">${doi}</a></li>"
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
        multiqc: referenceHtml('Ewels, P., Magnusson, M., Lundin, S., & Käller, M.', 2016, 'MultiQC: summarize analysis results for multiple tools and samples in a single report', 'Bioinformatics, 32(19), 3047-3048', '10.1093/bioinformatics/btw354')
    ]
    def minfi_norm_refs = [ noob: ['noob'], funnorm: ['funnorm'], quantile: ['quantile'], swan: ['swan', 'noob'] ]
    def keys = [
        sesame: ['sesame', 'noob'],
        minfi: ['minfi', 'minfi_epic'] + minfi_norm_refs[params.minfi_norm],
        enmix: ['minfi', 'enmix', 'relic', 'rcp'],
        combined: ['minfi', 'enmix', 'relic', 'rcp', 'sesame']
    ][params.preprocessing] ?: []
    return (keys + ['multiqc']).unique().collect { key -> refs[key] }.join(' ')
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

    // Tool references
    meta["tool_citations"] = toolCitationText()
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
