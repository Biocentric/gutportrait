//
// Subworkflow with functionality specific to the Biocentric/gutportrait pipeline
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
    min_panel_samples // integer: Minimum number of panel samples the samplesheet must contain
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

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --ampliseq_results <AMPLISEQ_OUTDIR> --outdir <OUTDIR>"

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
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Validate the samplesheet. gutportrait takes counts from the ampliseq results
    // directory and metadata from this sheet, so the analysis modules need the sheet
    // itself rather than a channel of per-sample files. Each row is checked against
    // assets/schema_input.json; the rules that span rows are checked here.
    //
    def rows = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
    def roles = validateSamplesheetRows(rows, min_panel_samples)
    log.info "Samplesheet: ${roles.personal} personal and ${roles.panel} panel sample(s) across ${roles.studies} studies, ${roles.repeats} repeat sample(s)"
    def ch_samplesheet = channel.value(file(input, checkIfExists: true))

    emit:
    samplesheet = ch_samplesheet
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
                multiqc_reports.getVal(),
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
// Checks that span rows, which a per-row JSON schema cannot express
//
def validateSamplesheetRows(rows, min_panel_samples) {
    // samplesheetToList puts the meta map first; every column here is a meta field
    def metas = rows.collect { row -> row instanceof List ? row[0] : row }

    def duplicated = metas.collect { meta -> meta.id }
        .countBy { id -> id }
        .findAll { _id, n -> n > 1 }
        .keySet()
    if (duplicated) {
        error("Samplesheet: sample names must be unique. Duplicated: ${duplicated.join(', ')}")
    }

    def n_personal = metas.count { meta -> meta.role == 'personal' }
    def n_panel    = metas.count { meta -> meta.role == 'panel' }
    def n_studies  = metas.findAll { meta -> meta.role == 'panel' }.collect { meta -> meta.study }.unique().size()

    if (n_personal < 1) {
        error("Samplesheet: no rows have role 'personal'. There is nothing to position against the panel.")
    }
    if (n_panel < min_panel_samples) {
        error("Samplesheet: ${n_panel} panel sample(s), but --min_panel_samples is ${min_panel_samples}. A percentile against fewer samples is too coarse to report.")
    }
    def mixed = metas.findAll { meta -> meta.role == 'personal' }.collect { meta -> meta.study }
        .intersect(metas.findAll { meta -> meta.role == 'panel' }.collect { meta -> meta.study })
    if (mixed) {
        error("Samplesheet: study ${mixed.unique().join(', ')} contains both personal and panel samples. Give the personal samples their own study, so they are never mistaken for a reference batch.")
    }

    // A set of bands ("30s|40s") describes a study's inclusion criteria, never a person
    // whose bands are known: the personal samples need single values.
    def ranged = metas.findAll { meta -> meta.role == 'personal' && [meta.age_cat, meta.sex, meta.bmi_cat].any { v -> v && v.toString().contains('|') } }
    if (ranged) {
        error("Samplesheet: personal sample(s) ${ranged.collect { meta -> meta.id }.join(', ')} give a set of bands (with '|'). Give your own age, sex and BMI band as single values; sets are for panel studies that publish only inclusion criteria.")
    }

    // A repeat is a second sample of a panel person: it must name that person in `subject`
    // (the subject of a panel row, which defaults to its sample name), in the same study.
    def panel_subjects = metas.findAll { meta -> meta.role == 'panel' }
        .collect { meta -> "${meta.study}\t${meta.subject ?: meta.id}".toString() } as Set
    def repeats = metas.findAll { meta -> meta.role == 'repeat' }
    def orphans = repeats.findAll { meta -> !meta.subject || !panel_subjects.contains("${meta.study}\t${meta.subject}".toString()) }
    if (orphans) {
        error("Samplesheet: repeat sample(s) ${orphans.collect { meta -> meta.id }.join(', ')} must name, in 'subject', a panel person of the same study (the panel row's subject, or its sample name if subject is empty).")
    }

    return [personal: n_personal, panel: n_panel, studies: n_studies, repeats: repeats.size()]
}
//
// Generate methods description for MultiQC
//
//
// ampliseq is a separate pipeline run, but its output is this pipeline's input, so the tools
// that produced it belong in the methods text. SILVA is cited only when it was the database.
//
def toolCitationText() {
    def taxonomy = params.taxonomy_db == 'silva' ? " and classified against SILVA (Quast et al. 2013)" : ""
    def citation_text = [
            "Amplicon sequence variants were inferred upstream with nf-core/ampliseq (Straub et al. 2020) using DADA2 (Callahan et al. 2016)${taxonomy}.",
            "Genus tables, centred-log-ratio percentiles, rarefied alpha diversity and deviation flags were computed with pandas (McKinney 2010) and NumPy (Harris et al. 2020), treating the data as compositional (Gloor et al. 2017).",
            "Between-study variance was partitioned by PERMANOVA (Anderson 2001) on Aitchison distance with vegan (Oksanen et al.), and reference studies were harmonised with MMUPHin (Ma et al. 2022), checked by an injected-signal positive control.",
            "Results were summarised with MultiQC (Ewels et al. 2016)."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            '<li>Straub, D., Blackwell, N., Langarica-Fuentes, A., Peltzer, A., Nahnsen, S., & Kleindienst, S. (2020). Interpretations of environmental microbial community studies are biased by the selected 16S rRNA (gene) amplicon sequencing pipeline. Frontiers in Microbiology, 11, 550420. doi: <a href="https://doi.org/10.3389/fmicb.2020.550420">10.3389/fmicb.2020.550420</a></li>',
            '<li>Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods, 13(7), 581–583. doi: <a href="https://doi.org/10.1038/nmeth.3869">10.1038/nmeth.3869</a></li>',
            params.taxonomy_db == 'silva' ? '<li>Quast, C., Pruesse, E., Yilmaz, P., Gerken, J., Schweer, T., Yarza, P., Peplies, J., & Glöckner, F. O. (2013). The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Research, 41(Database issue), D590–D596. doi: <a href="https://doi.org/10.1093/nar/gks1219">10.1093/nar/gks1219</a></li>' : '',
            '<li>McKinney, W. (2010). Data structures for statistical computing in Python. Proceedings of the 9th Python in Science Conference, 56–61. doi: <a href="https://doi.org/10.25080/Majora-92bf1922-00a">10.25080/Majora-92bf1922-00a</a></li>',
            '<li>Harris, C. R., Millman, K. J., van der Walt, S. J., et al. (2020). Array programming with NumPy. Nature, 585(7825), 357–362. doi: <a href="https://doi.org/10.1038/s41586-020-2649-2">10.1038/s41586-020-2649-2</a></li>',
            '<li>Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J. (2017). Microbiome datasets are compositional: and this is not optional. Frontiers in Microbiology, 8, 2224. doi: <a href="https://doi.org/10.3389/fmicb.2017.02224">10.3389/fmicb.2017.02224</a></li>',
            '<li>Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. Austral Ecology, 26(1), 32–46. doi: <a href="https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x">10.1111/j.1442-9993.2001.01070.pp.x</a></li>',
            '<li>Ma, S., Shungin, D., Mallick, H., Schirmer, M., Nguyen, L. H., Kolde, R., Franzosa, E., Vlamakis, H., Xavier, R., & Huttenhower, C. (2022). Population structure discovery in meta-analyzed microbial communities and inflammatory bowel disease using MMUPHin. Genome Biology, 23, 208. doi: <a href="https://doi.org/10.1186/s13059-022-02753-4">10.1186/s13059-022-02753-4</a></li>',
            '<li>Oksanen, J., Simpson, G. L., Blanchet, F. G., et al. vegan: Community Ecology Package. R package version 2.7-5. <a href="https://CRAN.R-project.org/package=vegan">https://CRAN.R-project.org/package=vegan</a></li>',
            '<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. doi: <a href="https://doi.org/10.1093/bioinformatics/btw354">10.1093/bioinformatics/btw354</a></li>'
        ].join(' ').trim()

    return reference_text
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
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
