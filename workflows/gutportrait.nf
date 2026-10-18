/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { COLLAPSE_GENUS         } from '../modules/local/collapse_genus/main'
include { TAXON_PERCENTILES      } from '../modules/local/taxon_percentiles/main'
include { ALPHA_DIVERSITY        } from '../modules/local/alpha_diversity/main'
include { BATCH_DIAGNOSE         } from '../modules/local/batch_diagnose/main'
include { DEVIATION_FLAGS        } from '../modules/local/deviation_flags/main'
include { BATCH_CORRECT          } from '../modules/local/batch_correct/main'
include { POSITIVE_CONTROL       } from '../modules/local/positive_control/main'
include { BATCH_DOMINANCE        } from '../modules/local/batch_dominance/main'
include { PERCENTILE_STABILITY   } from '../modules/local/percentile_stability/main'
include { TAXON_PERCENTILES as TAXON_PERCENTILES_CORRECTED } from '../modules/local/taxon_percentiles/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_gutportrait_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GUTPORTRAIT {

    take:
    ch_samplesheet   // channel: [ path(samplesheet) ], validated, value channel
    ampliseq_results // string:  path to an nf-core/ampliseq results directory
    taxonomy_db      // string:  reference database ampliseq classified against, e.g. 'silva'
    arm              // string:  'native' (V3-V4) or 'v4projected'
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // Read the ampliseq handoff. Only named-rank taxonomy is used: the positional
    // rel-table-N.tsv files change meaning with the reference database (under
    // SBDI-GTDB rel-table-6 is Family and no genus table exists at all).
    //
    def asv_table = file("${ampliseq_results}/dada2/ASV_table.tsv", checkIfExists: true)
    def taxonomy  = file("${ampliseq_results}/dada2/ASV_tax_species.${taxonomy_db}.tsv", checkIfExists: true)

    //
    // STAGE 2 - HARMONISE: collapse to genus, keeping unassigned reads as an explicit row
    //
    COLLAPSE_GENUS (
        channel.of([ [ id: arm ], asv_table, taxonomy ])
    )

    //
    // STAGE 4 - POSITION: per-genus CLR percentiles of each personal sample vs the panel
    //
    TAXON_PERCENTILES (
        COLLAPSE_GENUS.out.counts,
        ch_samplesheet
    )

    //
    // STAGE 3 - INTEGRATE, diagnosis only: how much of the panel's variation is study, and
    // how much is each metadata variable. Correction and the two gates come with BATCH_CORRECT.
    //
    BATCH_DIAGNOSE (
        COLLAPSE_GENUS.out.counts,
        ch_samplesheet
    )

    //
    // STAGE 3 - INTEGRATE, correction: harmonise the PANEL studies only (personal samples are
    // never corrected), check with a planted signal that correction preserves real effects
    // (the run stops if it does not), and recompute every percentile on the corrected panel.
    //
    BATCH_CORRECT (
        COLLAPSE_GENUS.out.counts,
        ch_samplesheet
    )
    POSITIVE_CONTROL (
        COLLAPSE_GENUS.out.counts,
        ch_samplesheet
    )
    TAXON_PERCENTILES_CORRECTED (
        BATCH_CORRECT.out.counts,
        ch_samplesheet
    )
    PERCENTILE_STABILITY (
        TAXON_PERCENTILES.out.percentiles.join(TAXON_PERCENTILES_CORRECTED.out.percentiles)
    )

    //
    // GATES - batch dominance: which contrasts and which personal deviations may be read as
    // biology at all, given how much variation is between studies
    //
    BATCH_DOMINANCE (
        BATCH_DIAGNOSE.out.permanova
            .join(BATCH_DIAGNOSE.out.summary)
            .join(COLLAPSE_GENUS.out.counts),
        ch_samplesheet
    )

    //
    // STAGE 4 - POSITION: alpha diversity at one fixed depth for every sample
    //
    ALPHA_DIVERSITY (
        channel.of([ [ id: arm ], asv_table ]),
        ch_samplesheet
    )

    //
    // STAGE 5 - DEVIATE: extremes-only flags. Provisional until both arms agree.
    //
    DEVIATION_FLAGS (
        channel.of([ [ id: arm ], asv_table, taxonomy ]).join(ALPHA_DIVERSITY.out.alpha),
        ch_samplesheet
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'gutportrait_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'gutportrait'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
