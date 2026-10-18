process BATCH_DIAGNOSE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/r-vegan:2.7.5'
        : 'biocontainers/r-vegan:2.7.5'}"

    input:
    tuple val(meta), path(genus_counts)
    path samplesheet

    output:
    tuple val(meta), path("*.permanova.tsv"), emit: permanova
    tuple val(meta), path("*.batch_summary.tsv"), emit: summary
    tuple val("${task.process}"), val('r-base'), eval('Rscript -e "cat(paste(R.version[c(\'major\', \'minor\')], collapse = \'.\'))"'), emit: versions_r_base, topic: versions
    tuple val("${task.process}"), val('vegan'), eval('Rscript -e "cat(as.character(packageVersion(\'vegan\')))"'), emit: versions_vegan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    batch_diagnose.R \\
        --counts ${genus_counts} \\
        --samplesheet ${samplesheet} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.permanova.tsv
    touch ${prefix}.batch_summary.tsv
    """
}
