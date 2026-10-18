process BATCH_CORRECT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioconductor-mmuphin:1.16.0--r43hdfd78af_0'
        : 'biocontainers/bioconductor-mmuphin:1.16.0--r43hdfd78af_0'}"

    input:
    tuple val(meta), path(genus_counts)
    path samplesheet

    output:
    tuple val(meta), path("*.genus_counts_corrected.tsv"), emit: counts
    tuple val(meta), path("*.correction_summary.tsv"), emit: summary
    tuple val("${task.process}"), val('r-base'), eval('Rscript -e "cat(paste(R.version[c(\'major\', \'minor\')], collapse = \'.\'))"'), emit: versions_r_base, topic: versions
    tuple val("${task.process}"), val('mmuphin'), eval('Rscript -e "cat(as.character(packageVersion(\'MMUPHin\')))"'), emit: versions_mmuphin, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    batch_correct.R \\
        --counts ${genus_counts} \\
        --samplesheet ${samplesheet} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genus_counts_corrected.tsv
    touch ${prefix}.correction_summary.tsv
    """
}
