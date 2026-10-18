process WITHIN_PERSON {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:2.2.1'
        : 'biocontainers/pandas:2.2.1'}"

    input:
    tuple val(meta), path(genus_counts), path(alpha)
    path samplesheet

    output:
    tuple val(meta), path("*.within_person_reference.tsv"), emit: reference
    tuple val(meta), path("*.within_person_summary.tsv")  , emit: summary
    tuple val(meta), path("*.personal_changes.tsv")       , emit: changes
    tuple val(meta), path("*.personal_genus_changes.tsv") , emit: genus_changes
    tuple val(meta), path("*.identity.tsv")               , emit: identity
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python3 -c "import pandas; print(pandas.__version__)"'), emit: versions_pandas, topic: versions
    tuple val("${task.process}"), val('numpy'), eval('python3 -c "import numpy; print(numpy.__version__)"'), emit: versions_numpy, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    within_person.py \\
        --counts ${genus_counts} \\
        --alpha ${alpha} \\
        --samplesheet ${samplesheet} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.within_person_reference.tsv
    touch ${prefix}.within_person_summary.tsv
    touch ${prefix}.personal_changes.tsv
    touch ${prefix}.personal_genus_changes.tsv
    touch ${prefix}.identity.tsv
    """
}
