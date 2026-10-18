process ALPHA_DIVERSITY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:2.2.1'
        : 'biocontainers/pandas:2.2.1'}"

    input:
    tuple val(meta), path(asv_table)
    path samplesheet

    output:
    tuple val(meta), path("*.alpha_diversity.tsv"), emit: alpha
    tuple val(meta), path("*.alpha_position.tsv"), emit: position
    tuple val(meta), path("*.alpha_params.tsv"), emit: params
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python3 -c "import pandas; print(pandas.__version__)"'), emit: versions_pandas, topic: versions
    tuple val("${task.process}"), val('numpy'), eval('python3 -c "import numpy; print(numpy.__version__)"'), emit: versions_numpy, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    alpha_diversity.py \\
        --asv-table ${asv_table} \\
        --samplesheet ${samplesheet} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alpha_diversity.tsv
    touch ${prefix}.alpha_position.tsv
    touch ${prefix}.alpha_params.tsv
    """
}
