process REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seaborn:0.13.2'
        : 'biocontainers/seaborn:0.13.2'}"

    input:
    tuple val(meta), path(tables, stageAs: 'tables/*')
    path samplesheet
    val pipeline_version

    output:
    tuple val(meta), path("*.report.html"), emit: html
    tuple val(meta), path("*.report.pdf") , emit: pdf
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('matplotlib'), eval('python3 -c "import matplotlib; print(matplotlib.__version__)"'), emit: versions_matplotlib, topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python3 -c "import pandas; print(pandas.__version__)"'), emit: versions_pandas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    make_report.py \\
        --tables tables \\
        --samplesheet ${samplesheet} \\
        --arm ${meta.id} \\
        --pipeline-version '${pipeline_version}' \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.report.html
    touch ${prefix}.report.pdf
    """
}
