process NANOREPERTOIRE_REPORT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plotly:4.14.3--py_0' :
        'community.wave.seqera.io/library/pip_plotly:d3e69b73a405fb83' }"

    input:
    path clustercounts
    path cdrcounts
    path cdrhists
    path clusterbig
    path fastaseq

    output:
    path "*.html", emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "nanorepertoire_report"
    """
    nanorepertoire_report.py \\
        --clustercounts ${clustercounts} \\
        --cdrcounts ${cdrcounts} \\
        --cdrhists ${cdrhists} \\
        --clusterbig ${clusterbig} \\
        --fastaseq ${fastaseq} \\
        --identity-threshold ${params.cdhit_identity * 100} \\
        --word-size ${params.cdhit_word_size} \\
        --output ${prefix}.html \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        plotly: \$(python3 -c "import plotly; print(plotly.__version__)")
    END_VERSIONS
    """
}
