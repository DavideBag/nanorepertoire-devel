process READCDHIT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78':
        'ghcr.io/nibscbioinformatics/biopython:v1.78' }"

    input:
    tuple val(meta), path(clusters)


    output:
    tuple val(meta), path("*.summary"), emit: summary
    path '*.summary'                  , emit: summaryonly
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    readcdout.py \
    ${clusters} \
    ${meta.id}_clusters.summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
