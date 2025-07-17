process NANOTRANSLATE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'ghcr.io/nibscbioinformatics/biopython:v1.78' }"

    input:
    tuple val(meta), path(mergedfq)

    output:
    tuple val(meta), path('*.fasta'), emit: fasta
    tuple val(meta), path('*.log')  , emit: log
    path 'versions.yml'             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -cd ${mergedfq} >sequences.fastq
    translate.py \
    sequences.fastq \
    ${meta.id}_AA_translated.fasta \
    >${meta.id}_AA_translated.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}

