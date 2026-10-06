process RENAME {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3'
        : 'biocontainers/seqtk:1.5--h577a1d6_1' }"

    input:
    tuple val(meta), path(reads)
    val single_end

    output:
    tuple val(meta), path('*.fastq.gz'), emit: renamed
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (single_end) {
        """
        seqtk rename ${reads[0]} | gzip > ${prefix}_renamed.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    }
    else {
        """
        seqtk rename ${reads[0]} | gzip > ${prefix}_R1_renamed.fastq.gz
        seqtk rename ${reads[1]} | gzip > ${prefix}_R2_renamed.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
    """
    }
}
