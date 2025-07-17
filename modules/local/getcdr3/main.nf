process GETCDR3 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78':
        'ghcr.io/nibscbioinformatics/biopython:v1.78' }"

    input:
    tuple val(meta), path(translated)


    output:
    tuple val("${meta.id}"), val("${meta.immunisation}"), val("${meta.boost}"), path("*.fasta"), emit: fasta
    tuple val(meta), path("*.hist")                                                            , emit: hist
    tuple val(meta), path("*.tsv")                                                             , emit: tsv
    path '*.hist'                                                                              , emit: histonly
    path '*.tsv'                                                                               , emit: tsvonly
    val meta                                                                                   , emit: metaonly
    path "versions.yml"                                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    getcdr3.py \
    -i ${translated} \
    -c ${meta.id}_cdr3.fasta \
    -o ${meta.id}_cdr3.hist \
    -t ${meta.id}_cdr3.tsv \
    -s ${meta.id} \
    -m ${meta.immunisation} \
    -b ${meta.boost}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

}
