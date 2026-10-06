process NANOCDRX {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/nanocdr-x:1.0.0--146125dcde3e4498'
        : 'community.wave.seqera.io/library/nanocdr-x:1.0.0--539f1c75ff746d30'}"

    input:
    tuple val(meta), path(translated)

    output:
    tuple val("${meta.id}"), val("${meta.immunisation}"), val("${meta.boost}"), path("*_cdr3.fasta"), emit: fasta
    tuple val(meta), path("*_cdr3.hist"), emit: hist
    tuple val(meta), path("*_cdr3.tsv"), emit: tsv
    path '*_cdr3.hist', emit: histonly
    path '*_cdr3.tsv', emit: tsvonly
    val meta, emit: metaonly
    tuple val(meta), path("*_cdr3_summary.tsv"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_nanocdrx.py \\
        -i ${translated} \\
        -p ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanocdr-x: \$(python -c "import importlib.metadata as m; print(m.version('nanocdr-x'))")
    END_VERSIONS
    """
}
