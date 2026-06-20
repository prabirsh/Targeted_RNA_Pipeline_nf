process FASTQC {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.html'), emit: html
    tuple val(meta), path('*.zip'),  emit: zip

    script:
    """
    fastqc --threads ${task.cpus} --quiet ${reads.join(' ')}
    """
}
