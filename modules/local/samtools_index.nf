process SAMTOOLS_INDEX {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.bai'),     emit: bai
    tuple val(meta), path('*.flagstat'), emit: flagstat

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    samtools flagstat -@ ${task.cpus} ${bam} > ${bam.baseName}.flagstat
    """
}
