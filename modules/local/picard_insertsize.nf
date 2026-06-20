process PICARD_INSERTSIZE {
    label 'process_medium'
    tag "${meta.id}"

    // Insert-size (fragment length) distribution. For FFPE this is short; a
    // collapsed insert-size distribution flags over-fragmented input.
    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.insert_size_metrics.txt"), emit: metrics
    tuple val(meta), path("*.insert_size_histogram.pdf"), emit: histogram, optional: true

    script:
    """
    picard CollectInsertSizeMetrics \\
        I=${bam} \\
        O=${meta.id}.insert_size_metrics.txt \\
        H=${meta.id}.insert_size_histogram.pdf \\
        M=0.5
    """
}
