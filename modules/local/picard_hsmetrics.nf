process PICARD_HSMETRICS {
    label 'process_medium'
    tag "${meta.id}"

    // On-target / panel capture efficiency. For a targeted panel the bait and
    // target intervals are the same panel definition (derived from the panel BED).
    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  panel_intervals

    output:
    tuple val(meta), path('*.hs_metrics.txt'),          emit: metrics
    tuple val(meta), path('*.per_target_coverage.txt'), emit: per_target, optional: true

    script:
    """
    picard CollectHsMetrics \\
        I=${bam} \\
        O=${meta.id}.hs_metrics.txt \\
        R=${fasta} \\
        BAIT_INTERVALS=${panel_intervals} \\
        TARGET_INTERVALS=${panel_intervals} \\
        PER_TARGET_COVERAGE=${meta.id}.per_target_coverage.txt
    """
}
