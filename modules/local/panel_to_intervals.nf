process PANEL_TO_INTERVALS {
    label 'process_low'
    tag "${panel_bed.baseName}"

    // Converts the selected panel BED into a Picard interval_list so the same
    // panel definition drives CollectHsMetrics. baits == targets for a panel.
    input:
    path panel_bed
    path dict

    output:
    path "*.interval_list", emit: intervals

    script:
    """
    picard BedToIntervalList \\
        I=${panel_bed} \\
        O=${panel_bed.baseName}.interval_list \\
        SD=${dict}
    """
}
