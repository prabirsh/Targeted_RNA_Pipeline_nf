process COVERAGE_SUMMARY {
    label 'process_low'
    tag "${meta.id}"

    // Rolls mosdepth + flagstat into one tidy clinical QC row:
    // total reads, mean target coverage, and % target bases >= 50/100/200/500/1000x.
    input:
    tuple val(meta), path(thresholds), path(summary), path(flagstat)

    output:
    tuple val(meta), path("*.coverage_summary.tsv"),     emit: tsv
    tuple val(meta), path("*.coverage_mqc.txt"),         emit: mqc

    script:
    """
    coverage_summary.py \\
        --sample ${meta.id} \\
        --thresholds ${thresholds} \\
        --summary ${summary} \\
        --flagstat ${flagstat} \\
        --out ${meta.id}.coverage_summary.tsv \\
        --mqc ${meta.id}.coverage_mqc.txt
    """
}
