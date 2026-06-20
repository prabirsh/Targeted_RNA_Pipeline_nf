process PICARD_RNASEQMETRICS {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)
    path  refflat
    path  rrna_intervals

    output:
    tuple val(meta), path('*.rna_metrics.txt'), emit: metrics

    script:
    def strand = meta.strandedness == 'reverse' ? 'SECOND_READ_TRANSCRIPTION_STRAND' :
                 meta.strandedness == 'forward' ? 'FIRST_READ_TRANSCRIPTION_STRAND'  : 'NONE'
    def rrna   = rrna_intervals ? "RIBOSOMAL_INTERVALS=${rrna_intervals}" : ''
    """
    picard CollectRnaSeqMetrics \\
        I=${bam} \\
        O=${meta.id}.rna_metrics.txt \\
        REF_FLAT=${refflat} \\
        STRAND_SPECIFICITY=${strand} \\
        ${rrna}
    """
}
