process FEATURECOUNTS {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam)
    path  gtf

    output:
    tuple val(meta), path('*.featureCounts.txt'),         emit: counts
    tuple val(meta), path('*.featureCounts.txt.summary'), emit: summary

    script:
    def strand = meta.strandedness == 'reverse' ? '2' :
                 meta.strandedness == 'forward' ? '1' : '0'
    def paired = meta.single_end ? '' : '-p --countReadPairs'
    """
    featureCounts \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -s ${strand} \\
        ${paired} \\
        -t exon -g gene_id \\
        -o ${meta.id}.featureCounts.txt \\
        ${bam}
    """
}
