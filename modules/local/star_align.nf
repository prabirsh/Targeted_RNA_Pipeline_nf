process STAR_ALIGN {
    label 'process_high_memory'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf

    output:
    tuple val(meta), path('*.Aligned.sortedByCoord.out.bam'),  emit: bam
    tuple val(meta), path('*.Chimeric.out.junction'),          emit: chimeric_junction, optional: true
    tuple val(meta), path('*.Log.final.out'),                  emit: log_final
    tuple val(meta), path('*.SJ.out.tab'),                     emit: junctions
    tuple val(meta), path('*.ReadsPerGene.out.tab'),           emit: gene_counts, optional: true

    script:
    def prefix    = "${meta.id}."
    def reads_cmd = meta.single_end ? "${reads[0]}" : "${reads[0]} ${reads[1]}"
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --sjdbGTFfile ${gtf} \\
        --readFilesIn ${reads_cmd} \\
        --readFilesCommand zcat \\
        --twopassMode Basic \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --outSAMstrandField intronMotif \\
        --quantMode GeneCounts \\
        --chimSegmentMin 12 \\
        --chimOutType WithinBAM Junctions \\
        --chimJunctionOverhangMin 10 \\
        --chimScoreMin 1 \\
        --chimScoreDropMax 30 \\
        --chimScoreJunctionNonGTAG 0 \\
        --chimScoreSeparation 1 \\
        --chimMultimapNmax 50 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --chimSegmentReadGapMax 3 \\
        --outFileNamePrefix ${prefix}
    """
}
