process MOSDEPTH {
    label 'process_medium'
    tag "${meta.id}"

    // Per-target coverage over the panel BED, plus the fraction of target bases
    // reaching 50x / 100x / 200x / 500x / 1000x (the clinical depth thresholds).
    input:
    tuple val(meta), path(bam), path(bai)
    path  panel_bed

    output:
    tuple val(meta), path("*.mosdepth.summary.txt"),        emit: summary
    tuple val(meta), path("*.thresholds.bed.gz"),           emit: thresholds
    tuple val(meta), path("*.regions.bed.gz"),              emit: regions
    tuple val(meta), path("*.mosdepth.region.dist.txt"),    emit: dist, optional: true
    tuple val(meta), path("*.mosdepth.global.dist.txt"),    emit: global_dist

    script:
    """
    mosdepth \\
        --threads ${task.cpus} \\
        --by ${panel_bed} \\
        --thresholds 50,100,200,500,1000 \\
        --no-per-base \\
        ${meta.id} \\
        ${bam}
    """
}
