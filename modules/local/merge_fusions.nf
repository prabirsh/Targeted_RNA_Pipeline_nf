process MERGE_FUSIONS {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(arriba), path(starfusion)
    val   has_starfusion

    output:
    tuple val(meta), path('*.fusions.merged.tsv'),     emit: merged
    tuple val(meta), path('*.fusions.reported.tsv'),   emit: reported

    script:
    def sf = has_starfusion ? "--starfusion ${starfusion}" : ''
    """
    merge_fusions.py \\
        --sample ${meta.id} \\
        --arriba ${arriba} \\
        ${sf} \\
        --merged ${meta.id}.fusions.merged.tsv \\
        --reported ${meta.id}.fusions.reported.tsv
    """
}
