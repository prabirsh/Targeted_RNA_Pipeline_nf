process STARFUSION {
    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(chimeric_junction)
    path  ctat_lib

    output:
    tuple val(meta), path('*.starfusion.fusions.tsv'), emit: fusions

    script:
    """
    STAR-Fusion \\
        --genome_lib_dir ${ctat_lib} \\
        --chimeric_junction ${chimeric_junction} \\
        --CPU ${task.cpus} \\
        --output_dir starfusion_out

    cp starfusion_out/star-fusion.fusion_predictions.tsv ${meta.id}.starfusion.fusions.tsv
    """
}
