process ARRIBA {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam)
    path  fasta
    path  gtf
    path  blacklist
    path  known_fusions
    path  protein_domains

    output:
    tuple val(meta), path('*.arriba.fusions.tsv'),           emit: fusions
    tuple val(meta), path('*.arriba.fusions.discarded.tsv'), emit: discarded, optional: true

    script:
    // STAR ran with --chimOutType WithinBAM, so chimeric reads are inside the main BAM.
    def bl = blacklist       ? "-b ${blacklist}"       : ""
    def kn = known_fusions   ? "-k ${known_fusions}"   : ""
    def pd = protein_domains ? "-p ${protein_domains}" : ""
    """
    arriba \\
        -x ${bam} \\
        -a ${fasta} \\
        -g ${gtf} \\
        ${bl} ${kn} ${pd} \\
        -o ${meta.id}.arriba.fusions.tsv \\
        -O ${meta.id}.arriba.fusions.discarded.tsv
    """
}
