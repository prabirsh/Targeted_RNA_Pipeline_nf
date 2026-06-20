process FASTP {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.fastp.json'),    emit: json
    tuple val(meta), path('*.fastp.html'),    emit: html
    tuple val(meta), path('*.fastp.log'),     emit: log

    script:
    def prefix = "${meta.id}"
    if (meta.single_end) {
        """
        fastp \\
            -i ${reads[0]} \\
            -o ${prefix}_R1.trim.fastq.gz \\
            --trim_poly_g \\
            --qualified_quality_phred 15 \\
            --length_required 30 \\
            --thread ${task.cpus} \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            2> ${prefix}.fastp.log
        """
    } else {
        """
        fastp \\
            -i ${reads[0]} -I ${reads[1]} \\
            -o ${prefix}_R1.trim.fastq.gz -O ${prefix}_R2.trim.fastq.gz \\
            --detect_adapter_for_pe \\
            --trim_poly_g \\
            --qualified_quality_phred 15 \\
            --length_required 30 \\
            --thread ${task.cpus} \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            2> ${prefix}.fastp.log
        """
    }
}
