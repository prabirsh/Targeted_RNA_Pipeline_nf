process MULTIQC {
    label 'process_low'

    input:
    path  qc_files
    path  config

    output:
    path '*.html',           emit: report
    path 'multiqc_data',     emit: data

    script:
    def cfg = config.name != 'NO_FILE' ? "--config ${config}" : ''
    """
    multiqc -f ${cfg} .
    """
}
