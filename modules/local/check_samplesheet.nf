process CHECK_SAMPLESHEET {
    label 'process_low'
    tag "$samplesheet"

    input:
    path samplesheet

    output:
    path 'samplesheet.valid.csv', emit: csv

    script:
    """
    check_samplesheet.py $samplesheet samplesheet.valid.csv
    """
}
