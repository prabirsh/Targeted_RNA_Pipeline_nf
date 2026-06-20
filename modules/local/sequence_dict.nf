process SEQUENCE_DICT {
    label 'process_low'
    tag "${fasta.baseName}"

    // Builds a Picard sequence dictionary from the genome FASTA, used to convert
    // panel BED files into interval_lists. Skipped if --dict is supplied.
    input:
    path fasta

    output:
    path "${fasta.baseName}.dict", emit: dict

    script:
    """
    picard CreateSequenceDictionary R=${fasta} O=${fasta.baseName}.dict
    """
}
