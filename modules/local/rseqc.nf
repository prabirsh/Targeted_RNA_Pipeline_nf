process RSEQC {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)
    path  panel_bed

    output:
    tuple val(meta), path('*.{txt,r,pdf,xls}'), emit: reports

    script:
    """
    # Read distribution over panel target regions
    read_distribution.py -i ${bam} -r ${panel_bed} > ${meta.id}.read_distribution.txt

    # Gene-body coverage (5'-3' bias; FFPE degradation surfaces here)
    geneBody_coverage.py -i ${bam} -r ${panel_bed} -o ${meta.id} || true

    # Inner distance (paired-end insert size profile)
    inner_distance.py -i ${bam} -r ${panel_bed} -o ${meta.id} || true
    """
}
