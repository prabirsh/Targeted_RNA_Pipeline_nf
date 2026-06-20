#!/usr/bin/env nextflow
/*
========================================================================================
    RNA Panel Pipeline  —  targeted-panel transcriptome analysis (DRAGEN-RNA-style)
    Clinical FFPE RNA, exome / targeted panel, fusion + expression + coverage QC
========================================================================================
*/

nextflow.enable.dsl = 2

include { CHECK_SAMPLESHEET    } from './modules/local/check_samplesheet.nf'
include { FASTQC               } from './modules/local/fastqc.nf'
include { FASTP                } from './modules/local/fastp.nf'
include { STAR_ALIGN           } from './modules/local/star_align.nf'
include { SAMTOOLS_INDEX       } from './modules/local/samtools_index.nf'
include { SEQUENCE_DICT        } from './modules/local/sequence_dict.nf'
include { PANEL_TO_INTERVALS   } from './modules/local/panel_to_intervals.nf'
include { PICARD_RNASEQMETRICS } from './modules/local/picard_rnaseqmetrics.nf'
include { PICARD_HSMETRICS     } from './modules/local/picard_hsmetrics.nf'
include { PICARD_INSERTSIZE    } from './modules/local/picard_insertsize.nf'
include { MOSDEPTH             } from './modules/local/mosdepth.nf'
include { COVERAGE_SUMMARY     } from './modules/local/coverage_summary.nf'
include { RSEQC                } from './modules/local/rseqc.nf'
include { ARRIBA               } from './modules/local/arriba.nf'
include { STARFUSION           } from './modules/local/starfusion.nf'
include { MERGE_FUSIONS        } from './modules/local/merge_fusions.nf'
include { FEATURECOUNTS        } from './modules/local/featurecounts.nf'
include { MULTIQC              } from './modules/local/multiqc.nf'

def helpMessage() {
    log.info """
    ============================================================
     RNA PANEL PIPELINE  v1.1  (clinical FFPE RNA targeted panel)
    ============================================================
    Usage:
      nextflow run main.nf -profile docker \\
        --input samplesheet.csv --outdir results \\
        --panel ST8 \\
        --star_index /refs/star_index --gtf /refs/annotation.gtf --fasta /refs/genome.fa \\
        --refflat /refs/refFlat.txt --rrna_intervals /refs/rRNA.interval_list \\
        --ctat_lib /refs/ctat_genome_lib \\
        --arriba_blacklist /refs/arriba/blacklist.tsv.gz

    Panel (required): --panel ST8 | CT     (or --panel_bed <file.bed> to override)
      ST8 = SureSelect XT V8 exome panel
      CT  = Indiegene 4baseCare RNA Fusion panel (GRCh37)

    Required references: --star_index --gtf --fasta
    Optional: --dict (else built from --fasta), --refflat, --rrna_intervals,
              --ctat_lib, --arriba_* , --gene_bed (RSeQC), --skip_fusion, --skip_rseqc
    """.stripIndent()
}

if (params.help) { helpMessage(); exit 0 }

// ---- fail fast --------------------------------------------------------------
def requireParam = { name, value ->
    if (!value) { log.error "ERROR: missing required parameter --${name}. Run with --help."; exit 1 }
}
requireParam('input',      params.input)
requireParam('star_index', params.star_index)
requireParam('gtf',        params.gtf)
requireParam('fasta',      params.fasta)

// ---- resolve the panel (the 'which panel?' decision) ------------------------
def panel_bed_path = params.panel_bed ?: (params.panel ? params.panels[params.panel] : null)
if (!panel_bed_path) {
    log.error "No panel selected. Use --panel <NAME> (one of: ${params.panels.keySet().join(', ')}) or --panel_bed <file.bed>."
    exit 1
}
def panel_name = params.panel ?: file(panel_bed_path).baseName
log.info "PANEL: ${panel_name}  ->  ${panel_bed_path}"

workflow {

    log.info "RNA PANEL PIPELINE v1.1  |  input=${params.input}  |  panel=${panel_name}  |  outdir=${params.outdir}"

    ch_panel_bed = file(panel_bed_path, checkIfExists: true)

    // Stage 0: samplesheet -> read channel
    CHECK_SAMPLESHEET(Channel.fromPath(params.input, checkIfExists: true))
    ch_reads = CHECK_SAMPLESHEET.out.csv
        .splitCsv(header: true)
        .map { row ->
            def meta = [ id: row.sample,
                         single_end: (row.fastq_2 == null || row.fastq_2 == ''),
                         strandedness: (row.strandedness ?: 'reverse'),
                         panel: panel_name ]
            def fastqs = meta.single_end ?
                [ file(row.fastq_1, checkIfExists: true) ] :
                [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ]
            tuple(meta, fastqs)
        }

    // Stage 1: raw QC + trim
    FASTQC(ch_reads)
    FASTP(ch_reads)

    // Stage 2: STAR 2-pass
    STAR_ALIGN(FASTP.out.reads, file(params.star_index), file(params.gtf))
    SAMTOOLS_INDEX(STAR_ALIGN.out.bam)
    ch_bam = STAR_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.bai)

    // Panel -> interval_list (needs a sequence dictionary)
    ch_dict = params.dict ? Channel.value(file(params.dict)) : SEQUENCE_DICT(file(params.fasta)).dict
    PANEL_TO_INTERVALS(ch_panel_bed, ch_dict)

    // Stage 3: post-alignment + panel QC
    PICARD_RNASEQMETRICS(ch_bam, file(params.refflat), file(params.rrna_intervals))
    PICARD_HSMETRICS(ch_bam, file(params.fasta), PANEL_TO_INTERVALS.out.intervals)
    PICARD_INSERTSIZE(ch_bam)

    // Coverage @50/100/200/500/1000x over the panel
    MOSDEPTH(ch_bam, ch_panel_bed)
    ch_cov_in = MOSDEPTH.out.thresholds
        .join(MOSDEPTH.out.summary)
        .join(SAMTOOLS_INDEX.out.flagstat)
    COVERAGE_SUMMARY(ch_cov_in)

    ch_rseqc = Channel.empty()
    if (!params.skip_rseqc && params.gene_bed) {
        RSEQC(ch_bam, file(params.gene_bed))
        ch_rseqc = RSEQC.out.reports
    }

    // Stage 4: expression
    FEATURECOUNTS(STAR_ALIGN.out.bam, file(params.gtf))

    // Stage 5: ensemble fusion (breakpoints + read support + read names)
    if (!params.skip_fusion) {
        ARRIBA(
            STAR_ALIGN.out.bam, file(params.fasta), file(params.gtf),
            params.arriba_blacklist ? file(params.arriba_blacklist) : [],
            params.arriba_known     ? file(params.arriba_known)     : [],
            params.arriba_protein   ? file(params.arriba_protein)   : []
        )
        if (params.ctat_lib) {
            STARFUSION(STAR_ALIGN.out.chimeric_junction, file(params.ctat_lib))
            MERGE_FUSIONS(ARRIBA.out.fusions.join(STARFUSION.out.fusions), true)
        } else {
            MERGE_FUSIONS(ARRIBA.out.fusions.map { meta, a -> tuple(meta, a, []) }, false)
        }
    }

    // Stage 9: aggregate QC
    ch_multiqc = Channel.empty()
        .mix( FASTQC.out.zip.map { meta, f -> f } )
        .mix( FASTP.out.json.map { meta, f -> f } )
        .mix( STAR_ALIGN.out.log_final.map { meta, f -> f } )
        .mix( PICARD_RNASEQMETRICS.out.metrics.map { meta, f -> f } )
        .mix( PICARD_HSMETRICS.out.metrics.map { meta, f -> f } )
        .mix( PICARD_INSERTSIZE.out.metrics.map { meta, f -> f } )
        .mix( MOSDEPTH.out.summary.map { meta, f -> f } )
        .mix( COVERAGE_SUMMARY.out.mqc.map { meta, f -> f } )
        .mix( ch_rseqc.map { meta, f -> f }.flatten() )
        .mix( FEATURECOUNTS.out.summary.map { meta, f -> f } )
        .collect()

    MULTIQC(ch_multiqc, file("${projectDir}/assets/multiqc_config.yml"))
}

workflow.onComplete {
    log.info "============================================================"
    log.info " Pipeline complete: ${ workflow.success ? 'SUCCESS' : 'FAILED' }"
    log.info " Panel:   ${panel_name}"
    log.info " Results: ${params.outdir}"
    log.info " Report:  ${params.outdir}/multiqc/multiqc_report.html"
    log.info "============================================================"
}
