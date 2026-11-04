#!/usr/bin/env nextflow
/*
========================================================================================
    RNA Panel Pipeline: targeted-panel RNA analysis
    QC, alignment, panel coverage, fusion detection and expression.
    Research use. Not validated for clinical decisions.
========================================================================================
*/

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
    return """
    RNA Panel Pipeline v${workflow.manifest.version}

    Usage:
      nextflow run main.nf -profile docker \\
        --input samplesheet.csv --outdir results --panel panel.bed \\
        --star_index /refs/star_index --gtf /refs/annotation.gtf --fasta /refs/genome.fa

    Required:
      --input             Samplesheet CSV: sample,fastq_1,fastq_2,strandedness
      --panel             BED file of the panel's target regions
      --star_index        STAR index directory
      --gtf               Annotation GTF
      --fasta             Genome FASTA, with .fai and .dict beside it

    Optional:
      --outdir            Results directory (default: results)
      --refflat           refFlat file; without it, CollectRnaSeqMetrics is skipped
      --rrna_intervals    Ribosomal interval_list for CollectRnaSeqMetrics
      --gene_bed          BED12 gene model; without it, RSeQC is skipped
      --ctat_lib          CTAT library; without it, STAR-Fusion is skipped
      --arriba_blacklist, --arriba_known, --arriba_protein    Arriba reference files
      --skip_fusion, --skip_rseqc

    Profiles: docker, singularity, conda, slurm, awsbatch, test
    """.stripIndent()
}

def requireParam(name, value) {
    if (!value) {
        error "Missing required parameter --${name}. Run with --help for usage."
    }
}

workflow {
    main:
    if (params.help) {
        log.info helpMessage()
        exit(0)
    }

    requireParam('input',      params.input)
    requireParam('panel',      params.panel)
    requireParam('star_index', params.star_index)
    requireParam('gtf',        params.gtf)
    requireParam('fasta',      params.fasta)

    def panel_bed  = file(params.panel, checkIfExists: true)
    def panel_name = panel_bed.baseName
    def fasta      = file(params.fasta, checkIfExists: true)
    def gtf        = file(params.gtf, checkIfExists: true)
    def fasta_fai  = file("${params.fasta}.fai", checkIfExists: true)
    def fasta_dict = params.dict
        ? file(params.dict, checkIfExists: true)
        : file("${fasta.parent}/${fasta.baseName}.dict")

    log.info "RNA Panel Pipeline v${workflow.manifest.version} | input=${params.input} | panel=${panel_name} | outdir=${params.outdir}"

    // Samplesheet -> read channel
    CHECK_SAMPLESHEET(channel.fromPath(params.input, checkIfExists: true))
    ch_reads = CHECK_SAMPLESHEET.out.csv
        .splitCsv(header: true)
        .map { row ->
            def single_end = !row.fastq_2
            def meta = [
                id: row.sample,
                single_end: single_end,
                strandedness: row.strandedness ?: 'reverse',
                panel: panel_name,
            ]
            def fastqs = single_end
                ? [file(row.fastq_1, checkIfExists: true)]
                : [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
            tuple(meta, fastqs)
        }

    // Raw QC and trimming
    FASTQC(ch_reads)
    FASTP(ch_reads)

    // Alignment
    STAR_ALIGN(FASTP.out.reads, file(params.star_index, checkIfExists: true), gtf)
    SAMTOOLS_INDEX(STAR_ALIGN.out.bam)
    ch_bam = STAR_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.bai)

    // The dictionary is built from the FASTA if it is not already present.
    ch_dict = fasta_dict.exists()
        ? channel.value(fasta_dict)
        : SEQUENCE_DICT(fasta).dict
    PANEL_TO_INTERVALS(panel_bed, ch_dict)

    // Post-alignment and panel QC
    ch_rnaseq_metrics = channel.empty()
    if (params.refflat) {
        PICARD_RNASEQMETRICS(
            ch_bam,
            file(params.refflat, checkIfExists: true),
            params.rrna_intervals ? file(params.rrna_intervals, checkIfExists: true) : []
        )
        ch_rnaseq_metrics = PICARD_RNASEQMETRICS.out.metrics
    }
    PICARD_HSMETRICS(ch_bam, fasta, fasta_fai, ch_dict, PANEL_TO_INTERVALS.out.intervals)
    PICARD_INSERTSIZE(ch_bam.filter { meta, _bam, _bai -> !meta.single_end })

    MOSDEPTH(ch_bam, panel_bed)
    COVERAGE_SUMMARY(
        MOSDEPTH.out.thresholds
            .join(MOSDEPTH.out.summary)
            .join(SAMTOOLS_INDEX.out.flagstat)
    )

    ch_rseqc = channel.empty()
    if (!params.skip_rseqc && params.gene_bed) {
        RSEQC(ch_bam, file(params.gene_bed, checkIfExists: true))
        ch_rseqc = RSEQC.out.reports
    }

    // Expression
    FEATURECOUNTS(STAR_ALIGN.out.bam, gtf)

    // Fusions
    if (!params.skip_fusion) {
        ARRIBA(
            STAR_ALIGN.out.bam, fasta, gtf,
            params.arriba_blacklist ? file(params.arriba_blacklist, checkIfExists: true) : [],
            params.arriba_known     ? file(params.arriba_known, checkIfExists: true)     : [],
            params.arriba_protein   ? file(params.arriba_protein, checkIfExists: true)   : []
        )
        if (params.ctat_lib) {
            STARFUSION(STAR_ALIGN.out.chimeric_junction, file(params.ctat_lib, checkIfExists: true))
            MERGE_FUSIONS(ARRIBA.out.fusions.join(STARFUSION.out.fusions), true)
        }
        else {
            MERGE_FUSIONS(ARRIBA.out.fusions.map { meta, arriba -> tuple(meta, arriba, []) }, false)
        }
    }

    // Aggregate QC
    ch_multiqc = channel.empty()
        .mix(FASTQC.out.zip.map { _meta, f -> f }.flatten())
        .mix(FASTP.out.json.map { _meta, f -> f })
        .mix(STAR_ALIGN.out.log_final.map { _meta, f -> f })
        .mix(ch_rnaseq_metrics.map { _meta, f -> f })
        .mix(PICARD_HSMETRICS.out.metrics.map { _meta, f -> f })
        .mix(PICARD_INSERTSIZE.out.metrics.map { _meta, f -> f })
        .mix(MOSDEPTH.out.summary.map { _meta, f -> f })
        .mix(MOSDEPTH.out.global_dist.map { _meta, f -> f })
        .mix(MOSDEPTH.out.dist.map { _meta, f -> f })
        .mix(COVERAGE_SUMMARY.out.mqc.map { _meta, f -> f })
        .mix(ch_rseqc.map { _meta, f -> f }.flatten())
        .mix(FEATURECOUNTS.out.summary.map { _meta, f -> f })
        .collect()

    MULTIQC(ch_multiqc, file("${projectDir}/assets/multiqc_config.yml"))

    workflow.onComplete = {
        log.info "Pipeline ${workflow.success ? 'succeeded' : 'FAILED'}. Results: ${params.outdir}"
    }
}
