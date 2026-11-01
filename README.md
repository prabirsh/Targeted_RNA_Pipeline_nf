# rna-panel-pipeline (v1.1)

Single-command, integrated **targeted-panel RNA-seq pipeline** for **clinical FFPE oncology**
samples (DRAGEN-RNA-style). Built on Nextflow DSL2. You pick a panel, point it at FASTQs +
references, and it runs sample QC, alignment, **panel coverage QC**, **gene-fusion detection**
(with breakpoints + read names), expression, and one aggregated QC report.



## What it produces

- **Coverage QC** (`qc/coverage/*.coverage_summary.tsv`): total reads, mean target coverage,
  and % of target bases at **50x / 100x / 200x / 500x / 1000x** (mosdepth).
- **Insert size** (Picard) and **on-target %** (Picard HsMetrics) for the chosen panel.
- **Fusions** (`fusion/merged/*.fusions.reported.tsv`): gene pair, **left/right breakpoints**,
  split + spanning **read support**, and supporting **read names** (Arriba + STAR-Fusion, reconciled).
- **Expression** (`expression/*.featureCounts.txt`).
- **One MultiQC report** (`multiqc/multiqc_report.html`).

## Panels (pick one per run)

| `--panel` | BED | Use |
|-----------|-----|-----|
| `panel2` | `panel2.bed` |  exome |
| `panel1`  | `panel1.bed` |  RNA Fusion (GRCh37) |

Pass `--panel panel1|panel2` (or `--panel_bed <file>`).
Run with no panel → the pipeline stops and lists the options (never guesses).

## Install (any system)

```bash
conda env create -f environment.yml
conda activate rna-panel-pipeline
```
Nextflow then pulls per-process containers (or use the conda tools directly).

## Run

```bash
nextflow run main.nf -profile docker \
  --input samplesheet.csv --outdir results --panel panel1 \
  --star_index /refs/star_index --gtf /refs/annotation.gtf --fasta /refs/genome.fa \
  --refflat /refs/refFlat.txt --rrna_intervals /refs/rRNA.interval_list \
  --ctat_lib /refs/ctat_genome_lib --arriba_blacklist /refs/arriba/blacklist.tsv.gz
```
`nextflow run main.nf --help` lists all options. Profiles: `docker`, `singularity`, `slurm`,
`awsbatch`, `test`.

## Reference bundle

| Param | What |
|-------|------|
| `--star_index` `--gtf` `--fasta` | genome + annotation (same build as the panel BED) |
| `--dict` | sequence dictionary (built from `--fasta` if omitted) |
| `--refflat` `--rrna_intervals` | Picard RnaSeqMetrics |
| `--ctat_lib` | STAR-Fusion CTAT library (match genome build) |
| `--arriba_blacklist` `--arriba_known` `--arriba_protein` | Arriba references |
| `--gene_bed` | BED12 gene model for RSeQC (optional) |

The Panel 1 panel is **GRCh37** (chr-prefixed). For GRCh37, use the GRCh37 CTAT lib and GRCh37 Arriba
files, and keep contig naming consistent across all references. No code changes needed for build.

## Samplesheet

```csv
sample,fastq_1,fastq_2,strandedness
TUMOR_001,/data/TUMOR_001_R1.fastq.gz,/data/TUMOR_001_R2.fastq.gz,reverse
```

## Output layout

See `example_output/` for a complete dummy run (synthetic, mirrors the real `results/` tree).

```
results/
├── qc/fastqc · qc/fastp · qc/picard · qc/coverage · qc/rseqc
├── star/<sample>/        BAM + index + logs + junction
├── expression/           featureCounts matrix
├── fusion/{arriba,starfusion,merged}/
├── multiqc/              multiqc_report.html
└── pipeline_info/        trace / timeline / DAG / software versions
```

## Clinical note

Research-grade; intended to progress to formal analytical validation (accuracy, LoD,
reproducibility, FFPE robustness). Pin container digests and lock QC thresholds before clinical
sign-out. The pipeline produces evidence; interpretation remains lab-director-owned.
