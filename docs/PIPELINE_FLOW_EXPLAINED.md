# RNA Panel Pipeline — Full Documentation (v1.1)

Clinical, targeted-panel RNA-seq for FFPE oncology samples. This document explains the
**whole pipeline**: for every tool, **why** it is here, **what** data it generates, and
**how** that data is used next. Read it with the visual `pipeline_flowchart.html`.

One sample (FASTQ) goes in; out come **gene fusions** (with breakpoints + read names),
**gene expression**, a **coverage QC table** (depth at 50/100/200/500/1000x), and a single
**MultiQC** report telling you whether to trust the run.

---

## 0. Panel selection (which panel?)

This is a clinical, multi-panel lab, so the run starts by choosing a panel:

| Name | BED file | Use |
|------|----------|-----|
| **ST8** | `SureSelectXT_V8_Covered.bed` | SureSelect XT V8 exome panel (~215k regions) |
| **CT**  | `Indiegene_Target_2109PD008-V1_4BaseCare_RNA_Fusion_GRCh37.bed` | 4baseCare RNA Fusion panel, GRCh37 (~5.5k regions) |

- `run.sh` **asks** at launch (menu: 1 = ST8, 2 = CT).
- Or pass `--panel ST8` / `--panel CT` directly; or `--panel_bed <file>` to override.
- If neither is given the pipeline **stops** with a message listing the panels — it never
  guesses a panel.
- The chosen BED is converted to a Picard `interval_list` (`PANEL_TO_INTERVALS`, using a
  sequence dictionary from `SEQUENCE_DICT`) and drives the two panel-specific QC steps.

---

## At a glance: tool → file → who uses it

| Stage | Tool | Generates | Used by |
|-------|------|-----------|---------|
| 0 | check_samplesheet.py | `samplesheet.valid.csv` | the workflow |
| 0 | SEQUENCE_DICT / PANEL_TO_INTERVALS (Picard) | `genome.dict`, `panel.interval_list` | Picard HsMetrics |
| 1 | FastQC | `*_fastqc.zip` | MultiQC |
| 1 | fastp | `*.trim.fastq.gz`, `*.fastp.json` | STAR; MultiQC |
| 2 | STAR (2-pass) | `Aligned.bam`, `Chimeric.out.junction`, `Log.final.out` | everything downstream |
| 2 | samtools | `.bai`, `.flagstat` | Picard/mosdepth; coverage_summary |
| 3 | Picard CollectRnaSeqMetrics | `rna_metrics.txt` | MultiQC |
| 3 | Picard CollectHsMetrics ★ | `hs_metrics.txt` | MultiQC |
| 3 | Picard CollectInsertSizeMetrics | `insert_size_metrics.txt` | MultiQC |
| 3 | mosdepth ★ | `*.thresholds.bed.gz`, `*.mosdepth.summary.txt` | coverage_summary; MultiQC |
| 3 | coverage_summary.py | `*.coverage_summary.tsv` | clinical QC |
| 4 | featureCounts | `counts.txt` (+ `.summary`) | expression deliverable; MultiQC |
| 5 | Arriba | `arriba.fusions.tsv` | merge_fusions |
| 5 | STAR-Fusion | `starfusion.fusions.tsv` | merge_fusions |
| 5 | merge_fusions.py | `fusions.reported.tsv` | **clinical deliverable** |
| 9 | MultiQC | `multiqc_report.html` | **QC sign-off** |

★ = panel-specific.

---

## 1. Intake & references

### check_samplesheet.py
**Why.** A clinical pipeline must reject bad input before spending compute.
**Generates.** `samplesheet.valid.csv` (unique IDs, real FASTQ paths, valid strandedness).
**Used next.** Builds the per-sample read channel.

### SEQUENCE_DICT + PANEL_TO_INTERVALS (Picard)
**Why.** Picard CollectHsMetrics needs the panel as an `interval_list`, which needs a genome
sequence dictionary. This lets the lab supply only a **BED per panel** and the pipeline builds
the rest.
**Generates.** `genome.dict` (if `--dict` not given) and `panel.interval_list`.
**Used next.** The interval_list is the bait+target for on-target metrics.

---

## 2. Raw QC & trimming — FastQC + fastp

**Why.** FFPE RNA is degraded and adapter-heavy; you must see raw quality and clean the reads.
**FastQC generates.** Per-FASTQ quality report (`*_fastqc.zip/.html`) — does not change reads.
**fastp generates.** Trimmed reads `*.trim.fastq.gz` (adapter + polyG removed, short reads
dropped) and `*.fastp.json` (incl. **total input reads** before/after).
**Used next.** Trimmed FASTQs feed STAR; the zip + json go to MultiQC. Trimming is gentle so
fusion breakpoints are not lost.

---

## 3. Alignment — STAR (2-pass) + samtools

**Why STAR.** RNA reads cross exon-exon junctions; a splice-aware aligner is required, and STAR
is also what finds the chimeric (fusion) reads.
**Generates.** `Aligned.sortedByCoord.out.bam` (the hub file, read by six tools),
`Chimeric.out.junction` (fusion evidence for STAR-Fusion), `Log.final.out` (mapping stats),
and per-gene counts. 2-pass improves junction + breakpoint sensitivity.
**samtools.** `index` makes the `.bai` (random access for Picard/mosdepth); `flagstat` gives
mapped/paired counts and the **total read count** used in the coverage summary.

---

## 4. Panel & coverage QC (the metrics you asked for)

### Picard CollectRnaSeqMetrics
**What/why.** RNA composition (exonic/intronic/intergenic) and 5′→3′ bias. High intronic =
DNA contamination; strong 3′ bias = degradation. → `rna_metrics.txt` → MultiQC.

### Picard CollectHsMetrics ★ — on-target / target coverage
**What/why.** The metric that makes this a *panel* assay: **% of reads on the panel** and
**mean target coverage**. → `hs_metrics.txt` → MultiQC.

### Picard CollectInsertSizeMetrics — insert length
**What/why.** Fragment-length (insert size) distribution. A collapsed insert size flags
over-fragmented FFPE input. → `insert_size_metrics.txt` → MultiQC.

### mosdepth ★ — coverage at 50/100/200/500/1000x
**What/why.** Depth over the panel BED. With `--thresholds 50,100,200,500,1000` it reports,
per target, the number of bases reaching each depth. → `*.thresholds.bed.gz` +
`*.mosdepth.summary.txt`.

### coverage_summary.py — one clinical QC row
**What/why.** Rolls mosdepth + flagstat into a single tidy line: **total reads**,
**target bases**, **mean target coverage**, and **% of target bases at ≥50x / 100x / 200x /
500x / 1000x**. → `*.coverage_summary.tsv` and a MultiQC block. This is the coverage table a
reviewer reads first.

---

## 5. Expression — featureCounts

**Why.** Gene-level counts for outlier-overexpression and signatures.
**Generates.** `counts.txt` (a deliverable) and `.summary` (assignment stats → MultiQC).

---

## 6. Fusion detection — Arriba + STAR-Fusion → merge_fusions.py

Fusions drive targeted therapy, so two orthogonal callers run and are reconciled.

### Arriba
**Why.** Fast, FFPE-tolerant; reads STAR's chimeric reads from the **main BAM**. Its output
includes a `read_identifiers` column — the **names of the reads** supporting each fusion.
→ `arriba.fusions.tsv`.

### STAR-Fusion
**Why.** Orthogonal method on `Chimeric.out.junction` + CTAT reference; the **full** prediction
file lists `JunctionReads` / `SpanningFrags` (read names) and left/right breakpoints.
→ `starfusion.fusions.tsv`.

### merge_fusions.py — the clinical fusion table
**What/why.** Normalises each gene pair to an unordered key (so `EML4-ALK` = `ALK--EML4`),
unions the callers, and emits one ranked table with exactly what you asked for:

| Column | Meaning |
|--------|---------|
| `fusion` | gene pair (left -- right) |
| `left_breakpoint` / `right_breakpoint` | genomic breakpoints |
| `split_reads` / `spanning_reads` / `total_support` | read support |
| `callers` / `n_callers` | which tools, how many |
| `arriba_confidence`, `type`, `starfusion_ffpm` | annotation |
| `n_read_names` / `read_names` | supporting **read names** (capped at 50) |

`*.fusions.reported.tsv` = promoted calls (≥2 callers, or a high-confidence single Arriba call);
`*.fusions.merged.tsv` = everything for review. See `example_output/` for a worked ALK--EML4.

---

## 9. Aggregate QC — MultiQC

**Why.** Eight tools emit QC fragments; a reviewer should open one file.
**Generates.** `multiqc_report.html` combining FastQC, fastp, STAR, all three Picard modules,
mosdepth coverage, coverage_summary, RSeQC (if enabled) and featureCounts.
**Used next.** The QC sign-off gate — if on-target %, coverage @100x, mapping rate or insert
size fail, the sample is flagged before fusion/expression results are trusted.

---

## QC cheat-sheet (tune at validation)

| Metric | Source | Watch for |
|--------|--------|-----------|
| Total input reads | fastp / flagstat | below target depth |
| Insert size (median) | Picard InsertSize | very short = over-fragmented FFPE |
| % uniquely mapped | STAR Log.final | low = contamination/degradation |
| On-target % | Picard HsMetrics ★ | low = capture failure |
| Mean target coverage | mosdepth / HsMetrics | below assay minimum |
| % bases ≥100x / ≥500x | coverage_summary ★ | low = under-covered targets |
| 5′-3′ bias | Picard RnaSeqMetrics | high = degradation |
| Positive-control fusion | merge_fusions | not recovered = batch invalid |

---

## Reference genome note (GRCh37 / GRCh38)

The pipeline logic is genome-agnostic — all references are paths. The bundled CT panel is
**GRCh37, chr-prefixed**; ST8 is also chr-prefixed. Supply a self-consistent bundle: same
contig naming across STAR index, GTF, FASTA, dict, panel BED, and fusion references; the
GRCh37 CTAT library for STAR-Fusion; and the GRCh37 Arriba files. See `README.md`.

## Tool citations

Every method is cited with PubMed PMIDs and DOIs in `docs/REFERENCES.md`.
