# Targeted RNA Panel Pipeline: Design

- **Date:** 2026-11-01
- **Status:** Awaiting user review of this written spec
- **Base:** `rna-panel-pipeline` v1.1.0 (commit `d1db7fb`), rebuilt in this repository

## 1. Goal

A Nextflow pipeline that analyses data from any targeted RNA panel the user supplies and produces gene expression and RNA fusion results with fewer false positives, measured evidence of performance, and full documentation.

### Success criteria

1. A new user can run the bundled test data from a fresh clone with one documented command.
2. A user can analyse a panel by supplying only a samplesheet, a panel BED and references. No panel name is built into the code.
3. On held-out truth samples, the validation report compares this pipeline directly with plain Arriba and plain STAR-Fusion. "Fewer false positives" may be claimed only if that report shows it, with sensitivity reported alongside.
4. Every reported fusion carries its evidence, and every removed or downgraded call keeps its reason.
5. Each negative fusion result states whether the gene could be assessed in that sample.
6. Every run records the versions, digests, checksums and parameters needed to reproduce it.

### Claims policy

- "Validated" means measured against truth data with the lab's acceptance limits. The software reports measurements and never declares an assay validated. The lab director signs off.
- A feature may be described as new only after a documented check against current tools. A search during this design found no equivalent for the candidates in section 4, but it read summaries only and is not proof.
- Each output type carries a status. It stays `research-only` until measured against truth data for that assay type.

## 2. Scope

### In scope, in build order

| # | Sub-project | Outcome |
|---|---|---|
| 1 | Reliable base | Confirmed defects fixed, working install, tests and CI, current tool versions |
| 2 | Any-panel support | User-supplied panels, input checks, panel-aware QC and expression |
| 3 | Fusion evidence | Breakpoint merging, FusionInspector, normal filter, carryover check, evidence levels, calibration, assessability |
| 4 | Validation and documentation | `validate` mode, reports, full documentation |
| 5 | Splicing | CTAT-Splicing with a normal-sample filter |
| 6 | Small variants | RNA small variant calling, `research-only` until measured |

Each sub-project gets its own implementation plan and is finished, tested and reviewed before the next starts. The first plan covers sub-project 1 only.

### Out of scope

- Clinical interpretation and clinical report generation.
- Direct expression comparison between different panels.
- Trained machine learning models. Rules come first, and a model may be considered later only with enough truth samples held back for independent testing.
- Any built-in panel definitions. Panel designs stay outside the repository.

## 3. Architecture

### Inputs

| Input | Required | Notes |
|---|---|---|
| `--input` samplesheet | Yes | `sample,fastq_1,fastq_2,strandedness`; strandedness defaults to `auto` |
| `--panel` BED | Yes | Path to the user's target BED |
| `--baits` BED | No | Defaults to the panel BED |
| `--assay_type` | Yes | `capture` or `amplicon` |
| References | Yes | FASTA with `.fai` and `.dict`, GTF, STAR index, refFlat, rRNA intervals, CTAT library, Arriba files; `--genome` is planned as a shortcut for standard bundles |
| `--pon` | No | Normal filter built by `build_pon` |
| `--calibration` | No | Thresholds built by `calibrate` |
| `--known_fusions` | No | Known-driver list |
| `--umi` | No | Opt-in UMI handling |

### Modes

`--mode run` (default), `build_pon`, `calibrate`, `validate`.

### Data flow

```
samplesheet + panel BED + references
        |
PREPARE_PANEL   check BED against FASTA (contigs, coordinates) and GTF (build, gene overlap);
                stop with a clear error on mismatch; build intervals and the panel gene list
        |
READ_QC         FastQC, fastp, optional UMI extraction
        |
ALIGN           STAR, index, duplicate marking, strandedness detection
        |
   +----+----------------+------------------+
PANEL_QC            EXPRESSION           FUSION
   +----+----------------+------------------+
        |
RUN_REPORT      MultiQC, results report, pass/fail gates, run record
        |
VALIDATE        optional, scores a run against truth data
```

### Code structure

- `main.nf` stays thin, and calls `workflows/`, then `subworkflows/` and modules.
- Standard tools use pinned nf-core modules under `modules/nf-core/`. Project logic goes under `modules/local/` and tested Python in `bin/`.
- nf-schema validates parameters before any task runs. There are no interactive prompts.
- Docker, Singularity and Conda profiles are provided. Each tool has its own environment, replacing the single `environment.yml`.

### Components removed

- Hard-coded `ST8` and `CT` panels and their documentation.
- `run.sh`, the interactive launcher.
- `scripts/git_init.sh`, which deletes `.git`.
- The nine research PDFs in `docs/`, replaced by the DOI links in `docs/REFERENCES.md`, unless the user asks to keep them.

### Tool versions

Move to current stable releases during sub-project 1 and pin container digests. Bioconda latest versions seen on 2026-11-01:

| Tool | Current pin | Latest seen |
|---|---|---|
| STAR | 2.7.11b | 2.7.11b |
| fastp | 0.23.4 | 1.4.0 |
| FastQC | 0.12.1 | 0.13.0 |
| samtools | 1.20 | 1.25 |
| Picard | 3.1.1 | 3.5.0 |
| mosdepth | 0.3.8 | 0.3.14 |
| RSeQC | 5.0.3 | 5.0.6 |
| Subread | 2.0.6 | 2.1.1 |
| Arriba | 2.4.0 | 2.5.1 |
| STAR-Fusion | 1.13.0 | 1.15.1 |
| MultiQC | 1.21 | 1.35 |
| Nextflow | 23.10.1 | 26.10.4 |

A version change is accepted only if tests pass and saved reference outputs are reviewed.

## 4. Fusion evidence

### Merging

- A fusion event is identified by ordered 5' gene, 3' gene, breakpoints and strand.
- Each distinct breakpoint keeps its own row. Events are also grouped by gene pair for reading.
- Callers agree only when breakpoints match within 10 bp (`--fusion_breakpoint_window`).
- Support is counted from unique read names across callers. Per-caller counts are kept alongside.
- Duplicates are excluded from support when UMIs or duplicate marks are available.

### Evidence recorded for each event

| Group | Fields |
|---|---|
| Reads | Unique split reads, unique spanning reads, anchor length, mismatches |
| Biology | Splice-site status, reading frame, partner gene sequence similarity |
| Confirmation | Caller agreement, FusionInspector result |
| Panel | Whether each gene is on the panel, breakpoint coverage |
| Sample quality | 5'/3' imbalance, degradation level |
| Artifact checks | Frequency in normals, same breakpoint elsewhere in the run |

### Normal filter

- `build_pon` uses fusion-negative normal samples from the same panel and records breakpoint frequency.
- The filter file records panel checksum, genome build, assay type and number of normals. A mismatch stops the run.
- A breakpoint found in at least 2 normals is flagged as a recurrent artifact (provisional default).
- A fusion on the known-driver list is never removed by this filter, only flagged.

### Same-run carryover check

A breakpoint seen in more than one sample of a run is flagged, with extra weight when one of them is a strong positive control. This is a flag only, because real recurrent fusions often share exon-exon breakpoints.

### Evidence levels

These are evidence levels and not AMP/ASCO/CAP clinical tiers. Read thresholds below are provisional defaults and results are stamped `uncalibrated` until replaced by calibration.

| Level | Rule |
|---|---|
| **Strong** | At least 5 unique reads and no artifact flag, with either two agreeing callers or one caller plus FusionInspector confirmation |
| **Moderate** | One high-confidence caller with at least 5 unique reads and no artifact flag, or a known driver with at least 3 unique reads |
| **Review** | A known driver with 1 or 2 unique reads, or a call otherwise meeting Strong or Moderate that has a conflicting flag |
| **Filtered** | Anything else; stays in the full table with reasons |

### Calibration

`calibrate` uses the calibration sample set to choose thresholds per panel and assay type. It writes a versioned file, and results are stamped with that version.

### Assessability

Each panel gene receives a label for each sample, using duplicate-free coverage over its targets:

| Label | Provisional rule |
|---|---|
| Assessable | At least 80% of target bases at 50x or more, and the sample passes gates |
| Limited | At least 20% of target bases at 50x or more |
| Not assessable | Anything lower, or a failed sample |

A negative result is reported as "no fusion detected, assessable" or "cannot rule out", with the reason. Dilution data replaces these rules during calibration.

### Assay types

`capture` and `amplicon` are calibrated and validated separately. Amplicon data can't use read start positions to identify duplicates, so duplicate logic depends on assay type and UMIs. An assay type with no truth data is marked unvalidated in its results.

## 5. Expression

- featureCounts gene-level counts, with strandedness resolved by detection plus the samplesheet.
- A contradiction between detected and declared strandedness produces a warning and is recorded in the report.
- Outputs: per-sample counts, a cohort count matrix, and panel-restricted counts per million.
- Each gene is marked on-panel or off-panel, so capture bias isn't mistaken for biology.

## 6. Panel QC and run gates

- Reported metrics: on-target rate, coverage at 50x, 100x, 200x, 500x and 1000x with and without duplicates, insert size, mapping, rRNA, 5'/3' bias.
- `CollectHsMetrics` is configured so its coverage cap doesn't limit deep panel data.
- Default gates, all marked provisional until the lab sets its own:

| Metric | Warn | Fail |
|---|---|---|
| Uniquely mapped reads | Below 70% | Below 50% |
| Total read pairs | Below 5 million | Below 1 million |
| On-target rate | Below 40% | Below 20% |
| Target bases at 50x or more | Below 80% | Below 50% |

These are engineering starting points, not evidence-based limits.

## 7. Error handling

- Invalid input or a reference mismatch stops the run before alignment, and the message says how to fix it.
- A sample failing a gate is flagged and stays visible in every report.
- Tool failures are not hidden. The existing `|| true` handling in RSeQC is removed, and only genuinely optional outputs are marked optional.
- Java tools receive memory based on the task's allocation.

## 8. Validation

`--mode validate` takes a run's results, a truth file and acceptance limits.

### Truth file

One row per expected result: sample, set (`calibration` or `validation`), expected 5' gene, 3' gene, breakpoints, known level, and explicit negative samples.

### Sample separation

Calibration and validation sets are declared up front. Scoring refuses any sample used for calibration.

### Fusion report

- Sensitivity, specificity and positive predictive value, always with 95% confidence intervals.
- Scored at gene-pair level and at breakpoint level.
- Broken down by evidence level and caller.
- Limit of detection from dilution series.
- Repeatability within a run and reproducibility between runs.
- Head-to-head comparison with plain Arriba and plain STAR-Fusion on the same samples.
- A failure review listing each missed or false call with its evidence.

### Expression report

- Replicate correlation.
- Agreement with an independent measurement such as qPCR when supplied.
- On-panel and off-panel results shown separately.

### Truth data available

The user confirmed access to commercial reference material and orthogonally confirmed samples, including negatives. Sub-project 3 calibration and sub-project 4 scoring depend on these samples being sequenced and supplied.

## 9. Testing

| Layer | Tool | Covers |
|---|---|---|
| Unit | pytest | All Python logic, with one regression test per confirmed defect |
| Module and workflow | nf-test | Each step, and full runs with stubbed tools |
| End to end | Small synthetic genome and reads in the repository | Real tools on planted fusions, single-end and paired-end |
| Automation | GitHub Actions | Linting and tests on every push, install and fresh-clone checks |

Saved reference outputs detect unintended result changes after updates.

## 10. Run record

Each run writes tool versions, container digests, reference and panel checksums, parameters, calibration version, normal filter version and evidence-level status.

## 11. Documentation

- Quick start, including a first run on bundled test data.
- Preparing a panel, references, normals and calibration samples.
- Tool and parameter rationale, with citations.
- Definition of every output column.
- Known limits.
- A step-by-step validation guide.

## 12. Known defects to fix in sub-project 1

Reproduced by running code:

| Location | Defect |
|---|---|
| `bin/merge_fusions.py:53`, `:78` | The last row for a gene pair overwrites earlier rows |
| `bin/merge_fusions.py:22` | Alphabetical sorting loses 5'/3' order and merges reciprocal fusions |
| `bin/merge_fusions.py:112` | Support is double counted across callers |
| `bin/merge_fusions.py:43` | Long read-name fields exceed the CSV field limit and crash |
| `bin/check_samplesheet.py:38` | A row with an extra field crashes with an unclear error |
| `environment.yml` | `star-fusion=1.13.0` and `samtools=1.20` can't be installed together |

Confirmed by reading code or file modes:

| Location | Defect |
|---|---|
| `scripts/git_init.sh:11` | Deletes `.git` when it exists, and the README recommends running it |
| `bin/*.py`, `run.sh` | Committed as non-executable (mode 100644) |
| `nextflow.config:26` | Panel files are absent, and README panel names don't match the code |
| `main.nf:109` | `--refflat` and `--rrna_intervals` are documented as optional but used unconditionally |
| `bin/coverage_summary.py:68` | Total reads includes secondary and supplementary alignments |
| `assets/multiqc_config.yml:8` | Report says v1.0 |
| `README.md`, `CHANGELOG.md` | `example_output/` is referenced but missing |
| `nextflow.config:48` | `read_length` is never used |

Suspected, and to be tested first:

| Location | Suspicion |
|---|---|
| `modules/local/arriba.nf:23` | Arriba needs `-I` to output read names |
| `modules/local/star_align.nf:33`, `:39` | STAR may reject these chimeric settings together |
| `modules/local/picard_hsmetrics.nf:21` | `.fai` and `.dict` aren't staged, which may fail every sample |
| `bin/coverage_summary.py:101` | The MultiQC file layout is transposed and likely creates phantom samples |
| `modules/local/fastp.nf:19`, `fastqc.nf:14`, `star_align.nf:19` | Single-end input may render a wrong path |
| `modules/local/fastqc.nf` | Samples sharing a FASTQ filename may collide in MultiQC |
| `assets/multiqc_config.yml:11` | `module_order` may exclude mosdepth from the report |

The current pipeline has never been run end to end in this review, because its install fails.

## 13. Risks

| Risk | Mitigation |
|---|---|
| Thresholds tuned to too few samples | Separate validation set, and confidence intervals always shown |
| A contaminated normal hides a real driver | Known drivers are never removed by the normal filter |
| Filters remove true low-level fusions | Sensitivity is reported beside false positives, and filtered calls stay visible with reasons |
| Scope grows too large to finish | Six independent sub-projects, each usable when complete |
| A tool update silently changes results | Pinned digests and saved reference outputs |
| RNA variant artifacts from editing and FFPE damage | Small variants stay `research-only` until measured |
