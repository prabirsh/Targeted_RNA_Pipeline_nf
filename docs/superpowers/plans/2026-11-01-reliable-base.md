# Reliable Base: Implementation Plan

- **Spec:** `docs/superpowers/specs/2026-11-01-targeted-rna-pipeline-design.md`, sub-project 1
- **Branch:** `feature/reliable-base`, created from `design/targeted-rna-pipeline`
- **Goal:** a pipeline that installs, runs end to end on synthetic data, has no known defects from spec section 12, and is covered by tests and CI.
- **Not in this plan:** `--assay_type`, panel checks, evidence levels, normal filter, nf-core module migration. Those belong to later sub-projects.

## Rules for every task

- Write the failing test first, see it fail, then fix.
- A suspected defect is reproduced before it is fixed. If it doesn't reproduce, record that and make no change.
- Commit after each task. Never commit `graphify-out/` or `PROJECT_CHANGES.md`.
- Docker isn't available on the development machine, so local runs use a new `conda` profile. Container runs are verified in CI, and any untested container path is reported as untested.

## Task 1: Safety and repository clean-up

Files: `scripts/git_init.sh`, `run.sh`, `docs/*.pdf`, `.gitignore`, `bin/*.py`

1. Delete `scripts/git_init.sh` and remove its README section.
2. Delete `run.sh` and remove its README references.
3. `git rm` the nine PDFs. Confirm `docs/REFERENCES.md` has a DOI for each.
4. Add `graphify-out/` and `PROJECT_CHANGES.md` to `.gitignore`.
5. Mark executable: `git update-index --chmod=+x bin/*.py`.
6. Check: `git ls-files -s bin` shows mode 100755.

## Task 2: Development environment

Files: `envs/dev.yml`, `envs/tools.yml`, `envs/star-fusion.yml`, delete `environment.yml`, `nextflow.config`

1. `envs/dev.yml`: Python, pytest, ruff, nextflow, nf-test.
2. `envs/tools.yml`: FastQC, fastp, STAR, samtools, Picard, mosdepth, RSeQC, Subread, Arriba, MultiQC, at the versions in spec section 3 where they solve together.
3. `envs/star-fusion.yml`: STAR-Fusion alone, which isolates its dependency conflict.
4. Add a `conda` profile that applies each environment by process.
5. Check: all three environments create successfully. Record any version that differs from the spec's table and why.

## Task 3: Samplesheet checks

Files: `tests/python/test_check_samplesheet.py`, `bin/check_samplesheet.py`

Failing tests, then fixes:

| Input | Expected |
|---|---|
| Row with an extra field | Clear error naming the line, no traceback |
| Byte-order mark in the header | Accepted |
| `fastq_1` equal to `fastq_2` | Error |
| Sample id containing `/` or other unsafe characters | Error |
| Valid paired and single-end rows | Accepted, unchanged output |

## Task 4: Fusion merging

Files: `tests/python/test_merge_fusions.py`, `bin/merge_fusions.py`

Failing tests, one per defect, then a rewrite of parsing and merging:

| Test | Expected |
|---|---|
| Two rows for one gene pair | Two events; high-confidence row survives and is reported |
| `TMPRSS2--ERG` | Order preserved |
| `TMPRSS2--ERG` and `ERG--TMPRSS2` | Separate events |
| Both callers, same breakpoint | `n_callers=2` |
| Both callers, breakpoints 11 bp or more apart | Separate events, each `n_callers=1` |
| Breakpoints within 10 bp | Treated as agreement |
| STAR-Fusion `chr1:100:-` and Arriba `chr1:100` | Same breakpoint |
| Shared read names | Support is the unique count, never the sum |
| Read names missing | Support is the larger caller's count, marked `support_source=caller_max` |
| 6000 read names in a field | No crash |
| Empty input from either caller | Valid empty outputs |
| Intergenic gene names | Kept verbatim, and don't crash |

Output columns change: per-caller split and spanning counts replace the summed `split_reads` and `spanning_reads`, and `total_support` holds the unique count. The promotion rule stays as it is, because evidence levels arrive in sub-project 3.

## Task 5: Coverage summary

Files: `tests/python/test_coverage_summary.py`, `bin/coverage_summary.py`

1. Failing test: flagstat with secondary and supplementary lines reports primary reads, and pairs where applicable.
2. Failing test: a real MultiQC run on the `_mqc` output shows one row per real sample and no phantom samples. Fix the file layout if this reproduces.

## Task 6: Configuration defects

Files: `nextflow.config`, `main.nf`, `conf/test.config`, `assets/multiqc_config.yml`

1. Remove the `panels` map and `read_length`.
2. `--panel` becomes a required BED path, and `--panel_bed` is removed.
3. Declare `aws_queue` and `aws_region`.
4. Make `--refflat` and `--rrna_intervals` truly optional: skip `PICARD_RNASEQMETRICS` without `refflat`, and omit ribosomal intervals without `rrna_intervals`.
5. Read the report's version from the manifest, not a hard-coded string.

## Task 7: Synthetic end-to-end data

Files: `tests/data/make_test_data.py`, `tests/data/README.md`, `conf/test.config`

1. A seeded generator writes a small FASTA, GTF, refFlat, rRNA intervals, BED, and paired-end and single-end FASTQs with planted fusions.
2. Generated output is deterministic, and small enough to keep in the repository or regenerate in CI.
3. STAR-Fusion's real CTAT library is too large to bundle, so the end-to-end test uses Arriba. STAR-Fusion gets parser tests with fixture files, plus a stub run. This limit is stated in the docs.
4. Check: `nextflow run main.nf -profile conda,test` starts and reaches its first real failure.

## Task 8: Resolve suspected module defects

Using the Task 7 run, test each suspicion from spec section 12, then fix only what reproduces:

| Module | Check | Fix if confirmed |
|---|---|---|
| `star_align.nf` | STAR accepts the chimeric settings | Adjust settings to Arriba and STAR-Fusion documentation |
| `arriba.nf` | `read_identifiers` populated | Add `-I` |
| `picard_hsmetrics.nf` | Runs with FASTA only | Stage `.fai` and `.dict` |
| `fastqc.nf`, `fastp.nf`, `star_align.nf` | Single-end sample runs | Normalise `reads` to a list |
| `fastqc.nf` | Two samples sharing a filename | Stage reads under sample-specific names |
| `multiqc_config.yml` | mosdepth section appears | Fix `module_order` |
| `picard_insertsize.nf` | Single-end sample completes | Skip the step for single-end |
| `rseqc.nf` | Failures are visible | Remove `|| true`, rename the misleading input |

Also apply these, which need no reproduction: `--limitBAMsortRAM` from the task's memory, `-Xmx` for Picard from the task's memory, and a `COVERAGE_CAP` high enough for deep panels.

Check: `nextflow run main.nf -profile conda,test` succeeds, with planted fusions in the merged table and read names present.

## Task 9: Workflow tests

Files: `tests/nf-test/`, `nf-test.config`, stub blocks in `modules/local/*.nf`

1. Add a `stub:` block to every process.
2. nf-test covers each module and the whole workflow with stubs.
3. An end-to-end nf-test uses the Task 7 data, and asserts planted fusions are found and outputs exist.
4. Save reference outputs for the end-to-end run.

## Task 10: Continuous integration

Files: `.github/workflows/ci.yml`

1. Ruff and pytest on every push.
2. nf-test stub runs on every push.
3. End-to-end run on every push, using containers.
4. A fresh-clone job, which proves that executable bits and paths work outside WSL.
5. Pin container digests in `conf/modules.config` and verify them here.

## Task 11: Documentation and release

Files: `README.md`, `CHANGELOG.md`, `VERSION`, `nextflow.config`, `docs/PIPELINE_FLOW_EXPLAINED.md`, `docs/pipeline_*.html`, `assets/panel/README.md`

1. Correct panel usage, the install and the output layout. Remove `example_output/` and `software versions` claims unless they're true.
2. Add an honest status line: research use, not validated.
3. Bring the flowcharts in line with the code, or mark them out of date.
4. Move version to 1.2.0 in all three version locations, and add a changelog entry that lists output column changes.

## Done when

- [ ] Every defect in spec section 12 is fixed, or recorded as not reproduced with evidence
- [ ] `pytest` passes
- [ ] nf-test passes
- [ ] `nextflow run main.nf -profile conda,test` succeeds from a fresh clone
- [ ] CI passes, once the user chooses to push the branch
- [ ] The user runs `/code-review ultra`, and its findings are handled
