# Changelog

All notable changes to the RNA panel pipeline are documented here.
This project follows [Semantic Versioning](https://semver.org).

## [1.1.0] - 2026-06-20
### Added
- Panel selection: `--panel ST8|CT` (or `--panel_bed`), resolved from a panel map.
  Interactive `run.sh` prompts for the panel at launch.
- Coverage QC via mosdepth: target coverage and % of target bases at
  50x / 100x / 200x / 500x / 1000x, summarised by `coverage_summary.py`.
- Insert-size metrics (Picard CollectInsertSizeMetrics).
- Panel BED -> Picard interval_list conversion (`PANEL_TO_INTERVALS`) and
  on-the-fly sequence dictionary (`SEQUENCE_DICT`).
- Fusion output enriched with LEFT/RIGHT breakpoints, split/spanning read support,
  and supporting READ NAMES (Arriba `read_identifiers` + STAR-Fusion read columns).
- `environment.yml` for one-command install on any system.
- `REFERENCES.md` (PubMed-sourced citations with DOIs).
- Redesigned data-flow flowchart and full pipeline documentation.
- `example_output/` dummy run so the output layout is self-explanatory.

### Changed
- `CollectHsMetrics` now derives bait/target intervals from the selected panel.
- MultiQC intake extended with coverage, insert-size and HsMetrics.
- Version bumped to 1.1.0.

## [1.0.0] - 2026-06-19
### Added
- Initial targeted-panel RNA-seq pipeline: samplesheet validation, FastQC, fastp,
  STAR 2-pass, samtools, Picard RnaSeqMetrics + HsMetrics, RSeQC, featureCounts,
  Arriba + STAR-Fusion ensemble with reconciliation, MultiQC.
