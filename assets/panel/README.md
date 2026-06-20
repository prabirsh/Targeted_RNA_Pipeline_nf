# Panel asset files

Place your targeted-panel definition files here. These drive the on-target QC and
the reportable range. All files must use the **same reference build** as the STAR index
and GTF (GRCh38 recommended).

Required for on-target metrics (Picard CollectHsMetrics):

| File | Format | How to make it |
|------|--------|----------------|
| `panel_baits.interval_list`   | Picard interval_list | `picard BedToIntervalList I=baits.bed O=panel_baits.interval_list SD=genome.dict` |
| `panel_targets.interval_list` | Picard interval_list | `picard BedToIntervalList I=targets.bed O=panel_targets.interval_list SD=genome.dict` |
| `panel_targets.bed`           | BED (RSeQC)          | exported from your panel design |

Notes
- `baits` = the probe/capture footprint; `targets` = the regions you actually report on.
  For many panels these are the same; keep both for correct HsMetrics.
- `genome.dict` is the Picard sequence dictionary: `picard CreateSequenceDictionary R=genome.fa`.
- The BED for RSeQC must be 12-column (BED12) for `geneBody_coverage.py` / `inner_distance.py`.

This folder ships empty on purpose — panel designs are assay-specific and often proprietary.
