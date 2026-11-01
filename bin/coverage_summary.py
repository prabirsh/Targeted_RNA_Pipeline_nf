#!/usr/bin/env python3
"""
Summarise targeted-panel coverage into one clinical QC row.

Inputs:
  --thresholds  mosdepth *.thresholds.bed.gz  (per-region bases >= each cutoff)
  --summary     mosdepth *.mosdepth.summary.txt (mean depth per region)
  --flagstat    samtools flagstat (total reads)

Outputs:
  --out  per-sample TSV: total reads, target bases, mean target coverage,
         and % of target bases at >= 50x / 100x / 200x / 500x / 1000x
  --mqc  a MultiQC-friendly key:value block of the same numbers
"""
import argparse
import gzip
import sys

CUTOFFS = ["50", "100", "200", "500", "1000"]


def read_thresholds(path):
    """mosdepth thresholds.bed[.gz]: chrom start end region 50X 100X 200X 500X 1000X"""
    opener = gzip.open if path.endswith(".gz") else open
    total_bases = 0
    ge = {c: 0 for c in CUTOFFS}
    header = None
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if line.startswith("#") or cols[0] == "#chrom" or cols[1] == "start":
                header = cols
                continue
            start, end = int(cols[1]), int(cols[2])
            total_bases += end - start
            # columns 4.. correspond to the requested cutoffs in order
            for idx, c in enumerate(CUTOFFS):
                col = 4 + idx
                if col < len(cols):
                    ge[c] += int(cols[col])
    return total_bases, ge


def read_mean_cov(path):
    """mosdepth summary: use the 'total_region' row mean if present, else 'total'."""
    mean_region = None
    mean_total = None
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            if cols[0] == "total_region":
                mean_region = cols[3]
            elif cols[0] == "total":
                mean_total = cols[3]
    return mean_region or mean_total or "NA"


def read_total_reads(path):
    with open(path) as fh:
        first = fh.readline()
    # "<N> + <M> in total (QC-passed reads + QC-failed reads)"
    try:
        return int(first.split()[0])
    except (ValueError, IndexError):
        return "NA"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--thresholds", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--flagstat", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mqc", required=True)
    args = ap.parse_args()

    total_bases, ge = read_thresholds(args.thresholds)
    mean_cov = read_mean_cov(args.summary)
    total_reads = read_total_reads(args.flagstat)

    pct = {}
    for c in CUTOFFS:
        pct[c] = round(100.0 * ge[c] / total_bases, 2) if total_bases else 0.0

    headers = ["sample", "total_reads", "target_bases", "mean_target_coverage"] + \
              [f"pct_bases_ge_{c}x" for c in CUTOFFS]
    values = [args.sample, total_reads, total_bases, mean_cov] + [pct[c] for c in CUTOFFS]

    with open(args.out, "w") as out:
        out.write("\t".join(headers) + "\n")
        out.write("\t".join(str(v) for v in values) + "\n")

    # MultiQC custom-content (key-value) block
    with open(args.mqc, "w") as out:
        out.write("# plot_type: 'generalstats'\n")
        out.write(f"Sample\t{args.sample}\n")
        out.write(f"total_reads\t{total_reads}\n")
        out.write(f"mean_target_coverage\t{mean_cov}\n")
        for c in CUTOFFS:
            out.write(f"pct_bases_ge_{c}x\t{pct[c]}\n")

    sys.stderr.write(
        f"[coverage_summary] {args.sample}: mean={mean_cov}x, "
        f">=100x={pct['100']}%, >=500x={pct['500']}%\n"
    )


if __name__ == "__main__":
    main()
