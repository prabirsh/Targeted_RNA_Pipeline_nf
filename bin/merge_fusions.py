#!/usr/bin/env python3
"""
Reconcile fusion calls from Arriba and STAR-Fusion into one ranked, clinical table.

For every fusion we report:
  - gene pair (left -- right)
  - LEFT breakpoint and RIGHT breakpoint (genomic coordinates)
  - split-read and spanning-read support
  - which callers found it and how many
  - the supporting READ NAMES (so a reviewer can inspect the exact reads)

Strategy: normalise each gene pair to an unordered key, union the calls, attach
evidence, rank by caller count then read support. 'reported' = calls promoted by
the rule (>=2 callers, or a single high-confidence Arriba call).
"""
import argparse
import csv
import sys


def norm_pair(a, b):
    return "--".join(sorted([a.strip(), b.strip()]))


def to_int(x, default=0):
    try:
        return int(float(x))
    except (ValueError, TypeError):
        return default


def split_reads(field):
    if not field or field in (".", "NA"):
        return set()
    return {r for r in field.replace(";", ",").split(",") if r and r != "."}


def parse_arriba(path):
    out = {}
    if not path:
        return out
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            g1 = row.get("#gene1") or row.get("gene1")
            g2 = row.get("gene2")
            if not g1 or not g2:
                continue
            key = norm_pair(g1, g2)
            split = to_int(row.get("split_reads1")) + to_int(row.get("split_reads2"))
            span = to_int(row.get("discordant_mates"))
            reads = split_reads(row.get("read_identifiers", ""))
            out[key] = {
                "left": row.get("breakpoint1", ""),
                "right": row.get("breakpoint2", ""),
                "arriba_split": split,
                "arriba_span": span,
                "confidence": row.get("confidence", ""),
                "type": row.get("type", ""),
                "reads": reads,
            }
    return out


def parse_starfusion(path):
    out = {}
    if not path:
        return out
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            name = row.get("#FusionName") or row.get("FusionName")
            if not name or "--" not in name:
                continue
            g1, g2 = name.split("--")[:2]
            key = norm_pair(g1, g2)
            reads = split_reads(row.get("JunctionReads", "")) | split_reads(row.get("SpanningFrags", ""))
            out[key] = {
                "left": row.get("LeftBreakpoint", ""),
                "right": row.get("RightBreakpoint", ""),
                "sf_split": to_int(row.get("JunctionReadCount")),
                "sf_span": to_int(row.get("SpanningFragCount")),
                "ffpm": row.get("FFPM", ""),
                "reads": reads,
            }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--arriba", required=True)
    ap.add_argument("--starfusion", default=None)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--reported", required=True)
    ap.add_argument("--min-callers", type=int, default=2)
    ap.add_argument("--max-readnames", type=int, default=50)
    args = ap.parse_args()

    arriba = parse_arriba(args.arriba)
    starf = parse_starfusion(args.starfusion)

    records = []
    for key in set(arriba) | set(starf):
        a = arriba.get(key, {})
        s = starf.get(key, {})
        callers = []
        if a:
            callers.append("arriba")
        if s:
            callers.append("starfusion")
        split = a.get("arriba_split", 0) + s.get("sf_split", 0)
        span = a.get("arriba_span", 0) + s.get("sf_span", 0)
        reads = sorted(a.get("reads", set()) | s.get("reads", set()))
        reads_capped = reads[:args.max_readnames]
        rec = {
            "sample": args.sample,
            "fusion": key,
            "left_breakpoint": a.get("left") or s.get("left", ""),
            "right_breakpoint": a.get("right") or s.get("right", ""),
            "split_reads": split,
            "spanning_reads": span,
            "total_support": split + span,
            "callers": ",".join(callers),
            "n_callers": len(callers),
            "arriba_confidence": a.get("confidence", ""),
            "type": a.get("type", ""),
            "starfusion_ffpm": s.get("ffpm", ""),
            "n_read_names": len(reads),
            "read_names": ",".join(reads_capped) if reads_capped else ".",
        }
        records.append(rec)

    records.sort(key=lambda r: (r["n_callers"], r["total_support"]), reverse=True)

    fields = ["sample", "fusion", "left_breakpoint", "right_breakpoint",
              "split_reads", "spanning_reads", "total_support", "callers",
              "n_callers", "arriba_confidence", "type", "starfusion_ffpm",
              "n_read_names", "read_names"]

    def write(path, rows):
        with open(path, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
            w.writeheader()
            for r in rows:
                w.writerow(r)

    def promote(r):
        if r["n_callers"] >= args.min_callers:
            return True
        if "arriba" in r["callers"] and r["arriba_confidence"].lower() == "high":
            return True
        return False

    write(args.merged, records)
    write(args.reported, [r for r in records if promote(r)])

    sys.stderr.write(
        "[merge_fusions] %s: %d pairs, %d reported\n"
        % (args.sample, len(records), sum(1 for r in records if promote(r)))
    )


if __name__ == "__main__":
    main()
