#!/usr/bin/env python3
"""
Reconcile Arriba and STAR-Fusion calls into one ranked table of fusion events.

An event is one breakpoint pair for an ordered 5'--3' gene pair. Every breakpoint
keeps its own row. Two callers agree only when the gene order matches and both
breakpoints fall within --breakpoint-window bases.

Both callers read the same alignments, so support is never summed across them.
total_support is the number of unique read names. If a caller reports support
without read names, it is the larger caller's count instead, and support_source
records which method was used. Per-caller counts are always kept.

'reported' holds events found by at least --min-callers callers, or a single
high-confidence Arriba call.
"""
import argparse
import csv
import re
import sys

FIELDS = [
    "sample", "fusion", "left_breakpoint", "right_breakpoint", "callers", "n_callers",
    "total_support", "support_source", "arriba_split", "arriba_span", "sf_split", "sf_span",
    "arriba_confidence", "type", "starfusion_ffpm", "n_read_names", "read_names",
]


def to_int(value, default=0):
    try:
        return int(float(value))
    except (ValueError, TypeError):
        return default


def parse_breakpoint(value):
    """'chr1:100' or 'chr1:100:+' -> ('chr1', 100). Returns None if unparseable."""
    parts = (value or "").strip().split(":")
    if len(parts) < 2 or not parts[1].isdigit():
        return None
    return parts[0], int(parts[1])


def split_reads(field):
    """Read names from a comma or semicolon list, in a form comparable across callers."""
    if not field or field in (".", "NA"):
        return set()
    names = set()
    for name in re.split(r"[,;]", field):
        name = re.sub(r"/[12]$", "", name.strip().lstrip("@"))
        if name and name != ".":
            names.add(name)
    return names


def read_rows(path):
    if not path:
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))


def parse_arriba(path):
    events = []
    for row in read_rows(path):
        gene1 = row.get("#gene1") or row.get("gene1")
        gene2 = row.get("gene2")
        if not gene1 or not gene2:
            continue
        events.append({
            "gene5": gene1.strip(),
            "gene3": gene2.strip(),
            "left": row.get("breakpoint1", ""),
            "right": row.get("breakpoint2", ""),
            "split": to_int(row.get("split_reads1")) + to_int(row.get("split_reads2")),
            "span": to_int(row.get("discordant_mates")),
            "confidence": row.get("confidence", ""),
            "type": row.get("type", ""),
            "reads": split_reads(row.get("read_identifiers", "")),
        })
    return events


def parse_starfusion(path):
    events = []
    for row in read_rows(path):
        name = row.get("#FusionName") or row.get("FusionName")
        if not name or "--" not in name:
            continue
        gene5, gene3 = name.split("--")[:2]
        events.append({
            "gene5": gene5.strip(),
            "gene3": gene3.strip(),
            "left": row.get("LeftBreakpoint", ""),
            "right": row.get("RightBreakpoint", ""),
            "split": to_int(row.get("JunctionReadCount")),
            "span": to_int(row.get("SpanningFragCount")),
            "ffpm": row.get("FFPM", ""),
            "reads": split_reads(row.get("JunctionReads", ""))
            | split_reads(row.get("SpanningFrags", "")),
        })
    return events


def breakpoint_distance(a, b):
    """Largest distance between matching breakpoints, or None if the events can't agree."""
    if (a["gene5"], a["gene3"]) != (b["gene5"], b["gene3"]):
        return None
    a_left, a_right = parse_breakpoint(a["left"]), parse_breakpoint(a["right"])
    b_left, b_right = parse_breakpoint(b["left"]), parse_breakpoint(b["right"])
    if not all((a_left, a_right, b_left, b_right)):
        return None
    if a_left[0] != b_left[0] or a_right[0] != b_right[0]:
        return None
    return max(abs(a_left[1] - b_left[1]), abs(a_right[1] - b_right[1]))


def pair_events(arriba, starfusion, window):
    """Pair each Arriba event with at most one STAR-Fusion event, closest first."""
    candidates = []
    for a_index, a in enumerate(arriba):
        for s_index, s in enumerate(starfusion):
            distance = breakpoint_distance(a, s)
            if distance is not None and distance <= window:
                candidates.append((distance, a_index, s_index))

    a_used, s_used, pairs = set(), set(), []
    for _, a_index, s_index in sorted(candidates):
        if a_index not in a_used and s_index not in s_used:
            a_used.add(a_index)
            s_used.add(s_index)
            pairs.append((arriba[a_index], starfusion[s_index]))

    pairs += [(a, None) for i, a in enumerate(arriba) if i not in a_used]
    pairs += [(None, s) for i, s in enumerate(starfusion) if i not in s_used]
    return pairs


def caller_total(event):
    return event["split"] + event["span"]


def support(a, s):
    """Return (total_support, support_source) without counting shared reads twice."""
    callers = [event for event in (a, s) if event]
    if all(event["reads"] or caller_total(event) == 0 for event in callers):
        return len(set().union(*(event["reads"] for event in callers))), "unique_read_names"
    return max(caller_total(event) for event in callers), "caller_max"


def build_record(sample, a, s, max_readnames):
    event = a or s
    callers = [name for name, found in (("arriba", a), ("starfusion", s)) if found]
    total, source = support(a, s)
    reads = sorted(set().union(*(e["reads"] for e in (a, s) if e)))
    return {
        "sample": sample,
        "fusion": f"{event['gene5']}--{event['gene3']}",
        "left_breakpoint": event["left"],
        "right_breakpoint": event["right"],
        "callers": ",".join(callers),
        "n_callers": len(callers),
        "total_support": total,
        "support_source": source,
        "arriba_split": a["split"] if a else 0,
        "arriba_span": a["span"] if a else 0,
        "sf_split": s["split"] if s else 0,
        "sf_span": s["span"] if s else 0,
        "arriba_confidence": a["confidence"] if a else "",
        "type": a["type"] if a else "",
        "starfusion_ffpm": s["ffpm"] if s else "",
        "n_read_names": len(reads),
        "read_names": ",".join(reads[:max_readnames]) if reads else ".",
    }


def write_table(path, rows):
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=FIELDS, delimiter="\t",
                                quoting=csv.QUOTE_NONE, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--arriba", required=True)
    ap.add_argument("--starfusion", default=None)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--reported", required=True)
    ap.add_argument("--min-callers", type=int, default=2)
    ap.add_argument("--max-readnames", type=int, default=50)
    ap.add_argument("--breakpoint-window", type=int, default=10)
    args = ap.parse_args()

    # Read-name columns grow with support and can exceed the csv module's default limit.
    csv.field_size_limit(sys.maxsize)

    pairs = pair_events(parse_arriba(args.arriba), parse_starfusion(args.starfusion),
                        args.breakpoint_window)
    records = [build_record(args.sample, a, s, args.max_readnames) for a, s in pairs]
    records.sort(key=lambda r: (-r["n_callers"], -r["total_support"], r["fusion"],
                                r["left_breakpoint"], r["right_breakpoint"]))

    def promoted(record):
        return (record["n_callers"] >= args.min_callers
                or record["arriba_confidence"].lower() == "high")

    reported = [r for r in records if promoted(r)]
    write_table(args.merged, records)
    write_table(args.reported, reported)

    sys.stderr.write(
        f"[merge_fusions] {args.sample}: {len(records)} events, {len(reported)} reported\n"
    )


if __name__ == "__main__":
    main()
