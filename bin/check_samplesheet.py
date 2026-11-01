#!/usr/bin/env python3
"""
Validate and normalise the RNA panel samplesheet.

Expected columns: sample,fastq_1,fastq_2,strandedness
  - sample       : unique sample identifier (letters, digits, '.', '_', '-')
  - fastq_1      : path to R1 FASTQ (.fastq.gz / .fq.gz)
  - fastq_2      : path to R2 FASTQ (blank for single-end)
  - strandedness : reverse | forward | unstranded  (default: reverse)

Emits a cleaned CSV on success; exits non-zero with a clear message on failure.
"""
import csv
import re
import sys

VALID_STRAND = {"reverse", "forward", "unstranded"}
REQUIRED = ["sample", "fastq_1"]
SAFE_SAMPLE_ID = re.compile(r"[A-Za-z0-9_.-]+")


def die(msg):
    sys.stderr.write(f"[samplesheet ERROR] {msg}\n")
    sys.exit(1)


def check_sample_id(sample, n, seen):
    if not sample:
        die(f"line {n}: empty sample id")
    # Sample ids become file names, so they must not contain path or shell characters.
    if not SAFE_SAMPLE_ID.fullmatch(sample) or sample.startswith("."):
        die(f"line {n}: sample id '{sample}' may only use letters, digits, '.', '_' and '-', "
            "and must not start with '.'")
    if sample in seen:
        die(f"line {n}: duplicate sample id '{sample}'")


def check_reads(row, n):
    if not row.get("fastq_1"):
        die(f"line {n}: fastq_1 is required")
    for key in ("fastq_1", "fastq_2"):
        value = row.get(key, "")
        if value and not value.endswith((".fastq.gz", ".fq.gz")):
            die(f"line {n}: {key} '{value}' must end in .fastq.gz or .fq.gz")
    if row.get("fastq_1") == row.get("fastq_2"):
        die(f"line {n}: fastq_1 and fastq_2 are the same file")


def main(infile, outfile):
    with open(infile, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            die("samplesheet is empty")
        header = [h.strip() for h in reader.fieldnames]
        for col in REQUIRED:
            if col not in header:
                die(f"missing required column: '{col}'. Found: {header}")

        seen = set()
        rows = []
        for n, raw in enumerate(reader, start=2):
            if None in raw:
                die(f"line {n}: more fields than the header's {len(header)} columns")
            row = {k.strip(): (v.strip() if v else "") for k, v in raw.items()}
            sample = row.get("sample", "")
            check_sample_id(sample, n, seen)
            seen.add(sample)
            check_reads(row, n)

            strand = row.get("strandedness", "") or "reverse"
            if strand not in VALID_STRAND:
                die(f"line {n}: strandedness '{strand}' not in {sorted(VALID_STRAND)}")
            row["strandedness"] = strand
            row.setdefault("fastq_2", "")
            rows.append(row)

        if not rows:
            die("no data rows found in samplesheet")

    with open(outfile, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=["sample", "fastq_1", "fastq_2", "strandedness"],
                                lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in writer.fieldnames})

    sys.stderr.write(f"[samplesheet OK] {len(rows)} sample(s) validated\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        die("usage: check_samplesheet.py <input.csv> <output.csv>")
    main(sys.argv[1], sys.argv[2])
