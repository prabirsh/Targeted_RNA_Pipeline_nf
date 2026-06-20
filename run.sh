#!/usr/bin/env bash
# Interactive launcher for the RNA panel pipeline.
# It ASKS which panel to use, then runs Nextflow.
# Usage:  ./run.sh [samplesheet.csv] [outdir]
set -euo pipefail

SAMPLESHEET="${1:-assets/samplesheet.example.csv}"
OUTDIR="${2:-results}"

# ---- reference bundle (edit to your environment) -----------------------------
REFS=/refs/GRCh37
STAR_INDEX=$REFS/star_index
GTF=$REFS/annotation.gtf
FASTA=$REFS/genome.fa
REFFLAT=$REFS/refFlat.txt
RRNA=$REFS/rRNA.interval_list
CTAT=$REFS/ctat_genome_lib
ARRIBA_BL=$REFS/arriba/blacklist.tsv.gz
ARRIBA_KN=$REFS/arriba/known_fusions.tsv.gz
ARRIBA_PD=$REFS/arriba/protein_domains.gff3

# ---- ASK which panel ---------------------------------------------------------
echo "Select the targeted panel for this run:"
echo "  1) ST8  - SureSelect XT V8 exome panel"
echo "  2) CT   - Indiegene 4baseCare RNA Fusion panel (GRCh37)"
read -rp "Enter 1 or 2 [1]: " choice
case "${choice:-1}" in
  1) PANEL=ST8 ;;
  2) PANEL=CT  ;;
  *) echo "Invalid choice"; exit 1 ;;
esac
echo "Using panel: $PANEL"

nextflow run main.nf \
    -profile docker \
    --input "$SAMPLESHEET" \
    --outdir "$OUTDIR" \
    --panel "$PANEL" \
    --star_index "$STAR_INDEX" \
    --gtf "$GTF" \
    --fasta "$FASTA" \
    --refflat "$REFFLAT" \
    --rrna_intervals "$RRNA" \
    --ctat_lib "$CTAT" \
    --arriba_blacklist "$ARRIBA_BL" \
    --arriba_known "$ARRIBA_KN" \
    --arriba_protein "$ARRIBA_PD" \
    -resume
