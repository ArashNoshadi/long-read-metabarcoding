#!/usr/bin/env bash
# cdhit_dedup.sh
# Usage: bash cdhit_dedup.sh input.fasta out_prefix identity_percentage
# example: bash cdhit_dedup.sh data/cox1.fasta results/cox1 100

IN="$1"
OUT_PREFIX="$2"
IDENTITY="$3"   # e.g., 100 for exact duplicates
if [ -z "$IN" ] || [ -z "$OUT_PREFIX" ] || [ -z "$IDENTITY" ]; then
  echo "Usage: bash cdhit_dedup.sh input.fasta out_prefix identity_percentage"
  exit 1
fi

OUT_DIR=$(dirname "$OUT_PREFIX")
mkdir -p "$OUT_DIR"

cd-hit -i "$IN" -o "${OUT_PREFIX}_cdhit.fasta" -c $(echo "scale=2; $IDENTITY/100" | bc) -n 5 -d 0
echo "CD-HIT output: ${OUT_PREFIX}_cdhit.fasta"
