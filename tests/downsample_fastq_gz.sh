#!/usr/bin/env bash
# Deterministically downsample a gzipped FASTQ.
#
# Two modes:
#  - seqtk (default): random sample with fixed seed, preserves FASTQ record structure
#  - head: take the first N reads (fastest; deterministic but may be biased if input is ordered)
#
# Usage:
#   ./downsample_fastq_gz.sh --in in.fastq.gz --out out.fastq.gz --reads 10000000 [--seed 1] [--mode seqtk|head]
#
# Notes:
# - Assumes standard 4-line FASTQ records.
# - For SLAM-seq dev fixtures, this is usually good enough. For perfect statistical sampling,
#   prefer seqtk mode.
set -euo pipefail

IN=""
OUT=""
READS=""
SEED="1"
MODE="seqtk"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --in) IN="${2:?}"; shift 2 ;;
    --out) OUT="${2:?}"; shift 2 ;;
    --reads) READS="${2:?}"; shift 2 ;;
    --seed) SEED="${2:?}"; shift 2 ;;
    --mode) MODE="${2:?}"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 --in in.fastq.gz --out out.fastq.gz --reads N [--seed 1] [--mode seqtk|head]"
      exit 0
      ;;
    *)
      echo "Unknown arg: $1" >&2
      exit 2
      ;;
  esac
done

if [[ -z "$IN" || -z "$OUT" || -z "$READS" ]]; then
  echo "Missing required args. Use --help." >&2
  exit 2
fi
if [[ ! -f "$IN" ]]; then
  echo "Input FASTQ not found: $IN" >&2
  exit 2
fi

mkdir -p "$(dirname "$OUT")"

case "$MODE" in
  seqtk)
    if ! command -v seqtk >/dev/null 2>&1; then
      echo "seqtk not found; install seqtk or use --mode head" >&2
      exit 2
    fi
    # seqtk sample outputs uncompressed FASTQ; pipe to pigz for fast compression
    if command -v pigz >/dev/null 2>&1; then
      seqtk sample -s"$SEED" "$IN" "$READS" | pigz -p 8 > "$OUT"
    else
      seqtk sample -s"$SEED" "$IN" "$READS" | gzip -c > "$OUT"
    fi
    ;;
  head)
    # Take first N reads = first 4*N lines
    LINES=$((READS * 4))
    if command -v pigz >/dev/null 2>&1; then
      zcat "$IN" | head -n "$LINES" | pigz -p 8 > "$OUT"
    else
      zcat "$IN" | head -n "$LINES" | gzip -c > "$OUT"
    fi
    ;;
  *)
    echo "Unknown --mode: $MODE (expected seqtk|head)" >&2
    exit 2
    ;;
esac

echo "Wrote: $OUT"
