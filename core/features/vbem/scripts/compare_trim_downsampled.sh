#!/usr/bin/env bash
set -euo pipefail

# Compare STAR-Flex external trimming against pre-trimmed FASTQs on a downsample.
# Downsamples by line count (default 100000 lines = 25000 reads).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TRIMVALIDATE_BIN="${TRIMVALIDATE_BIN:-$SCRIPT_DIR/../tools/trimvalidate/trimvalidate}"

RAW_R1="${RAW_R1:-}"
RAW_R2="${RAW_R2:-}"
TRIMMED_R1="${TRIMMED_R1:-}"
TRIMMED_R2="${TRIMMED_R2:-}"
OUT_DIR="${OUT_DIR:-/storage/production/trim_compare_downsampled}"
LINE_COUNT="${LINE_COUNT:-100000}"
QUALITY_CUTOFF="${QUALITY_CUTOFF:-20}"
MIN_LENGTH="${MIN_LENGTH:-20}"
ADAPTER_R1="${ADAPTER_R1:-}"
ADAPTER_R2="${ADAPTER_R2:-}"
COMPAT_MODE="${COMPAT_MODE:-}"

usage() {
    cat <<'USAGE'
Usage:
  compare_trim_downsampled.sh \
    --raw-r1 RAW_R1 --raw-r2 RAW_R2 \
    --trimmed-r1 TRIM_R1 --trimmed-r2 TRIM_R2 \
    --out-dir OUT_DIR [--lines N]

Options:
  --raw-r1 PATH        Raw R1 FASTQ(.gz)
  --raw-r2 PATH        Raw R2 FASTQ(.gz)
  --trimmed-r1 PATH    Pre-trimmed R1 FASTQ(.gz) (Trim Galore output)
  --trimmed-r2 PATH    Pre-trimmed R2 FASTQ(.gz)
  --out-dir DIR        Output directory
  --lines N            Number of lines to keep (default: 100000)
  --quality Q          Quality cutoff (default: 20)
  --min-length L       Minimum length after trimming (default: 20)
  --adapter-r1 SEQ     Adapter sequence for R1 (default: trimvalidate default)
  --adapter-r2 SEQ     Adapter sequence for R2 (default: trimvalidate default)
  --compat MODE        Compatibility mode: Off (default) or Cutadapt3

Env vars:
  RAW_R1, RAW_R2, TRIMMED_R1, TRIMMED_R2, OUT_DIR, LINE_COUNT, TRIMVALIDATE_BIN,
  QUALITY_CUTOFF, MIN_LENGTH, ADAPTER_R1, ADAPTER_R2, COMPAT_MODE
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --raw-r1) RAW_R1="$2"; shift 2 ;;
        --raw-r2) RAW_R2="$2"; shift 2 ;;
        --trimmed-r1) TRIMMED_R1="$2"; shift 2 ;;
        --trimmed-r2) TRIMMED_R2="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --lines) LINE_COUNT="$2"; shift 2 ;;
        --quality) QUALITY_CUTOFF="$2"; shift 2 ;;
        --min-length) MIN_LENGTH="$2"; shift 2 ;;
        --adapter-r1) ADAPTER_R1="$2"; shift 2 ;;
        --adapter-r2) ADAPTER_R2="$2"; shift 2 ;;
        --compat) COMPAT_MODE="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ -z "$RAW_R1" || -z "$RAW_R2" || -z "$TRIMMED_R1" || -z "$TRIMMED_R2" ]]; then
    echo "ERROR: raw and trimmed FASTQs are required." >&2
    usage
    exit 1
fi
if [[ ! -x "$TRIMVALIDATE_BIN" ]]; then
    echo "ERROR: trimvalidate not found: $TRIMVALIDATE_BIN" >&2
    exit 1
fi
if [[ ! "$LINE_COUNT" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --lines must be a non-negative integer." >&2
    exit 1
fi
if [[ ! "$QUALITY_CUTOFF" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --quality must be a non-negative integer." >&2
    exit 1
fi
if [[ ! "$MIN_LENGTH" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --min-length must be a non-negative integer." >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

downsample_fastq() {
    local src="$1"
    local dst="$2"
    if [[ "$src" == *.gz ]]; then
        if (( LINE_COUNT > 0 )); then
            ( gzip -dc -- "$src" || true ) | head -n "$LINE_COUNT" | gzip -c > "$dst"
        else
            : | gzip -c > "$dst"
        fi
    else
        if (( LINE_COUNT > 0 )); then
            head -n "$LINE_COUNT" -- "$src" > "$dst"
        else
            : > "$dst"
        fi
    fi
}

RAW_DS_R1="$OUT_DIR/raw_R1.downsampled.fq.gz"
RAW_DS_R2="$OUT_DIR/raw_R2.downsampled.fq.gz"
TRIM_DS_R1="$OUT_DIR/trimmed_R1.downsampled.fq.gz"
TRIM_DS_R2="$OUT_DIR/trimmed_R2.downsampled.fq.gz"

echo "Downsampling to $LINE_COUNT lines..."
downsample_fastq "$RAW_R1" "$RAW_DS_R1"
downsample_fastq "$RAW_R2" "$RAW_DS_R2"
downsample_fastq "$TRIMMED_R1" "$TRIM_DS_R1"
downsample_fastq "$TRIMMED_R2" "$TRIM_DS_R2"

RAW_DS_R1_PLAIN="$OUT_DIR/raw_R1.downsampled.fastq"
RAW_DS_R2_PLAIN="$OUT_DIR/raw_R2.downsampled.fastq"
TRIM_DS_R1_PLAIN="$OUT_DIR/trimmed_R1.downsampled.fastq"
TRIM_DS_R2_PLAIN="$OUT_DIR/trimmed_R2.downsampled.fastq"

gzip -dc "$RAW_DS_R1" > "$RAW_DS_R1_PLAIN"
gzip -dc "$RAW_DS_R2" > "$RAW_DS_R2_PLAIN"
gzip -dc "$TRIM_DS_R1" > "$TRIM_DS_R1_PLAIN"
gzip -dc "$TRIM_DS_R2" > "$TRIM_DS_R2_PLAIN"

echo "Running trimvalidate on raw downsample..."
trim_args=(--quality "$QUALITY_CUTOFF" --length "$MIN_LENGTH")
if [[ -n "$ADAPTER_R1" ]]; then
    trim_args+=(--adapter-r1 "$ADAPTER_R1")
fi
if [[ -n "$ADAPTER_R2" ]]; then
    trim_args+=(--adapter-r2 "$ADAPTER_R2")
fi
if [[ -n "$COMPAT_MODE" ]]; then
    trim_args+=(--compat "$COMPAT_MODE")
fi
"$TRIMVALIDATE_BIN" -1 "$RAW_DS_R1_PLAIN" -2 "$RAW_DS_R2_PLAIN" \
    -o1 "$OUT_DIR/raw_trimmed_R1.fastq" \
    -o2 "$OUT_DIR/raw_trimmed_R2.fastq" \
    "${trim_args[@]}" > "$OUT_DIR/trimvalidate_raw.log" 2>&1

echo "Running trimvalidate on pre-trimmed downsample..."
"$TRIMVALIDATE_BIN" -1 "$TRIM_DS_R1_PLAIN" -2 "$TRIM_DS_R2_PLAIN" \
    -o1 "$OUT_DIR/pretrim_trimmed_R1.fastq" \
    -o2 "$OUT_DIR/pretrim_trimmed_R2.fastq" \
    "${trim_args[@]}" > "$OUT_DIR/trimvalidate_pretrim.log" 2>&1

echo "Comparing outputs..."
diff -u "$OUT_DIR/pretrim_trimmed_R1.fastq" "$OUT_DIR/raw_trimmed_R1.fastq" > "$OUT_DIR/trim_R1.diff" || true
diff -u "$OUT_DIR/pretrim_trimmed_R2.fastq" "$OUT_DIR/raw_trimmed_R2.fastq" > "$OUT_DIR/trim_R2.diff" || true

echo "R1 diff: $OUT_DIR/trim_R1.diff"
echo "R2 diff: $OUT_DIR/trim_R2.diff"
echo "Done. Outputs in: $OUT_DIR"
