#!/usr/bin/env bash
# Run union parity (GEDI goal ∪ STAR mask) and add a random/background panel
# to estimate "are we calling everything?".
#
# Random panel is sampled from the BAM itself (read starts) using samtools -s,
# then reservoir-sampled to avoid chr1 bias from coordinate order.
#
# Usage:
#   ./test_union_plus_random_parity.sh <bam> <ref_fa> <star_mask.bed(.gz)> <gedi_goal.bed> <out_dir> \
#       [--randomN 30000] [--seed 1] [--samtoolsFrac 0.02]
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

UNION_SCRIPT="$SCRIPT_DIR/test_union_parity.sh"
PILEUP="$SCRIPT_DIR/pileup_snp"
SAMPLE_PY="$ROOT_DIR/tests/sample_random_bg_from_bam_stream.py"

if [[ $# -lt 5 ]]; then
  echo "Usage: $0 <bam> <ref_fa> <star_mask.bed(.gz)> <gedi_goal.bed> <out_dir> [--randomN 30000] [--seed 1] [--samtoolsFrac 0.02]" >&2
  exit 2
fi

BAM="$1"
REF_FA="$2"
STAR_BED="$3"
GEDI_BED="$4"
OUT_DIR="$5"
shift 5

RANDOM_N=30000
SEED=1
SAMTOOLS_FRAC=0.02

while [[ $# -gt 0 ]]; do
  case "$1" in
    --randomN) RANDOM_N="${2:?}"; shift 2 ;;
    --seed) SEED="${2:?}"; shift 2 ;;
    --samtoolsFrac) SAMTOOLS_FRAC="${2:?}"; shift 2 ;;
    *)
      echo "Unknown arg: $1" >&2
      exit 2
      ;;
  esac
done

mkdir -p "$OUT_DIR"

if [[ ! -x "$UNION_SCRIPT" ]]; then
  echo "ERROR: missing $UNION_SCRIPT" >&2
  exit 2
fi
if [[ ! -x "$PILEUP" ]]; then
  echo "ERROR: missing $PILEUP (run: cd tools/pileup_snp && make)" >&2
  exit 2
fi
if [[ ! -f "$SAMPLE_PY" ]]; then
  echo "ERROR: missing $SAMPLE_PY" >&2
  exit 2
fi

echo "=== Union + Random Parity ==="
echo "BAM: $BAM"
echo "REF: $REF_FA"
echo "STAR_BED: $STAR_BED"
echo "GEDI_BED: $GEDI_BED"
echo "OUT_DIR: $OUT_DIR"
echo "randomN=$RANDOM_N seed=$SEED samtoolsFrac=$SAMTOOLS_FRAC"
echo

echo "[1/3] Union parity..."
"$UNION_SCRIPT" "$BAM" "$REF_FA" "$STAR_BED" "$GEDI_BED" "$OUT_DIR" | tee "$OUT_DIR/union_run.log"

UNION="$OUT_DIR/union_3col.bed"
if [[ ! -f "$UNION" ]]; then
  echo "ERROR: union BED not found after union parity: $UNION" >&2
  exit 2
fi

echo
echo "[2/3] Random background panel (exclude UNION)..."
RAND_BED="$OUT_DIR/random_bg_${RANDOM_N}.seed${SEED}.bed"
SARG="${SEED}.${SAMTOOLS_FRAC#0.}"
echo "samtools view -F 4 -s $SARG"

samtools view -F 4 -s "$SARG" "$BAM" \
  | python3 "$SAMPLE_PY" --exclude "$UNION" --n "$RANDOM_N" --seed "$SEED" \
  > "$RAND_BED" 2> "$OUT_DIR/random_bg.stats"

echo "random_bed_lines=$(wc -l < "$RAND_BED")"
echo "random_stats=$(cat "$OUT_DIR/random_bg.stats")"

echo
echo "[3/3] Call random panel with pileup_snp (GEDI-compat settings)..."
"$PILEUP" \
  --bam "$BAM" \
  --bed "$RAND_BED" \
  --ref "$REF_FA" \
  --output "$OUT_DIR/pileup_gedi_compat_random.bed.gz" \
  --debug-tsv "$OUT_DIR/pileup_gedi_compat_random_debug.tsv" \
  --kMode any \
  --includeSecondary 1 \
  --minMapQ 0 \
  --minBaseQ 0 \
  --pErr 0.001 \
  --pval 0.001 \
  --minTcRatio 0.3 \
  --minCov 6 \
  --minAlt 1 \
  > "$OUT_DIR/pileup_gedi_compat_random.log" 2>&1

python3 - "$OUT_DIR/pileup_gedi_compat_random_debug.tsv" > "$OUT_DIR/random_parity_summary.txt" <<'PY'
import sys
import csv
path = sys.argv[1]
tot=cov6=called=called_cov6=0
with open(path) as f:
    r=csv.DictReader(f, delimiter="\t")
    for row in r:
        tot += 1
        n=int(row["n"])
        c=(row["called_any"]=="1")
        if n>=6: cov6 += 1
        if c: called += 1
        if c and n>=6: called_cov6 += 1
print("random_total", tot)
print("random_cov>=6", cov6)
print("random_called_any", called)
print("random_called_any_among_cov>=6", called_cov6)
print("call_rate_among_cov>=6", (called_cov6/cov6) if cov6 else "NA")
PY

cat "$OUT_DIR/random_parity_summary.txt"

echo
echo "DONE. See:"
echo "  - $OUT_DIR/union_parity_summary.txt"
echo "  - $OUT_DIR/random_parity_summary.txt"

