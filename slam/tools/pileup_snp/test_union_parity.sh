#!/usr/bin/env bash
# Union-parity harness:
#   - Build UNION loci = GEDI_goal_bed (3col) ∪ STAR_mask_bed(.gz) (first 3 cols, skip header)
#   - Run pileup_snp on UNION loci in two modes (star-compat and gedi-compat)
#   - Report TP/FP/FN vs GEDI and vs STAR, all restricted to UNION universe
#
# Usage:
#   ./test_union_parity.sh <bam> <ref_fa> <star_mask.bed.gz|bed> <gedi_goal.bed> <out_dir>
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOL="$SCRIPT_DIR/pileup_snp"

if [[ $# -lt 5 ]]; then
  echo "Usage: $0 <bam> <ref_fa> <star_mask.bed(.gz)> <gedi_goal.bed> <out_dir>" >&2
  exit 2
fi

BAM="$1"
REF_FA="$2"
STAR_BED="$3"
GEDI_BED="$4"
OUT_DIR="$5"

mkdir -p "$OUT_DIR"

if [[ ! -x "$TOOL" ]]; then
  echo "ERROR: pileup_snp binary not found at $TOOL (run: cd tools/pileup_snp && make)" >&2
  exit 2
fi
if [[ ! -f "$BAM" ]]; then
  echo "ERROR: BAM not found: $BAM" >&2
  exit 2
fi
if [[ ! -f "$REF_FA" ]]; then
  echo "ERROR: reference FASTA not found: $REF_FA" >&2
  exit 2
fi
if [[ ! -f "$GEDI_BED" ]]; then
  echo "ERROR: GEDI goal BED not found: $GEDI_BED" >&2
  exit 2
fi
if [[ ! -f "$STAR_BED" ]]; then
  echo "ERROR: STAR mask BED not found: $STAR_BED" >&2
  exit 2
fi

echo "=== Union Parity: pileup_snp vs GEDI goal + STAR mask (within union) ==="
echo "BAM: $BAM"
echo "REF: $REF_FA"
echo "STAR_BED: $STAR_BED"
echo "GEDI_BED: $GEDI_BED"
echo "OUT_DIR: $OUT_DIR"
echo

STAR_3="$OUT_DIR/star_3col.bed"
GEDI_3="$OUT_DIR/gedi_3col.bed"
UNION="$OUT_DIR/union_3col.bed"

echo "[1/5] Normalize inputs to 3-column BED + sort/uniq..."
if [[ "$STAR_BED" == *.gz ]]; then
  zcat "$STAR_BED" | grep -v '^#' | cut -f1-3 | sort -k1,1 -k2,2n | uniq > "$STAR_3"
else
  grep -v '^#' "$STAR_BED" | cut -f1-3 | sort -k1,1 -k2,2n | uniq > "$STAR_3"
fi

cat "$GEDI_BED" | grep -v '^#' | cut -f1-3 | sort -k1,1 -k2,2n | uniq > "$GEDI_3"

cat "$STAR_3" "$GEDI_3" | sort -k1,1 -k2,2n | uniq > "$UNION"

N_STAR=$(wc -l < "$STAR_3" || echo 0)
N_GEDI=$(wc -l < "$GEDI_3" || echo 0)
N_UNION=$(wc -l < "$UNION" || echo 0)
echo "STAR loci (uniq): $N_STAR"
echo "GEDI loci (uniq): $N_GEDI"
echo "UNION loci:       $N_UNION"
echo

echo "[2/5] Run pileup_snp on UNION in STAR-compat mode..."
"$TOOL" \
  --bam "$BAM" \
  --bed "$UNION" \
  --ref "$REF_FA" \
  --output "$OUT_DIR/pileup_star_compat_union.bed.gz" \
  --debug-tsv "$OUT_DIR/pileup_star_compat_union_debug.tsv" \
  --kMode conv \
  --includeSecondary 0 \
  --minMapQ 0 \
  --minBaseQ 0 \
  --pErr 0.001 \
  --pval 0.001 \
  --minTcRatio 0.3 \
  --minCov 6 \
  --minAlt 1 \
  > "$OUT_DIR/pileup_star_compat_union.log" 2>&1

echo "[3/5] Run pileup_snp on UNION in GEDI-compat mode..."
"$TOOL" \
  --bam "$BAM" \
  --bed "$UNION" \
  --ref "$REF_FA" \
  --output "$OUT_DIR/pileup_gedi_compat_union.bed.gz" \
  --debug-tsv "$OUT_DIR/pileup_gedi_compat_union_debug.tsv" \
  --kMode any \
  --includeSecondary 1 \
  --minMapQ 0 \
  --minBaseQ 0 \
  --pErr 0.001 \
  --pval 0.001 \
  --minTcRatio 0.3 \
  --minCov 6 \
  --minAlt 1 \
  > "$OUT_DIR/pileup_gedi_compat_union.log" 2>&1

P_STAR="$OUT_DIR/pileup_star_compat_union.3col.bed"
P_GEDI="$OUT_DIR/pileup_gedi_compat_union.3col.bed"
zcat "$OUT_DIR/pileup_star_compat_union.bed.gz" | grep -v '^#' | cut -f1-3 | sort -k1,1 -k2,2n | uniq > "$P_STAR"
zcat "$OUT_DIR/pileup_gedi_compat_union.bed.gz" | grep -v '^#' | cut -f1-3 | sort -k1,1 -k2,2n | uniq > "$P_GEDI"

N_PSTAR=$(wc -l < "$P_STAR" || echo 0)
N_PGEDI=$(wc -l < "$P_GEDI" || echo 0)
echo "pileup_snp STAR-compat calls (within union): $N_PSTAR"
echo "pileup_snp GEDI-compat calls (within union): $N_PGEDI"
echo

echo "[4/5] Compute confusion tables (within UNION)..."

tp() { bedtools intersect -a "$1" -b "$2" -u | wc -l; }
diffcnt() { bedtools intersect -a "$1" -b "$2" -v | wc -l; }

TP_GEDI_STAR=$(tp "$P_STAR" "$GEDI_3")
FP_GEDI_STAR=$(diffcnt "$P_STAR" "$GEDI_3")
FN_GEDI_STAR=$(diffcnt "$GEDI_3" "$P_STAR")

TP_GEDI_GEDI=$(tp "$P_GEDI" "$GEDI_3")
FP_GEDI_GEDI=$(diffcnt "$P_GEDI" "$GEDI_3")
FN_GEDI_GEDI=$(diffcnt "$GEDI_3" "$P_GEDI")

TP_STAR_STAR=$(tp "$P_STAR" "$STAR_3")
FP_STAR_STAR=$(diffcnt "$P_STAR" "$STAR_3")
FN_STAR_STAR=$(diffcnt "$STAR_3" "$P_STAR")

TP_STAR_GEDI=$(tp "$P_GEDI" "$STAR_3")
FP_STAR_GEDI=$(diffcnt "$P_GEDI" "$STAR_3")
FN_STAR_GEDI=$(diffcnt "$STAR_3" "$P_GEDI")

SUMMARY="$OUT_DIR/union_parity_summary.txt"
{
  echo "=== UNION Parity Summary ==="
  echo "UNION=$N_UNION  GEDI=$N_GEDI  STAR=$N_STAR"
  echo
  echo "pileup_snp STAR-compat (kMode=conv, secondary=0): calls=$N_PSTAR"
  echo "  vs GEDI: TP=$TP_GEDI_STAR  FP=$FP_GEDI_STAR  FN=$FN_GEDI_STAR"
  echo "  vs STAR: TP=$TP_STAR_STAR  FP=$FP_STAR_STAR  FN=$FN_STAR_STAR"
  echo
  echo "pileup_snp GEDI-compat (kMode=any, secondary=1): calls=$N_PGEDI"
  echo "  vs GEDI: TP=$TP_GEDI_GEDI  FP=$FP_GEDI_GEDI  FN=$FN_GEDI_GEDI"
  echo "  vs STAR: TP=$TP_STAR_GEDI  FP=$FP_STAR_GEDI  FN=$FN_STAR_GEDI"
} | tee "$SUMMARY"

echo
echo "[5/5] Done. Outputs in: $OUT_DIR"

