#!/usr/bin/env bash
# Sweep (minAlt, p_err) for GEDI-compat-like calling on:
#   - UNION loci (GEDI goal ∪ STAR mask) to compute TP/FP/FN vs GEDI
#   - RANDOM background loci to estimate over-calling rate
#
# This script is intentionally simple and file-based for reproducibility.
#
# Usage:
#   ./run_knob_sweep_union_random.sh \
#     --bam <bam> --ref <ref.fa> \
#     --union <union_3col.bed> --gedi <gedi_3col.bed> --random <random_3col.bed> \
#     --out <out.tsv> \
#     [--kMode any|conv] [--includeSecondary 0|1] [--minTcRatio 0.3] [--minCov 6] [--pval 0.001] \
#     [--minAltList 1,3,4,5,6] [--pErrList 0.001,0.05,0.1,0.15]
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PILEUP="$SCRIPT_DIR/pileup_snp"

BAM=""
REF=""
UNION=""
GEDI=""
RANDOM_BED=""
OUT=""

KMODE="any"
INCLUDE_SECONDARY=1
MIN_TCRATIO="0.3"
MIN_COV="6"
PVAL="0.001"

MINALT_LIST="1,3,4,5,6"
PERR_LIST="0.001,0.05,0.1,0.15"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="${2:?}"; shift 2 ;;
    --ref) REF="${2:?}"; shift 2 ;;
    --union) UNION="${2:?}"; shift 2 ;;
    --gedi) GEDI="${2:?}"; shift 2 ;;
    --random) RANDOM_BED="${2:?}"; shift 2 ;;
    --out) OUT="${2:?}"; shift 2 ;;

    --kMode) KMODE="${2:?}"; shift 2 ;;
    --includeSecondary) INCLUDE_SECONDARY="${2:?}"; shift 2 ;;
    --minTcRatio) MIN_TCRATIO="${2:?}"; shift 2 ;;
    --minCov) MIN_COV="${2:?}"; shift 2 ;;
    --pval) PVAL="${2:?}"; shift 2 ;;

    --minAltList) MINALT_LIST="${2:?}"; shift 2 ;;
    --pErrList) PERR_LIST="${2:?}"; shift 2 ;;
    -h|--help)
      echo "See header for usage." >&2
      exit 0
      ;;
    *)
      echo "Unknown arg: $1" >&2
      exit 2
      ;;
  esac
done

if [[ -z "$BAM" || -z "$REF" || -z "$UNION" || -z "$GEDI" || -z "$RANDOM_BED" || -z "$OUT" ]]; then
  echo "ERROR: missing required args (use --help)." >&2
  exit 2
fi
if [[ ! -x "$PILEUP" ]]; then
  echo "ERROR: pileup_snp not found at $PILEUP (run: cd tools/pileup_snp && make)" >&2
  exit 2
fi
for f in "$BAM" "$REF" "$UNION" "$GEDI" "$RANDOM_BED"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: file not found: $f" >&2
    exit 2
  fi
done

TMP_DIR="$(mktemp -d)"
cleanup() { rm -rf "$TMP_DIR"; }
trap cleanup EXIT

# Precompute random cov>=minCov once (independent of p_err/minAlt because minCov is fixed)
RAND_COV_TSV="$TMP_DIR/random_depth.tsv"
samtools depth -a -b "$RANDOM_BED" "$BAM" > "$RAND_COV_TSV"
RAND_COV_GE_MINCOV=$(awk -v mc="$MIN_COV" '$3>=mc{c++} END{print c+0}' "$RAND_COV_TSV")
RAND_TOTAL=$(wc -l < "$RANDOM_BED" | tr -d ' ')

IFS=',' read -r -a MINALTS <<< "$MINALT_LIST"
IFS=',' read -r -a PERRS <<< "$PERR_LIST"

mkdir -p "$(dirname "$OUT")"
{
  echo -e "minAlt\tp_err\tkMode\tincludeSecondary\tminCov\tminTcRatio\tpval\tunion_calls\tTP_vs_GEDI\tFP_vs_GEDI\tFN_vs_GEDI\trandom_total\trandom_cov_ge_minCov\trandom_calls\trandom_call_rate_among_cov_ge_minCov"

  for minAlt in "${MINALTS[@]}"; do
    for pErr in "${PERRS[@]}"; do
      # UNION call
      U_BEDGZ="$TMP_DIR/union_${minAlt}_${pErr}.bed.gz"
      R_BEDGZ="$TMP_DIR/random_${minAlt}_${pErr}.bed.gz"
      U_3="$TMP_DIR/union_${minAlt}_${pErr}.3col.bed"
      R_3="$TMP_DIR/random_${minAlt}_${pErr}.3col.bed"

      "$PILEUP" \
        --bam "$BAM" --bed "$UNION" --ref "$REF" \
        --output "$U_BEDGZ" \
        --kMode "$KMODE" --includeSecondary "$INCLUDE_SECONDARY" \
        --minMapQ 0 --minBaseQ 0 \
        --pErr "$pErr" --pval "$PVAL" --minTcRatio "$MIN_TCRATIO" --minCov "$MIN_COV" --minAlt "$minAlt" \
        > /dev/null 2>&1

      zcat "$U_BEDGZ" | grep -v '^#' | cut -f1-3 | sort -k1,1 -k2,2n | uniq > "$U_3"
      UNION_CALLS=$(wc -l < "$U_3" | tr -d ' ')
      TP=$(bedtools intersect -a "$U_3" -b "$GEDI" -u | wc -l | tr -d ' ')
      FP=$(bedtools intersect -a "$U_3" -b "$GEDI" -v | wc -l | tr -d ' ')
      FN=$(bedtools intersect -a "$GEDI" -b "$U_3" -v | wc -l | tr -d ' ')

      # RANDOM call
      "$PILEUP" \
        --bam "$BAM" --bed "$RANDOM_BED" --ref "$REF" \
        --output "$R_BEDGZ" \
        --kMode "$KMODE" --includeSecondary "$INCLUDE_SECONDARY" \
        --minMapQ 0 --minBaseQ 0 \
        --pErr "$pErr" --pval "$PVAL" --minTcRatio "$MIN_TCRATIO" --minCov "$MIN_COV" --minAlt "$minAlt" \
        > /dev/null 2>&1

      zcat "$R_BEDGZ" | grep -v '^#' | cut -f1-3 | sort -k1,1 -k2,2n | uniq > "$R_3"
      RANDOM_CALLS=$(wc -l < "$R_3" | tr -d ' ')
      if [[ "$RAND_COV_GE_MINCOV" -gt 0 ]]; then
        RATE=$(python3 - <<PY
rc=$RANDOM_CALLS
den=$RAND_COV_GE_MINCOV
print(rc/den)
PY
)
      else
        RATE="NA"
      fi

      echo -e "${minAlt}\t${pErr}\t${KMODE}\t${INCLUDE_SECONDARY}\t${MIN_COV}\t${MIN_TCRATIO}\t${PVAL}\t${UNION_CALLS}\t${TP}\t${FP}\t${FN}\t${RAND_TOTAL}\t${RAND_COV_GE_MINCOV}\t${RANDOM_CALLS}\t${RATE}"
    done
  done
} > "$OUT"

echo "Wrote sweep table: $OUT"

