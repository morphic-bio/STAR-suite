#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

STAR_M1_DIR="${STAR_M1_DIR:-/tmp/ucsf_pf_trace_foxd3_nxtbc_20260219_004502/iPSC2_1_AALG2}"
STAR_M0_DIR="${STAR_M0_DIR:-/tmp/ucsf_pf_trace_foxd3_nxtbc_m0_20260219_010611/iPSC2_1_AALG2}"
CR_MEX_DIR="${CR_MEX_DIR:-/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260217_200813/outs/raw_feature_bc_matrix}"
CR_CALLS_CSV="${CR_CALLS_CSV:-/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260217_200813/outs/crispr_analysis/protospacer_calls_per_cell.csv}"
TRANSLATION_FILE="${TRANSLATION_FILE:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt.gz}"
OUT_DIR="${OUT_DIR:-/tmp/ucsf_star1h_vs_cr_$(date +%Y%m%d_%H%M%S)}"

# Default examples from discrepancy deep-dive
INSPECT_BARCODES_TRU_DEFAULT="AGGTAGGAGAAGATCT,ACTATCTTCAACCTTT,CGGGTCAAGGCATCTT,CAACAACCAACAACAA"
INSPECT_BARCODES_TRU="${INSPECT_BARCODES_TRU:-$INSPECT_BARCODES_TRU_DEFAULT}"

usage() {
  cat <<'EOF'
Usage: scripts/run_ucsf_star1h_cr_analysis.sh [options]

Options:
  --star-m1-dir DIR
  --star-m0-dir DIR
  --cr-mex-dir DIR
  --cr-calls-csv FILE
  --translation-file FILE
  --out-dir DIR
  --inspect-barcodes-tru CSV   Comma-separated TRU barcodes for deep-dive inspection
  -h, --help

Environment overrides:
  STAR_M1_DIR STAR_M0_DIR CR_MEX_DIR CR_CALLS_CSV TRANSLATION_FILE OUT_DIR
  INSPECT_BARCODES_TRU
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --star-m1-dir) STAR_M1_DIR="$2"; shift 2 ;;
    --star-m0-dir) STAR_M0_DIR="$2"; shift 2 ;;
    --cr-mex-dir) CR_MEX_DIR="$2"; shift 2 ;;
    --cr-calls-csv) CR_CALLS_CSV="$2"; shift 2 ;;
    --translation-file) TRANSLATION_FILE="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --inspect-barcodes-tru) INSPECT_BARCODES_TRU="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

mkdir -p "$OUT_DIR"

BUILD_SCRIPT="$REPO_ROOT/scripts/ucsf_parity/build_star_m1_delta_vs_cr.py"
CLASSIFY_SCRIPT="$REPO_ROOT/scripts/ucsf_parity/classify_star_m1_cr_misses.py"
INSPECT_SCRIPT="$REPO_ROOT/scripts/ucsf_parity/inspect_barcode_feature_totals.py"
EXACT_SCRIPT="$REPO_ROOT/scripts/ucsf_parity/build_star_exact_vs_cr.py"

for f in "$BUILD_SCRIPT" "$CLASSIFY_SCRIPT" "$INSPECT_SCRIPT" "$EXACT_SCRIPT"; do
  if [[ ! -x "$f" ]]; then
    echo "ERROR: missing executable script: $f" >&2
    exit 1
  fi
done

python3 "$BUILD_SCRIPT" \
  --star-m1-dir "$STAR_M1_DIR" \
  --star-m0-dir "$STAR_M0_DIR" \
  --cr-mex-dir "$CR_MEX_DIR" \
  --cr-calls-csv "$CR_CALLS_CSV" \
  --translation-file "$TRANSLATION_FILE" \
  --out-dir "$OUT_DIR"

python3 "$CLASSIFY_SCRIPT" \
  --delta-table "$OUT_DIR/STAR_M1_DELTA_VS_CR.tsv" \
  --star-m1-dir "$STAR_M1_DIR" \
  --cr-mex-dir "$CR_MEX_DIR" \
  --out-dir "$OUT_DIR"

python3 "$INSPECT_SCRIPT" \
  --barcodes-tru "$INSPECT_BARCODES_TRU" \
  --star-m1-dir "$STAR_M1_DIR" \
  --cr-mex-dir "$CR_MEX_DIR" \
  --translation-file "$TRANSLATION_FILE" \
  --out-dir "$OUT_DIR"

python3 "$EXACT_SCRIPT" \
  --star-dir "$STAR_M0_DIR" \
  --cr-mex-dir "$CR_MEX_DIR" \
  --cr-calls-csv "$CR_CALLS_CSV" \
  --translation-file "$TRANSLATION_FILE" \
  --out-dir "$OUT_DIR"

echo "[done] Output directory: $OUT_DIR"
echo "[done] Main tables:"
echo "  $OUT_DIR/STAR_M1_DELTA_VS_CR.tsv"
echo "  $OUT_DIR/STAR_M1_DELTA_CR_MISSING_WITH_PARTNER_COUNTS.tsv"
echo "  $OUT_DIR/BARCODE_INSPECTION_SUMMARY.txt"
echo "  $OUT_DIR/STAR_EXACT_VS_CR.tsv"
