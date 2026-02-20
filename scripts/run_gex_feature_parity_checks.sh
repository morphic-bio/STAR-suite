#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)

STAR_RUN_DEFAULT="/storage/ucsf-full/bench_20260218_dynamic_first/runs/star_dynamic_off_full_20260218_203459"
CR_RUN_DEFAULT="/storage/ucsf-full/bench_20260218_dynamic_first/cellranger_runs/cr_full_iPSC2_1_AALG2_crstar32_20260218_205804"
BARCODE_TRANSLATION_DEFAULT="/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt.gz"

STAR_RUN="${STAR_RUN:-$STAR_RUN_DEFAULT}"
CR_RUN="${CR_RUN:-$CR_RUN_DEFAULT}"
OUT_DIR="${OUT_DIR:-/tmp/gex_feature_parity_$(date +%Y%m%d_%H%M%S)}"
CR_FILTERED_BARCODES="${CR_FILTERED_BARCODES:-}"
NORMALIZE_BARCODES="${NORMALIZE_BARCODES:-1}"
BARCODE_TRANSLATION="${BARCODE_TRANSLATION:-$BARCODE_TRANSLATION_DEFAULT}"
TRANSLATION_DIRECTION="${TRANSLATION_DIRECTION:-left-to-right}"
TRANSLATE_SIDE="${TRANSLATE_SIDE:-both}"
GENE_CORR_MIN_COUNTS="${GENE_CORR_MIN_COUNTS:-20}"
GENE_CORR_MIN_CELLS_PCT="${GENE_CORR_MIN_CELLS_PCT:-0.01}"
CR_RAW_MATRIX_BASENAME="${CR_RAW_MATRIX_BASENAME:-matrix.mtx}"
CR_FILTERED_MATRIX_BASENAME="${CR_FILTERED_MATRIX_BASENAME:-matrix.mtx}"
STAR_RAW_MATRIX_BASENAME="${STAR_RAW_MATRIX_BASENAME:-matrix.mtx}"
STAR_FILTERED_MATRIX_BASENAME="${STAR_FILTERED_MATRIX_BASENAME:-matrix.mtx}"

COMPARE_MEX_SCRIPT="${COMPARE_MEX_SCRIPT:-${REPO_ROOT}/tests/compare_feature_mex.py}"
COMPARE_BC_SCRIPT="${COMPARE_BC_SCRIPT:-${REPO_ROOT}/scripts/compare_barcode_sets.py}"
ADDITIONAL_METRICS_SCRIPT="${ADDITIONAL_METRICS_SCRIPT:-${REPO_ROOT}/scripts/report_additional_parity_metrics.py}"

usage() {
  cat <<EOF
Usage: $(basename "$0") [--star-run DIR] [--cr-run DIR] [--out-dir DIR] [options]

Runs GEX + CRISPR feature parity checks:
1) Raw (all available barcodes)
2) Raw restricted to CR filtered barcodes
3) Filtered-vs-filtered comparison

Options:
  --cr-filtered-barcodes FILE   Explicit CR filtered barcode file (txt/csv/.gz)
  --barcode-translation FILE    Two-column translation file for barcode normalization
  --translation-direction DIR   left-to-right or right-to-left (default: left-to-right)
  --translate-side SIDE         star|cr|both|none (default: both)
  --gene-corr-min-counts INT    Gene-level corr min total counts per gene (default: 20)
  --gene-corr-min-cells-pct F   Gene-level corr min expressing-cell fraction (default: 0.01)
  --cr-raw-matrix-basename N    CR raw matrix basename (default: matrix.mtx)
  --cr-filtered-matrix-basename N
                                CR filtered matrix basename (default: matrix.mtx)
  --star-raw-matrix-basename N  STAR raw matrix basename (default: matrix.mtx)
  --star-filtered-matrix-basename N
                                STAR filtered matrix basename (default: matrix.mtx)
  --no-barcode-normalization    Disable translation-based normalization

Defaults:
  STAR_RUN=${STAR_RUN_DEFAULT}
  CR_RUN=${CR_RUN_DEFAULT}
  OUT_DIR=/tmp/gex_feature_parity_<timestamp>
  BARCODE_TRANSLATION=${BARCODE_TRANSLATION_DEFAULT}
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --star-run)
      STAR_RUN="$2"
      shift 2
      ;;
    --cr-run)
      CR_RUN="$2"
      shift 2
      ;;
    --out-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --cr-filtered-barcodes)
      CR_FILTERED_BARCODES="$2"
      shift 2
      ;;
    --barcode-translation)
      BARCODE_TRANSLATION="$2"
      shift 2
      ;;
    --translation-direction)
      TRANSLATION_DIRECTION="$2"
      shift 2
      ;;
    --translate-side)
      TRANSLATE_SIDE="$2"
      shift 2
      ;;
    --gene-corr-min-counts)
      GENE_CORR_MIN_COUNTS="$2"
      shift 2
      ;;
    --gene-corr-min-cells-pct)
      GENE_CORR_MIN_CELLS_PCT="$2"
      shift 2
      ;;
    --cr-raw-matrix-basename)
      CR_RAW_MATRIX_BASENAME="$2"
      shift 2
      ;;
    --cr-filtered-matrix-basename)
      CR_FILTERED_MATRIX_BASENAME="$2"
      shift 2
      ;;
    --star-raw-matrix-basename)
      STAR_RAW_MATRIX_BASENAME="$2"
      shift 2
      ;;
    --star-filtered-matrix-basename)
      STAR_FILTERED_MATRIX_BASENAME="$2"
      shift 2
      ;;
    --no-barcode-normalization)
      NORMALIZE_BARCODES="0"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown arg: $1" >&2
      usage
      exit 2
      ;;
  esac
done

case "${TRANSLATION_DIRECTION}" in
  left-to-right|right-to-left) ;;
  *)
    echo "ERROR: invalid --translation-direction: ${TRANSLATION_DIRECTION}" >&2
    exit 2
    ;;
esac

case "${TRANSLATE_SIDE}" in
  star|cr|both|none) ;;
  *)
    echo "ERROR: invalid --translate-side: ${TRANSLATE_SIDE}" >&2
    exit 2
    ;;
esac

if ! [[ "${GENE_CORR_MIN_COUNTS}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: --gene-corr-min-counts must be a non-negative integer" >&2
  exit 2
fi

if ! [[ "${GENE_CORR_MIN_CELLS_PCT}" =~ ^([0-9]+([.][0-9]+)?|[.][0-9]+)$ ]]; then
  echo "ERROR: --gene-corr-min-cells-pct must be a numeric value in [0,1]" >&2
  exit 2
fi
if awk "BEGIN { exit !(${GENE_CORR_MIN_CELLS_PCT} >= 0 && ${GENE_CORR_MIN_CELLS_PCT} <= 1) }"; then
  :
else
  echo "ERROR: --gene-corr-min-cells-pct must be within [0,1]" >&2
  exit 2
fi

for required in \
  "${STAR_RUN}/Solo.out/GeneFull/raw" \
  "${STAR_RUN}/Solo.out/GeneFull/filtered" \
  "${STAR_RUN}/outs/raw_feature_bc_matrix" \
  "${STAR_RUN}/outs/filtered_feature_bc_matrix" \
  "${CR_RUN}/outs/raw_feature_bc_matrix" \
  "${CR_RUN}/outs/filtered_feature_bc_matrix" \
  "${COMPARE_MEX_SCRIPT}" \
  "${COMPARE_BC_SCRIPT}" \
  "${ADDITIONAL_METRICS_SCRIPT}"; do
  [[ -e "${required}" ]] || { echo "Missing required path: ${required}" >&2; exit 2; }
done

if [[ "${NORMALIZE_BARCODES}" == "1" ]]; then
  [[ -f "${BARCODE_TRANSLATION}" ]] || {
    echo "Missing barcode translation file: ${BARCODE_TRANSLATION}" >&2
    exit 2
  }
fi

if [[ -z "${CR_FILTERED_BARCODES}" && -d "${CR_RUN}/outs/per_sample_outs" ]]; then
  mapfile -t SAMPLE_FILTERS < <(
    find "${CR_RUN}/outs/per_sample_outs" -type f -name "sample_filtered_barcodes.csv" | sort
  )
  if [[ ${#SAMPLE_FILTERS[@]} -eq 1 ]]; then
    CR_FILTERED_BARCODES="${SAMPLE_FILTERS[0]}"
  elif [[ ${#SAMPLE_FILTERS[@]} -gt 1 ]]; then
    echo "ERROR: multiple sample_filtered_barcodes.csv files found under per_sample_outs." >&2
    echo "Provide one explicitly via --cr-filtered-barcodes." >&2
    exit 2
  fi
fi

if [[ -z "${CR_FILTERED_BARCODES}" ]]; then
  CR_FILTERED_BARCODES="${CR_RUN}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
fi
[[ -e "${CR_FILTERED_BARCODES}" ]] || {
  echo "Missing CR filtered barcode file: ${CR_FILTERED_BARCODES}" >&2
  exit 2
}

MEX_TRANSLATE_ARGS=()
BC_TRANSLATE_ARGS=()
if [[ "${NORMALIZE_BARCODES}" == "1" ]]; then
  MEX_TRANSLATE_ARGS+=(
    --barcode-translation "${BARCODE_TRANSLATION}"
    --translation-direction "${TRANSLATION_DIRECTION}"
    --translate "${TRANSLATE_SIDE}"
  )
  BC_TRANSLATE_ARGS+=(
    --barcode-translation "${BARCODE_TRANSLATION}"
    --translation-direction "${TRANSLATION_DIRECTION}"
  )
fi

mkdir -p "${OUT_DIR}"
REPORT_MAIN="${OUT_DIR}/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt"
REPORT_FILTERED_SETS="${OUT_DIR}/FILTERED_BARCODE_SET_OVERLAP.txt"

{
  echo "GEX + Feature Parity Report"
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "star_run=${STAR_RUN}"
  echo "cr_run=${CR_RUN}"
  echo "cr_filtered_barcodes=${CR_FILTERED_BARCODES}"
  echo "barcode_normalization=${NORMALIZE_BARCODES}"
  if [[ "${NORMALIZE_BARCODES}" == "1" ]]; then
    echo "barcode_translation=${BARCODE_TRANSLATION}"
    echo "translation_direction=${TRANSLATION_DIRECTION}"
    echo "translate_side=${TRANSLATE_SIDE}"
  fi
  echo "gene_corr_min_counts=${GENE_CORR_MIN_COUNTS}"
  echo "gene_corr_min_cells_pct=${GENE_CORR_MIN_CELLS_PCT}"
  echo "cr_raw_matrix_basename=${CR_RAW_MATRIX_BASENAME}"
  echo "cr_filtered_matrix_basename=${CR_FILTERED_MATRIX_BASENAME}"
  echo "star_raw_matrix_basename=${STAR_RAW_MATRIX_BASENAME}"
  echo "star_filtered_matrix_basename=${STAR_FILTERED_MATRIX_BASENAME}"
  echo

  echo "=== 1) GEX RAW (all available barcodes) ==="
  python3 "${COMPARE_MEX_SCRIPT}" \
    "${CR_RUN}/outs/raw_feature_bc_matrix" \
    "${STAR_RUN}/Solo.out/GeneFull/raw" \
    --matrix-basename-a "${CR_RAW_MATRIX_BASENAME}" \
    --matrix-basename-b "${STAR_RAW_MATRIX_BASENAME}" \
    "${MEX_TRANSLATE_ARGS[@]}" \
    --feature-type "Gene Expression"
  echo

  echo "=== 2) GEX RAW restricted to CR filtered cells ==="
  python3 "${COMPARE_MEX_SCRIPT}" \
    "${CR_RUN}/outs/raw_feature_bc_matrix" \
    "${STAR_RUN}/Solo.out/GeneFull/raw" \
    --matrix-basename-a "${CR_RAW_MATRIX_BASENAME}" \
    --matrix-basename-b "${STAR_RAW_MATRIX_BASENAME}" \
    "${MEX_TRANSLATE_ARGS[@]}" \
    --feature-type "Gene Expression" \
    --barcode-filter "${CR_FILTERED_BARCODES}" \
    --barcode-filter-side both
  echo

  echo "=== 3) GEX FILTERED vs FILTERED ==="
  python3 "${COMPARE_MEX_SCRIPT}" \
    "${CR_RUN}/outs/filtered_feature_bc_matrix" \
    "${STAR_RUN}/Solo.out/GeneFull/filtered" \
    --matrix-basename-a "${CR_FILTERED_MATRIX_BASENAME}" \
    --matrix-basename-b "${STAR_FILTERED_MATRIX_BASENAME}" \
    "${MEX_TRANSLATE_ARGS[@]}" \
    --feature-type "Gene Expression"
  echo

  echo "=== 4) FEATURE RAW (CRISPR Guide Capture; all available barcodes) ==="
  python3 "${COMPARE_MEX_SCRIPT}" \
    "${CR_RUN}/outs/raw_feature_bc_matrix" \
    "${STAR_RUN}/outs/raw_feature_bc_matrix" \
    "${MEX_TRANSLATE_ARGS[@]}" \
    --feature-type "CRISPR Guide Capture"
  echo

  echo "=== 5) FEATURE RAW restricted to CR filtered cells ==="
  python3 "${COMPARE_MEX_SCRIPT}" \
    "${CR_RUN}/outs/raw_feature_bc_matrix" \
    "${STAR_RUN}/outs/raw_feature_bc_matrix" \
    "${MEX_TRANSLATE_ARGS[@]}" \
    --feature-type "CRISPR Guide Capture" \
    --barcode-filter "${CR_FILTERED_BARCODES}" \
    --barcode-filter-side both
  echo

  echo "=== 6) FEATURE FILTERED vs FILTERED ==="
  python3 "${COMPARE_MEX_SCRIPT}" \
    "${CR_RUN}/outs/filtered_feature_bc_matrix" \
    "${STAR_RUN}/outs/filtered_feature_bc_matrix" \
    "${MEX_TRANSLATE_ARGS[@]}" \
    --feature-type "CRISPR Guide Capture"
  echo

  echo "=== 7) ADDITIONAL METRICS (Feature Calls + GEX Correlations) ==="
  ADDITIONAL_ARGS=()
  if [[ "${NORMALIZE_BARCODES}" == "1" ]]; then
    ADDITIONAL_ARGS+=(
      --barcode-translation "${BARCODE_TRANSLATION}"
      --translation-direction "${TRANSLATION_DIRECTION}"
      --translate "${TRANSLATE_SIDE}"
    )
  fi
  python3 "${ADDITIONAL_METRICS_SCRIPT}" \
    --cr-run "${CR_RUN}" \
    --star-run "${STAR_RUN}" \
    --cr-filtered-barcodes "${CR_FILTERED_BARCODES}" \
    --gene-corr-min-counts "${GENE_CORR_MIN_COUNTS}" \
    --gene-corr-min-cells-pct "${GENE_CORR_MIN_CELLS_PCT}" \
    --cr-raw-matrix-basename "${CR_RAW_MATRIX_BASENAME}" \
    --cr-filtered-matrix-basename "${CR_FILTERED_MATRIX_BASENAME}" \
    --star-raw-matrix-basename "${STAR_RAW_MATRIX_BASENAME}" \
    --star-filtered-matrix-basename "${STAR_FILTERED_MATRIX_BASENAME}" \
    "${ADDITIONAL_ARGS[@]}"
} | tee "${REPORT_MAIN}"

{
  echo "Filtered Barcode-Set Overlap"
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "star_run=${STAR_RUN}"
  echo "cr_run=${CR_RUN}"
  echo

  echo "=== GEX filtered set overlap (normalized barcodes) ==="
  python3 "${COMPARE_BC_SCRIPT}" \
    "${STAR_RUN}/Solo.out/GeneFull/filtered/barcodes.tsv" \
    "${CR_FILTERED_BARCODES}" \
    --label-a "star_gex_filtered" \
    --label-b "cr_filtered" \
    --strip-suffix \
    "${BC_TRANSLATE_ARGS[@]}"
  echo

  echo "=== Feature filtered set overlap (normalized barcodes) ==="
  python3 "${COMPARE_BC_SCRIPT}" \
    "${STAR_RUN}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" \
    "${CR_FILTERED_BARCODES}" \
    --label-a "star_feature_filtered" \
    --label-b "cr_filtered" \
    --strip-suffix \
    "${BC_TRANSLATE_ARGS[@]}"
} | tee "${REPORT_FILTERED_SETS}"

cat <<EOF
Reports written:
  ${REPORT_MAIN}
  ${REPORT_FILTERED_SETS}
EOF
