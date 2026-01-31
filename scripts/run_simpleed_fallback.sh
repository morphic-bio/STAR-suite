#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_simpleed_fallback.sh --raw-mex DIR --filtered-barcodes PATH [options]

Required:
  --raw-mex DIR             Solo raw MEX dir (matrix.mtx + barcodes.tsv)
  --filtered-barcodes PATH  Output barcodes.tsv path (may not exist)

Options:
  --min-cells N             Minimum filtered cells before fallback (default: 10)
  --force                   Force fallback run regardless of existing output
  --mode simple|full         Backend mode for libscrna (default: simple)
  --sim-n N                 Monte Carlo simulations for full ED
  --ed-retain-count N       Limit retained cells for ED (0 = all)
  --use-fdr-gate            Use FDR gate instead of raw p-value (default)
  --use-raw-pvalue          Use raw p-value gate
  --fdr X                   FDR threshold
  --raw-pvalue X            Raw p-value threshold
  --include-zero-umis       Include zero-UMI barcodes (default: exclude)
  --simpleed-bin PATH       Path to scrna_simpleed binary
  --out-dir DIR             Write detailed EmptyDrops outputs to DIR
EOF
}

RAW_MEX=""
FILTERED_BARCODES=""
MIN_CELLS=10
FORCE=0
SIMPLEED_BIN=""
OUT_DIR=""
MODE="simple"
SIM_N=""
ED_RETAIN=""
USE_FDR=1
FDR=""
RAW_PVALUE=""
INCLUDE_ZERO=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --raw-mex) RAW_MEX="$2"; shift 2;;
    --filtered-barcodes) FILTERED_BARCODES="$2"; shift 2;;
    --min-cells) MIN_CELLS="$2"; shift 2;;
    --force) FORCE=1; shift 1;;
    --mode) MODE="$2"; shift 2;;
    --sim-n) SIM_N="$2"; shift 2;;
    --ed-retain-count) ED_RETAIN="$2"; shift 2;;
    --use-fdr-gate) USE_FDR=1; shift 1;;
    --use-raw-pvalue) USE_FDR=0; shift 1;;
    --fdr) FDR="$2"; shift 2;;
    --raw-pvalue) RAW_PVALUE="$2"; shift 2;;
    --include-zero-umis) INCLUDE_ZERO=1; shift 1;;
    --simpleed-bin) SIMPLEED_BIN="$2"; shift 2;;
    --out-dir) OUT_DIR="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1" >&2; usage; exit 2;;
  esac
done

if [[ -z "${RAW_MEX}" || -z "${FILTERED_BARCODES}" ]]; then
  usage
  exit 2
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -z "${SIMPLEED_BIN}" ]]; then
  SIMPLEED_BIN="${REPO_ROOT}/core/features/libscrna/bin/scrna_simpleed"
fi

if [[ ! -x "${SIMPLEED_BIN}" ]]; then
  echo "scrna_simpleed not found at ${SIMPLEED_BIN}; building..." >&2
  make -C "${REPO_ROOT}/core/features/libscrna" tools
fi

matrix_path=""
barcodes_path=""
tmp_matrix=""
tmp_barcodes=""

if [[ -f "${RAW_MEX}/matrix.mtx" ]]; then
  matrix_path="${RAW_MEX}/matrix.mtx"
elif [[ -f "${RAW_MEX}/matrix.mtx.gz" ]]; then
  tmp_matrix="$(mktemp /tmp/simpleed_matrix.XXXXXX.mtx)"
  zcat "${RAW_MEX}/matrix.mtx.gz" > "${tmp_matrix}"
  matrix_path="${tmp_matrix}"
else
  echo "Missing matrix.mtx(.gz) in ${RAW_MEX}" >&2
  exit 1
fi

if [[ -f "${RAW_MEX}/barcodes.tsv" ]]; then
  barcodes_path="${RAW_MEX}/barcodes.tsv"
elif [[ -f "${RAW_MEX}/barcodes.tsv.gz" ]]; then
  tmp_barcodes="$(mktemp /tmp/simpleed_barcodes.XXXXXX.tsv)"
  zcat "${RAW_MEX}/barcodes.tsv.gz" > "${tmp_barcodes}"
  barcodes_path="${tmp_barcodes}"
else
  echo "Missing barcodes.tsv(.gz) in ${RAW_MEX}" >&2
  exit 1
fi

cleanup() {
  [[ -n "${tmp_matrix}" ]] && rm -f "${tmp_matrix}"
  [[ -n "${tmp_barcodes}" ]] && rm -f "${tmp_barcodes}"
  return 0
}
trap cleanup EXIT

should_fallback=0
if [[ "${FORCE}" -eq 1 ]]; then
  should_fallback=1
else
  if [[ ! -f "${FILTERED_BARCODES}" ]]; then
    should_fallback=1
  else
    count="$(wc -l < "${FILTERED_BARCODES}" || echo 0)"
    if [[ "${count}" -lt "${MIN_CELLS}" ]]; then
      should_fallback=1
    fi
  fi
fi

if [[ "${should_fallback}" -eq 0 ]]; then
  echo "SimpleED fallback not triggered (filtered barcodes OK)." >&2
  exit 0
fi

mkdir -p "$(dirname "${FILTERED_BARCODES}")"

echo "Running SimpleED fallback..." >&2
args=(--matrix "${matrix_path}" --barcodes "${barcodes_path}" --out "${FILTERED_BARCODES}" --mode "${MODE}")
if [[ -n "${SIM_N}" ]]; then
  args+=(--sim-n "${SIM_N}")
fi
if [[ -n "${ED_RETAIN}" ]]; then
  args+=(--ed-retain-count "${ED_RETAIN}")
fi
if [[ -n "${FDR}" ]]; then
  args+=(--fdr "${FDR}")
fi
if [[ -n "${RAW_PVALUE}" ]]; then
  args+=(--raw-pvalue "${RAW_PVALUE}")
fi
if [[ "${USE_FDR}" -eq 1 ]]; then
  args+=(--use-fdr-gate)
fi
if [[ "${INCLUDE_ZERO}" -eq 1 ]]; then
  args+=(--include-zero-umis)
fi
if [[ -n "${OUT_DIR}" ]]; then
  args+=(--out-dir "${OUT_DIR}")
fi

"${SIMPLEED_BIN}" "${args[@]}"

marker_dir="$(dirname "${FILTERED_BARCODES}")"
touch "${marker_dir}/USED_LIBSCRNA_SIMPLEED"

echo "SimpleED fallback complete: ${FILTERED_BARCODES}" >&2
