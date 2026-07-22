#!/usr/bin/env bash
set -euo pipefail

INPUT_H5AD=""
OUTPUT_DIR=""
EXPECTED_CELLS=""
TOTAL_DROPLETS=""
EPOCHS=""
MODEL="full"
CPU_THREADS="${CELLBENDER_CPU_CORES:-4}"
LOW_COUNT_THRESHOLD=""
FORCE_CELL_UMI_PRIOR=""
FORCE_EMPTY_UMI_PRIOR=""
ESTIMATOR=""

usage() {
  cat <<'EOF'
Usage: run_ucsf_cellbender_stage.sh --input-h5ad FILE --output-dir DIR [options]

Run only CellBender remove-background. Layer integration is a separate CPU stage.

Options:
  --expected-cells N
  --total-droplets-included N
  --epochs N
  --model NAME
  --cpu-threads N
  --low-count-threshold N
  --force-cell-umi-prior N
  --force-empty-umi-prior N
  --estimator NAME
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-h5ad) INPUT_H5AD="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --expected-cells) EXPECTED_CELLS="$2"; shift 2 ;;
    --total-droplets-included) TOTAL_DROPLETS="$2"; shift 2 ;;
    --epochs) EPOCHS="$2"; shift 2 ;;
    --model) MODEL="$2"; shift 2 ;;
    --cpu-threads) CPU_THREADS="$2"; shift 2 ;;
    --low-count-threshold) LOW_COUNT_THRESHOLD="$2"; shift 2 ;;
    --force-cell-umi-prior) FORCE_CELL_UMI_PRIOR="$2"; shift 2 ;;
    --force-empty-umi-prior) FORCE_EMPTY_UMI_PRIOR="$2"; shift 2 ;;
    --estimator) ESTIMATOR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -f "${INPUT_H5AD}" ]] || { echo "ERROR: missing --input-h5ad ${INPUT_H5AD}" >&2; exit 2; }
[[ -n "${OUTPUT_DIR}" ]] || { echo "ERROR: --output-dir is required" >&2; exit 2; }

INPUT_H5AD="$(realpath "${INPUT_H5AD}")"
OUTPUT_DIR="$(realpath -m "${OUTPUT_DIR}")"
mkdir -p \
  "${OUTPUT_DIR}" \
  "${OUTPUT_DIR}/.numba" \
  "${OUTPUT_DIR}/.matplotlib" \
  "${OUTPUT_DIR}/.tmp"
export NUMBA_CACHE_DIR="${OUTPUT_DIR}/.numba"
export MPLCONFIGDIR="${OUTPUT_DIR}/.matplotlib"
export TMPDIR="${OUTPUT_DIR}/.tmp"

OUTPUT_H5="${OUTPUT_DIR}/cellbender_counts.h5"
LOG="${OUTPUT_DIR}/cellbender_counts.log"
rm -f "${OUTPUT_H5}" "${OUTPUT_DIR}/CELLBENDER_FAILED.txt"

CMD=(
  cellbender remove-background
  --input "${INPUT_H5AD}"
  --output "${OUTPUT_H5}"
  --cuda
  --cpu-threads "${CPU_THREADS}"
  --model "${MODEL}"
)
[[ -n "${EXPECTED_CELLS}" ]] && CMD+=(--expected-cells "${EXPECTED_CELLS}")
[[ -n "${TOTAL_DROPLETS}" ]] && CMD+=(--total-droplets-included "${TOTAL_DROPLETS}")
[[ -n "${EPOCHS}" ]] && CMD+=(--epochs "${EPOCHS}")
[[ -n "${LOW_COUNT_THRESHOLD}" ]] && CMD+=(--low-count-threshold "${LOW_COUNT_THRESHOLD}")
[[ -n "${FORCE_CELL_UMI_PRIOR}" ]] && CMD+=(--force-cell-umi-prior "${FORCE_CELL_UMI_PRIOR}")
[[ -n "${FORCE_EMPTY_UMI_PRIOR}" ]] && CMD+=(--force-empty-umi-prior "${FORCE_EMPTY_UMI_PRIOR}")
[[ -n "${ESTIMATOR}" ]] && CMD+=(--estimator "${ESTIMATOR}")

(
  cd "${OUTPUT_DIR}"
  "${CMD[@]}"
) 2>&1 | tee "${LOG}"

[[ -s "${OUTPUT_H5}" ]] || {
  printf 'CellBender produced no output\ninput_h5ad=%s\nlog=%s\n' \
    "${INPUT_H5AD}" "${LOG}" > "${OUTPUT_DIR}/CELLBENDER_FAILED.txt"
  exit 1
}

python3 - "${INPUT_H5AD}" "${OUTPUT_H5}" "${LOG}" "${OUTPUT_DIR}/CELLBENDER_COMPLETE.json" <<'PY'
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

input_h5ad, output_h5, log_path, manifest_path = map(Path, sys.argv[1:])
payload = {
    "schema": "star_suite.ucsf_cellbender/v1",
    "status": "complete",
    "completed_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
    "input_h5ad": str(input_h5ad),
    "output_h5": str(output_h5),
    "output_bytes": output_h5.stat().st_size,
    "log": str(log_path),
}
manifest_path.write_text(json.dumps(payload, indent=2) + "\n")
PY

echo "PASS: UCSF CellBender GPU stage: ${OUTPUT_H5}"
