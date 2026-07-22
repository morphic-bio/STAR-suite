#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RUNNER="${REPO_ROOT}/scripts/run_remote_scrna_downstream_rsync.sh"

REMOTE_HOST="${REMOTE_HOST:-128.208.252.232}"
REMOTE_ROOT="${REMOTE_ROOT:-/home/lhhung/ucsf_remote_downstream_smoke}"
SAMPLE_DIR="${SAMPLE_DIR:-${REPO_ROOT}/tests/ucsf_corrected_production_100k_output_20260329_171919/samples/EBs2_2}"
CELLBENDER_CPU_CORES="${CELLBENDER_CPU_CORES:-8}"
DOCKER_IMAGE="${DOCKER_IMAGE:-biodepot/scrna-matrices:20260330-ce24c46cbaf9}"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUTPUT_NAME="downstream_remote_100k_smoke_${STAMP}"
LOG_FILE="${REPO_ROOT}/tests/remote_ucsf_downstream_100k_smoke_${STAMP}.log"

[[ -x "${RUNNER}" ]] || { echo "ERROR: missing runner ${RUNNER}" >&2; exit 1; }
[[ -d "${SAMPLE_DIR}/run" ]] || { echo "ERROR: missing sample run dir ${SAMPLE_DIR}/run" >&2; exit 1; }

echo "=== Remote UCSF 100K downstream smoke ==="
echo "Sample dir: ${SAMPLE_DIR}"
echo "Remote host: ${REMOTE_HOST}"
echo "Remote root: ${REMOTE_ROOT}"
echo "Output name: ${OUTPUT_NAME}"
echo "Docker image: ${DOCKER_IMAGE}"
echo "Log file: ${LOG_FILE}"

"${RUNNER}" \
  --sample-dir "${SAMPLE_DIR}" \
  --remote-host "${REMOTE_HOST}" \
  --remote-root "${REMOTE_ROOT}" \
  --output-name "${OUTPUT_NAME}" \
  --docker-image "${DOCKER_IMAGE}" \
  --adaptive-filter \
  --run-cellbender \
  --cellbender-cpu-cores "${CELLBENDER_CPU_CORES}" \
  --local-log "${LOG_FILE}"

OUTPUT_DIR="${SAMPLE_DIR}/${OUTPUT_NAME}"
[[ -f "${OUTPUT_DIR}/summary.txt" ]] || { echo "ERROR: missing ${OUTPUT_DIR}/summary.txt" >&2; exit 1; }
if [[ ! -f "${OUTPUT_DIR}/final_counts.h5ad" && ! -f "${OUTPUT_DIR}/cellbender/CELLBENDER_FAILED.txt" ]]; then
  echo "ERROR: expected final_counts.h5ad or CELLBENDER_FAILED.txt in ${OUTPUT_DIR}" >&2
  exit 1
fi

echo "PASS: remote UCSF 100K downstream smoke"
echo "Output dir: ${OUTPUT_DIR}"
echo "Log: ${LOG_FILE}"
