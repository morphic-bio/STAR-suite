#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RUNNER="${REPO_ROOT}/scripts/run_remote_cellbender_rsync.sh"

REMOTE_HOST="${REMOTE_HOST:-128.208.252.232}"
REMOTE_ROOT="${REMOTE_ROOT:-/home/lhhung/ucsf_remote_cellbender_smoke}"
SOURCE_DOWNSTREAM_DIR="${SOURCE_DOWNSTREAM_DIR:-${REPO_ROOT}/tests/ucsf_corrected_production_100k_output_20260329_171919/samples/EBs2_2/downstream_genefull_velocyto_cellbender}"
CELLBENDER_IMAGE="${CELLBENDER_IMAGE:-biodepot/cellbender:0.3.2}"
CELLBENDER_CPU_CORES="${CELLBENDER_CPU_CORES:-8}"
NO_SYNC_IMAGE="${NO_SYNC_IMAGE:-1}"
STAMP="$(date +%Y%m%d_%H%M%S)"
SMOKE_DIR="${REPO_ROOT}/tests/remote_cellbender_100k_smoke_${STAMP}"
LOG_FILE="${REPO_ROOT}/tests/remote_cellbender_100k_smoke_${STAMP}.log"

[[ -x "${RUNNER}" ]] || { echo "ERROR: missing runner ${RUNNER}" >&2; exit 1; }
[[ -d "${SOURCE_DOWNSTREAM_DIR}" ]] || { echo "ERROR: missing source downstream dir ${SOURCE_DOWNSTREAM_DIR}" >&2; exit 1; }

mkdir -p "${SMOKE_DIR}"
for rel in \
  counts.h5ad \
  unfiltered_counts.h5ad \
  filtered_counts.h5ad \
  default_singlet_filtered_counts.h5ad
do
  cp -f "${SOURCE_DOWNSTREAM_DIR}/${rel}" "${SMOKE_DIR}/${rel}"
done

echo "=== Remote CellBender 100K smoke ==="
echo "Source downstream dir: ${SOURCE_DOWNSTREAM_DIR}"
echo "Smoke dir: ${SMOKE_DIR}"
echo "Remote host: ${REMOTE_HOST}"
echo "Remote root: ${REMOTE_ROOT}"
echo "CellBender image: ${CELLBENDER_IMAGE}"
echo "No sync image: ${NO_SYNC_IMAGE}"
echo "Log file: ${LOG_FILE}"

RUNNER_ARGS=(
  --downstream-dir "${SMOKE_DIR}"
  --remote-host "${REMOTE_HOST}"
  --remote-root "${REMOTE_ROOT}"
  --cellbender-image "${CELLBENDER_IMAGE}"
  --cellbender-cpu-cores "${CELLBENDER_CPU_CORES}"
  --local-log "${LOG_FILE}"
)
if [[ "${NO_SYNC_IMAGE}" == "1" ]]; then
  RUNNER_ARGS+=(--no-sync-image)
fi

"${RUNNER}" \
  "${RUNNER_ARGS[@]}"

if [[ ! -f "${SMOKE_DIR}/cellbender/cellbender_counts.h5" && ! -f "${SMOKE_DIR}/cellbender/CELLBENDER_FAILED.txt" ]]; then
  echo "ERROR: expected remote CellBender output or failure note in ${SMOKE_DIR}/cellbender" >&2
  exit 1
fi

[[ -f "${SMOKE_DIR}/final_counts.h5ad" ]] || { echo "ERROR: missing ${SMOKE_DIR}/final_counts.h5ad" >&2; exit 1; }
[[ -f "${SMOKE_DIR}/summary.txt" ]] || { echo "ERROR: missing ${SMOKE_DIR}/summary.txt" >&2; exit 1; }

echo "PASS: remote CellBender 100K smoke"
echo "Smoke dir: ${SMOKE_DIR}"
echo "Log: ${LOG_FILE}"
