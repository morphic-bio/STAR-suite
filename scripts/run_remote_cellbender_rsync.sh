#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

REMOTE_HOST=""
REMOTE_ROOT=""
DOWNSTREAM_DIR=""
KEEP_REMOTE="0"
SYNC_IMAGE="1"
CELLBENDER_CPU_CORES="8"
CELLBENDER_IMAGE="biodepot/cellbender:0.3.2"
CELLBENDER_LAYER="cellbender"
CELLBENDER_FLAGS=""
LOCAL_LOG=""
EXTRA_REMOTE_ENV=()

usage() {
  cat <<'EOF'
Usage:
  run_remote_cellbender_rsync.sh --downstream-dir PATH --remote-host HOST --remote-root PATH [options]

This stages an existing downstream output directory to a remote host, runs only
the CellBender remove-background step there, copies the CellBender outputs back,
and locally propagates the denoised layer into the h5ad outputs.

Options:
  --downstream-dir PATH       Local downstream output directory containing
                              counts.h5ad, unfiltered_counts.h5ad, filtered_counts.h5ad,
                              and default_singlet_filtered_counts.h5ad
  --remote-host HOST          SSH target
  --remote-root PATH          Remote staging root on local disk (not NFS)
  --cellbender-image IMG      CellBender image (default: biodepot/cellbender:0.3.2)
  --cellbender-cpu-cores N    CPU cores passed through to remove_noise.sh
  --cellbender-layer NAME     Output layer name (default: cellbender)
  --cellbender-flags STR      Extra CellBender flags
  --local-log PATH            Local log path
  --no-sync-image            Do not sync the local CellBender image to remote
  --keep-remote               Leave remote staging behind
  --help                      Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --downstream-dir)
      DOWNSTREAM_DIR="$2"
      shift 2
      ;;
    --remote-host)
      REMOTE_HOST="$2"
      shift 2
      ;;
    --remote-root)
      REMOTE_ROOT="$2"
      shift 2
      ;;
    --cellbender-image)
      CELLBENDER_IMAGE="$2"
      shift 2
      ;;
    --cellbender-cpu-cores)
      CELLBENDER_CPU_CORES="$2"
      shift 2
      ;;
    --cellbender-layer)
      CELLBENDER_LAYER="$2"
      shift 2
      ;;
    --cellbender-flags)
      CELLBENDER_FLAGS="$2"
      shift 2
      ;;
    --local-log)
      LOCAL_LOG="$2"
      shift 2
      ;;
    --no-sync-image)
      SYNC_IMAGE="0"
      shift
      ;;
    --keep-remote)
      KEEP_REMOTE="1"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

[[ -n "${DOWNSTREAM_DIR}" ]] || { echo "ERROR: --downstream-dir is required" >&2; exit 1; }
[[ -n "${REMOTE_HOST}" ]] || { echo "ERROR: --remote-host is required" >&2; exit 1; }
[[ -n "${REMOTE_ROOT}" ]] || { echo "ERROR: --remote-root is required" >&2; exit 1; }

DOWNSTREAM_DIR="$(realpath "${DOWNSTREAM_DIR}")"
[[ -d "${DOWNSTREAM_DIR}" ]] || { echo "ERROR: missing downstream dir ${DOWNSTREAM_DIR}" >&2; exit 1; }

COUNTS_H5AD="${DOWNSTREAM_DIR}/counts.h5ad"
UNFILTERED_H5AD="${DOWNSTREAM_DIR}/unfiltered_counts.h5ad"
FILTERED_H5AD="${DOWNSTREAM_DIR}/filtered_counts.h5ad"
DEFAULT_SINGLET_H5AD="${DOWNSTREAM_DIR}/default_singlet_filtered_counts.h5ad"
FINAL_H5AD="${DOWNSTREAM_DIR}/final_counts.h5ad"
PROPAGATE_LAYER="${REPO_ROOT}/scripts/propagate_anndata_layer.py"
INSPECT_ANNDATA="${REPO_ROOT}/../scRNA-seq/utilities/inspect_anndata.py"

for required in \
  "${COUNTS_H5AD}" \
  "${UNFILTERED_H5AD}" \
  "${FILTERED_H5AD}" \
  "${DEFAULT_SINGLET_H5AD}" \
  "${PROPAGATE_LAYER}" \
  "${INSPECT_ANNDATA}"
do
  [[ -f "${required}" ]] || { echo "ERROR: missing required file ${required}" >&2; exit 1; }
done

command -v rsync >/dev/null 2>&1 || { echo "ERROR: rsync is required" >&2; exit 1; }
command -v ssh >/dev/null 2>&1 || { echo "ERROR: ssh is required" >&2; exit 1; }
command -v docker >/dev/null 2>&1 || { echo "ERROR: local docker is required" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: local python3 is required" >&2; exit 1; }

STAMP="$(date +%Y%m%d_%H%M%S)"
SAMPLE_NAME="$(basename "$(dirname "${DOWNSTREAM_DIR}")")"
REMOTE_JOB_ROOT="${REMOTE_ROOT%/}/${SAMPLE_NAME}_cellbender_${STAMP}"
REMOTE_WORK_DIR="${REMOTE_JOB_ROOT}/work"
REMOTE_OUTPUT_DIR="${REMOTE_WORK_DIR}/cellbender"
LOCAL_CELLBENDER_DIR="${DOWNSTREAM_DIR}/cellbender"
LOCAL_LOG="${LOCAL_LOG:-${DOWNSTREAM_DIR}/remote_cellbender_${STAMP}.log}"
CELLBENDER_CB_FILE="${LOCAL_CELLBENDER_DIR}/cellbender_counts.h5"
CELLBENDER_FAILURE_NOTE="${LOCAL_CELLBENDER_DIR}/CELLBENDER_FAILED.txt"

echo "=== Remote CellBender rsync runner ==="
echo "Downstream dir: ${DOWNSTREAM_DIR}"
echo "Remote host: ${REMOTE_HOST}"
echo "Remote job root: ${REMOTE_JOB_ROOT}"
echo "Local log: ${LOCAL_LOG}"

sync_remote_image() {
  local image="$1"
  local local_id remote_id
  local_id="$(docker image inspect "${image}" --format '{{.Id}}' 2>/dev/null || true)"
  [[ -n "${local_id}" ]] || { echo "ERROR: local Docker image not found: ${image}" >&2; exit 1; }
  remote_id="$(ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new "${REMOTE_HOST}" \
    "docker image inspect '${image}' --format '{{.Id}}' 2>/dev/null || true")"
  if [[ -n "${remote_id}" && "${remote_id}" == "${local_id}" ]]; then
    echo "Remote CellBender image already matches local: ${image} (${local_id})"
    return 0
  fi
  echo "Syncing CellBender image to remote: ${image}"
  docker save "${image}" | ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new "${REMOTE_HOST}" docker load >/dev/null
  remote_id="$(ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new "${REMOTE_HOST}" \
    "docker image inspect '${image}' --format '{{.Id}}'")"
  [[ "${remote_id}" == "${local_id}" ]] || {
    echo "ERROR: remote Docker image mismatch after sync for ${image}" >&2
    echo "  local:  ${local_id}" >&2
    echo "  remote: ${remote_id}" >&2
    exit 1
  }
}

if [[ "${SYNC_IMAGE}" == "1" ]]; then
  sync_remote_image "${CELLBENDER_IMAGE}"
fi

ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new "${REMOTE_HOST}" \
  "mkdir -p '${REMOTE_WORK_DIR}' '${REMOTE_OUTPUT_DIR}'"

rsync -az "${UNFILTERED_H5AD}" "${REMOTE_HOST}:${REMOTE_WORK_DIR}/unfiltered_counts.h5ad"

printf 'Remote CellBender image: %s\n' "${CELLBENDER_IMAGE}" > "${LOCAL_LOG}"
printf 'Remote work dir: %s\n' "${REMOTE_WORK_DIR}" >> "${LOCAL_LOG}"

REMOTE_STATUS=0
set +e
ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new "${REMOTE_HOST}" bash -s -- \
  "${REMOTE_WORK_DIR}" \
  "${REMOTE_OUTPUT_DIR}" \
  "${CELLBENDER_IMAGE}" \
  "${CELLBENDER_LAYER}" \
  "${CELLBENDER_CPU_CORES}" \
  "${CELLBENDER_FLAGS}" <<'EOF' >> "${LOCAL_LOG}" 2>&1
set -euo pipefail
REMOTE_WORK_DIR="$1"
REMOTE_OUTPUT_DIR="$2"
CELLBENDER_IMAGE="$3"
CELLBENDER_LAYER="$4"
CELLBENDER_CPU_CORES="$5"
CELLBENDER_FLAGS="${6-}"

mkdir -p "${REMOTE_WORK_DIR}/.numba" "${REMOTE_WORK_DIR}/.matplotlib" "${REMOTE_OUTPUT_DIR}"

DOCKER_ARGS=(
  run --rm
  --user "$(id -u):$(id -g)"
  -v "${REMOTE_WORK_DIR}:${REMOTE_WORK_DIR}"
  -w "${REMOTE_WORK_DIR}"
  -e "NUMBA_CACHE_DIR=${REMOTE_WORK_DIR}/.numba"
  -e "MPLCONFIGDIR=${REMOTE_WORK_DIR}/.matplotlib"
)
CELLBENDER_CMD=(
  cellbender
  remove-background
  --input "${REMOTE_WORK_DIR}/unfiltered_counts.h5ad"
  --output "${REMOTE_OUTPUT_DIR}/cellbender_counts.h5"
  --cpu-threads "${CELLBENDER_CPU_CORES}"
)
if [[ -n "${CELLBENDER_FLAGS}" ]]; then
  # shellcheck disable=SC2206
  EXTRA_FLAGS=( ${CELLBENDER_FLAGS} )
  CELLBENDER_CMD+=("${EXTRA_FLAGS[@]}")
fi

docker "${DOCKER_ARGS[@]}" "${CELLBENDER_IMAGE}" "${CELLBENDER_CMD[@]}"
EOF
REMOTE_STATUS=$?
set -e

mkdir -p "${LOCAL_CELLBENDER_DIR}"
rsync -az "${REMOTE_HOST}:${REMOTE_OUTPUT_DIR}/" "${LOCAL_CELLBENDER_DIR}/"

if [[ -f "${CELLBENDER_CB_FILE}" ]]; then
  docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "${DOWNSTREAM_DIR}:${DOWNSTREAM_DIR}" \
    -w "${DOWNSTREAM_DIR}" \
    -e "NUMBA_CACHE_DIR=${DOWNSTREAM_DIR}/.numba" \
    -e "MPLCONFIGDIR=${DOWNSTREAM_DIR}/.matplotlib" \
    -e "overwrite_layer=1" \
    "${CELLBENDER_IMAGE}" \
    addCounts.py \
    -c "${CELLBENDER_CB_FILE}" \
    -l "${CELLBENDER_LAYER}" \
    -i "${COUNTS_H5AD}" \
    -o "${COUNTS_H5AD}"

  docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "${DOWNSTREAM_DIR}:${DOWNSTREAM_DIR}" \
    -w "${DOWNSTREAM_DIR}" \
    -e "NUMBA_CACHE_DIR=${DOWNSTREAM_DIR}/.numba" \
    -e "MPLCONFIGDIR=${DOWNSTREAM_DIR}/.matplotlib" \
    -e "overwrite_layer=1" \
    "${CELLBENDER_IMAGE}" \
    addCounts.py \
    -c "${CELLBENDER_CB_FILE}" \
    -l "${CELLBENDER_LAYER}" \
    -i "${UNFILTERED_H5AD}" \
    -o "${UNFILTERED_H5AD}"

  python3 "${PROPAGATE_LAYER}" \
    --source-h5ad "${UNFILTERED_H5AD}" \
    --target-h5ad "${FILTERED_H5AD}" \
    --output-h5ad "${FILTERED_H5AD}" \
    --layer-name "${CELLBENDER_LAYER}"

  python3 "${PROPAGATE_LAYER}" \
    --source-h5ad "${UNFILTERED_H5AD}" \
    --target-h5ad "${DEFAULT_SINGLET_H5AD}" \
    --output-h5ad "${DEFAULT_SINGLET_H5AD}" \
    --layer-name "${CELLBENDER_LAYER}"

  cp -f "${UNFILTERED_H5AD}" "${FINAL_H5AD}"
  python3 "${INSPECT_ANNDATA}" "${FINAL_H5AD}" > "${DOWNSTREAM_DIR}/final_counts.summary.txt"
  python3 "${INSPECT_ANNDATA}" "${FINAL_H5AD}" > "${DOWNSTREAM_DIR}/summary.txt"
  rm -f "${CELLBENDER_FAILURE_NOTE}"
else
  {
    echo "CellBender did not produce ${CELLBENDER_CB_FILE}"
    echo "remote_exit_code=${REMOTE_STATUS}"
    echo "input_h5ad=${UNFILTERED_H5AD}"
    echo "counts_h5ad=${COUNTS_H5AD}"
    echo "fallback_h5ad=${FINAL_H5AD}"
    echo "reason=sparse_or_prefiltered_input_can_fail_prior_estimation"
  } > "${CELLBENDER_FAILURE_NOTE}"
  cp -f "${UNFILTERED_H5AD}" "${FINAL_H5AD}"
  python3 "${INSPECT_ANNDATA}" "${FINAL_H5AD}" > "${DOWNSTREAM_DIR}/final_counts.summary.txt"
  python3 "${INSPECT_ANNDATA}" "${FINAL_H5AD}" > "${DOWNSTREAM_DIR}/summary.txt"
fi

if [[ "${KEEP_REMOTE}" != "1" ]]; then
  ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new "${REMOTE_HOST}" \
    "rm -rf '${REMOTE_JOB_ROOT}'"
fi

echo "PASS: remote CellBender step complete"
echo "Downstream dir: ${DOWNSTREAM_DIR}"
echo "Log: ${LOCAL_LOG}"
