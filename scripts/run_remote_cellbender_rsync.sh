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
CELLBENDER_LAYER="denoised"
CELLBENDER_FLAGS=""
CELLBENDER_USE_GPU="0"
CELLBENDER_GPU_DEVICE=""
FEATURE_GATHER_IMAGE="${FEATURE_GATHER_IMAGE:-biodepot/gather_features:latest}"
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
  --cellbender-gpu            Run CellBender with --gpus all on the remote host
  --cellbender-gpu-device ID  Pin CellBender to a specific remote GPU device
  --cellbender-cpu-cores N    CPU cores passed through to remove-background
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
    --cellbender-gpu)
      CELLBENDER_USE_GPU="1"
      shift
      ;;
    --cellbender-gpu-device)
      CELLBENDER_USE_GPU="1"
      CELLBENDER_GPU_DEVICE="$2"
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
RUN_DIR="$(realpath "$(dirname "${DOWNSTREAM_DIR}")/run")"
PROPAGATE_LAYER="${REPO_ROOT}/scripts/propagate_anndata_layer.py"
ADD_CELLBENDER_LAYER="${REPO_ROOT}/scripts/add_cellbender_layer_from_h5.py"
FEATURE_GATHER_SCRIPT="${REPO_ROOT}/scripts/integrate_feature_library.py"
INSPECT_ANNDATA="${REPO_ROOT}/../scRNA-seq/utilities/inspect_anndata.py"

for required in \
  "${COUNTS_H5AD}" \
  "${UNFILTERED_H5AD}" \
  "${FILTERED_H5AD}" \
  "${DEFAULT_SINGLET_H5AD}" \
  "${PROPAGATE_LAYER}" \
  "${ADD_CELLBENDER_LAYER}" \
  "${FEATURE_GATHER_SCRIPT}" \
  "${INSPECT_ANNDATA}"
do
  [[ -f "${required}" ]] || { echo "ERROR: missing required file ${required}" >&2; exit 1; }
done
[[ -d "${RUN_DIR}" ]] || { echo "ERROR: missing run dir ${RUN_DIR}" >&2; exit 1; }

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
echo "CellBender GPU: ${CELLBENDER_USE_GPU}"
if [[ -n "${CELLBENDER_GPU_DEVICE}" ]]; then
  echo "CellBender GPU device: ${CELLBENDER_GPU_DEVICE}"
fi

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
  "${CELLBENDER_USE_GPU}" \
  "${CELLBENDER_GPU_DEVICE:-__EMPTY__}" \
  "${CELLBENDER_LAYER}" \
  "${CELLBENDER_CPU_CORES}" \
  "${CELLBENDER_FLAGS:-__EMPTY__}" <<'EOF' >> "${LOCAL_LOG}" 2>&1
set -euo pipefail
REMOTE_WORK_DIR="$1"
REMOTE_OUTPUT_DIR="$2"
CELLBENDER_IMAGE="$3"
CELLBENDER_USE_GPU="$4"
CELLBENDER_GPU_DEVICE="$5"
CELLBENDER_LAYER="$6"
CELLBENDER_CPU_CORES="$7"
CELLBENDER_FLAGS="${8-}"
if [[ "${CELLBENDER_GPU_DEVICE}" == "__EMPTY__" ]]; then
  CELLBENDER_GPU_DEVICE=""
fi
if [[ "${CELLBENDER_FLAGS}" == "__EMPTY__" ]]; then
  CELLBENDER_FLAGS=""
fi

mkdir -p "${REMOTE_WORK_DIR}/.numba" "${REMOTE_WORK_DIR}/.matplotlib" "${REMOTE_OUTPUT_DIR}"

DOCKER_ARGS=(
  run --rm
  --user "$(id -u):$(id -g)"
  -v "${REMOTE_WORK_DIR}:${REMOTE_WORK_DIR}"
  -w "${REMOTE_WORK_DIR}"
  -e "NUMBA_CACHE_DIR=${REMOTE_WORK_DIR}/.numba"
  -e "MPLCONFIGDIR=${REMOTE_WORK_DIR}/.matplotlib"
)
if [[ "${CELLBENDER_USE_GPU}" == "1" ]]; then
  if [[ -n "${CELLBENDER_GPU_DEVICE}" ]]; then
    DOCKER_ARGS+=(--gpus "device=${CELLBENDER_GPU_DEVICE}")
  else
    DOCKER_ARGS+=(--gpus all)
  fi
fi
CELLBENDER_CMD=(
  cellbender
  remove-background
  --input "${REMOTE_WORK_DIR}/unfiltered_counts.h5ad"
  --output "${REMOTE_OUTPUT_DIR}/cellbender_counts.h5"
  --cpu-threads "${CELLBENDER_CPU_CORES}"
)
if [[ "${CELLBENDER_USE_GPU}" == "1" ]]; then
  CELLBENDER_CMD+=(--cuda)
fi
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
    -v "${REPO_ROOT}:${REPO_ROOT}:ro" \
    -w "${DOWNSTREAM_DIR}" \
    -e "NUMBA_CACHE_DIR=${DOWNSTREAM_DIR}/.numba" \
    -e "MPLCONFIGDIR=${DOWNSTREAM_DIR}/.matplotlib" \
    "${CELLBENDER_IMAGE}" \
    python \
    "${ADD_CELLBENDER_LAYER}" \
    --cellbender-h5 "${CELLBENDER_CB_FILE}" \
    --layer-name "${CELLBENDER_LAYER}" \
    --input-h5ad "${COUNTS_H5AD}" \
    --output-h5ad "${COUNTS_H5AD}"

  docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "${DOWNSTREAM_DIR}:${DOWNSTREAM_DIR}" \
    -v "${REPO_ROOT}:${REPO_ROOT}:ro" \
    -w "${DOWNSTREAM_DIR}" \
    -e "NUMBA_CACHE_DIR=${DOWNSTREAM_DIR}/.numba" \
    -e "MPLCONFIGDIR=${DOWNSTREAM_DIR}/.matplotlib" \
    "${CELLBENDER_IMAGE}" \
    python \
    "${ADD_CELLBENDER_LAYER}" \
    --cellbender-h5 "${CELLBENDER_CB_FILE}" \
    --layer-name "${CELLBENDER_LAYER}" \
    --input-h5ad "${UNFILTERED_H5AD}" \
    --output-h5ad "${UNFILTERED_H5AD}"

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

FEATURE_OUTPUT_ROOT="${DOWNSTREAM_DIR}/feature_libraries"
mapfile -t FEATURE_LIBRARY_DIRS < <(find "${RUN_DIR}/cr_assign" -type f -name 'pf_library_provenance.tsv' -print 2>/dev/null | sed 's#/pf_library_provenance.tsv$##' | sort)
if (( ${#FEATURE_LIBRARY_DIRS[@]} > 0 )); then
  CRISPR_CALLS_CSV="${RUN_DIR}/outs/crispr_analysis/protospacer_calls_per_cell.csv"
  CRISPR_LIBRARY_COUNT=0
  for feature_library in "${FEATURE_LIBRARY_DIRS[@]}"; do
    feature_type_dir="$(basename "$(dirname "$(dirname "${feature_library}")")")"
    if [[ "${feature_type_dir}" == "CRISPR_Guide_Capture" ]]; then
      ((CRISPR_LIBRARY_COUNT += 1))
    fi
  done

  if (( CRISPR_LIBRARY_COUNT > 1 )) && [[ -f "${CRISPR_CALLS_CSV}" ]]; then
    echo "ERROR: Found ${CRISPR_LIBRARY_COUNT} CRISPR feature libraries but only one global call file at ${CRISPR_CALLS_CSV}" >&2
    exit 1
  fi

  COUNTS_TARGETS=("${COUNTS_H5AD}" "${UNFILTERED_H5AD}" "${FILTERED_H5AD}" "${DEFAULT_SINGLET_H5AD}")
  if [[ -f "${FINAL_H5AD}" ]]; then
    COUNTS_TARGETS+=("${FINAL_H5AD}")
  fi

  for feature_library in "${FEATURE_LIBRARY_DIRS[@]}"; do
    FEATURE_GATHER_ARGS=(
      run --rm
      --user "$(id -u):$(id -g)"
      -v "${RUN_DIR}:${RUN_DIR}:ro"
      -v "${DOWNSTREAM_DIR}:${DOWNSTREAM_DIR}"
      -v "${REPO_ROOT}:${REPO_ROOT}:ro"
      "${FEATURE_GATHER_IMAGE}"
      python3 "${FEATURE_GATHER_SCRIPT}"
      --library-dir "${feature_library}"
      --feature-output-root "${FEATURE_OUTPUT_ROOT}"
    )
    for counts_target in "${COUNTS_TARGETS[@]}"; do
      FEATURE_GATHER_ARGS+=(--counts-h5ad "${counts_target}")
    done

    feature_type_dir="$(basename "$(dirname "$(dirname "${feature_library}")")")"
    if [[ "${feature_type_dir}" == "CRISPR_Guide_Capture" ]] && [[ -f "${CRISPR_CALLS_CSV}" ]]; then
      FEATURE_GATHER_ARGS+=(--calls-csv "${CRISPR_CALLS_CSV}")
      if (( ${#FEATURE_LIBRARY_DIRS[@]} == 1 )); then
        FEATURE_GATHER_ARGS+=(--set-generic-aliases)
      fi
    fi

    docker "${FEATURE_GATHER_ARGS[@]}"
  done
fi

if [[ "${KEEP_REMOTE}" != "1" ]]; then
  ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new "${REMOTE_HOST}" \
    "rm -rf '${REMOTE_JOB_ROOT}'"
fi

echo "PASS: remote CellBender step complete"
echo "Downstream dir: ${DOWNSTREAM_DIR}"
echo "Log: ${LOCAL_LOG}"
