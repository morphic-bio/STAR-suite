#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

REMOTE_HOST=""
REMOTE_ROOT=""
SAMPLE_DIR=""
OUTPUT_NAME="downstream_genefull_velocyto_cellbender_remote"
RUN_CELLBENDER="0"
ADAPTIVE_FILTER="0"
KEEP_REMOTE="0"
SYNC_IMAGES="1"
CELLBENDER_CPU_CORES=""
LOCAL_LOG=""
DOCKER_IMAGE=""
CELLBENDER_IMAGE=""
FEATURE_GATHER_IMAGE=""
EXTRA_ARGS=()

usage() {
  cat <<'EOF'
Usage:
  run_remote_scrna_downstream_rsync.sh --sample-dir PATH --remote-host HOST --remote-root PATH [options]

This stages only the UCSF downstream inputs to a remote host via rsync,
executes the existing downstream wrapper there, and rsyncs the output dir back.

Options:
  --sample-dir PATH          Local sample directory containing run/
  --remote-host HOST         SSH target
  --remote-root PATH         Remote staging root on local disk (not NFS)
  --output-name NAME         Output directory name copied back under sample-dir
                             (default: downstream_genefull_velocyto_cellbender_remote)
  --run-cellbender           Enable CellBender
  --adaptive-filter          Enable adaptive n_genes and MT percentage filtering
  --cellbender-cpu-cores N   Pass through to downstream wrapper
  --docker-image IMAGE       Pass through to downstream wrapper
  --cellbender-image IMAGE   Pass through to downstream wrapper
  --feature-gather-image IMG Pass through to downstream wrapper
  --local-log PATH           Local log file for remote execution
  --no-sync-images          Do not sync local Docker images to remote first
  --keep-remote              Leave staged remote directory in place
  --help                     Show this help

Any extra arguments are passed through to run_scrna_downstream_gene_full_velocyto.sh.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-dir)
      SAMPLE_DIR="$2"
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
    --output-name)
      OUTPUT_NAME="$2"
      shift 2
      ;;
    --run-cellbender)
      RUN_CELLBENDER="1"
      shift
      ;;
    --adaptive-filter)
      ADAPTIVE_FILTER="1"
      shift
      ;;
    --cellbender-cpu-cores)
      CELLBENDER_CPU_CORES="$2"
      shift 2
      ;;
    --docker-image)
      DOCKER_IMAGE="$2"
      shift 2
      ;;
    --cellbender-image)
      CELLBENDER_IMAGE="$2"
      shift 2
      ;;
    --feature-gather-image)
      FEATURE_GATHER_IMAGE="$2"
      shift 2
      ;;
    --local-log)
      LOCAL_LOG="$2"
      shift 2
      ;;
    --no-sync-images)
      SYNC_IMAGES="0"
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
      EXTRA_ARGS+=("$1")
      shift
      ;;
  esac
done

[[ -n "${SAMPLE_DIR}" ]] || { echo "ERROR: --sample-dir is required" >&2; exit 1; }
[[ -n "${REMOTE_HOST}" ]] || { echo "ERROR: --remote-host is required" >&2; exit 1; }
[[ -n "${REMOTE_ROOT}" ]] || { echo "ERROR: --remote-root is required" >&2; exit 1; }

SAMPLE_DIR="$(realpath "${SAMPLE_DIR}")"
RUN_DIR="${SAMPLE_DIR}/run"
[[ -d "${RUN_DIR}" ]] || { echo "ERROR: missing run dir ${RUN_DIR}" >&2; exit 1; }

for required in \
  "${RUN_DIR}/outs/filtered_feature_bc_matrix" \
  "${RUN_DIR}/outs/raw_feature_bc_matrix" \
  "${RUN_DIR}/outs/raw_velocyto_feature_bc_matrix"
do
  [[ -d "${required}" ]] || { echo "ERROR: missing required input ${required}" >&2; exit 1; }
done

INSPECT_ANNDATA_LOCAL="/mnt/pikachu/scRNA-seq/utilities/inspect_anndata.py"
[[ -f "${INSPECT_ANNDATA_LOCAL}" ]] || { echo "ERROR: missing ${INSPECT_ANNDATA_LOCAL}" >&2; exit 1; }
command -v rsync >/dev/null 2>&1 || { echo "ERROR: rsync is required" >&2; exit 1; }
command -v ssh >/dev/null 2>&1 || { echo "ERROR: ssh is required" >&2; exit 1; }
command -v docker >/dev/null 2>&1 || { echo "ERROR: local docker is required" >&2; exit 1; }

STAMP="$(date +%Y%m%d_%H%M%S)"
SAMPLE_NAME="$(basename "${SAMPLE_DIR}")"
REMOTE_JOB_ROOT="${REMOTE_ROOT%/}/${SAMPLE_NAME}_${STAMP}"
REMOTE_REPO_ROOT="${REMOTE_JOB_ROOT}/repo"
REMOTE_SCRNA_ROOT="${REMOTE_JOB_ROOT}/scRNA-seq"
REMOTE_SAMPLE_ROOT="${REMOTE_JOB_ROOT}/sample"
REMOTE_RUN_DIR="${REMOTE_SAMPLE_ROOT}/run"
REMOTE_OUTPUT_DIR="${REMOTE_SAMPLE_ROOT}/${OUTPUT_NAME}"
LOCAL_OUTPUT_DIR="${SAMPLE_DIR}/${OUTPUT_NAME}"
LOCAL_LOG="${LOCAL_LOG:-${SAMPLE_DIR}/${OUTPUT_NAME}.remote.log}"
EFFECTIVE_DOCKER_IMAGE="${DOCKER_IMAGE:-biodepot/scrna-matrices:latest}"
EFFECTIVE_CELLBENDER_IMAGE="${CELLBENDER_IMAGE:-biodepot/cellbender:0.3.2}"
EFFECTIVE_FEATURE_GATHER_IMAGE="${FEATURE_GATHER_IMAGE:-biodepot/gather_features:latest}"

SCRIPTS_TO_STAGE=(
  "scripts/run_scrna_downstream_gene_full_velocyto.sh"
  "scripts/build_gene_full_velocyto_h5ad.py"
  "scripts/run_star_cell_doublets.R"
  "scripts/integrate_feature_library.py"
  "scripts/postprocess_downstream_filters.py"
  "scripts/compute_adaptive_qc_threshold.py"
  "scripts/scrna_mt_adaptive.py"
  "scripts/apply_adaptive_mt_filter.py"
  "scripts/generate_qc_histogram_mt_adaptive.py"
  "scripts/propagate_anndata_layer.py"
  "scripts/add_cellbender_layer_from_h5.py"
)

echo "=== Remote UCSF downstream rsync runner ==="
echo "Sample dir: ${SAMPLE_DIR}"
echo "Remote host: ${REMOTE_HOST}"
echo "Remote job root: ${REMOTE_JOB_ROOT}"
echo "Output dir: ${LOCAL_OUTPUT_DIR}"
echo "Local log: ${LOCAL_LOG}"

sync_remote_image() {
  local image="$1"
  local local_id remote_id
  local_id="$(docker image inspect "${image}" --format '{{.Id}}' 2>/dev/null || true)"
  [[ -n "${local_id}" ]] || { echo "ERROR: local Docker image not found: ${image}" >&2; exit 1; }
  remote_id="$(ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new "${REMOTE_HOST}" \
    "docker image inspect '${image}' --format '{{.Id}}' 2>/dev/null || true")"
  if [[ -n "${remote_id}" && "${remote_id}" == "${local_id}" ]]; then
    echo "Remote image already matches local: ${image} (${local_id})"
    return 0
  fi
  echo "Syncing Docker image to remote: ${image}"
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

if [[ "${SYNC_IMAGES}" == "1" ]]; then
  sync_remote_image "${EFFECTIVE_DOCKER_IMAGE}"
  sync_remote_image "${EFFECTIVE_FEATURE_GATHER_IMAGE}"
  if [[ "${RUN_CELLBENDER}" == "1" ]]; then
    sync_remote_image "${EFFECTIVE_CELLBENDER_IMAGE}"
  fi
fi

ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new "${REMOTE_HOST}" \
  "mkdir -p '${REMOTE_REPO_ROOT}/scripts' '${REMOTE_SCRNA_ROOT}/utilities' '${REMOTE_RUN_DIR}/outs' '${REMOTE_SAMPLE_ROOT}'"

for rel in "${SCRIPTS_TO_STAGE[@]}"; do
  rsync -az "${REPO_ROOT}/${rel}" "${REMOTE_HOST}:${REMOTE_REPO_ROOT}/${rel}"
done
rsync -az "${INSPECT_ANNDATA_LOCAL}" "${REMOTE_HOST}:${REMOTE_SCRNA_ROOT}/utilities/inspect_anndata.py"

for rel in \
  "outs/filtered_feature_bc_matrix" \
  "outs/raw_feature_bc_matrix" \
  "outs/raw_velocyto_feature_bc_matrix"
do
  rsync -az "${RUN_DIR}/${rel}/" "${REMOTE_HOST}:${REMOTE_RUN_DIR}/${rel}/"
done

if [[ -d "${RUN_DIR}/outs/crispr_analysis" ]]; then
  rsync -az "${RUN_DIR}/outs/crispr_analysis/" "${REMOTE_HOST}:${REMOTE_RUN_DIR}/outs/crispr_analysis/"
fi
if [[ -d "${RUN_DIR}/cr_assign" ]]; then
  rsync -az "${RUN_DIR}/cr_assign/" "${REMOTE_HOST}:${REMOTE_RUN_DIR}/cr_assign/"
fi

REMOTE_CMD=(
  "${REMOTE_REPO_ROOT}/scripts/run_scrna_downstream_gene_full_velocyto.sh"
  --run-dir "${REMOTE_RUN_DIR}"
  --output-dir "${REMOTE_OUTPUT_DIR}"
)
if [[ "${RUN_CELLBENDER}" == "1" ]]; then
  REMOTE_CMD+=(--run-cellbender)
fi
if [[ "${ADAPTIVE_FILTER}" == "1" ]]; then
  REMOTE_CMD+=(--adaptive-filter)
fi
if [[ -n "${CELLBENDER_CPU_CORES}" ]]; then
  REMOTE_CMD+=(--cellbender-cpu-cores "${CELLBENDER_CPU_CORES}")
fi
if [[ -n "${DOCKER_IMAGE}" ]]; then
  REMOTE_CMD+=(--docker-image "${DOCKER_IMAGE}")
fi
if [[ -n "${CELLBENDER_IMAGE}" ]]; then
  REMOTE_CMD+=(--cellbender-image "${CELLBENDER_IMAGE}")
fi
if [[ -n "${FEATURE_GATHER_IMAGE}" ]]; then
  REMOTE_CMD+=(--feature-gather-image "${FEATURE_GATHER_IMAGE}")
fi
if (( ${#EXTRA_ARGS[@]} > 0 )); then
  REMOTE_CMD+=("${EXTRA_ARGS[@]}")
fi

printf 'Remote command:' > "${LOCAL_LOG}"
for arg in "${REMOTE_CMD[@]}"; do
  printf ' %q' "${arg}" >> "${LOCAL_LOG}"
done
printf '\n' >> "${LOCAL_LOG}"

ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new "${REMOTE_HOST}" bash -s -- \
  "${REMOTE_REPO_ROOT}" \
  "${REMOTE_SCRNA_ROOT}" \
  "${LOCAL_LOG}" \
  "${REMOTE_CMD[@]}" <<'EOF' >> "${LOCAL_LOG}" 2>&1
set -euo pipefail
REMOTE_REPO_ROOT="$1"
REMOTE_SCRNA_ROOT="$2"
shift 3
export SC_RNA_SEQ_ROOT="${REMOTE_SCRNA_ROOT}"
export INSPECT_ANNDATA="${REMOTE_SCRNA_ROOT}/utilities/inspect_anndata.py"
"$@"
EOF

rm -rf "${LOCAL_OUTPUT_DIR}"
rsync -az "${REMOTE_HOST}:${REMOTE_OUTPUT_DIR}/" "${LOCAL_OUTPUT_DIR}/"

if [[ "${KEEP_REMOTE}" != "1" ]]; then
  ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new "${REMOTE_HOST}" \
    "rm -rf '${REMOTE_JOB_ROOT}'"
fi

echo "PASS: remote downstream complete"
echo "Output dir: ${LOCAL_OUTPUT_DIR}"
echo "Log: ${LOCAL_LOG}"
