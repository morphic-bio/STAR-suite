#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

REMOTE_HOST=""
REMOTE_ROOT=""
SAMPLE_DIR=""
OUTPUT_NAME="downstream_genefull_velocyto_cellbender"
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
  run_remote_multiome_post_mex_rsync.sh --sample-dir PATH --remote-host HOST --remote-root PATH [options]

Stages only post-MEX multiome inputs to a remote host: STAR GeneFull/Velocyto
MEX plus locally materialized ATAC peak MEX/metrics. The remote host runs RNA
downstream h5ad/CellBender and MuData construction, then syncs final downstream
and mudata outputs back. It does not need STAR, Chromap, or libchromap.

Options:
  --sample-dir PATH          Local star_sample directory containing run/
  --remote-host HOST         SSH target
  --remote-root PATH         Remote staging root on local disk
  --output-name NAME         RNA downstream output dir name under sample-dir
  --run-cellbender
  --adaptive-filter          Enable adaptive n_genes and MT percentage filtering
  --cellbender-cpu-cores N
  --docker-image IMAGE
  --cellbender-image IMAGE
  --feature-gather-image IMG
  --local-log PATH
  --no-sync-images
  --keep-remote
  --help

Any extra arguments are passed through to run_scrna_downstream_gene_full_velocyto.sh.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-dir) SAMPLE_DIR="$2"; shift 2 ;;
    --remote-host) REMOTE_HOST="$2"; shift 2 ;;
    --remote-root) REMOTE_ROOT="$2"; shift 2 ;;
    --output-name) OUTPUT_NAME="$2"; shift 2 ;;
    --run-cellbender) RUN_CELLBENDER="1"; shift ;;
    --adaptive-filter) ADAPTIVE_FILTER="1"; shift ;;
    --cellbender-cpu-cores) CELLBENDER_CPU_CORES="$2"; shift 2 ;;
    --docker-image) DOCKER_IMAGE="$2"; shift 2 ;;
    --cellbender-image) CELLBENDER_IMAGE="$2"; shift 2 ;;
    --feature-gather-image) FEATURE_GATHER_IMAGE="$2"; shift 2 ;;
    --local-log) LOCAL_LOG="$2"; shift 2 ;;
    --no-sync-images) SYNC_IMAGES="0"; shift ;;
    --keep-remote) KEEP_REMOTE="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) EXTRA_ARGS+=("$1"); shift ;;
  esac
done

[[ -n "${SAMPLE_DIR}" ]] || { echo "ERROR: --sample-dir is required" >&2; exit 1; }
[[ -n "${REMOTE_HOST}" ]] || { echo "ERROR: --remote-host is required" >&2; exit 1; }
[[ -n "${REMOTE_ROOT}" ]] || { echo "ERROR: --remote-root is required" >&2; exit 1; }

SAMPLE_DIR="$(realpath "${SAMPLE_DIR}")"
OUT_DIR="$(dirname "${SAMPLE_DIR}")"
RUN_DIR="${SAMPLE_DIR}/run"
ATAC_DIR="${OUT_DIR}/atac"
ATAC_MEX="${ATAC_DIR}/peak_mex"
ATAC_METRICS="${ATAC_DIR}/atac_metrics.tsv"
ATAC_SIDECAR="${RUN_DIR}/atac_fragments.bin"
ATAC_PEAKS="${RUN_DIR}/atac_peaks.narrowPeak"
LOCAL_OUTPUT_DIR="${SAMPLE_DIR}/${OUTPUT_NAME}"
LOCAL_MUDATA_DIR="${OUT_DIR}/mudata"
LOCAL_LOG="${LOCAL_LOG:-${OUT_DIR}/logs/remote_post_mex.log}"

[[ -d "${RUN_DIR}" ]] || { echo "ERROR: missing run dir ${RUN_DIR}" >&2; exit 1; }
for required in \
  "${RUN_DIR}/outs/filtered_feature_bc_matrix" \
  "${RUN_DIR}/outs/raw_feature_bc_matrix" \
  "${RUN_DIR}/outs/raw_velocyto_feature_bc_matrix" \
  "${ATAC_MEX}/matrix.mtx.gz" \
  "${ATAC_MEX}/barcodes.tsv.gz" \
  "${ATAC_MEX}/features.tsv.gz" \
  "${ATAC_METRICS}"
do
  [[ -e "${required}" ]] || { echo "ERROR: missing required post-MEX input ${required}" >&2; exit 1; }
done

INSPECT_ANNDATA_LOCAL="/mnt/pikachu/scRNA-seq/utilities/inspect_anndata.py"
[[ -f "${INSPECT_ANNDATA_LOCAL}" ]] || { echo "ERROR: missing ${INSPECT_ANNDATA_LOCAL}" >&2; exit 1; }
command -v rsync >/dev/null 2>&1 || { echo "ERROR: rsync is required" >&2; exit 1; }
command -v ssh >/dev/null 2>&1 || { echo "ERROR: ssh is required" >&2; exit 1; }
command -v docker >/dev/null 2>&1 || { echo "ERROR: local docker is required" >&2; exit 1; }

STAMP="$(date +%Y%m%d_%H%M%S)"
SAMPLE_NAME="$(basename "${OUT_DIR}")"
REMOTE_JOB_ROOT="${REMOTE_ROOT%/}/${SAMPLE_NAME}_${STAMP}"
REMOTE_REPO_ROOT="${REMOTE_JOB_ROOT}/repo"
REMOTE_SCRNA_ROOT="${REMOTE_JOB_ROOT}/scRNA-seq"
REMOTE_SAMPLE_ROOT="${REMOTE_JOB_ROOT}/sample"
REMOTE_RUN_DIR="${REMOTE_SAMPLE_ROOT}/run"
REMOTE_ATAC_DIR="${REMOTE_SAMPLE_ROOT}/atac"
REMOTE_ATAC_MEX="${REMOTE_ATAC_DIR}/peak_mex"
REMOTE_ATAC_METRICS="${REMOTE_ATAC_DIR}/atac_metrics.tsv"
REMOTE_OUTPUT_DIR="${REMOTE_SAMPLE_ROOT}/${OUTPUT_NAME}"
REMOTE_MUDATA_DIR="${REMOTE_SAMPLE_ROOT}/mudata"
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
  "scripts/build_multiome_mudata.py"
)

echo "=== Remote multiome post-MEX rsync runner ==="
echo "Sample dir: ${SAMPLE_DIR}"
echo "Remote host: ${REMOTE_HOST}"
echo "Remote job root: ${REMOTE_JOB_ROOT}"
echo "RNA output dir: ${LOCAL_OUTPUT_DIR}"
echo "MuData dir: ${LOCAL_MUDATA_DIR}"
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
  "mkdir -p '${REMOTE_REPO_ROOT}/scripts' '${REMOTE_SCRNA_ROOT}/utilities' '${REMOTE_RUN_DIR}/outs' '${REMOTE_ATAC_MEX}' '${REMOTE_MUDATA_DIR}'"

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

if [[ -d "${RUN_DIR}/outs/filtered_velocyto_feature_bc_matrix" ]]; then
  rsync -az "${RUN_DIR}/outs/filtered_velocyto_feature_bc_matrix/" \
    "${REMOTE_HOST}:${REMOTE_RUN_DIR}/outs/filtered_velocyto_feature_bc_matrix/"
fi
rsync -az "${ATAC_MEX}/" "${REMOTE_HOST}:${REMOTE_ATAC_MEX}/"
rsync -az "${ATAC_METRICS}" "${REMOTE_HOST}:${REMOTE_ATAC_METRICS}"
[[ -f "${ATAC_PEAKS}" ]] && rsync -az "${ATAC_PEAKS}" "${REMOTE_HOST}:${REMOTE_SAMPLE_ROOT}/run/atac_peaks.narrowPeak"

REMOTE_DOWNSTREAM_CMD=(
  "${REMOTE_REPO_ROOT}/scripts/run_scrna_downstream_gene_full_velocyto.sh"
  --run-dir "${REMOTE_RUN_DIR}"
  --output-dir "${REMOTE_OUTPUT_DIR}"
)
[[ "${RUN_CELLBENDER}" == "1" ]] && REMOTE_DOWNSTREAM_CMD+=(--run-cellbender)
[[ "${ADAPTIVE_FILTER}" == "1" ]] && REMOTE_DOWNSTREAM_CMD+=(--adaptive-filter)
[[ -n "${CELLBENDER_CPU_CORES}" ]] && REMOTE_DOWNSTREAM_CMD+=(--cellbender-cpu-cores "${CELLBENDER_CPU_CORES}")
[[ -n "${DOCKER_IMAGE}" ]] && REMOTE_DOWNSTREAM_CMD+=(--docker-image "${DOCKER_IMAGE}")
[[ -n "${CELLBENDER_IMAGE}" ]] && REMOTE_DOWNSTREAM_CMD+=(--cellbender-image "${CELLBENDER_IMAGE}")
[[ -n "${FEATURE_GATHER_IMAGE}" ]] && REMOTE_DOWNSTREAM_CMD+=(--feature-gather-image "${FEATURE_GATHER_IMAGE}")
if (( ${#EXTRA_ARGS[@]} > 0 )); then
  REMOTE_DOWNSTREAM_CMD+=("${EXTRA_ARGS[@]}")
fi

printf 'Remote downstream command:' > "${LOCAL_LOG}"
for arg in "${REMOTE_DOWNSTREAM_CMD[@]}"; do
  printf ' %q' "${arg}" >> "${LOCAL_LOG}"
done
printf '\n' >> "${LOCAL_LOG}"

ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new "${REMOTE_HOST}" bash -s -- \
  "${REMOTE_REPO_ROOT}" \
  "${REMOTE_SCRNA_ROOT}" \
  "${REMOTE_OUTPUT_DIR}" \
  "${REMOTE_ATAC_MEX}" \
  "${REMOTE_ATAC_METRICS}" \
  "${REMOTE_MUDATA_DIR}" \
  "${LOCAL_OUTPUT_DIR}/final_counts.h5ad" \
  "${LOCAL_OUTPUT_DIR}/filtered_counts.h5ad" \
  "${ATAC_MEX}" \
  "${ATAC_SIDECAR}" \
  "${ATAC_PEAKS}" \
  "${ATAC_METRICS}" \
  "${REMOTE_DOWNSTREAM_CMD[@]}" <<'EOF' >> "${LOCAL_LOG}" 2>&1
set -euo pipefail
REMOTE_REPO_ROOT="$1"
REMOTE_SCRNA_ROOT="$2"
REMOTE_OUTPUT_DIR="$3"
REMOTE_ATAC_MEX="$4"
REMOTE_ATAC_METRICS="$5"
REMOTE_MUDATA_DIR="$6"
LOCAL_RNA_FINAL_SOURCE="$7"
LOCAL_RNA_FILTERED_SOURCE="$8"
LOCAL_ATAC_MEX_SOURCE="$9"
LOCAL_ATAC_SIDECAR_SOURCE="${10}"
LOCAL_ATAC_PEAKS_SOURCE="${11}"
LOCAL_ATAC_METRICS_SOURCE="${12}"
shift 12

export SC_RNA_SEQ_ROOT="${REMOTE_SCRNA_ROOT}"
export INSPECT_ANNDATA="${REMOTE_SCRNA_ROOT}/utilities/inspect_anndata.py"
"$@"

RNA_UNFILTERED="${REMOTE_OUTPUT_DIR}/final_counts.h5ad"
if [[ ! -f "${RNA_UNFILTERED}" ]]; then
  RNA_UNFILTERED="${REMOTE_OUTPUT_DIR}/unfiltered_counts.h5ad"
fi
RNA_FILTERED="${REMOTE_OUTPUT_DIR}/filtered_counts.h5ad"
[[ -f "${RNA_UNFILTERED}" ]] || { echo "ERROR: missing remote RNA h5ad ${RNA_UNFILTERED}" >&2; exit 1; }
[[ -f "${RNA_FILTERED}" ]] || { echo "ERROR: missing remote filtered RNA h5ad ${RNA_FILTERED}" >&2; exit 1; }

MUDATA_PYTHON=python3
if ! python3 - <<'PY' >/dev/null 2>&1
import mudata
PY
then
  VENV="${REMOTE_MUDATA_DIR}/mudata_venv"
  python3 -m venv --system-site-packages "${VENV}"
  "${VENV}/bin/python" -m pip install --quiet mudata
  MUDATA_PYTHON="${VENV}/bin/python"
fi

"${MUDATA_PYTHON}" "${REMOTE_REPO_ROOT}/scripts/build_multiome_mudata.py" \
  --rna-h5ad "${RNA_UNFILTERED}" \
  --atac-mex-dir "${REMOTE_ATAC_MEX}" \
  --per-barcode-metrics "${REMOTE_ATAC_METRICS}" \
  --metrics-barcode-column barcode \
  --require-rna-velocyto-layers \
  --cell-call-source star_downstream_h5ad_chromap_atac \
  --rna-source "${LOCAL_RNA_FINAL_SOURCE}" \
  --atac-source "${LOCAL_ATAC_MEX_SOURCE}" \
  --fragments-source "${LOCAL_ATAC_SIDECAR_SOURCE}" \
  --peaks-source "${LOCAL_ATAC_PEAKS_SOURCE}" \
  --evidence-source "${LOCAL_ATAC_METRICS_SOURCE}" \
  --y-removal-enabled true \
  --output-h5mu "${REMOTE_MUDATA_DIR}/star_chromap_unfiltered_multiome.h5mu"

"${MUDATA_PYTHON}" "${REMOTE_REPO_ROOT}/scripts/build_multiome_mudata.py" \
  --rna-h5ad "${RNA_FILTERED}" \
  --atac-mex-dir "${REMOTE_ATAC_MEX}" \
  --per-barcode-metrics "${REMOTE_ATAC_METRICS}" \
  --metrics-barcode-column barcode \
  --all-barcodes-are-cells \
  --allow-empty-barcode-intersection \
  --require-rna-velocyto-layers \
  --cell-call-source star_downstream_filtered_h5ad_chromap_atac \
  --rna-source "${LOCAL_RNA_FILTERED_SOURCE}" \
  --atac-source "${LOCAL_ATAC_MEX_SOURCE}" \
  --fragments-source "${LOCAL_ATAC_SIDECAR_SOURCE}" \
  --peaks-source "${LOCAL_ATAC_PEAKS_SOURCE}" \
  --evidence-source "${LOCAL_ATAC_METRICS_SOURCE}" \
  --y-removal-enabled true \
  --output-h5mu "${REMOTE_MUDATA_DIR}/star_chromap_filtered_multiome.h5mu"

"${MUDATA_PYTHON}" - <<PY
import mudata as md
for path in ["${REMOTE_MUDATA_DIR}/star_chromap_unfiltered_multiome.h5mu", "${REMOTE_MUDATA_DIR}/star_chromap_filtered_multiome.h5mu"]:
    m = md.read_h5mu(path)
    rna = m.mod["rna"]
    atac = m.mod["atac"]
    missing = sorted({"counts", "spliced", "unspliced", "ambiguous"} - set(rna.layers))
    if missing:
        raise SystemExit(f"{path}: missing RNA layers {missing}")
    if "counts" not in atac.layers:
        raise SystemExit(f"{path}: missing ATAC counts layer")
    print(path)
    print(f"  obs={m.n_obs} rna_vars={rna.n_vars} atac_vars={atac.n_vars}")
    if m.n_obs == 0:
        print(f"WARNING: {path} has zero observations after RNA/ATAC barcode intersection; acceptable for sparse smoke tests.")
PY
EOF

{
  echo "Local sync-back start: $(date -Is)"
  echo "  remote downstream: ${REMOTE_HOST}:${REMOTE_OUTPUT_DIR}/"
  echo "  local downstream:  ${LOCAL_OUTPUT_DIR}/"
  echo "  remote mudata:     ${REMOTE_HOST}:${REMOTE_MUDATA_DIR}/"
  echo "  local mudata:      ${LOCAL_MUDATA_DIR}/"
  rm -rf "${LOCAL_OUTPUT_DIR}" "${LOCAL_MUDATA_DIR}"
  mkdir -p "${LOCAL_OUTPUT_DIR}" "${LOCAL_MUDATA_DIR}"
  rsync -az --info=stats2 "${REMOTE_HOST}:${REMOTE_OUTPUT_DIR}/" "${LOCAL_OUTPUT_DIR}/"
  rsync -az --exclude 'mudata_venv/' --info=stats2 "${REMOTE_HOST}:${REMOTE_MUDATA_DIR}/" "${LOCAL_MUDATA_DIR}/"
  echo "Local sync-back complete: $(date -Is)"
} >> "${LOCAL_LOG}" 2>&1

if [[ "${KEEP_REMOTE}" != "1" ]]; then
  ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new "${REMOTE_HOST}" \
    "rm -rf '${REMOTE_JOB_ROOT}'"
fi

[[ -f "${LOCAL_MUDATA_DIR}/star_chromap_unfiltered_multiome.h5mu" ]] || {
  echo "ERROR: missing synced unfiltered MuData output" >&2
  exit 1
}
[[ -f "${LOCAL_MUDATA_DIR}/star_chromap_filtered_multiome.h5mu" ]] || {
  echo "ERROR: missing synced filtered MuData output" >&2
  exit 1
}

cat > "${OUT_DIR}/REMOTE_POST_MEX_READY.txt" <<EOF
remote_host=${REMOTE_HOST}
remote_job_root=${REMOTE_JOB_ROOT}
local_downstream_dir=${LOCAL_OUTPUT_DIR}
local_mudata_dir=${LOCAL_MUDATA_DIR}
completed_at=$(date -Is)
EOF

echo "PASS: remote multiome post-MEX complete" | tee -a "${LOCAL_LOG}"
