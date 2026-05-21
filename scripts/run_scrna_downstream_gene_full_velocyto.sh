#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

SC_RNA_SEQ_ROOT="${SC_RNA_SEQ_ROOT:-/mnt/pikachu/scRNA-seq}"
BUILD_COUNTS="${REPO_ROOT}/scripts/build_gene_full_velocyto_h5ad.py"
DOUBLET_SCRIPT="${REPO_ROOT}/scripts/run_star_cell_doublets.R"
FEATURE_GATHER_SCRIPT="${REPO_ROOT}/scripts/integrate_feature_library.py"
POSTPROCESS_FILTERS="${REPO_ROOT}/scripts/postprocess_downstream_filters.py"
COMPUTE_ADAPTIVE_QC="${REPO_ROOT}/scripts/compute_adaptive_qc_threshold.py"
APPLY_ADAPTIVE_MT="${REPO_ROOT}/scripts/apply_adaptive_mt_filter.py"
GENERATE_QC_HISTOGRAM_MT="${REPO_ROOT}/scripts/generate_qc_histogram_mt_adaptive.py"
PROPAGATE_LAYER="${REPO_ROOT}/scripts/propagate_anndata_layer.py"
ADD_CELLBENDER_LAYER="${REPO_ROOT}/scripts/add_cellbender_layer_from_h5.py"
INSPECT_ANNDATA="${INSPECT_ANNDATA:-${SC_RNA_SEQ_ROOT}/utilities/inspect_anndata.py}"
DOCKER_IMAGE="${SCRNA_DOWNSTREAM_IMAGE:-biodepot/scrna-matrices:latest}"
CELLBENDER_IMAGE="${CELLBENDER_IMAGE:-biodepot/cellbender:0.3.2}"
FEATURE_GATHER_IMAGE="${FEATURE_GATHER_IMAGE:-biodepot/gather_features:latest}"

RUN_DIR=""
OUTPUT_DIR=""
MITO_GENES=""
MIN_GENES="${MIN_GENES:-200}"
MAX_GENES="${MAX_GENES:-2500}"
MT_PCT_CUTOFF="${MT_PCT_CUTOFF:-5}"
ADAPTIVE_FILTER="1"
N_MAD="${N_MAD:-3}"
RUN_CELLBENDER="0"
REUSE_CELLBENDER="0"
CELLBENDER_USE_GPU="0"
CELLBENDER_CPU_CORES="${CELLBENDER_CPU_CORES:-8}"
CELLBENDER_LAYER="${CELLBENDER_LAYER:-denoised}"
CELLBENDER_FLAGS="${CELLBENDER_FLAGS:-}"
PYTHON_BACKEND="${SCRNA_DOWNSTREAM_PYTHON_BACKEND:-auto}"

usage() {
  cat <<'EOF'
Usage:
  run_scrna_downstream_gene_full_velocyto.sh --run-dir <run-dir> [options]

Options:
  --run-dir PATH         STAR/CR-compat run directory containing outs/
  --output-dir PATH      Output directory (default: <run-dir>/downstream_genefull_velocyto)
  --mito-genes PATH      Optional mito genes file for combineFilters.py
  --min-genes INT        Minimum genes cutoff (default: 200)
  --max-genes INT        Maximum genes cutoff (default: 2500)
  --mt-pct-cutoff FLOAT  MT percentage floor/fixed cutoff (default: 5)
  --adaptive-filter      Use adaptive n_genes and MT percentage thresholds
  --n-mad FLOAT          Number of MADs for adaptive QC thresholds (default: 3)
  --docker-image IMAGE   Downstream container image (default: biodepot/scrna-matrices:latest)
  --feature-gather-image IMAGE
                         Feature-library integration image
                         (default: biodepot/gather_features:latest)
  --run-cellbender       Run CellBender on raw-backed unfiltered_counts.h5ad and add denoised layer
  --reuse-cellbender     Skip CellBender denoising but integrate existing cellbender_counts.h5 layer
  --cellbender-image IMG CellBender image (default: biodepot/cellbender:0.3.2)
  --cellbender-gpu       Run CellBender with --gpus all instead of CPU mode
  --cellbender-cpu-cores INT
                         CPU cores passed to CellBender (default: 8)
  --cellbender-layer NAME
                         Layer name for CellBender output (default: denoised)
  --cellbender-flags STR Additional flags passed to cellbender remove-background
  --python-backend MODE  Python helper backend: auto, host, or docker (default: auto)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-dir)
      RUN_DIR="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --mito-genes)
      MITO_GENES="$2"
      shift 2
      ;;
    --min-genes)
      MIN_GENES="$2"
      shift 2
      ;;
    --max-genes)
      MAX_GENES="$2"
      shift 2
      ;;
    --mt-pct-cutoff)
      MT_PCT_CUTOFF="$2"
      shift 2
      ;;
    --adaptive-filter)
      ADAPTIVE_FILTER="1"
      shift
      ;;
    --n-mad)
      N_MAD="$2"
      shift 2
      ;;
    --docker-image)
      DOCKER_IMAGE="$2"
      shift 2
      ;;
    --feature-gather-image)
      FEATURE_GATHER_IMAGE="$2"
      shift 2
      ;;
    --run-cellbender)
      RUN_CELLBENDER="1"
      shift
      ;;
    --reuse-cellbender)
      REUSE_CELLBENDER="1"
      shift
      ;;
    --cellbender-image)
      CELLBENDER_IMAGE="$2"
      shift 2
      ;;
    --cellbender-gpu)
      CELLBENDER_USE_GPU="1"
      shift
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
    --python-backend)
      PYTHON_BACKEND="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "${RUN_DIR}" ]]; then
  echo "ERROR: --run-dir is required" >&2
  usage >&2
  exit 1
fi

RUN_DIR="$(realpath "${RUN_DIR}")"
OUTPUT_DIR="${OUTPUT_DIR:-${RUN_DIR}/downstream_genefull_velocyto}"
OUTPUT_DIR="$(realpath -m "${OUTPUT_DIR}")"
COUNTS_H5AD="${OUTPUT_DIR}/counts.h5ad"
FILTERED_H5AD="${OUTPUT_DIR}/filtered_counts.h5ad"
UNFILTERED_H5AD="${OUTPUT_DIR}/unfiltered_counts.h5ad"
DEFAULT_SINGLET_FILTERED_H5AD="${OUTPUT_DIR}/default_singlet_filtered_counts.h5ad"
FINAL_H5AD="${OUTPUT_DIR}/final_counts.h5ad"
PRIMARY_H5AD="${FILTERED_H5AD}"

[[ -f "${BUILD_COUNTS}" ]] || { echo "ERROR: Missing helper ${BUILD_COUNTS}" >&2; exit 1; }
[[ -f "${DOUBLET_SCRIPT}" ]] || { echo "ERROR: Missing helper ${DOUBLET_SCRIPT}" >&2; exit 1; }
[[ -f "${FEATURE_GATHER_SCRIPT}" ]] || { echo "ERROR: Missing helper ${FEATURE_GATHER_SCRIPT}" >&2; exit 1; }
[[ -f "${POSTPROCESS_FILTERS}" ]] || { echo "ERROR: Missing helper ${POSTPROCESS_FILTERS}" >&2; exit 1; }
[[ -f "${COMPUTE_ADAPTIVE_QC}" ]] || { echo "ERROR: Missing helper ${COMPUTE_ADAPTIVE_QC}" >&2; exit 1; }
[[ -f "${APPLY_ADAPTIVE_MT}" ]] || { echo "ERROR: Missing helper ${APPLY_ADAPTIVE_MT}" >&2; exit 1; }
[[ -f "${GENERATE_QC_HISTOGRAM_MT}" ]] || { echo "ERROR: Missing helper ${GENERATE_QC_HISTOGRAM_MT}" >&2; exit 1; }
[[ -f "${PROPAGATE_LAYER}" ]] || { echo "ERROR: Missing helper ${PROPAGATE_LAYER}" >&2; exit 1; }
[[ -f "${ADD_CELLBENDER_LAYER}" ]] || { echo "ERROR: Missing helper ${ADD_CELLBENDER_LAYER}" >&2; exit 1; }
[[ -f "${INSPECT_ANNDATA}" ]] || { echo "ERROR: Missing helper ${INSPECT_ANNDATA}" >&2; exit 1; }
[[ -d "${RUN_DIR}/outs/filtered_feature_bc_matrix" ]] || { echo "ERROR: Missing ${RUN_DIR}/outs/filtered_feature_bc_matrix" >&2; exit 1; }
[[ -d "${RUN_DIR}/outs/raw_velocyto_feature_bc_matrix" ]] || { echo "ERROR: Missing ${RUN_DIR}/outs/raw_velocyto_feature_bc_matrix" >&2; exit 1; }
command -v docker >/dev/null 2>&1 || { echo "ERROR: docker is required" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 is required" >&2; exit 1; }

mkdir -p "${OUTPUT_DIR}" "${OUTPUT_DIR}/.numba" "${OUTPUT_DIR}/.matplotlib"

PYTHON_BACKEND_EFFECTIVE="${PYTHON_BACKEND}"
case "${PYTHON_BACKEND_EFFECTIVE}" in
  auto)
    if python3 - <<'PY' >/dev/null 2>&1
import anndata
import scanpy
import numpy
import pandas
import scipy
PY
    then
      PYTHON_BACKEND_EFFECTIVE="host"
    else
      PYTHON_BACKEND_EFFECTIVE="docker"
    fi
    ;;
  host|docker)
    ;;
  *)
    echo "ERROR: --python-backend must be auto, host, or docker" >&2
    exit 1
    ;;
esac

PYTHON_DOCKER_ARGS=(
  run --rm
  --user "$(id -u):$(id -g)"
  -v "${RUN_DIR}:${RUN_DIR}"
  -v "${OUTPUT_DIR}:${OUTPUT_DIR}"
  -v "${REPO_ROOT}:${REPO_ROOT}:ro"
  -e "NUMBA_CACHE_DIR=${OUTPUT_DIR}/.numba"
  -e "NUMBA_DISABLE_JIT=1"
  -e "MPLCONFIGDIR=${OUTPUT_DIR}/.matplotlib"
)
if [[ -d "${SC_RNA_SEQ_ROOT}" ]]; then
  PYTHON_DOCKER_ARGS+=(-v "${SC_RNA_SEQ_ROOT}:${SC_RNA_SEQ_ROOT}:ro")
fi

if [[ "${PYTHON_BACKEND_EFFECTIVE}" == "docker" ]]; then
  docker image inspect "${DOCKER_IMAGE}" >/dev/null 2>&1 || {
    echo "ERROR: Python helper backend selected Docker, but image is missing: ${DOCKER_IMAGE}" >&2
    exit 1
  }
fi

run_py() {
  if [[ "${PYTHON_BACKEND_EFFECTIVE}" == "docker" ]]; then
    docker "${PYTHON_DOCKER_ARGS[@]}" "${DOCKER_IMAGE}" python3 "$@"
  else
    python3 "$@"
  fi
}

run_py_stdin() {
  if [[ "${PYTHON_BACKEND_EFFECTIVE}" == "docker" ]]; then
    docker "${PYTHON_DOCKER_ARGS[@]}" -i "${DOCKER_IMAGE}" python3 - "$@"
  else
    python3 - "$@"
  fi
}

echo "=== Downstream GeneFull + Velocyto ==="
echo "Run dir: ${RUN_DIR}"
echo "Output dir: ${OUTPUT_DIR}"
echo "Docker image: ${DOCKER_IMAGE}"
echo "Feature gather image: ${FEATURE_GATHER_IMAGE}"
echo "Cell calls: STAR filtered barcodes"
echo "Python helper backend: ${PYTHON_BACKEND_EFFECTIVE}"
if [[ "${RUN_CELLBENDER}" == "1" ]]; then
  echo "CellBender image: ${CELLBENDER_IMAGE}"
fi

run_py "${BUILD_COUNTS}" \
  --run-dir "${RUN_DIR}" \
  --output-h5ad "${COUNTS_H5AD}"

DOCKER_ARGS=(
  run --rm
  --user "$(id -u):$(id -g)"
  -v "${OUTPUT_DIR}:${OUTPUT_DIR}"
  -v "${REPO_ROOT}:${REPO_ROOT}:ro"
  -e "min_genes=${MIN_GENES}"
  -e "max_genes=${MAX_GENES}"
  -e "mt_pct_cutoff=${MT_PCT_CUTOFF}"
  -e "NUMBA_CACHE_DIR=${OUTPUT_DIR}/.numba"
  -e "NUMBA_DISABLE_JIT=1"
  -e "MPLCONFIGDIR=${OUTPUT_DIR}/.matplotlib"
)

if [[ -n "${MITO_GENES}" ]]; then
  MITO_GENES="$(realpath "${MITO_GENES}")"
  [[ -f "${MITO_GENES}" ]] || { echo "ERROR: Missing mito genes file ${MITO_GENES}" >&2; exit 1; }
  DOCKER_ARGS+=(-v "${MITO_GENES}:${MITO_GENES}:ro")
fi

docker "${DOCKER_ARGS[@]}" \
  "${DOCKER_IMAGE}" \
  Rscript "${DOUBLET_SCRIPT}" "${COUNTS_H5AD}"

COMBINE_ARGS=(--input_file "${COUNTS_H5AD}")
if [[ -n "${MITO_GENES}" ]]; then
  COMBINE_ARGS+=(--mito_genes "${MITO_GENES}")
fi

EFFECTIVE_MAX_GENES="${MAX_GENES}"
RAW_ADAPTIVE_MAX_GENES=""
if [[ "${ADAPTIVE_FILTER}" == "1" ]]; then
  ADAPTIVE_QC_JSON="${OUTPUT_DIR}/adaptive_qc_threshold.json"
  run_py "${COMPUTE_ADAPTIVE_QC}" \
    --counts-h5ad "${COUNTS_H5AD}" \
    --non-empty-barcodes "${OUTPUT_DIR}/non_empty_barcodes.txt" \
    --doublet-barcodes "${OUTPUT_DIR}/doublet_barcodes.txt" \
    --min-genes "${MIN_GENES}" \
    --n-mad "${N_MAD}" \
    --output-json "${ADAPTIVE_QC_JSON}" >/dev/null

  EFFECTIVE_MAX_GENES="$(python3 - <<'PY' "${ADAPTIVE_QC_JSON}"
import json, sys
with open(sys.argv[1], "r", encoding="utf-8") as handle:
    data = json.load(handle)
print(data["effective_max_genes"])
PY
)"
  RAW_ADAPTIVE_MAX_GENES="$(python3 - <<'PY' "${ADAPTIVE_QC_JSON}"
import json, sys
with open(sys.argv[1], "r", encoding="utf-8") as handle:
    data = json.load(handle)
print(data["raw_adaptive_max_genes"])
PY
)"
  if [[ "${EFFECTIVE_MAX_GENES}" != "${RAW_ADAPTIVE_MAX_GENES}" ]]; then
    echo "WARNING: adaptive max_genes (${RAW_ADAPTIVE_MAX_GENES}) fell below min_genes (${MIN_GENES}); clamping applied max_genes to ${EFFECTIVE_MAX_GENES}"
  else
    echo "Adaptive max_genes: ${EFFECTIVE_MAX_GENES}"
  fi
fi

DOCKER_ARGS+=(-e "max_genes=${EFFECTIVE_MAX_GENES}")

docker "${DOCKER_ARGS[@]}" \
  "${DOCKER_IMAGE}" \
  python3 /usr/local/bin/combineFilters.py "${COMBINE_ARGS[@]}"

if [[ "${ADAPTIVE_FILTER}" == "1" ]]; then
  run_py "${APPLY_ADAPTIVE_MT}" \
    --input-h5ad "${UNFILTERED_H5AD}" \
    --threshold-json "${ADAPTIVE_QC_JSON}" \
    --mt-floor "${MT_PCT_CUTOFF}" \
    --n-mad "${N_MAD}"

  docker "${DOCKER_ARGS[@]}" \
    "${DOCKER_IMAGE}" \
    python3 "${GENERATE_QC_HISTOGRAM_MT}" \
      --input-h5ad "${UNFILTERED_H5AD}" \
      --output-dir "${OUTPUT_DIR}" \
      --threshold-json "${ADAPTIVE_QC_JSON}"
fi

run_py "${POSTPROCESS_FILTERS}" \
  --unfiltered-h5ad "${UNFILTERED_H5AD}" \
  --qc-output-h5ad "${FILTERED_H5AD}" \
  --default-singlet-output-h5ad "${DEFAULT_SINGLET_FILTERED_H5AD}"

if [[ "${RUN_CELLBENDER}" == "1" || "${REUSE_CELLBENDER}" == "1" ]]; then
  CELLBENDER_CB_FILE="${OUTPUT_DIR}/cellbender/cellbender_counts.h5"
  CELLBENDER_FAILURE_NOTE="${OUTPUT_DIR}/cellbender/CELLBENDER_FAILED.txt"
  mkdir -p "$(dirname "${CELLBENDER_CB_FILE}")"
  rm -f "${CELLBENDER_FAILURE_NOTE}"

  if [[ "${REUSE_CELLBENDER}" == "1" ]]; then
    echo "Reusing existing CellBender output: ${CELLBENDER_CB_FILE}"
    [[ -f "${CELLBENDER_CB_FILE}" ]] || { echo "ERROR: --reuse-cellbender requires existing ${CELLBENDER_CB_FILE}" >&2; exit 1; }
  else
    CELLBENDER_ARGS=(
      run --rm
      --user "$(id -u):$(id -g)"
      -v "${OUTPUT_DIR}:${OUTPUT_DIR}"
      -w "${OUTPUT_DIR}"
      -e "NUMBA_CACHE_DIR=${OUTPUT_DIR}/.numba"
      -e "MPLCONFIGDIR=${OUTPUT_DIR}/.matplotlib"
    )
    CELLBENDER_CMD=(
      cellbender
      remove-background
      --input "${UNFILTERED_H5AD}"
      --output "${CELLBENDER_CB_FILE}"
      --cpu-threads "${CELLBENDER_CPU_CORES}"
    )
    if [[ "${CELLBENDER_USE_GPU}" == "1" ]]; then
      CELLBENDER_ARGS+=(--gpus all)
    else
      :
    fi
    if [[ "${CELLBENDER_USE_GPU}" == "1" ]]; then
      CELLBENDER_CMD+=(--cuda)
    fi
    if [[ -n "${CELLBENDER_FLAGS}" ]]; then
      # shellcheck disable=SC2206
      EXTRA_CB_FLAGS=( ${CELLBENDER_FLAGS} )
      CELLBENDER_CMD+=("${EXTRA_CB_FLAGS[@]}")
    fi

    if ! docker "${CELLBENDER_ARGS[@]}" \
      "${CELLBENDER_IMAGE}" \
      "${CELLBENDER_CMD[@]}"
    then
      echo "WARNING: CellBender remove-background failed; downstream will continue without ${CELLBENDER_LAYER} layer if no output was produced." >&2
    fi
  fi

  if [[ -f "${CELLBENDER_CB_FILE}" ]]; then
    CELLBENDER_ADD_ARGS=(
      run --rm
      --user "$(id -u):$(id -g)"
      -v "${OUTPUT_DIR}:${OUTPUT_DIR}"
      -v "${REPO_ROOT}:${REPO_ROOT}:ro"
      -w "${OUTPUT_DIR}"
      -e "NUMBA_CACHE_DIR=${OUTPUT_DIR}/.numba"
      -e "MPLCONFIGDIR=${OUTPUT_DIR}/.matplotlib"
      "${CELLBENDER_IMAGE}"
      python
      "${ADD_CELLBENDER_LAYER}"
      --cellbender-h5 "${CELLBENDER_CB_FILE}"
      --layer-name "${CELLBENDER_LAYER}"
    )

    for target_h5ad in "${COUNTS_H5AD}" "${UNFILTERED_H5AD}"; do
      docker "${CELLBENDER_ADD_ARGS[@]}" --input-h5ad "${target_h5ad}" --output-h5ad "${target_h5ad}"
    done

    run_py "${PROPAGATE_LAYER}" \
      --source-h5ad "${UNFILTERED_H5AD}" \
      --target-h5ad "${FILTERED_H5AD}" \
      --output-h5ad "${FILTERED_H5AD}" \
      --layer-name "${CELLBENDER_LAYER}"
    run_py "${PROPAGATE_LAYER}" \
      --source-h5ad "${UNFILTERED_H5AD}" \
      --target-h5ad "${DEFAULT_SINGLET_FILTERED_H5AD}" \
      --output-h5ad "${DEFAULT_SINGLET_FILTERED_H5AD}" \
      --layer-name "${CELLBENDER_LAYER}"

    cp -f "${UNFILTERED_H5AD}" "${FINAL_H5AD}"
    PRIMARY_H5AD="${FILTERED_H5AD}"
  else
    {
      echo "CellBender did not produce ${CELLBENDER_CB_FILE}"
      echo "input_h5ad=${UNFILTERED_H5AD}"
      echo "counts_h5ad=${COUNTS_H5AD}"
      echo "fallback_h5ad=${FINAL_H5AD}"
      echo "reason=sparse_or_prefiltered_input_can_fail_prior_estimation"
    } > "${CELLBENDER_FAILURE_NOTE}"
    echo "WARNING: CellBender did not produce ${CELLBENDER_CB_FILE}" >&2
    echo "WARNING: Continuing without denoised layer; wrote ${CELLBENDER_FAILURE_NOTE}" >&2
    cp -f "${UNFILTERED_H5AD}" "${FINAL_H5AD}"
    PRIMARY_H5AD="${FILTERED_H5AD}"
  fi
fi

FEATURE_OUTPUT_ROOT="${OUTPUT_DIR}/feature_libraries"
mapfile -t FEATURE_LIBRARY_DIRS < <(find "${RUN_DIR}/cr_assign" -type f -name 'pf_library_provenance.tsv' -print 2>/dev/null | sed 's#/pf_library_provenance.tsv$##' | sort)
if (( ${#FEATURE_LIBRARY_DIRS[@]} > 0 )); then
  echo "Feature libraries: ${#FEATURE_LIBRARY_DIRS[@]}"
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
    echo "ERROR: Explicit per-library call-file mapping is required before automatic integration is safe." >&2
    exit 1
  fi

  COUNTS_TARGETS=("${COUNTS_H5AD}" "${UNFILTERED_H5AD}" "${FILTERED_H5AD}" "${DEFAULT_SINGLET_FILTERED_H5AD}")
  if [[ -f "${FINAL_H5AD}" ]]; then
    COUNTS_TARGETS+=("${FINAL_H5AD}")
  fi

  for feature_library in "${FEATURE_LIBRARY_DIRS[@]}"; do
    echo "Integrating feature library: ${feature_library}"
    FEATURE_GATHER_ARGS=(
      run --rm
      --user "$(id -u):$(id -g)"
      -v "${RUN_DIR}:${RUN_DIR}:ro"
      -v "${OUTPUT_DIR}:${OUTPUT_DIR}"
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
else
  echo "Feature libraries: none detected"
fi

rm -f \
  "${OUTPUT_DIR}/counts.summary.txt" \
  "${OUTPUT_DIR}/unfiltered_counts.summary.txt" \
  "${OUTPUT_DIR}/filtered_counts.summary.txt" \
  "${OUTPUT_DIR}/final_counts.summary.txt"
SUMMARY_H5AD="$(run_py_stdin "${PRIMARY_H5AD}" "${FINAL_H5AD}" "${UNFILTERED_H5AD}" "${COUNTS_H5AD}" "${FILTERED_H5AD}" "${DEFAULT_SINGLET_FILTERED_H5AD}" <<'PY'
import sys
from pathlib import Path
import anndata as ad

for candidate in sys.argv[1:]:
    path = Path(candidate)
    if not path.exists():
        continue
    try:
        adata = ad.read_h5ad(path, backed="r")
        n_obs = adata.n_obs
        adata.file.close()
    except Exception:
        continue
    if n_obs > 0:
        print(path)
        break
else:
    print(sys.argv[1])
PY
)"
run_py "${INSPECT_ANNDATA}" "${SUMMARY_H5AD}" > "${OUTPUT_DIR}/summary.txt"

echo "PASS: downstream GeneFull + Velocyto"
echo "counts.h5ad: ${COUNTS_H5AD}"
echo "unfiltered_counts.h5ad: ${UNFILTERED_H5AD}"
echo "filtered_counts.h5ad: ${FILTERED_H5AD}"
echo "default_singlet_filtered_counts.h5ad: ${DEFAULT_SINGLET_FILTERED_H5AD}"
if [[ "${RUN_CELLBENDER}" == "1" || "${REUSE_CELLBENDER}" == "1" ]]; then
  echo "final_counts.h5ad: ${FINAL_H5AD}"
fi
if (( ${#FEATURE_LIBRARY_DIRS[@]} > 0 )); then
  echo "feature_libraries/: ${FEATURE_OUTPUT_ROOT}"
fi
echo "summary.txt: ${OUTPUT_DIR}/summary.txt"
echo "summary_source_h5ad: ${SUMMARY_H5AD}"
