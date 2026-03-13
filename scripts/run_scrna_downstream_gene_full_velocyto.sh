#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

SC_RNA_SEQ_ROOT="${SC_RNA_SEQ_ROOT:-/mnt/pikachu/scRNA-seq}"
BUILD_COUNTS="${REPO_ROOT}/scripts/build_gene_full_velocyto_h5ad.py"
DOUBLET_SCRIPT="${REPO_ROOT}/scripts/run_star_cell_doublets.R"
FEATURE_GATHER_SCRIPT="${REPO_ROOT}/scripts/integrate_feature_library.py"
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
ADAPTIVE_FILTER="0"
N_MAD="${N_MAD:-3}"
RUN_CELLBENDER="0"
CELLBENDER_USE_GPU="0"
CELLBENDER_CPU_CORES="${CELLBENDER_CPU_CORES:-8}"
CELLBENDER_LAYER="${CELLBENDER_LAYER:-denoised}"
CELLBENDER_FLAGS="${CELLBENDER_FLAGS:-}"

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
  --mt-pct-cutoff FLOAT  MT percentage cutoff (default: 5)
  --adaptive-filter      Use adaptive max_genes threshold
  --n-mad FLOAT          Number of MADs for adaptive filtering (default: 3)
  --docker-image IMAGE   Downstream container image (default: biodepot/scrna-matrices:latest)
  --feature-gather-image IMAGE
                         Feature-library integration image
                         (default: biodepot/gather_features:latest)
  --run-cellbender       Run CellBender on filtered_counts.h5ad and add denoised layer
  --cellbender-image IMG CellBender image (default: biodepot/cellbender:0.3.2)
  --cellbender-gpu       Run CellBender with --gpus all instead of CPU mode
  --cellbender-cpu-cores INT
                         CPU cores passed to CellBender (default: 8)
  --cellbender-layer NAME
                         Layer name for CellBender output (default: denoised)
  --cellbender-flags STR Additional flags passed to cellbender remove-background
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
FINAL_H5AD="${OUTPUT_DIR}/final_counts.h5ad"
PRIMARY_H5AD="${FILTERED_H5AD}"

[[ -f "${BUILD_COUNTS}" ]] || { echo "ERROR: Missing helper ${BUILD_COUNTS}" >&2; exit 1; }
[[ -f "${DOUBLET_SCRIPT}" ]] || { echo "ERROR: Missing helper ${DOUBLET_SCRIPT}" >&2; exit 1; }
[[ -f "${FEATURE_GATHER_SCRIPT}" ]] || { echo "ERROR: Missing helper ${FEATURE_GATHER_SCRIPT}" >&2; exit 1; }
[[ -f "${INSPECT_ANNDATA}" ]] || { echo "ERROR: Missing helper ${INSPECT_ANNDATA}" >&2; exit 1; }
[[ -d "${RUN_DIR}/outs/filtered_feature_bc_matrix" ]] || { echo "ERROR: Missing ${RUN_DIR}/outs/filtered_feature_bc_matrix" >&2; exit 1; }
[[ -d "${RUN_DIR}/outs/raw_velocyto_feature_bc_matrix" ]] || { echo "ERROR: Missing ${RUN_DIR}/outs/raw_velocyto_feature_bc_matrix" >&2; exit 1; }
command -v docker >/dev/null 2>&1 || { echo "ERROR: docker is required" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 is required" >&2; exit 1; }

mkdir -p "${OUTPUT_DIR}" "${OUTPUT_DIR}/.numba" "${OUTPUT_DIR}/.matplotlib"

echo "=== Downstream GeneFull + Velocyto ==="
echo "Run dir: ${RUN_DIR}"
echo "Output dir: ${OUTPUT_DIR}"
echo "Docker image: ${DOCKER_IMAGE}"
echo "Feature gather image: ${FEATURE_GATHER_IMAGE}"
echo "Cell calls: STAR filtered barcodes"
if [[ "${RUN_CELLBENDER}" == "1" ]]; then
  echo "CellBender image: ${CELLBENDER_IMAGE}"
fi

python3 "${BUILD_COUNTS}" \
  --run-dir "${RUN_DIR}" \
  --feature-raw-dir "${RUN_DIR}/outs/filtered_feature_bc_matrix" \
  --feature-filtered-dir "${RUN_DIR}/outs/filtered_feature_bc_matrix" \
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
if [[ "${ADAPTIVE_FILTER}" == "1" ]]; then
  COMBINE_ARGS+=(--adaptive_filter --n_mad "${N_MAD}")
fi

docker "${DOCKER_ARGS[@]}" \
  "${DOCKER_IMAGE}" \
  python3 /usr/local/bin/combineFilters.py "${COMBINE_ARGS[@]}"

if [[ "${RUN_CELLBENDER}" == "1" ]]; then
  CELLBENDER_ARGS=(
    run --rm
    --user "$(id -u):$(id -g)"
    -v "${OUTPUT_DIR}:${OUTPUT_DIR}"
    -e "alignsDir=${OUTPUT_DIR}"
    -e "input_pattern=filtered_counts.h5ad"
    -e "output_pattern=$(basename "${FINAL_H5AD}")"
    -e "cb_subdir=cellbender"
    -e "cb_file=cellbender_counts.h5"
    -e "layername=${CELLBENDER_LAYER}"
    -e "nThreads=1"
    -e "cpu_cores=${CELLBENDER_CPU_CORES}"
  )
  if [[ "${CELLBENDER_USE_GPU}" == "1" ]]; then
    CELLBENDER_ARGS+=(--gpus all)
  else
    CELLBENDER_ARGS+=(-e "usecpu=1")
  fi
  if [[ -n "${CELLBENDER_FLAGS}" ]]; then
    CELLBENDER_ARGS+=(-e "additional_flags=${CELLBENDER_FLAGS}")
  fi

  docker "${CELLBENDER_ARGS[@]}" \
    "${CELLBENDER_IMAGE}" \
    remove_noise.sh

  if [[ -f "${OUTPUT_DIR}/cellbender/$(basename "${FINAL_H5AD}")" ]]; then
    cp -f "${OUTPUT_DIR}/cellbender/$(basename "${FINAL_H5AD}")" "${FINAL_H5AD}"
    PRIMARY_H5AD="${FINAL_H5AD}"
    python3 "${INSPECT_ANNDATA}" "${PRIMARY_H5AD}" > "${OUTPUT_DIR}/summary.txt"
  else
    echo "ERROR: CellBender did not produce ${OUTPUT_DIR}/cellbender/$(basename "${FINAL_H5AD}")" >&2
    echo "ERROR: For this 100K smoke, filtered_counts.h5ad is too prefiltered/sparse for CellBender." >&2
    echo "ERROR: CellBender needs a raw unfiltered barcode matrix to infer priors." >&2
    exit 1
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

  COUNTS_TARGETS=("${COUNTS_H5AD}" "${UNFILTERED_H5AD}" "${FILTERED_H5AD}")
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
python3 "${INSPECT_ANNDATA}" "${PRIMARY_H5AD}" > "${OUTPUT_DIR}/summary.txt"

echo "PASS: downstream GeneFull + Velocyto"
echo "counts.h5ad: ${COUNTS_H5AD}"
echo "unfiltered_counts.h5ad: ${UNFILTERED_H5AD}"
echo "filtered_counts.h5ad: ${FILTERED_H5AD}"
if [[ "${RUN_CELLBENDER}" == "1" ]]; then
  echo "final_counts.h5ad: ${FINAL_H5AD}"
fi
if (( ${#FEATURE_LIBRARY_DIRS[@]} > 0 )); then
  echo "feature_libraries/: ${FEATURE_OUTPUT_ROOT}"
fi
echo "summary.txt: ${OUTPUT_DIR}/summary.txt"
