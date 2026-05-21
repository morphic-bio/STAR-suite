#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

EXTRACT_MEX="${REPO_ROOT}/scripts/extract_cr_feature_type_mex.py"
BUILD_MUDATA="${REPO_ROOT}/scripts/build_multiome_mudata.py"
PACKAGE_GENEFULL="${REPO_ROOT}/scripts/package_star_genefull_mex.py"
PREPARE_VELOCYTO="${REPO_ROOT}/scripts/prepare_velocyto_mex.py"
ALLOW_LEGACY_PREPARE_VELOCYTO="${ALLOW_LEGACY_PREPARE_VELOCYTO:-0}"
LOCAL_DOWNSTREAM="${REPO_ROOT}/scripts/run_scrna_downstream_gene_full_velocyto.sh"
REMOTE_DOWNSTREAM="${REPO_ROOT}/scripts/run_remote_scrna_downstream_rsync.sh"
REMOTE_CELLBENDER="${REPO_ROOT}/scripts/run_remote_cellbender_rsync.sh"

ARC_OUT=""
STAR_RUN_DIR=""
SAMPLE_DIR=""
DOWNSTREAM_DIR=""
OUT_DIR=""
REMOTE_HOST="${MULTIOME_REMOTE_HOST:-}"
REMOTE_ROOT="${MULTIOME_REMOTE_ROOT:-}"
REMOTE_OUTPUT_NAME="${MULTIOME_REMOTE_OUTPUT_NAME:-downstream_genefull_velocyto_cellbender}"
CELLBENDER_CPU_CORES="${MULTIOME_CELLBENDER_CPU_CORES:-}"
RUN_ARC_ONLY="1"
RUN_STAR_HYBRID="1"
FORCE="0"
ALLOW_LOCAL_DOWNSTREAM="0"
REMOTE_CELLBENDER_GPU="0"
NO_SYNC_IMAGES="0"
KEEP_REMOTE="0"

usage() {
  cat <<'EOF'
Usage:
  run_multiome_mudata_smoke.sh --arc-outs PATH [options]

Builds production-shaped PBMC/JAX-style multiome MuData smoke outputs. The ARC
path splits combined Cell Ranger MEX into RNA and ATAC peak MEX. When a STAR
run is provided, the wrapper verifies/regenerates packaged GeneFull and
Velocyto MEX, runs the GeneFull+Velocyto downstream h5ad path, offloads
CellBender remotely when remote parameters are supplied, then builds STAR-RNA +
ARC-ATAC h5mu outputs.

Required:
  --arc-outs PATH            Cell Ranger ARC outs directory with raw/filtered
                             feature_bc_matrix and per_barcode_metrics.csv

STAR RNA options:
  --star-run-dir PATH        STAR run directory containing outs/ and Solo.out/
  --sample-dir PATH          Sample directory containing run/; inferred when
                             --star-run-dir ends in /run
  --downstream-dir PATH      Existing downstream GeneFull+Velocyto h5ad dir
  --allow-local-downstream   Permit local downstream h5ad generation without
                             remote CellBender. Default requires remote args.

Remote CellBender/downstream options:
  --remote-host HOST         SSH target, e.g. 10.159.4.53
  --remote-root PATH         Remote local-disk staging root
  --remote-output-name NAME  Downstream output dir name under sample-dir
                             (default: downstream_genefull_velocyto_cellbender)
  --cellbender-cpu-cores N   CPU cores passed to remote CellBender
  --cellbender-gpu           Pass GPU mode through to remote downstream/CellBender
  --no-sync-images           Do not sync local Docker images to the remote host
  --keep-remote              Leave remote staging directory in place

Output/control:
  --out-dir PATH             Output directory (default:
                             tests/multiome_mudata_smoke_output_<timestamp>)
  --skip-arc-only            Do not build ARC-only h5mu outputs
  --skip-star-hybrid         Do not build STAR-RNA + ARC-ATAC h5mu outputs
  --force                    Regenerate split MEX and downstream h5mu outputs
  --help                     Show this help

Environment fallbacks:
  MULTIOME_REMOTE_HOST, MULTIOME_REMOTE_ROOT, MULTIOME_CELLBENDER_CPU_CORES
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --arc-outs)
      ARC_OUT="$2"
      shift 2
      ;;
    --star-run-dir)
      STAR_RUN_DIR="$2"
      shift 2
      ;;
    --sample-dir)
      SAMPLE_DIR="$2"
      shift 2
      ;;
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
    --remote-output-name)
      REMOTE_OUTPUT_NAME="$2"
      shift 2
      ;;
    --cellbender-cpu-cores)
      CELLBENDER_CPU_CORES="$2"
      shift 2
      ;;
    --cellbender-gpu)
      REMOTE_CELLBENDER_GPU="1"
      shift
      ;;
    --no-sync-images)
      NO_SYNC_IMAGES="1"
      shift
      ;;
    --keep-remote)
      KEEP_REMOTE="1"
      shift
      ;;
    --out-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --skip-arc-only)
      RUN_ARC_ONLY="0"
      shift
      ;;
    --skip-star-hybrid)
      RUN_STAR_HYBRID="0"
      shift
      ;;
    --force)
      FORCE="1"
      shift
      ;;
    --allow-local-downstream)
      ALLOW_LOCAL_DOWNSTREAM="1"
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

[[ -n "${ARC_OUT}" ]] || { echo "ERROR: --arc-outs is required" >&2; exit 1; }
ARC_OUT="$(realpath "${ARC_OUT}")"
[[ -d "${ARC_OUT}/raw_feature_bc_matrix" ]] || { echo "ERROR: missing ${ARC_OUT}/raw_feature_bc_matrix" >&2; exit 1; }
[[ -d "${ARC_OUT}/filtered_feature_bc_matrix" ]] || { echo "ERROR: missing ${ARC_OUT}/filtered_feature_bc_matrix" >&2; exit 1; }
[[ -f "${ARC_OUT}/per_barcode_metrics.csv" ]] || { echo "ERROR: missing ${ARC_OUT}/per_barcode_metrics.csv" >&2; exit 1; }

for helper in "${EXTRACT_MEX}" "${BUILD_MUDATA}" "${PACKAGE_GENEFULL}" "${LOCAL_DOWNSTREAM}"; do
  [[ -f "${helper}" ]] || { echo "ERROR: missing helper ${helper}" >&2; exit 1; }
done
if [[ "${ALLOW_LEGACY_PREPARE_VELOCYTO}" == "1" ]]; then
  [[ -f "${PREPARE_VELOCYTO}" ]] || { echo "ERROR: missing helper ${PREPARE_VELOCYTO}" >&2; exit 1; }
fi

OUT_DIR="${OUT_DIR:-${REPO_ROOT}/tests/multiome_mudata_smoke_output_$(date -u +%Y%m%dT%H%M%SZ)}"
OUT_DIR="$(realpath -m "${OUT_DIR}")"
SPLIT_DIR="${OUT_DIR}/split_mex"
mkdir -p "${SPLIT_DIR}" "${OUT_DIR}/logs"

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"
}

run_extract() {
  local input_mex="$1"
  local feature_type="$2"
  local out_dir="$3"
  if [[ "${FORCE}" == "1" || ! -f "${out_dir}/matrix.mtx.gz" || ! -f "${out_dir}/features.tsv.gz" || ! -f "${out_dir}/barcodes.tsv.gz" ]]; then
    rm -rf "${out_dir}"
    log "Extracting ${feature_type} from ${input_mex} -> ${out_dir}"
    python3 "${EXTRACT_MEX}" \
      --input-mex-dir "${input_mex}" \
      --feature-type "${feature_type}" \
      --out-dir "${out_dir}" > "${OUT_DIR}/logs/$(basename "${out_dir}").log"
  else
    log "Reusing split MEX ${out_dir}"
  fi
}

ensure_arc_splits() {
  RAW_GEX_MEX="${SPLIT_DIR}/raw_gex_mex"
  FILTERED_GEX_MEX="${SPLIT_DIR}/filtered_gex_mex"
  RAW_ATAC_MEX="${SPLIT_DIR}/raw_atac_peak_mex"
  FILTERED_ATAC_MEX="${SPLIT_DIR}/filtered_atac_peak_mex"

  run_extract "${ARC_OUT}/raw_feature_bc_matrix" "Gene Expression" "${RAW_GEX_MEX}"
  run_extract "${ARC_OUT}/filtered_feature_bc_matrix" "Gene Expression" "${FILTERED_GEX_MEX}"
  run_extract "${ARC_OUT}/raw_feature_bc_matrix" "Peaks" "${RAW_ATAC_MEX}"
  run_extract "${ARC_OUT}/filtered_feature_bc_matrix" "Peaks" "${FILTERED_ATAC_MEX}"
}

genefull_packaging_ok() {
  local run_dir="$1"
  local outs="${run_dir}/outs"
  for required in \
    "${outs}/raw_feature_bc_matrix/matrix.mtx.gz" \
    "${outs}/raw_feature_bc_matrix/features.tsv.gz" \
    "${outs}/raw_feature_bc_matrix/barcodes.tsv.gz" \
    "${outs}/filtered_feature_bc_matrix/matrix.mtx.gz" \
    "${outs}/filtered_feature_bc_matrix/features.tsv.gz" \
    "${outs}/filtered_feature_bc_matrix/barcodes.tsv.gz" \
    "${outs}/gene_full_feature_bc_matrix_manifest.json"
  do
    [[ -f "${required}" ]] || return 1
  done
  return 0
}

ensure_genefull_packaging() {
  local run_dir="$1"
  if genefull_packaging_ok "${run_dir}"; then
    log "GeneFull MEX packaging is present"
    return 0
  fi

  log "GeneFull MEX packaging missing or incomplete; regenerating from Solo.out/GeneFull"
  python3 "${PACKAGE_GENEFULL}" --run-dir "${run_dir}" > "${OUT_DIR}/logs/package_star_genefull_mex.log"
  if ! genefull_packaging_ok "${run_dir}"; then
    echo "ERROR: GeneFull MEX packaging is still incomplete after regeneration" >&2
    exit 1
  fi
}

velocyto_packaging_ok() {
  local run_dir="$1"
  local outs="${run_dir}/outs"
  for required in \
    "${outs}/raw_velocyto_feature_bc_matrix/matrix.mtx.gz" \
    "${outs}/raw_velocyto_feature_bc_matrix/spliced.mtx.gz" \
    "${outs}/raw_velocyto_feature_bc_matrix/unspliced.mtx.gz" \
    "${outs}/raw_velocyto_feature_bc_matrix/ambiguous.mtx.gz" \
    "${outs}/filtered_velocyto_feature_bc_matrix/matrix.mtx.gz" \
    "${outs}/velocyto_feature_bc_matrix_manifest.json"
  do
    [[ -f "${required}" ]] || return 1
  done
  return 0
}

ensure_velocyto_packaging() {
  local run_dir="$1"
  if velocyto_packaging_ok "${run_dir}"; then
    log "Native Velocyto MEX packaging is present"
    return 0
  fi

  if [[ "${ALLOW_LEGACY_PREPARE_VELOCYTO}" == "1" ]]; then
    log "WARNING: native Velocyto MEX missing; using legacy prepare_velocyto_mex.py fallback"
    python3 "${PREPARE_VELOCYTO}" --run-dir "${run_dir}" > "${OUT_DIR}/logs/prepare_velocyto_mex.log"
  else
    echo "ERROR: Native Velocyto MEX packaging is missing under ${run_dir}/outs" >&2
    exit 1
  fi
  if ! velocyto_packaging_ok "${run_dir}"; then
    echo "ERROR: Velocyto MEX packaging is still incomplete after regeneration" >&2
    exit 1
  fi
}

infer_sample_dir() {
  local run_dir="$1"
  if [[ -n "${SAMPLE_DIR}" ]]; then
    realpath "${SAMPLE_DIR}"
    return 0
  fi
  if [[ "$(basename "${run_dir}")" == "run" ]]; then
    dirname "${run_dir}"
    return 0
  fi
  return 1
}

run_or_reuse_downstream() {
  local run_dir="$1"
  local sample_dir=""

  if [[ -n "${DOWNSTREAM_DIR}" ]]; then
    DOWNSTREAM_DIR="$(realpath "${DOWNSTREAM_DIR}")"
    [[ -d "${DOWNSTREAM_DIR}" ]] || { echo "ERROR: missing downstream dir ${DOWNSTREAM_DIR}" >&2; exit 1; }
  else
    if sample_dir="$(infer_sample_dir "${run_dir}")"; then
      :
    else
      echo "ERROR: --sample-dir is required when --star-run-dir is not named run" >&2
      exit 1
    fi

    if [[ -n "${REMOTE_HOST}" && -n "${REMOTE_ROOT}" ]]; then
      [[ -x "${REMOTE_DOWNSTREAM}" ]] || { echo "ERROR: missing remote downstream helper ${REMOTE_DOWNSTREAM}" >&2; exit 1; }
      local remote_args=(
        "${REMOTE_DOWNSTREAM}"
        --sample-dir "${sample_dir}"
        --remote-host "${REMOTE_HOST}"
        --remote-root "${REMOTE_ROOT}"
        --output-name "${REMOTE_OUTPUT_NAME}"
        --run-cellbender
        --adaptive-filter
      )
      if [[ -n "${CELLBENDER_CPU_CORES}" ]]; then
        remote_args+=(--cellbender-cpu-cores "${CELLBENDER_CPU_CORES}")
      fi
      if [[ "${NO_SYNC_IMAGES}" == "1" ]]; then
        remote_args+=(--no-sync-images)
      fi
      if [[ "${KEEP_REMOTE}" == "1" ]]; then
        remote_args+=(--keep-remote)
      fi
      if [[ "${REMOTE_CELLBENDER_GPU}" == "1" ]]; then
        remote_args+=(--cellbender-gpu)
      fi
      log "Running remote GeneFull+Velocyto downstream with CellBender"
      "${remote_args[@]}" | tee "${OUT_DIR}/logs/remote_downstream.log"
      DOWNSTREAM_DIR="${sample_dir}/${REMOTE_OUTPUT_NAME}"
    elif [[ "${ALLOW_LOCAL_DOWNSTREAM}" == "1" ]]; then
      DOWNSTREAM_DIR="${sample_dir}/downstream_genefull_velocyto"
      log "Running local GeneFull+Velocyto downstream without CellBender"
      "${LOCAL_DOWNSTREAM}" \
        --run-dir "${run_dir}" \
        --output-dir "${DOWNSTREAM_DIR}" \
        --adaptive-filter | tee "${OUT_DIR}/logs/local_downstream.log"
    else
      echo "ERROR: provide --remote-host and --remote-root for remote CellBender, or pass --downstream-dir for completed downstream outputs" >&2
      exit 1
    fi
  fi

  if [[ ! -f "${DOWNSTREAM_DIR}/cellbender/cellbender_counts.h5" && -f "${DOWNSTREAM_DIR}/cellbender/CELLBENDER_FAILED.txt" ]]; then
    log "WARNING: CellBender failure note exists; keeping downstream fallback"
  elif [[ ! -f "${DOWNSTREAM_DIR}/cellbender/cellbender_counts.h5" && -n "${REMOTE_HOST}" && -n "${REMOTE_ROOT}" ]]; then
    [[ -x "${REMOTE_CELLBENDER}" ]] || { echo "ERROR: missing remote CellBender helper ${REMOTE_CELLBENDER}" >&2; exit 1; }
    local cb_args=(
      "${REMOTE_CELLBENDER}"
      --downstream-dir "${DOWNSTREAM_DIR}"
      --remote-host "${REMOTE_HOST}"
      --remote-root "${REMOTE_ROOT}"
    )
    if [[ -n "${CELLBENDER_CPU_CORES}" ]]; then
      cb_args+=(--cellbender-cpu-cores "${CELLBENDER_CPU_CORES}")
    fi
    if [[ "${REMOTE_CELLBENDER_GPU}" == "1" ]]; then
      cb_args+=(--cellbender-gpu)
    fi
    if [[ "${NO_SYNC_IMAGES}" == "1" ]]; then
      cb_args+=(--no-sync-image)
    fi
    if [[ "${KEEP_REMOTE}" == "1" ]]; then
      cb_args+=(--keep-remote)
    fi
    log "Running remote CellBender on existing downstream directory"
    "${cb_args[@]}" | tee "${OUT_DIR}/logs/remote_cellbender.log"
  elif [[ ! -f "${DOWNSTREAM_DIR}/cellbender/cellbender_counts.h5" && "${ALLOW_LOCAL_DOWNSTREAM}" != "1" ]]; then
    echo "ERROR: downstream output lacks cellbender/cellbender_counts.h5; provide --remote-host/--remote-root for remote CellBender or pass --allow-local-downstream to proceed without it" >&2
    exit 1
  elif [[ ! -f "${DOWNSTREAM_DIR}/cellbender/cellbender_counts.h5" ]]; then
    log "WARNING: proceeding without CellBender output because --allow-local-downstream was set"
  fi
}

choose_unfiltered_rna() {
  local downstream="$1"
  if [[ -f "${downstream}/final_counts.h5ad" ]]; then
    echo "${downstream}/final_counts.h5ad"
  elif [[ -f "${downstream}/unfiltered_counts.h5ad" ]]; then
    echo "${downstream}/unfiltered_counts.h5ad"
  else
    echo "${downstream}/counts.h5ad"
  fi
}

build_h5mu() {
  local output="$1"
  shift
  if [[ "${FORCE}" == "1" || ! -f "${output}" ]]; then
    log "Building ${output}"
    python3 "${BUILD_MUDATA}" "$@" --output-h5mu "${output}" | tee "${output%.h5mu}.build.log"
  else
    log "Reusing ${output}"
  fi
}

write_manifest() {
  local manifest="${OUT_DIR}/RUN_MANIFEST.txt"
  local y_removal_mode="unknown"
  if [[ -n "${STAR_RUN_DIR}" ]]; then
    y_removal_mode="true"
  fi
  {
    printf 'date_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'arc_out=%s\n' "${ARC_OUT}"
    printf 'star_run_dir=%s\n' "${STAR_RUN_DIR}"
    printf 'downstream_dir=%s\n' "${DOWNSTREAM_DIR}"
    printf 'remote_host=%s\n' "${REMOTE_HOST}"
    printf 'remote_root=%s\n' "${REMOTE_ROOT}"
    printf 'y_removal_enabled=%s\n' "${y_removal_mode}"
    printf 'gex_genome_dir=%s\n' "/storage/autoindex_110_44/bulk_index"
  } > "${manifest}"
  log "Wrote ${manifest}"
}

ensure_arc_splits

if [[ "${RUN_ARC_ONLY}" == "1" ]]; then
  build_h5mu "${OUT_DIR}/arc_unfiltered_multiome.h5mu" \
    --rna-mex-dir "${RAW_GEX_MEX}" \
    --atac-mex-dir "${RAW_ATAC_MEX}" \
    --per-barcode-metrics "${ARC_OUT}/per_barcode_metrics.csv" \
    --cell-call-source arc_per_barcode_metrics \
    --rna-source "${ARC_OUT}/raw_feature_bc_matrix" \
    --atac-source "${ARC_OUT}/raw_feature_bc_matrix" \
    --y-removal-enabled unknown

  build_h5mu "${OUT_DIR}/arc_filtered_multiome.h5mu" \
    --rna-mex-dir "${FILTERED_GEX_MEX}" \
    --atac-mex-dir "${FILTERED_ATAC_MEX}" \
    --per-barcode-metrics "${ARC_OUT}/per_barcode_metrics.csv" \
    --all-barcodes-are-cells \
    --cell-call-source arc_filtered_feature_bc_matrix \
    --rna-source "${ARC_OUT}/filtered_feature_bc_matrix" \
    --atac-source "${ARC_OUT}/filtered_feature_bc_matrix" \
    --y-removal-enabled unknown
fi

if [[ -n "${STAR_RUN_DIR}" && "${RUN_STAR_HYBRID}" == "1" ]]; then
  STAR_RUN_DIR="$(realpath "${STAR_RUN_DIR}")"
  [[ -d "${STAR_RUN_DIR}" ]] || { echo "ERROR: missing STAR run dir ${STAR_RUN_DIR}" >&2; exit 1; }
  ensure_genefull_packaging "${STAR_RUN_DIR}"
  ensure_velocyto_packaging "${STAR_RUN_DIR}"
  run_or_reuse_downstream "${STAR_RUN_DIR}"

  RNA_UNFILTERED_H5AD="$(choose_unfiltered_rna "${DOWNSTREAM_DIR}")"
  RNA_FILTERED_H5AD="${DOWNSTREAM_DIR}/filtered_counts.h5ad"
  [[ -f "${RNA_UNFILTERED_H5AD}" ]] || { echo "ERROR: missing RNA h5ad ${RNA_UNFILTERED_H5AD}" >&2; exit 1; }
  [[ -f "${RNA_FILTERED_H5AD}" ]] || { echo "ERROR: missing filtered RNA h5ad ${RNA_FILTERED_H5AD}" >&2; exit 1; }

  build_h5mu "${OUT_DIR}/star_arc_unfiltered_multiome.h5mu" \
    --rna-h5ad "${RNA_UNFILTERED_H5AD}" \
    --atac-mex-dir "${RAW_ATAC_MEX}" \
    --per-barcode-metrics "${ARC_OUT}/per_barcode_metrics.csv" \
    --require-rna-velocyto-layers \
    --cell-call-source star_downstream_h5ad \
    --rna-source "${RNA_UNFILTERED_H5AD}" \
    --atac-source "${ARC_OUT}/raw_feature_bc_matrix" \
    --strip-barcode-suffix \
    --y-removal-enabled true

  build_h5mu "${OUT_DIR}/star_arc_filtered_multiome.h5mu" \
    --rna-h5ad "${RNA_FILTERED_H5AD}" \
    --atac-mex-dir "${RAW_ATAC_MEX}" \
    --per-barcode-metrics "${ARC_OUT}/per_barcode_metrics.csv" \
    --all-barcodes-are-cells \
    --require-rna-velocyto-layers \
    --cell-call-source star_downstream_filtered_h5ad \
    --rna-source "${RNA_FILTERED_H5AD}" \
    --atac-source "${ARC_OUT}/raw_feature_bc_matrix" \
    --strip-barcode-suffix \
    --y-removal-enabled true
fi

write_manifest

echo "PASS: multiome MuData smoke pipeline complete"
echo "Output dir: ${OUT_DIR}"
