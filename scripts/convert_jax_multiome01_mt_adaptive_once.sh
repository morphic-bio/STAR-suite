#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

RUN_ROOT="${RUN_ROOT:-/mnt/pikachu/JAX_Multiome01_processed/star_multiome_prod_globus_20260517T183219Z}"
OUTPUT_NAME="${OUTPUT_NAME:-downstream_genefull_velocyto_cellbender}"
MT_FLOOR="${MT_FLOOR:-5}"
N_MAD="${N_MAD:-3}"
SAMPLES="${SAMPLES:-}"
FORCE="${FORCE:-0}"
MUDATA_PYTHON="${MUDATA_PYTHON:-}"

usage() {
  cat <<'EOF'
Usage:
  convert_jax_multiome01_mt_adaptive_once.sh [options]

One-time retrofit for JAX_Multiome01 production outputs. For each completed
sample, applies the adaptive mitochondrial-percentage QC guard to downstream RNA
h5ads, regenerates filtered h5ads/QC histograms, and rebuilds filtered and
unfiltered RNA+ATAC MuData outputs.

Options:
  --run-root PATH       Production root (default: JAX_Multiome01 20260517 run)
  --output-name NAME    Downstream output directory name
  --samples CSV         Restrict to comma-separated sample directory names
  --mt-floor FLOAT      MT percentage floor (default: 5)
  --n-mad FLOAT         MAD multiplier (default: 3)
  --mudata-python PATH  Python with mudata installed
  --force              Re-run even if a sample has a conversion marker
  --help

Environment variables with the same uppercase names are also honored.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-root) RUN_ROOT="$2"; shift 2 ;;
    --output-name) OUTPUT_NAME="$2"; shift 2 ;;
    --samples) SAMPLES="$2"; shift 2 ;;
    --mt-floor) MT_FLOOR="$2"; shift 2 ;;
    --n-mad) N_MAD="$2"; shift 2 ;;
    --mudata-python) MUDATA_PYTHON="$2"; shift 2 ;;
    --force) FORCE="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument $1" >&2; usage >&2; exit 1 ;;
  esac
done

RUN_ROOT="$(realpath "${RUN_ROOT}")"
SAMPLES_ROOT="${RUN_ROOT}/samples"
LOG_DIR="${RUN_ROOT}/logs/mt_adaptive_conversion"
CONVERT_MARKER_NAME="MT_ADAPTIVE_CONVERTED.txt"

APPLY_ADAPTIVE_MT="${REPO_ROOT}/scripts/apply_adaptive_mt_filter.py"
POSTPROCESS_FILTERS="${REPO_ROOT}/scripts/postprocess_downstream_filters.py"
GENERATE_QC_HISTOGRAM_MT="${REPO_ROOT}/scripts/generate_qc_histogram_mt_adaptive.py"
BUILD_MUDATA="${REPO_ROOT}/scripts/build_multiome_mudata.py"

for required in \
  "${APPLY_ADAPTIVE_MT}" \
  "${POSTPROCESS_FILTERS}" \
  "${GENERATE_QC_HISTOGRAM_MT}" \
  "${BUILD_MUDATA}"
do
  [[ -f "${required}" ]] || { echo "ERROR: missing helper ${required}" >&2; exit 1; }
done
[[ -d "${SAMPLES_ROOT}" ]] || { echo "ERROR: missing samples root ${SAMPLES_ROOT}" >&2; exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 is required" >&2; exit 1; }
mkdir -p "${LOG_DIR}"

python_has_mudata() {
  "$1" - <<'PY' >/dev/null 2>&1
import mudata
PY
}

if [[ -z "${MUDATA_PYTHON}" ]]; then
  if python_has_mudata python3; then
    MUDATA_PYTHON="python3"
  else
    for candidate in \
      "${REPO_ROOT}/tests/jax_multiome_lane_smoke_native_20260517T073556Z/mudata_venv/bin/python" \
      "${REPO_ROOT}/tests/jax_multiome_lane_smoke_20260517T052512Z/mudata_venv/bin/python" \
      "${REPO_ROOT}/tests/multiome_mudata_smoke_output_star_end2end_20260516T045156Z/mudata_venv/bin/python"
    do
      if [[ -x "${candidate}" ]] && python_has_mudata "${candidate}"; then
        MUDATA_PYTHON="${candidate}"
        break
      fi
    done
  fi
fi
if [[ -z "${MUDATA_PYTHON}" ]]; then
  VENV="${LOG_DIR}/mudata_venv"
  python3 -m venv --system-site-packages "${VENV}"
  "${VENV}/bin/python" -m pip install --quiet mudata
  MUDATA_PYTHON="${VENV}/bin/python"
fi
python_has_mudata "${MUDATA_PYTHON}" || {
  echo "ERROR: selected Python cannot import mudata: ${MUDATA_PYTHON}" >&2
  exit 1
}

RUN_LOG="${LOG_DIR}/conversion_$(date -u +%Y%m%dT%H%M%SZ).log"

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "${RUN_LOG}"
}

sample_selected() {
  local sample="$1"
  [[ -z "${SAMPLES}" ]] && return 0
  local item
  IFS=',' read -r -a wanted <<< "${SAMPLES}"
  for item in "${wanted[@]}"; do
    [[ "${item}" == "${sample}" ]] && return 0
  done
  return 1
}

backup_once() {
  local src="$1"
  [[ -e "${src}" ]] || return 0
  local dst="${src}.pre_mt_adaptive"
  if [[ -e "${dst}" ]]; then
    return 0
  fi
  cp -an --reflink=auto "${src}" "${dst}"
}

summarize_threshold() {
  local sample="$1"
  local json_path="$2"
  python3 - "${sample}" "${json_path}" <<'PY' | tee -a "${RUN_LOG}"
import json
import sys
sample, path = sys.argv[1:3]
with open(path, "r", encoding="utf-8") as handle:
    data = json.load(handle)
keys = [
    "mt_pct_median",
    "mt_pct_mad",
    "mt_pct_raw_threshold",
    "mt_pct_threshold",
    "mt_pct_threshold_was_floored",
    "filter_cells_strict_mt5",
    "filter_cells_mt_adaptive",
    "singlet_filtered_cells_strict_mt5",
    "singlet_filtered_cells_mt_adaptive",
    "mt_pct_flag",
]
print(sample + "\t" + "\t".join(f"{key}={data.get(key)}" for key in keys))
PY
}

mapfile -t sample_dirs < <(find "${SAMPLES_ROOT}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)

log "JAX_Multiome01 adaptive MT conversion"
log "Run root: ${RUN_ROOT}"
log "Output name: ${OUTPUT_NAME}"
log "MT floor: ${MT_FLOOR}; n_mad: ${N_MAD}"
log "MuData Python: ${MUDATA_PYTHON}"
log "Log: ${RUN_LOG}"

converted=0
skipped=0
for sample in "${sample_dirs[@]}"; do
  sample_dir="${SAMPLES_ROOT}/${sample}"
  sample_selected "${sample}" || continue

  downstream="${sample_dir}/star_sample/${OUTPUT_NAME}"
  mudata_dir="${sample_dir}/mudata"
  marker="${sample_dir}/${CONVERT_MARKER_NAME}"
  if [[ "${FORCE}" != "1" && -f "${marker}" ]]; then
    log "Skipping ${sample}; marker exists (${marker})"
    skipped=$((skipped + 1))
    continue
  fi

  threshold_json="${downstream}/adaptive_qc_threshold.json"
  unfiltered_h5ad="${downstream}/unfiltered_counts.h5ad"
  final_h5ad="${downstream}/final_counts.h5ad"
  filtered_h5ad="${downstream}/filtered_counts.h5ad"
  default_singlet_h5ad="${downstream}/default_singlet_filtered_counts.h5ad"
  atac_mex="${sample_dir}/atac/peak_mex"
  atac_metrics="${sample_dir}/atac/atac_metrics.tsv"
  atac_sidecar="${sample_dir}/star_sample/run/atac_fragments.bin"
  atac_peaks="${sample_dir}/star_sample/run/atac_peaks.narrowPeak"
  unfiltered_h5mu="${mudata_dir}/star_chromap_unfiltered_multiome.h5mu"
  filtered_h5mu="${mudata_dir}/star_chromap_filtered_multiome.h5mu"

  for required in \
    "${threshold_json}" \
    "${unfiltered_h5ad}" \
    "${final_h5ad}" \
    "${atac_mex}/matrix.mtx.gz" \
    "${atac_mex}/barcodes.tsv.gz" \
    "${atac_mex}/features.tsv.gz" \
    "${atac_metrics}"
  do
    [[ -e "${required}" ]] || { echo "ERROR: ${sample}: missing ${required}" >&2; exit 1; }
  done

  log "Converting ${sample}"
  mkdir -p "${mudata_dir}"
  for path in \
    "${threshold_json}" \
    "${downstream}/gene_quantile_histogram.html" \
    "${downstream}/gene_quantile_histogram.png" \
    "${unfiltered_h5ad}" \
    "${final_h5ad}" \
    "${filtered_h5ad}" \
    "${default_singlet_h5ad}" \
    "${unfiltered_h5mu}" \
    "${filtered_h5mu}"
  do
    backup_once "${path}"
  done

  python3 "${APPLY_ADAPTIVE_MT}" \
    --input-h5ad "${unfiltered_h5ad}" \
    --threshold-json "${threshold_json}" \
    --mt-floor "${MT_FLOOR}" \
    --n-mad "${N_MAD}" >> "${RUN_LOG}" 2>&1

  python3 "${APPLY_ADAPTIVE_MT}" \
    --input-h5ad "${final_h5ad}" \
    --threshold-json "${threshold_json}" \
    --mt-floor "${MT_FLOOR}" \
    --n-mad "${N_MAD}" >> "${RUN_LOG}" 2>&1

  python3 "${POSTPROCESS_FILTERS}" \
    --unfiltered-h5ad "${unfiltered_h5ad}" \
    --qc-output-h5ad "${filtered_h5ad}" \
    --default-singlet-output-h5ad "${default_singlet_h5ad}" >> "${RUN_LOG}" 2>&1

  python3 "${GENERATE_QC_HISTOGRAM_MT}" \
    --input-h5ad "${unfiltered_h5ad}" \
    --output-dir "${downstream}" \
    --threshold-json "${threshold_json}" >> "${RUN_LOG}" 2>&1

  "${MUDATA_PYTHON}" "${BUILD_MUDATA}" \
    --rna-h5ad "${final_h5ad}" \
    --atac-mex-dir "${atac_mex}" \
    --per-barcode-metrics "${atac_metrics}" \
    --metrics-barcode-column barcode \
    --require-rna-velocyto-layers \
    --cell-call-source star_downstream_h5ad_chromap_atac \
    --rna-source "${final_h5ad}" \
    --atac-source "${atac_mex}" \
    --fragments-source "${atac_sidecar}" \
    --peaks-source "${atac_peaks}" \
    --evidence-source "${atac_metrics}" \
    --y-removal-enabled true \
    --output-h5mu "${unfiltered_h5mu}" >> "${RUN_LOG}" 2>&1

  "${MUDATA_PYTHON}" "${BUILD_MUDATA}" \
    --rna-h5ad "${filtered_h5ad}" \
    --atac-mex-dir "${atac_mex}" \
    --per-barcode-metrics "${atac_metrics}" \
    --metrics-barcode-column barcode \
    --all-barcodes-are-cells \
    --allow-empty-barcode-intersection \
    --require-rna-velocyto-layers \
    --cell-call-source star_downstream_filtered_h5ad_chromap_atac \
    --rna-source "${filtered_h5ad}" \
    --atac-source "${atac_mex}" \
    --fragments-source "${atac_sidecar}" \
    --peaks-source "${atac_peaks}" \
    --evidence-source "${atac_metrics}" \
    --y-removal-enabled true \
    --output-h5mu "${filtered_h5mu}" >> "${RUN_LOG}" 2>&1

  summarize_threshold "${sample}" "${threshold_json}"
  {
    echo "sample=${sample}"
    echo "run_root=${RUN_ROOT}"
    echo "downstream=${downstream}"
    echo "mudata_dir=${mudata_dir}"
    echo "mt_floor=${MT_FLOOR}"
    echo "n_mad=${N_MAD}"
    echo "converted_at=$(date -Is)"
    echo "log=${RUN_LOG}"
  } > "${marker}"
  converted=$((converted + 1))
  log "Completed ${sample}"
done

log "Done. converted=${converted} skipped=${skipped}"
