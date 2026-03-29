#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
WORKFLOW="${REPO_ROOT}/scripts/paper/run_ucsf_corrected_production_workflow.sh"

SAMPLE="${UCSF_100K_SMOKE_SAMPLE:-EBs2_2}"
DATASET_ROOT="${UCSF_100K_SMOKE_DATASET_ROOT:-/mnt/pikachu/ucsf-perturb-seq-corrected}"
FEATURE_REF="${UCSF_100K_SMOKE_FEATURE_REF:-/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv}"
GENOME_DIR="${UCSF_100K_SMOKE_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
SOLO_CB_WHITELIST="${UCSF_100K_SMOKE_SOLO_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt}"
CR_WHITELIST="${UCSF_100K_SMOKE_CR_WHITELIST:-${SOLO_CB_WHITELIST}}"
STAR_BIN="${UCSF_100K_SMOKE_STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR.release}"
THREADS="${UCSF_100K_SMOKE_THREADS:-24}"
READS="${UCSF_100K_SMOKE_READS:-100000}"
SOURCE_SAMPLE_ROOT="${UCSF_100K_SMOKE_SOURCE_SAMPLE_ROOT:-${SCRIPT_DIR}/ucsf_corrected_production_100k_output_20260329_074313/samples/${SAMPLE}}"
OUT_ROOT="${UCSF_100K_SMOKE_OUT_ROOT:-${SCRIPT_DIR}/ucsf_corrected_production_100k_output_$(date +%Y%m%d_%H%M%S)}"
SUMMARY_FILE="${OUT_ROOT}/SUMMARY.txt"
REUSE_OUT_ROOT="${UCSF_100K_SMOKE_REUSE_OUT_ROOT:-0}"
TEMP_FIXTURE_ROOT="${UCSF_100K_SMOKE_TEMP_FIXTURE_ROOT:-/tmp/ucsf_corrected_production_fixture_${SAMPLE}_${READS}}"

die() {
  echo "ERROR: $*" >&2
  exit 1
}

[[ -x "${WORKFLOW}" ]] || die "Missing workflow wrapper: ${WORKFLOW}"
[[ -x "${STAR_BIN}" ]] || die "Missing release binary: ${STAR_BIN}"
[[ -d "${DATASET_ROOT}/${SAMPLE}/GEX" ]] || die "Missing corrected GEX dir for ${SAMPLE}"
[[ -d "${DATASET_ROOT}/${SAMPLE}/guides" ]] || die "Missing corrected guide dir for ${SAMPLE}"

if [[ "${REUSE_OUT_ROOT}" != "1" ]]; then
  rm -rf "${OUT_ROOT}"
fi
mkdir -p "${OUT_ROOT}"

USE_EXISTING_STAGED_FIXTURE=0
if [[ -d "${SOURCE_SAMPLE_ROOT}/staged_input" && -f "${SOURCE_SAMPLE_ROOT}/pf_multi_config.csv" ]]; then
  USE_EXISTING_STAGED_FIXTURE=1
fi

echo "=== UCSF corrected production 100K smoke ==="
echo "Sample: ${SAMPLE}"
echo "Output root: ${OUT_ROOT}"
echo "Threads: ${THREADS}"
echo "Downsample reads per library: ${READS}"
if [[ "${USE_EXISTING_STAGED_FIXTURE}" == "1" ]]; then
  echo "Reusing staged 100K fixture from: ${SOURCE_SAMPLE_ROOT}"
fi

if [[ "${USE_EXISTING_STAGED_FIXTURE}" == "1" ]]; then
  rm -rf "${TEMP_FIXTURE_ROOT}"
  mkdir -p "${TEMP_FIXTURE_ROOT}/${SAMPLE}"
  ln -s "${SOURCE_SAMPLE_ROOT}/staged_input/GEX/${SAMPLE}" "${TEMP_FIXTURE_ROOT}/${SAMPLE}/GEX"
  ln -s "${SOURCE_SAMPLE_ROOT}/staged_input/guides/${SAMPLE}" "${TEMP_FIXTURE_ROOT}/${SAMPLE}/guides"
  cp -f "${SOURCE_SAMPLE_ROOT}/pf_multi_config.csv" "${TEMP_FIXTURE_ROOT}/${SAMPLE}/pf_multi_config.csv"
  "${WORKFLOW}" \
    --samples "${SAMPLE}" \
    --dataset-root "${TEMP_FIXTURE_ROOT}" \
    --feature-ref "${FEATURE_REF}" \
    --genome-dir "${GENOME_DIR}" \
    --solo-cb-whitelist "${SOLO_CB_WHITELIST}" \
    --cr-whitelist "${CR_WHITELIST}" \
    --out-root "${OUT_ROOT}" \
    --threads "${THREADS}" \
    --star-bin "${STAR_BIN}" \
    --cellbender-cpu-cores "${THREADS}"
else
  "${WORKFLOW}" \
    --samples "${SAMPLE}" \
    --dataset-root "${DATASET_ROOT}" \
    --feature-ref "${FEATURE_REF}" \
    --genome-dir "${GENOME_DIR}" \
    --solo-cb-whitelist "${SOLO_CB_WHITELIST}" \
    --cr-whitelist "${CR_WHITELIST}" \
    --out-root "${OUT_ROOT}" \
    --threads "${THREADS}" \
    --star-bin "${STAR_BIN}" \
    --cellbender-cpu-cores "${THREADS}" \
    --downsample-reads "${READS}"
fi

SAMPLE_ROOT="${OUT_ROOT}/samples/${SAMPLE}"
RUN_DIR="${SAMPLE_ROOT}/run"
DOWNSTREAM_DIR="${SAMPLE_ROOT}/downstream_genefull_velocyto_cellbender"

for path in \
  "${OUT_ROOT}/RUN_WRAPPER_COMMAND.sh" \
  "${OUT_ROOT}/RUNS.tsv" \
  "${SAMPLE_ROOT}/RUN_COMMAND.sh" \
  "${SAMPLE_ROOT}/RUN_MANIFEST.txt" \
  "${SAMPLE_ROOT}/RUN_COMPLETE.ok" \
  "${SAMPLE_ROOT}/globus_batch.tsv" \
  "${SAMPLE_ROOT}/RUN_GLOBUS_TRANSFER.sh" \
  "${RUN_DIR}/Aligned.out_Y.bam" \
  "${RUN_DIR}/Aligned.out_noY.bam" \
  "${RUN_DIR}/outs/raw_feature_bc_matrix/matrix.mtx.gz" \
  "${RUN_DIR}/outs/filtered_feature_bc_matrix/matrix.mtx.gz" \
  "${RUN_DIR}/outs/raw_velocyto_feature_bc_matrix/matrix.mtx.gz" \
  "${RUN_DIR}/outs/filtered_velocyto_feature_bc_matrix/matrix.mtx.gz" \
  "${RUN_DIR}/outs/velocyto_feature_bc_matrix_manifest.json" \
  "${DOWNSTREAM_DIR}/counts.h5ad" \
  "${DOWNSTREAM_DIR}/unfiltered_counts.h5ad" \
  "${DOWNSTREAM_DIR}/filtered_counts.h5ad" \
  "${DOWNSTREAM_DIR}/default_singlet_filtered_counts.h5ad" \
  "${DOWNSTREAM_DIR}/final_counts.h5ad" \
  "${DOWNSTREAM_DIR}/summary.txt"; do
  [[ -f "${path}" ]] || die "Missing expected output: ${path}"
done

if [[ ! -f "${DOWNSTREAM_DIR}/cellbender/cellbender_counts.h5" ]] && [[ ! -f "${DOWNSTREAM_DIR}/cellbender/CELLBENDER_FAILED.txt" ]]; then
  die "Expected either CellBender output or failure note in ${DOWNSTREAM_DIR}/cellbender"
fi

[[ -d "${RUN_DIR}/y_separated" ]] || die "Missing y_separated directory"
find "${RUN_DIR}/y_separated" -maxdepth 1 -type f -name '*.fastq.gz' | grep -q . || die "No Y/noY FASTQs emitted"

readarray -t COUNTS < <(python3 - <<'PY' "${RUN_DIR}/outs/velocyto_feature_bc_matrix_manifest.json" "${DOWNSTREAM_DIR}/summary.txt"
import json
import sys
from pathlib import Path

manifest_path = Path(sys.argv[1])
summary_path = Path(sys.argv[2])
manifest = json.loads(manifest_path.read_text())

print(f"RAW_VELOCYTO_BARCODES={manifest['raw']['barcodes']}")
print(f"FILTERED_VELOCYTO_BARCODES={manifest['filtered']['barcodes']}")
print(f"RAW_VELOCYTO_NNZ={manifest['raw']['nnz_total']}")
print(f"FILTERED_VELOCYTO_NNZ={manifest['filtered']['nnz_total']}")

summary_lines = [line.strip() for line in summary_path.read_text().splitlines() if line.strip()]
print(f"SUMMARY_LINES={len(summary_lines)}")
print(f"CELLBENDER_STATUS={'produced_h5' if (summary_path.parent / 'cellbender' / 'cellbender_counts.h5').exists() else 'fallback_note'}")
PY
)

declare -A METRICS=()
for line in "${COUNTS[@]}"; do
  key="${line%%=*}"
  value="${line#*=}"
  METRICS["${key}"]="${value}"
done

{
  echo "UCSF corrected production 100K smoke: PASS"
  echo "Sample: ${SAMPLE}"
  echo "STAR binary: ${STAR_BIN}"
  echo "Threads: ${THREADS}"
  echo "Downsample reads per library: ${READS}"
  echo "Output root: ${OUT_ROOT}"
  echo "Raw velocyto barcodes: ${METRICS[RAW_VELOCYTO_BARCODES]}"
  echo "Filtered velocyto barcodes: ${METRICS[FILTERED_VELOCYTO_BARCODES]}"
  echo "Raw velocyto nnz: ${METRICS[RAW_VELOCYTO_NNZ]}"
  echo "Filtered velocyto nnz: ${METRICS[FILTERED_VELOCYTO_NNZ]}"
  echo "CellBender status: ${METRICS[CELLBENDER_STATUS]}"
  echo "Downstream summary lines: ${METRICS[SUMMARY_LINES]}"
} > "${SUMMARY_FILE}"

cat "${SUMMARY_FILE}"
