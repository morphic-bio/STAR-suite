#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
ASSIGN_BIN="${ASSIGN_BIN:-${REPO_ROOT}/core/features/process_features/assignBarcodes}"
OUT_DIR="${CATATAC_GUIDE_SMOKE_OUT:-${REPO_ROOT}/tests/catatac_guide_arm_smoke_output}"
FIXTURE_DIR="${CATATAC_GUIDE_FIXTURE:-/mnt/pikachu/catatac_gse288996/guide_redump/fixture}"
FEATURE_REF="${CATATAC_GUIDE_FEATURE_REF:-${REPO_ROOT}/core/features/process_features/feature_lists/catatac_crispri_guide_capture.csv}"
ATAC_WHITELIST="${CATATAC_ATAC_WHITELIST:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt}"
ATAC2GEX="${CATATAC_ATAC2GEX_MAP:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/atac2gex.tsv}"

rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"
ORACLE_OUT="${OUT_DIR}/oracle"
NATIVE_OUT="${OUT_DIR}/native"

if [[ ! -x "${ASSIGN_BIN}" ]]; then
  make -C "${REPO_ROOT}/core/features/process_features" assignBarcodes
fi
if [[ ! -x "${REPO_ROOT}/core/features/process_features/tests/test_catatac_split_read" ]]; then
  make -C "${REPO_ROOT}/core/features/process_features" tests/test_catatac_split_read
fi

export CATATAC_GUIDE_MAX_READS="${CATATAC_GUIDE_MAX_READS:-50000}"

python3 "${REPO_ROOT}/tests/catatac_guide_oracle.py" \
  --fixture-dir "${FIXTURE_DIR}" \
  --feature-ref "${FEATURE_REF}" \
  --whitelist "${ATAC_WHITELIST}" \
  --output-map "${ATAC2GEX}" \
  --oracle-out "${ORACLE_OUT}" \
  --native-out "${NATIVE_OUT}" \
  --limit "${CATATAC_GUIDE_MAX_READS}" \
  --gex-whitelist "${CATATAC_GEX_WHITELIST:-/mnt/pikachu/GEX_whitelist/737K-arc-v1.txt}"

bash "${SCRIPT_DIR}/test_catatac_srr_fastq_naming.sh"

echo "CAT-ATAC guide-arm smoke completed: ${OUT_DIR}"
