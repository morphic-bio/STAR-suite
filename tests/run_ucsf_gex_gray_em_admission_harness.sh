#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BASE_RUN="${UCSF_GRAY_EM_BASE_RUN:-/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_optimized_current_zcat_20260326}"
EM_RUN="${UCSF_GRAY_EM_REFERENCE_RUN:-/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_cellgeni_historical_7a7fb08_truwhitelist_genefullonly_20260326}"
CR_FILTERED="${UCSF_GRAY_EM_CR_FILTERED:-/storage/cr9_ebs2_2_benchmark_20260318/cr9_ebs2_2/outs/per_sample_outs/cr9_ebs2_2/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz}"
OUTDIR="${UCSF_GRAY_EM_OUTDIR:-tests/ucsf_gex_gray_em_harness_output_$(date +%Y%m%d_%H%M%S)}"
WINDOW_BEFORE="${UCSF_GRAY_EM_WINDOW_BEFORE:-2000}"
WINDOW_AFTER="${UCSF_GRAY_EM_WINDOW_AFTER:-2000}"

BASE_RAW_DIR="${UCSF_GRAY_EM_BASE_RAW_DIR:-${BASE_RUN}/Solo.out/GeneFull/raw}"
BASE_FILTERED="${UCSF_GRAY_EM_BASE_FILTERED:-${BASE_RUN}/Solo.out/GeneFull/filtered/barcodes.tsv}"
BASE_UMI_SORTED="${UCSF_GRAY_EM_BASE_UMI_SORTED:-${BASE_RUN}/Solo.out/GeneFull/UMIperCellSorted.txt}"
EM_FILTERED="${UCSF_GRAY_EM_REFERENCE_FILTERED:-${EM_RUN}/output/GeneFull/filtered/barcodes.tsv}"

mkdir -p "${OUTDIR}"

python3 "${REPO_ROOT}/scripts/evaluate_gray_em_admission.py" \
  --base-raw-dir "${BASE_RAW_DIR}" \
  --base-filtered-barcodes "${BASE_FILTERED}" \
  --base-umi-sorted "${BASE_UMI_SORTED}" \
  --em-filtered-barcodes "${EM_FILTERED}" \
  --cr-filtered-barcodes "${CR_FILTERED}" \
  --window-before "${WINDOW_BEFORE}" \
  --window-after "${WINDOW_AFTER}" \
  --outdir "${OUTDIR}"

echo "Output directory: ${OUTDIR}"
