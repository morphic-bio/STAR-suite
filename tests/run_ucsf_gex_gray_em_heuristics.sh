#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

HARNESS_OUTDIR="${UCSF_GRAY_EM_HARNESS_OUTDIR:-tests/ucsf_gex_gray_em_harness_output_20260326_175204}"
BASE_RUN="${UCSF_GRAY_EM_BASE_RUN:-/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_optimized_current_zcat_20260326}"
EM_RUN="${UCSF_GRAY_EM_REFERENCE_RUN:-/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_cellgeni_historical_7a7fb08_truwhitelist_genefullonly_20260326}"
CR_FILTERED="${UCSF_GRAY_EM_CR_FILTERED:-/storage/cr9_ebs2_2_benchmark_20260318/cr9_ebs2_2/outs/per_sample_outs/cr9_ebs2_2/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz}"
OUTDIR="${UCSF_GRAY_EM_HEUR_OUTDIR:-tests/ucsf_gex_gray_em_heuristics_output_$(date +%Y%m%d_%H%M%S)}"

python3 "${REPO_ROOT}/scripts/analyze_gray_admission_heuristics.py" \
  --base-raw-dir "${BASE_RUN}/Solo.out/GeneFull/raw" \
  --gray-harness-outdir "${HARNESS_OUTDIR}" \
  --base-filtered-barcodes "${BASE_RUN}/Solo.out/GeneFull/filtered/barcodes.tsv" \
  --em-filtered-barcodes "${EM_RUN}/output/GeneFull/filtered/barcodes.tsv" \
  --cr-filtered-barcodes "${CR_FILTERED}" \
  --outdir "${OUTDIR}"

echo "Output directory: ${OUTDIR}"
