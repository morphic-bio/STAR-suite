#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

HEUR_DIR="${UCSF_GRAY_HEUR_DIR:-tests/ucsf_gex_gray_em_heuristics_output_20260326_175723}"
HARNESS_DIR="${UCSF_GRAY_HARNESS_DIR:-tests/ucsf_gex_gray_em_harness_output_20260326_175204}"
BASE_RUN="${UCSF_GRAY_BASE_RUN:-/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_optimized_current_zcat_20260326}"
CR_FILTERED="${UCSF_GRAY_CR_FILTERED:-/storage/cr9_ebs2_2_benchmark_20260318/cr9_ebs2_2/outs/per_sample_outs/cr9_ebs2_2/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz}"
FOLDS="${UCSF_GRAY_CV_FOLDS:-5}"
SEED="${UCSF_GRAY_CV_SEED:-42}"
MODE="${UCSF_GRAY_CV_MODE:-rescue_only}"
OUTDIR="${UCSF_GRAY_CV_OUTDIR:-tests/ucsf_gex_gray_threshold_cv_output_$(date +%Y%m%d_%H%M%S)}"

python3 "${REPO_ROOT}/scripts/validate_gray_admission_thresholds.py" \
  --per-barcode-metrics "${HEUR_DIR}/per_barcode_metrics.tsv" \
  --gray-zone-barcodes "${HARNESS_DIR}/gray_zone_barcodes.tsv" \
  --base-filtered-barcodes "${BASE_RUN}/Solo.out/GeneFull/filtered/barcodes.tsv" \
  --cr-filtered-barcodes "${CR_FILTERED}" \
  --mode "${MODE}" \
  --folds "${FOLDS}" \
  --seed "${SEED}" \
  --outdir "${OUTDIR}"

echo "Output directory: ${OUTDIR}"
