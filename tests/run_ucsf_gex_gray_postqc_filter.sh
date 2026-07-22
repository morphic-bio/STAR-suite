#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BASE_RUN="${UCSF_GRAY_POSTQC_BASE_RUN:-/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_optimized_current_zcat_20260326}"
HARNESS_DIR="${UCSF_GRAY_POSTQC_HARNESS_DIR:-tests/ucsf_gex_gray_em_harness_output_20260326_175204}"
CR_FILTERED="${UCSF_GRAY_POSTQC_CR_FILTERED:-/storage/cr9_ebs2_2_benchmark_20260318/cr9_ebs2_2/outs/per_sample_outs/cr9_ebs2_2/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz}"
OUTDIR="${UCSF_GRAY_POSTQC_OUTDIR:-tests/ucsf_gex_gray_postqc_filter_output_$(date +%Y%m%d_%H%M%S)}"
MIN_GENES="${UCSF_GRAY_POSTQC_MIN_GENES:-200}"
MAX_GENES="${UCSF_GRAY_POSTQC_MAX_GENES:-2500}"
MT_PCT_CUTOFF="${UCSF_GRAY_POSTQC_MT_PCT_CUTOFF:-5}"
ADAPTIVE_MAX_GENES="${UCSF_GRAY_POSTQC_ADAPTIVE_MAX_GENES:-0}"
N_MAD="${UCSF_GRAY_POSTQC_N_MAD:-3}"

EXTRA_ARGS=()
if [[ "${ADAPTIVE_MAX_GENES}" == "1" ]]; then
  EXTRA_ARGS+=(--adaptive-max-genes --n-mad "${N_MAD}")
else
  EXTRA_ARGS+=(--max-genes "${MAX_GENES}")
fi

python3 "${REPO_ROOT}/scripts/evaluate_gray_postqc_filter.py" \
  --base-raw-dir "${BASE_RUN}/Solo.out/GeneFull/raw" \
  --base-filtered-barcodes "${BASE_RUN}/Solo.out/GeneFull/filtered/barcodes.tsv" \
  --gray-zone-barcodes "${HARNESS_DIR}/gray_zone_barcodes.tsv" \
  --cr-filtered-barcodes "${CR_FILTERED}" \
  --min-genes "${MIN_GENES}" \
  --mt-pct-cutoff "${MT_PCT_CUTOFF}" \
  "${EXTRA_ARGS[@]}" \
  --outdir "${OUTDIR}"

echo "Output directory: ${OUTDIR}"
