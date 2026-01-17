#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${WORK_DIR:-${SCRIPT_DIR}/tmp_vcf_mask_overlap}"
VCF_MASK="${VCF_MASK:-}"
MASK_A="${MASK_A:-/storage/gedi_common_mask_20260114_202352/gedi_mask.bed.gz}"
MASK_B="${MASK_B:-/storage/slam_e2e_kmode_baseline_maskON_20260114_190244/mask/wt0.mask.bed.gz}"
OUT_REPORT="${OUT_REPORT:-${WORK_DIR}/vcf_mask_overlap.tsv}"
MAX_SPAN="${MAX_SPAN:-1000}"

if [[ -z "${VCF_MASK}" ]]; then
  echo "ERROR: Set VCF_MASK to the VCF-derived BED path." >&2
  exit 1
fi

mkdir -p "${WORK_DIR}"

python3 "${SCRIPT_DIR}/compare_mask_overlap.py" \
  --vcf-mask "${VCF_MASK}" \
  --mask-a "${MASK_A}" \
  --mask-b "${MASK_B}" \
  --max-span "${MAX_SPAN}" \
  --out "${OUT_REPORT}"
