#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
RUN_SCRIPT="${REPO_ROOT}/scripts/run_ucsf_perturb_yremove_batch.sh"
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
CR_CONFIG="${CR_CONFIG:-/mnt/pikachu/ucsf-perturb-yremove-smoke_1m_integrated_20260311_221028/cellranger_iPSC2_1_AALG2_1m_match/config.csv}"
OUT_ROOT="${OUT_ROOT:-/tmp/ucsf_cr_config_1m_smoke_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-8}"

if [[ ! -f "${CR_CONFIG}" ]]; then
  echo "SKIP: missing UCSF 1M Cell Ranger config: ${CR_CONFIG}"
  exit 0
fi

"${RUN_SCRIPT}" \
  --star-bin "${STAR_BIN}" \
  --cr-config "${CR_CONFIG}" \
  --threads "${THREADS}" \
  --out-root "${OUT_ROOT}"

sample_root="${OUT_ROOT}/samples/iPSC2_1_AALG2"

[[ -f "${sample_root}/RUN_COMPLETE.ok" ]]
[[ -f "${sample_root}/run/Aligned.out_Y.bam" ]]
[[ -f "${sample_root}/run/Aligned.out_noY.bam" ]]
[[ -d "${sample_root}/run/y_separated" ]]
[[ "$(find "${sample_root}/run/y_separated" -maxdepth 1 -type f -name '*.fastq.gz' | wc -l)" == "4" ]]
[[ -f "${sample_root}/run/outs/filtered_feature_bc_matrix/matrix.mtx" || -f "${sample_root}/run/outs/filtered_feature_bc_matrix/matrix.mtx.gz" ]]
[[ -f "${sample_root}/run/outs/crispr_analysis/protospacer_calls_per_cell.csv" ]]

y_count="$(samtools view -c "${sample_root}/run/Aligned.out_Y.bam")"
noy_count="$(samtools view -c "${sample_root}/run/Aligned.out_noY.bam")"
(( y_count >= 0 ))
(( noy_count > 0 ))

grep -F "cr_config=${CR_CONFIG}" "${sample_root}/RUN_MANIFEST.txt" >/dev/null

echo "PASS: ${OUT_ROOT}"
