#!/usr/bin/env bash
set -euo pipefail

# Cheap simulation: evaluate whether switching STAR's SNP statistic from conversion-only -> mismatch-agnostic
# (GEDI-like) could plausibly close divergence, without rerunning STAR/GEDI.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

COMPARE_DIR="/storage/slam_e2e_arid1a_20260113/compare/fixed"
BAM="/storage/slam_e2e_arid1a_20260113/star/wt_Aligned.sortedByCoord.out.bam"
REF="/storage/autoindex_110_44/bulk_index/cellranger_ref/genome.fa"
REF_FAI="${REF}.fai"

STAR_ONLY="${COMPARE_DIR}/star_only.bed"
GEDI_ONLY="${COMPARE_DIR}/gedi_only.bed"

OUT_DIR="/storage/simulate_gedi_compat_knobs_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${OUT_DIR}"
OUT_PREFIX="${OUT_DIR}/sim"

echo "OUT_DIR=${OUT_DIR}"

python3 "${ROOT_DIR}/tests/simulate_gedi_compat_knobs.py" \
  --bam "${BAM}" \
  --ref "${REF}" \
  --ref-fai "${REF_FAI}" \
  --star-only-bed "${STAR_ONLY}" \
  --gedi-only-bed "${GEDI_ONLY}" \
  --sample-per-set 200 \
  --seed 1 \
  --p-tc 0.001 \
  --p-any -1 \
  --p-any-sample-n 2000 \
  --pval 0.001 \
  --min-ratio 0.3 \
  --min-cov 6 \
  --min-alt 1 \
  --out-prefix "${OUT_PREFIX}" | tee "${OUT_DIR}/summary.txt"

echo ""
echo "Wrote:"
echo "  ${OUT_PREFIX}.tsv"
echo "  ${OUT_PREFIX}.sampled_loci.bed"
echo "  ${OUT_DIR}/summary.txt"

