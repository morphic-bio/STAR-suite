#!/usr/bin/env bash
set -euo pipefail

# Create tiered downsample FASTQs for A375 (GEX + CRISPR) without overwriting existing outputs.
# Uses the process_features downsample script, staging via symlinks to avoid clobbering.

ROOT="/storage/A375"
FASTQ_ROOT="${ROOT}/fastqs/1k_CRISPR_5p_gemx_fastqs"
GEX_DIR="${FASTQ_ROOT}/gex"
CRISPR_DIR="${FASTQ_ROOT}/crispr"
DOWN_SCRIPT="/mnt/pikachu/process_features/scripts/downsample_fastq_directory.sh"

# Tiers are counts per FASTQ file (reads = 4 lines).
TIERS=(1000000 500000 100000 50000 10000)

if [[ ! -x "${DOWN_SCRIPT}" ]]; then
  echo "Missing downsample script: ${DOWN_SCRIPT}" >&2
  exit 1
fi

stage_and_downsample() {
  local src_dir="$1"
  local nreads="$2"
  local out_dir="${src_dir}/downsampled_${nreads}_v2"

  if [[ -d "${out_dir}" ]]; then
    echo "Skip: ${out_dir} already exists"
    return
  fi

  local tmp_dir
  tmp_dir="$(mktemp -d)"
  trap 'rm -rf "${tmp_dir}"' EXIT

  # Hardlink FASTQs into staging directory (downsample script ignores symlinks).
  shopt -s nullglob
  for fq in "${src_dir}"/*.fastq "${src_dir}"/*.fq "${src_dir}"/*.fastq.gz "${src_dir}"/*.fq.gz; do
    ln "${fq}" "${tmp_dir}/$(basename "${fq}")"
  done
  shopt -u nullglob

  if ! ls "${tmp_dir}"/* >/dev/null 2>&1; then
    echo "No FASTQs found in ${src_dir}" >&2
    rm -rf "${tmp_dir}"
    trap - EXIT
    return 1
  fi

  "${DOWN_SCRIPT}" "${nreads}" "${tmp_dir}"
  mv "${tmp_dir}/downsampled" "${out_dir}"

  rm -rf "${tmp_dir}"
  trap - EXIT

  echo "Wrote: ${out_dir}"
}

for nreads in "${TIERS[@]}"; do
  stage_and_downsample "${GEX_DIR}" "${nreads}"
  stage_and_downsample "${CRISPR_DIR}" "${nreads}"
done

echo "Done. Tiered outputs under:"
echo "  ${GEX_DIR}/downsampled_<N>_v2"
echo "  ${CRISPR_DIR}/downsampled_<N>_v2"
