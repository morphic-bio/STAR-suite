#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASSIGN_BIN="${ASSIGN_BIN:-${ROOT_DIR}/core/features/process_features/assignBarcodes}"
WHITELIST="${WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-february-2018_TRU.txt}"
FEATURE_CSV="${FEATURE_CSV:-/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv}"
GUIDE_FASTQ_DIR="${GUIDE_FASTQ_DIR:-/storage/ucsf-2M/guides/iPSC2_1_AALG2}"
# Optional. Default empty to benchmark matcher behavior without filtered-MEX side effects.
FILTERED_BARCODES="${FILTERED_BARCODES:-}"

MAX_READS="${MAX_READS:-100000}"
FEATURE_OFFSET="${FEATURE_OFFSET:-0}"
LIMIT_SEARCH="${LIMIT_SEARCH:-1}"
FEATURE_N="${FEATURE_N:-1}"
BARCODE_N="${BARCODE_N:-2}"
BARCODE_OFFSET="${BARCODE_OFFSET:-0}"
CONSUMER_THREADS="${CONSUMER_THREADS:-8}"
SEARCH_THREADS="${SEARCH_THREADS:-1}"

OUT_BASE="${OUT_BASE:-/tmp/pf_limited_fallback_ucsf2m_$(date +%Y%m%d_%H%M%S)}"
MODES="${MODES:-full simple}"
HAMMINGS="${HAMMINGS:-1 2 3 4 5}"

if [[ ! -x "${ASSIGN_BIN}" ]]; then
  echo "ERROR: assignBarcodes binary not found or not executable: ${ASSIGN_BIN}" >&2
  echo "Build first: make feature-barcodes-tools" >&2
  exit 1
fi
if [[ ! -f "${WHITELIST}" ]]; then
  echo "ERROR: whitelist not found: ${WHITELIST}" >&2
  exit 1
fi
if [[ ! -f "${FEATURE_CSV}" ]]; then
  echo "ERROR: feature CSV not found: ${FEATURE_CSV}" >&2
  exit 1
fi
if [[ ! -d "${GUIDE_FASTQ_DIR}" ]]; then
  echo "ERROR: guide FASTQ dir not found: ${GUIDE_FASTQ_DIR}" >&2
  exit 1
fi

mapfile -t R1_FILES < <(find "${GUIDE_FASTQ_DIR}" -maxdepth 1 -type f \( -name '*_R1_*.fastq.gz' -o -name '*_R1_*.fq.gz' \) | sort)
mapfile -t R2_FILES < <(find "${GUIDE_FASTQ_DIR}" -maxdepth 1 -type f \( -name '*_R2_*.fastq.gz' -o -name '*_R2_*.fq.gz' \) | sort)

if [[ ${#R1_FILES[@]} -eq 0 || ${#R2_FILES[@]} -eq 0 ]]; then
  echo "ERROR: could not find R1/R2 FASTQs in ${GUIDE_FASTQ_DIR}" >&2
  exit 1
fi
if [[ ${#R1_FILES[@]} -ne ${#R2_FILES[@]} ]]; then
  echo "ERROR: R1/R2 file count mismatch (${#R1_FILES[@]} vs ${#R2_FILES[@]})" >&2
  exit 1
fi

mkdir -p "${OUT_BASE}"
SUMMARY_TSV="${OUT_BASE}/summary.tsv"
{
  echo -e "mode\tmax_hamming\twall_sec\tmax_rss_kb\ttotal_feature_counts\ttotal_deduped_feature_counts\ttotal_barcodes\ttotal_unmatched_reads\tlimited_exact_checks\tlimited_exact_hits\tlimited_simple_fallback_calls\tlimited_simple_fallback_hits\tlimited_full_fallback_calls\tlimited_full_fallback_hits\tstats_file"
} > "${SUMMARY_TSV}"

extract_value() {
  local key="$1"
  local file="$2"
  awk -v k="$key" '$1==k {print $2; exit}' "$file"
}

for mode in ${MODES}; do
  for h in ${HAMMINGS}; do
    run_dir="${OUT_BASE}/${mode}_h${h}"
    mkdir -p "${run_dir}"

    cmd=(
      "${ASSIGN_BIN}"
      -w "${WHITELIST}"
      -f "${FEATURE_CSV}"
      -d "${run_dir}"
      -m "${h}"
      -o "${FEATURE_OFFSET}"
      --limit_search "${LIMIT_SEARCH}"
      --feature_limited_fallback "${mode}"
      --feature_n "${FEATURE_N}"
      --barcode_n "${BARCODE_N}"
      -B "${BARCODE_OFFSET}"
      -c "${CONSUMER_THREADS}"
      -S "${SEARCH_THREADS}"
      -t 1
      --max_reads "${MAX_READS}"
      "${GUIDE_FASTQ_DIR}"
    )

    if [[ -n "${FILTERED_BARCODES}" && -f "${FILTERED_BARCODES}" ]]; then
      cmd+=(--filtered_barcodes "${FILTERED_BARCODES}")
    fi

    printf '%q ' "${cmd[@]}" > "${run_dir}/cmd.sh"
    printf '\n' >> "${run_dir}/cmd.sh"

    echo "[run] mode=${mode} h=${h}"
    (
      /usr/bin/time -f "wall_sec=%e\nuser_sec=%U\nsys_sec=%S\nmax_rss_kb=%M" -o "${run_dir}/time.txt" "${cmd[@]}"
    ) > "${run_dir}/stdout.log" 2> "${run_dir}/stderr.log"

    stats_file="$(find "${run_dir}" -type f -name stats.txt | head -n 1 || true)"
    if [[ -z "${stats_file}" ]]; then
      echo "ERROR: stats.txt not found for mode=${mode} h=${h}" >&2
      exit 1
    fi

    wall_sec="$(awk -F= '/^wall_sec=/{print $2; exit}' "${run_dir}/time.txt")"
    max_rss_kb="$(awk -F= '/^max_rss_kb=/{print $2; exit}' "${run_dir}/time.txt")"
    total_feature_counts="$(awk '/^Total feature counts /{print $4; exit}' "${stats_file}")"
    total_deduped="$(awk '/^Total deduped feature counts /{print $5; exit}' "${stats_file}")"
    total_barcodes="$(awk '/^Total barcodes /{print $3; exit}' "${stats_file}")"
    total_unmatched="$(awk '/^Total_unmatched_reads /{print $2; exit}' "${stats_file}")"

    limited_exact_checks="$(extract_value Limited_exact_checks "${stats_file}")"
    limited_exact_hits="$(extract_value Limited_exact_hits "${stats_file}")"
    limited_simple_calls="$(extract_value Limited_simple_fallback_calls "${stats_file}")"
    limited_simple_hits="$(extract_value Limited_simple_fallback_hits "${stats_file}")"
    limited_full_calls="$(extract_value Limited_full_fallback_calls "${stats_file}")"
    limited_full_hits="$(extract_value Limited_full_fallback_hits "${stats_file}")"

    {
      echo -e "${mode}\t${h}\t${wall_sec}\t${max_rss_kb}\t${total_feature_counts}\t${total_deduped}\t${total_barcodes}\t${total_unmatched}\t${limited_exact_checks}\t${limited_exact_hits}\t${limited_simple_calls}\t${limited_simple_hits}\t${limited_full_calls}\t${limited_full_hits}\t${stats_file}"
    } >> "${SUMMARY_TSV}"
  done
done

echo "Done. Summary: ${SUMMARY_TSV}"
