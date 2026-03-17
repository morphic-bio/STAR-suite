#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
RUNNER="${TRIM_QC_FASTQ_RUNNER:-${SCRIPT_DIR}/run_trim_qc_fastq.sh}"

ROOT="${TRIM_QC_FASTQ_ROOT:-/mnt/pikachu/ucsf-perturb-seq/GEX}"
OUTPUT_ROOT="${TRIM_QC_FASTQ_OUTPUT_ROOT:-${ROOT}}"
REPORT_BASENAME="${TRIM_QC_FASTQ_BASENAME:-read_qc_r2}"
STAGE="${TRIM_QC_FASTQ_STAGE:-fastq_replay_r2}"
MAX_READS="${TRIM_QC_FASTQ_MAX_READS:-0}"
MATE_COUNT="${TRIM_QC_FASTQ_MATE_COUNT:-1}"
MATE_INDEX="${TRIM_QC_FASTQ_MATE_INDEX:-1}"
SAMPLE_FILTER="${TRIM_QC_FASTQ_SAMPLE_FILTER:-}"

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --root DIR             Sample root containing one directory per sample
  --output-root DIR      Output sample root for qc/QC directories
  --basename NAME        Report basename inside qc/QC dir (default: read_qc_r2)
  --stage LABEL          Stage label for trim-QC output (default: fastq_replay_r2)
  --max-reads N          Limit reads per sample (default: 0 = all)
  --mate-count N         Total mate count reported in JSON (default: 1)
  --mate-index N         Mate slot to populate (default: 1)
  --sample NAME          Restrict to one sample directory name
  -h, --help             Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root)
      ROOT="$2"
      shift 2
      ;;
    --basename)
      REPORT_BASENAME="$2"
      shift 2
      ;;
    --output-root)
      OUTPUT_ROOT="$2"
      shift 2
      ;;
    --stage)
      STAGE="$2"
      shift 2
      ;;
    --max-reads)
      MAX_READS="$2"
      shift 2
      ;;
    --mate-count)
      MATE_COUNT="$2"
      shift 2
      ;;
    --mate-index)
      MATE_INDEX="$2"
      shift 2
      ;;
    --sample)
      SAMPLE_FILTER="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ ! -d "${ROOT}" ]]; then
  echo "ERROR: sample root not found: ${ROOT}" >&2
  exit 1
fi
if [[ ! -d "${OUTPUT_ROOT}" ]]; then
  echo "ERROR: output sample root not found: ${OUTPUT_ROOT}" >&2
  exit 1
fi

if [[ ! -x "${RUNNER}" ]]; then
  echo "ERROR: missing trim-QC runner: ${RUNNER}" >&2
  exit 1
fi

run_id="$(date -u +%Y%m%d_%H%M%S)"
status_path="${OUTPUT_ROOT}/trim_qc_fastq_batch_${run_id}.tsv"
printf 'sample\tstatus\tqc_dir\treport_prefix\tr2_fastqs\tjson\thtml\tsummary\n' > "${status_path}"

mapfile -t sample_dirs < <(find "${ROOT}" -mindepth 1 -maxdepth 1 -type d | sort)

for sample_dir in "${sample_dirs[@]}"; do
  sample="$(basename "${sample_dir}")"
  if [[ -n "${SAMPLE_FILTER}" && "${sample}" != "${SAMPLE_FILTER}" ]]; then
    continue
  fi

  output_sample_dir="${OUTPUT_ROOT}/${sample}"
  if [[ ! -d "${output_sample_dir}" ]]; then
    printf '%s\tmissing_output_sample\t-\t-\t0\t-\t-\t-\n' \
      "${sample}" >> "${status_path}"
    continue
  fi

  qc_dir=""
  if [[ -d "${output_sample_dir}/qc" ]]; then
    qc_dir="${output_sample_dir}/qc"
  elif [[ -d "${output_sample_dir}/QC" ]]; then
    qc_dir="${output_sample_dir}/QC"
  else
    qc_dir="${output_sample_dir}/qc"
    mkdir -p "${qc_dir}"
  fi

  mapfile -t r2_fastqs < <(find "${sample_dir}" -maxdepth 1 -type f -name '*_R2_001.fastq.gz' | sort)
  if [[ "${#r2_fastqs[@]}" -eq 0 ]]; then
    printf '%s\tmissing_r2\t%s\t-\t0\t-\t-\t-\n' \
      "${sample}" "${qc_dir}" >> "${status_path}"
    continue
  fi

  report_prefix="${qc_dir}/${REPORT_BASENAME}"
  summary_path="${qc_dir}/${REPORT_BASENAME}.summary.tsv"
  raw_json_path="${report_prefix}.trim_qc.json"
  raw_html_path="${report_prefix}.trim_qc.html"
  json_path="${report_prefix}.json"
  html_path="${report_prefix}.html"

  input_arg="$(printf '%s,' "${r2_fastqs[@]}")"
  input_arg="${input_arg%,}"

  if "${RUNNER}" \
      --input "${input_arg}" \
      --report "${report_prefix}" \
      --stage "${STAGE}" \
      --mate-count "${MATE_COUNT}" \
      --mate-index "${MATE_INDEX}" \
      --max-reads "${MAX_READS}" \
      > "${summary_path}"; then
    mv -f "${raw_json_path}" "${json_path}"
    mv -f "${raw_html_path}" "${html_path}"
    printf '%s\tdone\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${sample}" "${qc_dir}" "${report_prefix}" "${#r2_fastqs[@]}" \
      "${json_path}" "${html_path}" "${summary_path}" >> "${status_path}"
  else
    printf '%s\tfailed\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${sample}" "${qc_dir}" "${report_prefix}" "${#r2_fastqs[@]}" \
      "${json_path}" "${html_path}" "${summary_path}" >> "${status_path}"
    exit 1
  fi
done

printf 'batch_status\tdone\t%s\n' "${status_path}"
