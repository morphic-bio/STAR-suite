#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
WRAPPER="${REPO_ROOT}/scripts/run_scrna_downstream_gene_full_velocyto.sh"

SAMPLES_ROOT=""
OUTPUT_NAME="downstream_genefull_velocyto_cellbender"
RUN_CELLBENDER="0"
ADAPTIVE_FILTER="1"
MIN_GENES=""
MAX_GENES=""
MT_PCT_CUTOFF=""
N_MAD=""
CELLBENDER_CPU_CORES=""
CELLBENDER_GPU="0"
SKIP_EXISTING_CELLBENDER="0"
SKIP_EXISTING_OUTPUT="0"
EXTRA_ARGS=()

usage() {
  cat <<'EOF'
Usage:
  run_ucsf_downstream_batch.sh --samples-root <samples-root> [options]

Options:
  --samples-root PATH       Root containing per-sample directories with run/
  --output-name NAME        Per-sample output dir name
                            (default: downstream_genefull_velocyto_cellbender)
  --adaptive-filter         Enable adaptive n_genes and MT percentage filtering
  --run-cellbender          Enable CellBender
  --min-genes INT           Pass through to downstream wrapper
  --max-genes INT           Pass through to downstream wrapper
  --mt-pct-cutoff FLOAT     Pass through to downstream wrapper
  --n-mad FLOAT             Pass through to downstream wrapper
  --cellbender-cpu-cores N  Pass through to downstream wrapper
  --cellbender-gpu          Pass through to downstream wrapper
  --skip-existing-cellbender
                            Skip samples that already have cellbender/cellbender_counts.h5
  --skip-existing-output
                            Skip samples that already have summary.txt in the output dir
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples-root)
      SAMPLES_ROOT="$2"
      shift 2
      ;;
    --output-name)
      OUTPUT_NAME="$2"
      shift 2
      ;;
    --adaptive-filter)
      ADAPTIVE_FILTER="1"
      shift
      ;;
    --run-cellbender)
      RUN_CELLBENDER="1"
      shift
      ;;
    --min-genes)
      MIN_GENES="$2"
      shift 2
      ;;
    --max-genes)
      MAX_GENES="$2"
      shift 2
      ;;
    --mt-pct-cutoff)
      MT_PCT_CUTOFF="$2"
      shift 2
      ;;
    --n-mad)
      N_MAD="$2"
      shift 2
      ;;
    --cellbender-cpu-cores)
      CELLBENDER_CPU_CORES="$2"
      shift 2
      ;;
    --cellbender-gpu)
      CELLBENDER_GPU="1"
      shift
      ;;
    --skip-existing-cellbender)
      SKIP_EXISTING_CELLBENDER="1"
      shift
      ;;
    --skip-existing-output)
      SKIP_EXISTING_OUTPUT="1"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      EXTRA_ARGS+=("$1")
      shift
      ;;
  esac
done

[[ -n "${SAMPLES_ROOT}" ]] || { echo "ERROR: --samples-root is required" >&2; exit 1; }
[[ -x "${WRAPPER}" ]] || { echo "ERROR: missing executable wrapper ${WRAPPER}" >&2; exit 1; }

SAMPLES_ROOT="$(realpath "${SAMPLES_ROOT}")"
RUNS_TSV="${SAMPLES_ROOT}/DOWNSTREAM_RUNS.tsv"
LOG_ROOT="${SAMPLES_ROOT}/downstream_logs_${OUTPUT_NAME}"
mkdir -p "${LOG_ROOT}"
printf "sample\tstatus\toutput_dir\tlog\n" > "${RUNS_TSV}"

mapfile -t SAMPLE_DIRS < <(find "${SAMPLES_ROOT}" -mindepth 1 -maxdepth 1 -type d | sort)

for sample_dir in "${SAMPLE_DIRS[@]}"; do
  sample="$(basename "${sample_dir}")"
  run_dir="${sample_dir}/run"
  output_dir="${sample_dir}/${OUTPUT_NAME}"
  log_file="${LOG_ROOT}/${sample}.log"
  cellbender_h5="${output_dir}/cellbender/cellbender_counts.h5"
  output_summary="${output_dir}/summary.txt"

  if [[ ! -d "${run_dir}" ]]; then
    printf "%s\tmissing_run_dir\t%s\t%s\n" "${sample}" "${output_dir}" "${log_file}" >> "${RUNS_TSV}"
    continue
  fi
  if [[ "${SKIP_EXISTING_CELLBENDER}" == "1" && -f "${cellbender_h5}" ]]; then
    printf "%s\tskipped_existing_cellbender\t%s\t%s\n" "${sample}" "${output_dir}" "${log_file}" >> "${RUNS_TSV}"
    continue
  fi
  if [[ "${SKIP_EXISTING_OUTPUT}" == "1" && -f "${output_summary}" ]]; then
    printf "%s\tskipped_existing_output\t%s\t%s\n" "${sample}" "${output_dir}" "${log_file}" >> "${RUNS_TSV}"
    continue
  fi

  args=(--run-dir "${run_dir}" --output-dir "${output_dir}")
  [[ "${ADAPTIVE_FILTER}" == "1" ]] && args+=(--adaptive-filter)
  [[ "${RUN_CELLBENDER}" == "1" ]] && args+=(--run-cellbender)
  [[ -n "${MIN_GENES}" ]] && args+=(--min-genes "${MIN_GENES}")
  [[ -n "${MAX_GENES}" ]] && args+=(--max-genes "${MAX_GENES}")
  [[ -n "${MT_PCT_CUTOFF}" ]] && args+=(--mt-pct-cutoff "${MT_PCT_CUTOFF}")
  [[ -n "${N_MAD}" ]] && args+=(--n-mad "${N_MAD}")
  [[ -n "${CELLBENDER_CPU_CORES}" ]] && args+=(--cellbender-cpu-cores "${CELLBENDER_CPU_CORES}")
  [[ "${CELLBENDER_GPU}" == "1" ]] && args+=(--cellbender-gpu)
  if (( ${#EXTRA_ARGS[@]} > 0 )); then
    args+=("${EXTRA_ARGS[@]}")
  fi

  if "${WRAPPER}" "${args[@]}" >"${log_file}" 2>&1; then
    printf "%s\tdone\t%s\t%s\n" "${sample}" "${output_dir}" "${log_file}" >> "${RUNS_TSV}"
  else
    status=$?
    printf "%s\tfailed(%s)\t%s\t%s\n" "${sample}" "${status}" "${output_dir}" "${log_file}" >> "${RUNS_TSV}"
  fi
done
