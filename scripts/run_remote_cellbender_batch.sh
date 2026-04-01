#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
HELPER="${REPO_ROOT}/scripts/run_remote_cellbender_rsync.sh"

SAMPLES_ROOT=""
REMOTE_HOST=""
REMOTE_ROOT=""
OUTPUT_NAME="downstream_genefull_velocyto_cellbender"
CELLBENDER_IMAGE="biodepot/cellbender:0.3.2"
CELLBENDER_CPU_CORES="8"
GPU_SLOTS="0,1,0,1"
NO_SYNC_IMAGE="1"

usage() {
  cat <<'EOF'
Usage:
  run_remote_cellbender_batch.sh --samples-root PATH --remote-host HOST --remote-root PATH [options]

Runs remote CellBender jobs in fixed-size batches, one job per GPU slot.
For example, GPU slots "0,1,0,1" runs 4 jobs at a time: 2 on GPU 0 and 2 on GPU 1.

Options:
  --samples-root PATH        Root containing per-sample directories
  --remote-host HOST         SSH target
  --remote-root PATH         Remote staging root
  --output-name NAME         Downstream output dir name
                             (default: downstream_genefull_velocyto_cellbender)
  --cellbender-image IMG     CellBender image (default: biodepot/cellbender:0.3.2)
  --cellbender-cpu-cores N   CPU threads per job (default: 8)
  --gpu-slots CSV            Comma-separated GPU assignment slots
                             (default: 0,1,0,1)
  --sync-image               Sync local image instead of using remote pull
  --help                     Show help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples-root)
      SAMPLES_ROOT="$2"
      shift 2
      ;;
    --remote-host)
      REMOTE_HOST="$2"
      shift 2
      ;;
    --remote-root)
      REMOTE_ROOT="$2"
      shift 2
      ;;
    --output-name)
      OUTPUT_NAME="$2"
      shift 2
      ;;
    --cellbender-image)
      CELLBENDER_IMAGE="$2"
      shift 2
      ;;
    --cellbender-cpu-cores)
      CELLBENDER_CPU_CORES="$2"
      shift 2
      ;;
    --gpu-slots)
      GPU_SLOTS="$2"
      shift 2
      ;;
    --sync-image)
      NO_SYNC_IMAGE="0"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

[[ -n "${SAMPLES_ROOT}" ]] || { echo "ERROR: --samples-root is required" >&2; exit 1; }
[[ -n "${REMOTE_HOST}" ]] || { echo "ERROR: --remote-host is required" >&2; exit 1; }
[[ -n "${REMOTE_ROOT}" ]] || { echo "ERROR: --remote-root is required" >&2; exit 1; }
[[ -x "${HELPER}" ]] || { echo "ERROR: missing helper ${HELPER}" >&2; exit 1; }

SAMPLES_ROOT="$(realpath "${SAMPLES_ROOT}")"
STAMP="$(date +%Y%m%d_%H%M%S)"
RUNS_TSV="${SAMPLES_ROOT}/REMOTE_CELLBENDER_BATCH_${STAMP}.tsv"
printf "sample\tgpu_device\tstatus\tlog\n" > "${RUNS_TSV}"

IFS=',' read -r -a GPU_SLOT_ARRAY <<< "${GPU_SLOTS}"
(( ${#GPU_SLOT_ARRAY[@]} > 0 )) || { echo "ERROR: --gpu-slots must not be empty" >&2; exit 1; }

mapfile -t SAMPLE_DIRS < <(find "${SAMPLES_ROOT}" -mindepth 1 -maxdepth 1 -type d | sort)

ELIGIBLE_SAMPLES=()
for sample_dir in "${SAMPLE_DIRS[@]}"; do
  sample="$(basename "${sample_dir}")"
  [[ "${sample}" == downstream_logs_* ]] && continue
  output_dir="${sample_dir}/${OUTPUT_NAME}"
  if [[ -f "${output_dir}/cellbender/cellbender_counts.h5" ]]; then
    printf "%s\t-\tskipped_existing\t-\n" "${sample}" >> "${RUNS_TSV}"
    continue
  fi
  if [[ ! -f "${output_dir}/unfiltered_counts.h5ad" ]]; then
    printf "%s\t-\tskipped_missing_unfiltered\t-\n" "${sample}" >> "${RUNS_TSV}"
    continue
  fi
  ELIGIBLE_SAMPLES+=("${sample}")
done

run_batch() {
  local -a batch_samples=("$@")
  local -a pids=()
  local -a pid_samples=()
  local -a pid_logs=()
  local idx=0

  for sample in "${batch_samples[@]}"; do
    local sample_dir="${SAMPLES_ROOT}/${sample}"
    local output_dir="${sample_dir}/${OUTPUT_NAME}"
    local gpu_device="${GPU_SLOT_ARRAY[$idx]}"
    local log_file="${output_dir}/remote_batch_${STAMP}_gpu${gpu_device}.log"
    local -a args=(
      "${HELPER}"
      --downstream-dir "${output_dir}"
      --remote-host "${REMOTE_HOST}"
      --remote-root "${REMOTE_ROOT}"
      --cellbender-image "${CELLBENDER_IMAGE}"
      --cellbender-cpu-cores "${CELLBENDER_CPU_CORES}"
      --cellbender-gpu
      --cellbender-gpu-device "${gpu_device}"
      --local-log "${log_file}"
    )
    if [[ "${NO_SYNC_IMAGE}" == "1" ]]; then
      args+=(--no-sync-image)
    fi

    (
      "${args[@]}"
    ) > "${log_file}.launcher" 2>&1 &
    pids+=("$!")
    pid_samples+=("${sample}")
    pid_logs+=("${log_file}")
    printf "%s\t%s\tlaunched\t%s\n" "${sample}" "${gpu_device}" "${log_file}" >> "${RUNS_TSV}"
    idx=$((idx + 1))
  done

  for i in "${!pids[@]}"; do
    if wait "${pids[$i]}"; then
      printf "%s\t-\tdone\t%s\n" "${pid_samples[$i]}" "${pid_logs[$i]}" >> "${RUNS_TSV}"
    else
      status=$?
      printf "%s\t-\tfailed(%s)\t%s\n" "${pid_samples[$i]}" "${status}" "${pid_logs[$i]}" >> "${RUNS_TSV}"
    fi
  done
}

batch_size="${#GPU_SLOT_ARRAY[@]}"
for (( start=0; start<${#ELIGIBLE_SAMPLES[@]}; start+=batch_size )); do
  batch=("${ELIGIBLE_SAMPLES[@]:start:batch_size}")
  run_batch "${batch[@]}"
done

echo "PASS: remote CellBender batch complete"
echo "Run ledger: ${RUNS_TSV}"
