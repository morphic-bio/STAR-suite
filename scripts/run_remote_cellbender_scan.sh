#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DOWNSTREAM_WRAPPER="${REPO_ROOT}/scripts/run_scrna_downstream_gene_full_velocyto.sh"
REMOTE_HELPER="${REPO_ROOT}/scripts/run_remote_cellbender_rsync.sh"

SAMPLES_ROOT=""
REMOTE_HOST=""
REMOTE_ROOT=""
OUTPUT_NAME="downstream_genefull_velocyto_cellbender"
CELLBENDER_IMAGE="biodepot/cellbender:0.3.2"
CELLBENDER_CPU_CORES="8"
CELLBENDER_LAYER="denoised"
CELLBENDER_FLAGS=""
GPU_SLOTS="0,1,0,1"
POLL_SECONDS="60"
EXPECTED_SAMPLES=""
NO_SYNC_IMAGE="1"
DRY_RUN="0"
ONCE="0"
ADAPTIVE_FILTER="1"
MIN_GENES=""
MAX_GENES=""
MT_PCT_CUTOFF=""
N_MAD=""

usage() {
  cat <<'EOF'
Usage:
  run_remote_cellbender_scan.sh --samples-root PATH --remote-host HOST --remote-root PATH [options]

Continuously scans a STAR samples root, prepares downstream h5ad outputs for
samples whose STAR stage has completed, and then launches remote CellBender
jobs for samples that do not already have completed CellBender outputs.

Options:
  --samples-root PATH        Root containing per-sample directories with run/
  --remote-host HOST         SSH target for remote CellBender
  --remote-root PATH         Remote local-disk staging root
  --output-name NAME         Downstream output dir name
                             (default: downstream_genefull_velocyto_cellbender)
  --cellbender-image IMG     CellBender image (default: biodepot/cellbender:0.3.2)
  --cellbender-cpu-cores N   CPU threads per remote CellBender job (default: 8)
  --cellbender-layer NAME    Denoised layer name (default: cellbender)
  --cellbender-flags STR     Extra CellBender flags
  --gpu-slots CSV            Comma-separated GPU slots for concurrent jobs
                             (default: 0,1,0,1)
  --poll-seconds N           Scan interval in seconds (default: 60)
  --expected-samples N       Exit once at least N sample dirs exist and all N
                             have completed CellBender
  --adaptive-filter          Pass through adaptive n_genes and MT percentage QC
  --min-genes INT            Pass through to local downstream prep
  --max-genes INT            Pass through to local downstream prep
  --mt-pct-cutoff FLOAT      Pass through to local downstream prep
  --n-mad FLOAT              Pass through to local downstream prep
  --sync-image               Sync local CellBender image instead of using remote pull
  --no-sync-image            Use the remote image as-is; do not sync local image
  --dry-run                  Report what would launch without launching
  --once                     Run one scan iteration and exit
  --help                     Show this help
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
    --cellbender-layer)
      CELLBENDER_LAYER="$2"
      shift 2
      ;;
    --cellbender-flags)
      CELLBENDER_FLAGS="$2"
      shift 2
      ;;
    --gpu-slots)
      GPU_SLOTS="$2"
      shift 2
      ;;
    --poll-seconds)
      POLL_SECONDS="$2"
      shift 2
      ;;
    --expected-samples)
      EXPECTED_SAMPLES="$2"
      shift 2
      ;;
    --adaptive-filter)
      ADAPTIVE_FILTER="1"
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
    --sync-image)
      NO_SYNC_IMAGE="0"
      shift
      ;;
    --no-sync-image)
      NO_SYNC_IMAGE="1"
      shift
      ;;
    --dry-run)
      DRY_RUN="1"
      shift
      ;;
    --once)
      ONCE="1"
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
[[ -x "${DOWNSTREAM_WRAPPER}" ]] || { echo "ERROR: missing downstream wrapper ${DOWNSTREAM_WRAPPER}" >&2; exit 1; }
[[ -x "${REMOTE_HELPER}" ]] || { echo "ERROR: missing remote helper ${REMOTE_HELPER}" >&2; exit 1; }
command -v ssh >/dev/null 2>&1 || { echo "ERROR: ssh is required" >&2; exit 1; }
command -v rsync >/dev/null 2>&1 || { echo "ERROR: rsync is required" >&2; exit 1; }

SAMPLES_ROOT="$(realpath "${SAMPLES_ROOT}")"
STAMP="$(date +%Y%m%d_%H%M%S)"
WATCH_LOG="${SAMPLES_ROOT}/REMOTE_CELLBENDER_SCAN_${STAMP}.log"
WATCH_TSV="${SAMPLES_ROOT}/REMOTE_CELLBENDER_SCAN_${STAMP}.tsv"
PREP_LOG_ROOT="${SAMPLES_ROOT}/downstream_prep_logs_${OUTPUT_NAME}"
mkdir -p "${PREP_LOG_ROOT}"
printf "timestamp\tsample\tstage\tstatus\tdetail\n" > "${WATCH_TSV}"

IFS=',' read -r -a GPU_SLOT_ARRAY <<< "${GPU_SLOTS}"
(( ${#GPU_SLOT_ARRAY[@]} > 0 )) || { echo "ERROR: --gpu-slots must not be empty" >&2; exit 1; }

log() {
  local ts
  ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "[${ts}] $*" | tee -a "${WATCH_LOG}"
}

record_tsv() {
  local sample="$1"
  local stage="$2"
  local status="$3"
  local detail="$4"
  printf "%s\t%s\t%s\t%s\t%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${sample}" "${stage}" "${status}" "${detail}" >> "${WATCH_TSV}"
}

sample_dirs() {
  find "${SAMPLES_ROOT}" -mindepth 1 -maxdepth 1 -type d ! -name 'downstream_*' | sort
}

star_complete() {
  local sample_dir="$1"
  local run_dir="${sample_dir}/run"
  [[ -f "${run_dir}/Solo.out/GeneFull/filtered/barcodes.tsv" ]] || return 1
  if [[ -f "${sample_dir}/pf_multi_config.csv" || -d "${run_dir}/cr_assign" ]]; then
    [[ -f "${run_dir}/outs/crispr_analysis/protospacer_calls_summary.csv" ]] || return 1
  fi
}

downstream_dir_for() {
  local sample_dir="$1"
  printf "%s/%s" "${sample_dir}" "${OUTPUT_NAME}"
}

downstream_prepared() {
  local sample_dir="$1"
  local output_dir
  output_dir="$(downstream_dir_for "${sample_dir}")"
  [[ -f "${output_dir}/counts.h5ad" ]] || return 1
  [[ -f "${output_dir}/unfiltered_counts.h5ad" ]] || return 1
  [[ -f "${output_dir}/filtered_counts.h5ad" ]] || return 1
  [[ -f "${output_dir}/default_singlet_filtered_counts.h5ad" ]] || return 1
  [[ -f "${output_dir}/summary.txt" ]] || return 1
}

cellbender_complete() {
  local sample_dir="$1"
  local output_dir
  output_dir="$(downstream_dir_for "${sample_dir}")"
  if [[ -f "${output_dir}/cellbender/cellbender_counts.h5" ]]; then
    return 0
  fi
  [[ -f "${output_dir}/cellbender/CELLBENDER_FAILED.txt" ]]
}

marker_pid_alive() {
  local marker="$1"
  [[ -f "${marker}" ]] || return 1
  local pid
  pid="$(head -n 1 "${marker}" 2>/dev/null || true)"
  [[ "${pid}" =~ ^[0-9]+$ ]] || return 1
  kill -0 "${pid}" 2>/dev/null
}

cleanup_stale_marker() {
  local marker="$1"
  if [[ -f "${marker}" ]] && ! marker_pid_alive "${marker}"; then
    rm -f "${marker}"
  fi
}

prep_marker_for() {
  local sample_dir="$1"
  local output_dir
  output_dir="$(downstream_dir_for "${sample_dir}")"
  printf "%s/.remote_scan_prepare.pid" "${output_dir}"
}

remote_marker_for() {
  local sample_dir="$1"
  local output_dir
  output_dir="$(downstream_dir_for "${sample_dir}")"
  printf "%s/.remote_scan_cellbender.pid" "${output_dir}"
}

prep_active() {
  local sample_dir="$1"
  local marker
  marker="$(prep_marker_for "${sample_dir}")"
  cleanup_stale_marker "${marker}"
  marker_pid_alive "${marker}"
}

remote_active() {
  local sample_dir="$1"
  local marker
  marker="$(remote_marker_for "${sample_dir}")"
  cleanup_stale_marker "${marker}"
  marker_pid_alive "${marker}"
}

declare -a SLOT_PIDS
declare -a SLOT_SAMPLES
declare -a SLOT_LOGS
declare -a SLOT_MARKERS
declare -a SLOT_DEVICES
for _ in "${GPU_SLOT_ARRAY[@]}"; do
  SLOT_PIDS+=("")
  SLOT_SAMPLES+=("")
  SLOT_LOGS+=("")
  SLOT_MARKERS+=("")
  SLOT_DEVICES+=("")
done
PREP_PID=""
PREP_SAMPLE=""
PREP_LOG=""
PREP_MARKER=""

refresh_prep_state() {
  [[ -n "${PREP_PID}" ]] || return 0
  if kill -0 "${PREP_PID}" 2>/dev/null; then
    return 0
  fi
  local status=0
  if ! wait "${PREP_PID}"; then
    status=$?
  fi
  rm -f "${PREP_MARKER}"
  if [[ "${status}" -eq 0 ]]; then
    log "Prepared downstream ${PREP_SAMPLE}"
    record_tsv "${PREP_SAMPLE}" "prepare" "done" "${PREP_LOG}"
  else
    log "Preparation failed for ${PREP_SAMPLE} (exit ${status})"
    record_tsv "${PREP_SAMPLE}" "prepare" "failed(${status})" "${PREP_LOG}"
  fi
  PREP_PID=""
  PREP_SAMPLE=""
  PREP_LOG=""
  PREP_MARKER=""
}

refresh_remote_state() {
  local i pid status sample marker log_file
  for i in "${!SLOT_PIDS[@]}"; do
    pid="${SLOT_PIDS[$i]}"
    [[ -n "${pid}" ]] || continue
    if kill -0 "${pid}" 2>/dev/null; then
      continue
    fi
    status=0
    if ! wait "${pid}"; then
      status=$?
    fi
    sample="${SLOT_SAMPLES[$i]}"
    marker="${SLOT_MARKERS[$i]}"
    log_file="${SLOT_LOGS[$i]}"
    rm -f "${marker}"
    if [[ "${status}" -eq 0 ]]; then
      log "Remote CellBender finished for ${sample} on GPU ${SLOT_DEVICES[$i]}"
      record_tsv "${sample}" "cellbender" "done" "${log_file}"
    else
      log "Remote CellBender failed for ${sample} on GPU ${SLOT_DEVICES[$i]} (exit ${status})"
      record_tsv "${sample}" "cellbender" "failed(${status})" "${log_file}"
    fi
    SLOT_PIDS[$i]=""
    SLOT_SAMPLES[$i]=""
    SLOT_LOGS[$i]=""
    SLOT_MARKERS[$i]=""
    SLOT_DEVICES[$i]=""
  done
}

launch_prep() {
  local sample="$1"
  local sample_dir="$2"
  local run_dir="${sample_dir}/run"
  local output_dir
  output_dir="$(downstream_dir_for "${sample_dir}")"
  mkdir -p "${output_dir}"
  local log_file="${PREP_LOG_ROOT}/${sample}.log"
  local marker
  marker="$(prep_marker_for "${sample_dir}")"
  local args=(--run-dir "${run_dir}" --output-dir "${output_dir}")
  [[ "${ADAPTIVE_FILTER}" == "1" ]] && args+=(--adaptive-filter)
  [[ -n "${MIN_GENES}" ]] && args+=(--min-genes "${MIN_GENES}")
  [[ -n "${MAX_GENES}" ]] && args+=(--max-genes "${MAX_GENES}")
  [[ -n "${MT_PCT_CUTOFF}" ]] && args+=(--mt-pct-cutoff "${MT_PCT_CUTOFF}")
  [[ -n "${N_MAD}" ]] && args+=(--n-mad "${N_MAD}")

  if [[ "${DRY_RUN}" == "1" ]]; then
    log "DRY RUN: would prepare downstream ${sample}"
    record_tsv "${sample}" "prepare" "dry_run" "${log_file}"
    return 0
  fi

  (
    "${DOWNSTREAM_WRAPPER}" "${args[@]}"
  ) > "${log_file}" 2>&1 &
  PREP_PID="$!"
  PREP_SAMPLE="${sample}"
  PREP_LOG="${log_file}"
  PREP_MARKER="${marker}"
  printf "%s\n" "${PREP_PID}" > "${marker}"
  log "Started downstream prep for ${sample}"
  record_tsv "${sample}" "prepare" "launched" "${log_file}"
}

launch_remote() {
  local slot_idx="$1"
  local sample="$2"
  local sample_dir="$3"
  local output_dir
  output_dir="$(downstream_dir_for "${sample_dir}")"
  local gpu_device="${GPU_SLOT_ARRAY[$slot_idx]}"
  local log_file="${output_dir}/remote_scan_${STAMP}_gpu${gpu_device}.log"
  local marker
  marker="$(remote_marker_for "${sample_dir}")"

  if [[ "${DRY_RUN}" == "1" ]]; then
    log "DRY RUN: would launch remote CellBender for ${sample} on GPU ${gpu_device}"
    record_tsv "${sample}" "cellbender" "dry_run" "${log_file}"
    return 0
  fi

  local -a args=(
    "${REMOTE_HELPER}"
    --downstream-dir "${output_dir}"
    --remote-host "${REMOTE_HOST}"
    --remote-root "${REMOTE_ROOT}"
    --cellbender-image "${CELLBENDER_IMAGE}"
    --cellbender-cpu-cores "${CELLBENDER_CPU_CORES}"
    --cellbender-layer "${CELLBENDER_LAYER}"
    --cellbender-gpu
    --cellbender-gpu-device "${gpu_device}"
    --local-log "${log_file}"
  )
  [[ -n "${CELLBENDER_FLAGS}" ]] && args+=(--cellbender-flags "${CELLBENDER_FLAGS}")
  [[ "${NO_SYNC_IMAGE}" == "1" ]] && args+=(--no-sync-image)

  (
    "${args[@]}"
  ) > "${log_file}.launcher" 2>&1 &
  SLOT_PIDS[$slot_idx]="$!"
  SLOT_SAMPLES[$slot_idx]="${sample}"
  SLOT_LOGS[$slot_idx]="${log_file}"
  SLOT_MARKERS[$slot_idx]="${marker}"
  SLOT_DEVICES[$slot_idx]="${gpu_device}"
  printf "%s\n" "${SLOT_PIDS[$slot_idx]}" > "${marker}"
  log "Started remote CellBender for ${sample} on GPU ${gpu_device}"
  record_tsv "${sample}" "cellbender" "launched" "${log_file}"
}

sample_count() {
  sample_dirs | wc -l | tr -d ' '
}

all_expected_complete() {
  local dirs=()
  mapfile -t dirs < <(sample_dirs)
  [[ -n "${EXPECTED_SAMPLES}" ]] || return 1
  (( ${#dirs[@]} >= EXPECTED_SAMPLES )) || return 1
  local complete=0
  local sample_dir
  for sample_dir in "${dirs[@]}"; do
    if cellbender_complete "${sample_dir}"; then
      ((complete += 1))
    fi
  done
  (( complete >= EXPECTED_SAMPLES ))
}

log "Starting remote CellBender scan"
log "Samples root: ${SAMPLES_ROOT}"
log "Remote host: ${REMOTE_HOST}"
log "Remote root: ${REMOTE_ROOT}"
log "Output name: ${OUTPUT_NAME}"
log "GPU slots: ${GPU_SLOTS}"
[[ -n "${EXPECTED_SAMPLES}" ]] && log "Expected samples: ${EXPECTED_SAMPLES}"
[[ "${DRY_RUN}" == "1" ]] && log "Dry run enabled"

while true; do
  refresh_prep_state
  refresh_remote_state

  mapfile -t SAMPLE_DIRS < <(sample_dirs)
  star_ready=0
  prepared=0
  cb_done=0
  for sample_dir in "${SAMPLE_DIRS[@]}"; do
    star_complete "${sample_dir}" && ((star_ready += 1)) || true
    downstream_prepared "${sample_dir}" && ((prepared += 1)) || true
    cellbender_complete "${sample_dir}" && ((cb_done += 1)) || true
  done
  log "Scan status: sample_dirs=${#SAMPLE_DIRS[@]} star_ready=${star_ready} downstream_prepared=${prepared} cellbender_done=${cb_done}"

  if [[ -z "${PREP_PID}" ]]; then
    for sample_dir in "${SAMPLE_DIRS[@]}"; do
      sample="$(basename "${sample_dir}")"
      star_complete "${sample_dir}" || continue
      downstream_prepared "${sample_dir}" && continue
      cellbender_complete "${sample_dir}" && continue
      prep_active "${sample_dir}" && continue
      remote_active "${sample_dir}" && continue
      launch_prep "${sample}" "${sample_dir}"
      break
    done
  fi

  for slot_idx in "${!GPU_SLOT_ARRAY[@]}"; do
    [[ -z "${SLOT_PIDS[$slot_idx]}" ]] || continue
    for sample_dir in "${SAMPLE_DIRS[@]}"; do
      sample="$(basename "${sample_dir}")"
      star_complete "${sample_dir}" || continue
      downstream_prepared "${sample_dir}" || continue
      cellbender_complete "${sample_dir}" && continue
      prep_active "${sample_dir}" && continue
      remote_active "${sample_dir}" && continue
      if [[ "${PREP_SAMPLE}" == "${sample}" ]]; then
        continue
      fi
      launch_remote "${slot_idx}" "${sample}" "${sample_dir}"
      break
    done
  done

  if [[ "${ONCE}" == "1" ]]; then
    break
  fi

  if all_expected_complete && [[ -z "${PREP_PID}" ]]; then
    active_remote=0
    for pid in "${SLOT_PIDS[@]}"; do
      [[ -n "${pid}" ]] && ((active_remote += 1))
    done
    if (( active_remote == 0 )); then
      log "All expected samples completed"
      break
    fi
  fi

  sleep "${POLL_SECONDS}"
done

refresh_prep_state
refresh_remote_state
log "Watcher exiting"
echo "Watch log: ${WATCH_LOG}"
echo "Watch ledger: ${WATCH_TSV}"
