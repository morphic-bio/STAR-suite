#!/usr/bin/env bash
set -euo pipefail

RUN_ROOT=""
SOURCE_ENDPOINT="${JAX_MULTIOME_GLOBUS_SOURCE_ENDPOINT:-07446cad-33b8-11f0-8c0c-0afffb017b7d}"
DEST_ENDPOINT="${JAX_MULTIOME_GLOBUS_DEST_ENDPOINT:-61fb8b9a-9b52-456e-928c-30c0fb0140bf}"
DEST_ROOT="${JAX_MULTIOME_GLOBUS_DEST_ROOT:-/JAX_Multiome01_processed/large_files}"
SAMPLES=""
WATCH="0"
INTERVAL_SECONDS="300"
TASK_POLL_SECONDS="60"
DRY_RUN="0"
FORCE="0"
DELETE_GENERATED_AFTER_TRANSFER="${JAX_MULTIOME_DELETE_GENERATED_AFTER_TRANSFER:-1}"
INCLUDE_INPUT_FASTQS="1"
INCLUDE_GENERATED_FASTQS="1"
INCLUDE_BAMS="1"
INCLUDE_FRAGMENTS="0"

usage() {
  cat <<'EOF'
Usage:
  upload_jax_multiome01_large_files_globus.sh --run-root PATH [options]

Submits checksum-synced Globus transfers for large JAX_Multiome01 sample files
after each sample has completed. By default this includes input FASTQs from the
production manifest, generated Y/noY FASTQs, and STAR BAMs. The uploader writes
one batch file and one Globus task per completed sample.

Options:
  --run-root PATH            Production output root containing metadata/sample_manifest.tsv
  --source-endpoint UUID     Local Globus endpoint
  --dest-endpoint UUID       Destination Globus endpoint
  --dest-root PATH           Destination root (default: /JAX_Multiome01_processed/large_files)
  --samples CSV             Restrict to sample labels or slugs
  --watch                   Poll until all selected samples are complete and submitted
  --interval-seconds N      Poll interval for --watch (default: 300)
  --task-poll-seconds N     Globus task wait polling interval before deletion (default: 60)
  --dry-run                 Write batch files and validate with Globus dry-run only
  --force                   Resubmit even if a sample is already in the upload state file
  --delete-generated-after-transfer
                            Wait for Globus success, then delete generated local BAM/FASTQ files (default)
  --no-delete-generated-after-transfer
                            Leave generated local BAM/FASTQ files in place after transfer
  --no-input-fastqs         Do not upload raw input FASTQs listed in the manifest
  --no-generated-fastqs     Do not upload STAR-generated Y/noY FASTQs
  --no-bams                 Do not upload BAMs
  --include-fragments       Also upload generated ATAC fragment sidecars/files
  --help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-root) RUN_ROOT="$2"; shift 2 ;;
    --source-endpoint) SOURCE_ENDPOINT="$2"; shift 2 ;;
    --dest-endpoint) DEST_ENDPOINT="$2"; shift 2 ;;
    --dest-root) DEST_ROOT="$2"; shift 2 ;;
    --samples) SAMPLES="$2"; shift 2 ;;
    --watch) WATCH="1"; shift ;;
    --interval-seconds) INTERVAL_SECONDS="$2"; shift 2 ;;
    --task-poll-seconds) TASK_POLL_SECONDS="$2"; shift 2 ;;
    --dry-run) DRY_RUN="1"; shift ;;
    --force) FORCE="1"; shift ;;
    --delete-generated-after-transfer) DELETE_GENERATED_AFTER_TRANSFER="1"; shift ;;
    --no-delete-generated-after-transfer) DELETE_GENERATED_AFTER_TRANSFER="0"; shift ;;
    --no-input-fastqs) INCLUDE_INPUT_FASTQS="0"; shift ;;
    --no-generated-fastqs) INCLUDE_GENERATED_FASTQS="0"; shift ;;
    --no-bams) INCLUDE_BAMS="0"; shift ;;
    --include-fragments) INCLUDE_FRAGMENTS="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "${RUN_ROOT}" ]] || { echo "ERROR: --run-root is required" >&2; exit 1; }
RUN_ROOT="$(realpath "${RUN_ROOT}")"
MANIFEST="${RUN_ROOT}/metadata/sample_manifest.tsv"
[[ -f "${MANIFEST}" ]] || { echo "ERROR: missing manifest ${MANIFEST}" >&2; exit 1; }
command -v globus >/dev/null || { echo "ERROR: globus CLI is not available" >&2; exit 1; }

if ! [[ "${INTERVAL_SECONDS}" =~ ^[0-9]+$ && "${INTERVAL_SECONDS}" -gt 0 ]]; then
  echo "ERROR: --interval-seconds must be a positive integer" >&2
  exit 1
fi
if ! [[ "${TASK_POLL_SECONDS}" =~ ^[0-9]+$ && "${TASK_POLL_SECONDS}" -gt 0 ]]; then
  echo "ERROR: --task-poll-seconds must be a positive integer" >&2
  exit 1
fi

RUN_LABEL="$(basename "${RUN_ROOT}")"
LOG_DIR="${RUN_ROOT}/logs/globus_large_files"
BATCH_DIR="${LOG_DIR}/batches"
STATE_FILE="${LOG_DIR}/upload_state.tsv"
mkdir -p "${BATCH_DIR}"
if [[ ! -f "${STATE_FILE}" ]]; then
  printf 'date_utc\tsample\tsample_slug\ttask_id\tfile_count\tbyte_count\tdelete_count\tstatus\tbatch_file\n' > "${STATE_FILE}"
fi

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "${LOG_DIR}/upload.log"
}

batch_quote() {
  python3 - "$1" <<'PY'
import shlex
import sys
print(shlex.quote(sys.argv[1]))
PY
}

sample_selected() {
  local sample="$1"
  local slug="$2"
  [[ -z "${SAMPLES}" ]] && return 0
  local item
  IFS=',' read -r -a wanted <<< "${SAMPLES}"
  for item in "${wanted[@]}"; do
    [[ "${item}" == "${sample}" || "${item}" == "${slug}" ]] && return 0
  done
  return 1
}

sample_already_submitted() {
  local slug="$1"
  [[ "${FORCE}" == "1" ]] && return 1
  awk -F'\t' -v slug="${slug}" 'NR > 1 && $3 == slug { found = 1 } END { exit found ? 0 : 1 }' "${STATE_FILE}"
}

append_pair() {
  local batch_file="$1"
  local inventory_file="$2"
  local role="$3"
  local source_path="$4"
  local dest_path="$5"
  [[ -f "${source_path}" ]] || return 0
  printf '%s %s\n' "$(batch_quote "${source_path}")" "$(batch_quote "${dest_path}")" >> "${batch_file}"
  printf '%s\t%s\t%s\t%s\n' "$(stat -c '%s' "${source_path}")" "${role}" "${source_path}" "${dest_path}" >> "${inventory_file}"
}

append_csv_fastqs() {
  local batch_file="$1"
  local inventory_file="$2"
  local csv="$3"
  local sample_slug="$4"
  local source_path
  IFS=',' read -r -a fastqs <<< "${csv}"
  for source_path in "${fastqs[@]}"; do
    [[ -n "${source_path}" ]] || continue
    append_pair \
      "${batch_file}" \
      "${inventory_file}" \
      "raw_fastq" \
      "${source_path}" \
      "${DEST_ROOT}/${RUN_LABEL}/raw/${sample_slug}/$(basename "${source_path}")"
  done
}

append_generated_files() {
  local batch_file="$1"
  local inventory_file="$2"
  local sample_dir="$3"
  local sample_slug="$4"
  local role="$5"
  shift 5
  local source_path rel_path
  while IFS= read -r -d '' source_path; do
    rel_path="${source_path#${sample_dir}/}"
    append_pair \
      "${batch_file}" \
      "${inventory_file}" \
      "${role}" \
      "${source_path}" \
      "${DEST_ROOT}/${RUN_LABEL}/samples/${sample_slug}/${rel_path}"
  done < <(find "${sample_dir}" -type f "$@" -print0 2>/dev/null)
}

delete_generated_files() {
  local inventory_file="$1"
  local delete_log="$2"
  local source_path role
  local deleted=0
  : > "${delete_log}"
  while IFS=$'\t' read -r bytes role source_path dest_path; do
    [[ "${bytes}" != "bytes" ]] || continue
    case "${role}" in
      generated_bam|generated_fastq|generated_fragment)
        if [[ -f "${source_path}" ]]; then
          rm -f -- "${source_path}"
          printf '%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${role}" "${source_path}" >> "${delete_log}"
          deleted=$((deleted + 1))
        fi
        ;;
    esac
  done < "${inventory_file}"
  echo "${deleted}"
}

submit_sample() {
  local sample="$1"
  local sample_slug="$2"
  local gex_r1="$3"
  local gex_r2="$4"
  local atac_r1="$5"
  local atac_barcode="$6"
  local atac_r2="$7"

  local sample_dir="${RUN_ROOT}/samples/${sample_slug}"
  local complete_marker="${sample_dir}/mudata/star_chromap_filtered_multiome.h5mu"
  if [[ ! -f "${complete_marker}" ]]; then
    return 2
  fi
  if sample_already_submitted "${sample_slug}"; then
    return 3
  fi

  local batch_file="${BATCH_DIR}/${sample_slug}.batch"
  local inventory_file="${BATCH_DIR}/${sample_slug}.inventory.tsv"
  : > "${batch_file}"
  printf 'bytes\trole\tsource\tdestination\n' > "${inventory_file}"

  if [[ "${INCLUDE_INPUT_FASTQS}" == "1" ]]; then
    append_csv_fastqs "${batch_file}" "${inventory_file}" "${gex_r1}" "${sample_slug}"
    append_csv_fastqs "${batch_file}" "${inventory_file}" "${gex_r2}" "${sample_slug}"
    append_csv_fastqs "${batch_file}" "${inventory_file}" "${atac_r1}" "${sample_slug}"
    append_csv_fastqs "${batch_file}" "${inventory_file}" "${atac_barcode}" "${sample_slug}"
    append_csv_fastqs "${batch_file}" "${inventory_file}" "${atac_r2}" "${sample_slug}"
  fi

  if [[ "${INCLUDE_BAMS}" == "1" ]]; then
    append_generated_files "${batch_file}" "${inventory_file}" "${sample_dir}" "${sample_slug}" "generated_bam" -name '*.bam'
  fi
  if [[ "${INCLUDE_GENERATED_FASTQS}" == "1" ]]; then
    append_generated_files "${batch_file}" "${inventory_file}" "${sample_dir}" "${sample_slug}" "generated_fastq" -path '*/y_separated/*' -name '*.fastq.gz'
  fi
  if [[ "${INCLUDE_FRAGMENTS}" == "1" ]]; then
    append_generated_files "${batch_file}" "${inventory_file}" "${sample_dir}" "${sample_slug}" "generated_fragment" \( -name 'atac_fragments.tsv.gz' -o -name 'atac_fragments.bin' -o -name 'atac_fragments.bin.chroms.tsv' \)
  fi

  local file_count
  file_count="$(grep -vc '^[[:space:]]*$' "${batch_file}" || true)"
  if [[ "${file_count}" == "0" ]]; then
    log "No transferable large files found for ${sample}"
    return 4
  fi
  local byte_count
  byte_count="$(awk -F'\t' 'NR > 1 { total += $1 } END { printf "%.0f", total }' "${inventory_file}")"
  local generated_count
  generated_count="$(awk -F'\t' 'NR > 1 && $2 != "raw_fastq" { total += 1 } END { printf "%.0f", total }' "${inventory_file}")"
  local label="JAX_Multiome01 ${RUN_LABEL} ${sample_slug} large files"

  if [[ "${DRY_RUN}" == "1" ]]; then
    globus transfer \
      "${SOURCE_ENDPOINT}" \
      "${DEST_ENDPOINT}" \
      --batch "${batch_file}" \
      --sync-level checksum \
      --label "${label}" \
      --dry-run > "${BATCH_DIR}/${sample_slug}.dry_run.txt"
    log "DRY-RUN ${sample}: ${file_count} files, ${byte_count} bytes; batch=${batch_file}"
    return 0
  fi

  local submit_json task_id delete_count status
  submit_json="$(globus transfer \
    "${SOURCE_ENDPOINT}" \
    "${DEST_ENDPOINT}" \
    --batch "${batch_file}" \
    --sync-level checksum \
    --label "${label}" \
    --format json)"
  printf '%s\n' "${submit_json}" > "${BATCH_DIR}/${sample_slug}.submit.json"
  task_id="$(python3 -c 'import json,sys; print(json.load(sys.stdin).get("task_id",""))' <<< "${submit_json}")"
  [[ -n "${task_id}" ]] || { echo "ERROR: Globus submission did not return a task_id for ${sample}" >&2; exit 1; }
  log "Submitted ${sample}: task=${task_id}, files=${file_count}, bytes=${byte_count}"
  delete_count=0
  status="submitted"
  if [[ "${DELETE_GENERATED_AFTER_TRANSFER}" == "1" && "${generated_count}" -gt 0 ]]; then
    log "Waiting for Globus task ${task_id} before deleting ${generated_count} generated local files for ${sample}"
    globus task wait --polling-interval "${TASK_POLL_SECONDS}" "${task_id}" > "${BATCH_DIR}/${sample_slug}.task_wait.log"
    delete_count="$(delete_generated_files "${inventory_file}" "${BATCH_DIR}/${sample_slug}.deleted.tsv")"
    status="transferred_deleted_generated"
    log "Deleted ${delete_count} generated local files for ${sample}; raw input FASTQs were preserved"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    "${sample}" \
    "${sample_slug}" \
    "${task_id}" \
    "${file_count}" \
    "${byte_count}" \
    "${delete_count}" \
    "${status}" \
    "${batch_file}" >> "${STATE_FILE}"
  return 0
}

scan_once() {
  local complete_count=0
  local selected_count=0
  local submitted_or_skipped=0
  local pending_count=0
  local rc
  while IFS=$'\t' read -r sample sample_slug atac_lib gex_lib gex_r1 gex_r2 atac_r1 atac_barcode atac_r2 gex_run_ids atac_run_ids; do
    [[ "${sample}" != "sample" ]] || continue
    if ! sample_selected "${sample}" "${sample_slug}"; then
      continue
    fi
    selected_count=$((selected_count + 1))
    set +e
    submit_sample "${sample}" "${sample_slug}" "${gex_r1}" "${gex_r2}" "${atac_r1}" "${atac_barcode}" "${atac_r2}"
    rc=$?
    set -e
    case "${rc}" in
      0|3)
        complete_count=$((complete_count + 1))
        submitted_or_skipped=$((submitted_or_skipped + 1))
        ;;
      2)
        pending_count=$((pending_count + 1))
        ;;
      4)
        complete_count=$((complete_count + 1))
        ;;
      *)
        return "${rc}"
        ;;
    esac
  done < "${MANIFEST}"
  log "Scan summary: selected=${selected_count}, complete=${complete_count}, submitted_or_previously_submitted=${submitted_or_skipped}, pending=${pending_count}"
  [[ "${selected_count}" -gt 0 ]] || { echo "ERROR: no samples selected" >&2; return 1; }
  [[ "${pending_count}" == "0" ]]
}

log "Starting Globus large-file uploader for ${RUN_ROOT}"
if [[ "${WATCH}" == "1" ]]; then
  while true; do
    if scan_once; then
      log "All selected samples are complete and submitted"
      exit 0
    fi
    sleep "${INTERVAL_SECONDS}"
  done
else
  scan_once || true
fi
