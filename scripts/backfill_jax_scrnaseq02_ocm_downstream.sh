#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BATCH_ROOT="${BATCH_ROOT:-/mnt/pikachu/JAX_scRNAseq02_processed/ocm_prod_batch_20260520T024749Z}"
FIRST_ROOT="${FIRST_ROOT:-/mnt/pikachu/JAX_scRNAseq02_processed/ocm_prod_25E32-L3_aggressive_lowmem_20260519T201025Z}"
REMOTE_HOST="${REMOTE_HOST:-10.159.4.53}"
REMOTE_ROOT="${REMOTE_ROOT:-/home/lhhung/jax_scrnaseq02_remote_downstream}"
OUTPUT_NAME="${OUTPUT_NAME:-downstream_genefull_velocyto_cellbender_remote}"
LOCK="${LOCK:-${BATCH_ROOT}/logs/post_phase.lock}"
SLEEP_SECONDS="${SLEEP_SECONDS:-300}"

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"
}

run_one() {
  local root="$1"
  local sample="$2"
  local sample_dir="${root}/samples/${sample}"
  local final_h5ad="${sample_dir}/${OUTPUT_NAME}/final_counts.h5ad"
  local local_log="${root}/logs/${sample}.remote_downstream.backfill.log"

  [[ -f "${final_h5ad}" ]] && return 0
  [[ -d "${sample_dir}/run/outs/raw_feature_bc_matrix" ]] || return 1
  [[ -d "${sample_dir}/run/outs/raw_velocyto_feature_bc_matrix" ]] || return 1

  mkdir -p "${root}/logs"
  log "BACKFILL_START sample=${sample} root=${root}"
  flock "${LOCK}" \
    "${REPO_ROOT}/scripts/run_remote_scrna_downstream_rsync.sh" \
      --sample-dir "${sample_dir}" \
      --remote-host "${REMOTE_HOST}" \
      --remote-root "${REMOTE_ROOT}" \
      --output-name "${OUTPUT_NAME}" \
      --run-cellbender \
      --adaptive-filter \
      --n-mad 3 \
      --cellbender-gpu \
      --cellbender-cpu-cores 16 \
      --no-sync-images \
      --local-log "${local_log}"

  [[ -f "${final_h5ad}" ]] || {
    log "BACKFILL_FAIL sample=${sample} missing final=${final_h5ad}"
    return 2
  }
  log "BACKFILL_PASS sample=${sample} final=${final_h5ad}"
}

check_lib() {
  local root="$1"
  shift
  local sample
  local pending=0
  for sample in "$@"; do
    if [[ ! -f "${root}/samples/${sample}/${OUTPUT_NAME}/final_counts.h5ad" ]]; then
      pending=1
      run_one "${root}" "${sample}" || true
    fi
  done
  return "${pending}"
}

main() {
  log "BACKFILL_WATCH_START batch_root=${BATCH_ROOT}"
  while :; do
    local pending=0
    check_lib "${FIRST_ROOT}" \
      GCM1-Day-4 GRHL1-Day-4 OVOL1-Day-4 WT-PrS-20pct-Day-4 || pending=1
    check_lib "${BATCH_ROOT}/25E32-L4" \
      EPAS1-Day-4 WT-PrS-3pct-Day-4 || pending=1
    check_lib "${BATCH_ROOT}/25E34-L3" \
      GCM1-Day-5 GRHL1-Day-5 OVOL1-Day-5 WT-PrS-20pct-Day-5 || pending=1
    check_lib "${BATCH_ROOT}/25E34-L4" \
      ISL1-Day-5 EPAS1-Day-5 WT-PrS-3pct-Day-5 WT-ExM-Day-5 || pending=1
    check_lib "${BATCH_ROOT}/25E35-L3" \
      GCM1-Day-6 GRHL1-Day-6 OVOL1-Day-6 WT-PrS-20pct-Day-6 || pending=1
    check_lib "${BATCH_ROOT}/25E35-L4" \
      ISL1-Day-6 EPAS1-Day-6 WT-PrS-3pct-Day-6 WT-ExM-Day-6 || pending=1

    if [[ "${pending}" == "0" ]]; then
      log "BACKFILL_COMPLETE all expected final h5ad files are present"
      return 0
    fi
    sleep "${SLEEP_SECONDS}"
  done
}

main "$@"
