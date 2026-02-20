#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

EXTERNAL_ENV="${REPO_ROOT}/tests/external_fixtures_env.sh"
if [[ -f "${EXTERNAL_ENV}" ]]; then
  # shellcheck disable=SC1091
  source "${EXTERNAL_ENV}"
fi

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
OUT_BASE="${PF_DYNAMIC_100K_OUT_BASE:-/tmp/pf_dynamic_permit_100k_$(date +%Y%m%d_%H%M%S)_$$}"

TIER="${PF_DYNAMIC_100K_TIER:-100000}"
THREADS="${PF_DYNAMIC_100K_THREADS:-8}"
MAP_PERMITS="${PF_DYNAMIC_100K_MAP_PERMITS:-4}"
READ_MAP_NUMBER="${PF_DYNAMIC_100K_READ_MAP_NUMBER:-}"
RUN_BASELINE="${PF_DYNAMIC_100K_RUN_BASELINE:-1}"
RUN_VARIABLE_CASES="${PF_DYNAMIC_100K_RUN_VARIABLE_CASES:-1}"
VAR_RETUNE_EVERY_ACQUIRES="${PF_DYNAMIC_100K_VAR_RETUNE_EVERY_ACQUIRES:-1}"
VAR_READ_MAP_NUMBER="${PF_DYNAMIC_100K_VAR_READ_MAP_NUMBER:-${READ_MAP_NUMBER}}"
MAX_FEATURE_WAIT_NS="${PF_DYNAMIC_100K_MAX_FEATURE_WAIT_NS:-20000000000}"
RUN_FORCED_EXIT_PROBE="${PF_DYNAMIC_100K_RUN_FORCED_EXIT_PROBE:-1}"
FORCED_EXIT_TIMEOUT_SEC="${PF_DYNAMIC_100K_FORCED_EXIT_TIMEOUT_SEC:-2}"
FORCED_EXIT_RECOVERY_READ_MAP_NUMBER="${PF_DYNAMIC_100K_FORCED_EXIT_RECOVERY_READ_MAP_NUMBER:-5000}"
RUN_PF_CONTROLLER_CASES="${PF_DYNAMIC_100K_RUN_PF_CONTROLLER_CASES:-1}"
PF_CONTROLLER_INTERVAL_MS="${PF_DYNAMIC_100K_PF_CONTROLLER_INTERVAL_MS:-5}"
PF_CONTROLLER_SEQUENCE="${PF_DYNAMIC_100K_PF_CONTROLLER_SEQUENCE:-2 4}"
PF_CONTROLLER_MAX_UPDATES="${PF_DYNAMIC_100K_PF_CONTROLLER_MAX_UPDATES:-64}"
PF_CONTROLLER_READ_MAP_NUMBER="${PF_DYNAMIC_100K_PF_CONTROLLER_READ_MAP_NUMBER:-${VAR_READ_MAP_NUMBER}}"

CASE_3_2_4_THREADS="${PF_DYNAMIC_100K_CASE_3_2_4_THREADS:-4}"
CASE_3_2_4_INIT_PERMITS="${PF_DYNAMIC_100K_CASE_3_2_4_INIT_PERMITS:-3}"
CASE_3_2_4_SEQUENCE="${PF_DYNAMIC_100K_CASE_3_2_4_SEQUENCE:-2 4}"
CASE_1_2_1_THREADS="${PF_DYNAMIC_100K_CASE_1_2_1_THREADS:-2}"
CASE_1_2_1_INIT_PERMITS="${PF_DYNAMIC_100K_CASE_1_2_1_INIT_PERMITS:-1}"
CASE_1_2_1_SEQUENCE="${PF_DYNAMIC_100K_CASE_1_2_1_SEQUENCE:-2 1}"

FASTQ_ROOT="${CR_MULTI_FASTQ_ROOT:-/storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs}"
GEX_BASE="${CR_MULTI_GEX_DIR:-${FASTQ_ROOT}/gex}"
CRISPR_BASE="${CR_MULTI_CRISPR_DIR:-${FASTQ_ROOT}/crispr}"
GEX_TIER_DIR="${PF_DYNAMIC_100K_GEX_DIR:-${GEX_BASE}/downsampled_${TIER}_v2}"
CRISPR_TIER_DIR="${PF_DYNAMIC_100K_CRISPR_DIR:-${CRISPR_BASE}/downsampled_${TIER}_v2}"

FEATURE_REF="${CR_MULTI_FEATURE_REF:-/storage/A375/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv}"
WHITELIST="${CR_MULTI_WHITELIST:-/storage/A375/3M-5pgex-jan-2023.txt}"
GENOME_DIR="${CR_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
LAST_API_RUN=""

log() {
  printf '%s\n' "$*"
}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

require_path() {
  local p="$1"
  local msg="$2"
  [[ -e "${p}" ]] || die "${msg}: ${p}"
}

get_kv_value() {
  local file="$1"
  local key="$2"
  awk -F= -v k="${key}" '$1==k{print $2}' "${file}" | tail -n 1
}

num_gt_zero() {
  local v="$1"
  awk -v x="${v:-0}" 'BEGIN{exit !(x+0 > 0)}'
}

num_eq_zero_or_empty() {
  local v="$1"
  if [[ -z "${v}" ]]; then
    return 0
  fi
  awk -v x="${v}" 'BEGIN{exit !(x+0 == 0)}'
}

num_le() {
  local value="$1"
  local max_allowed="$2"
  awk -v x="${value:-0}" -v y="${max_allowed:-0}" 'BEGIN{exit !(x+0 <= y+0)}'
}

hash_text_or_gz() {
  local file="$1"
  if [[ "${file}" == *.gz ]]; then
    zcat "${file}" | sha256sum | awk '{print $1}'
  else
    sha256sum "${file}" | awk '{print $1}'
  fi
}

assert_parity_file() {
  local off_file="$1"
  local on_file="$2"
  local label="$3"
  [[ -f "${off_file}" ]] || die "missing baseline file for ${label}: ${off_file}"
  [[ -f "${on_file}" ]] || die "missing dynamic file for ${label}: ${on_file}"
  local off_hash on_hash
  off_hash="$(hash_text_or_gz "${off_file}")"
  on_hash="$(hash_text_or_gz "${on_file}")"
  if [[ "${off_hash}" != "${on_hash}" ]]; then
    die "parity mismatch for ${label}: ${off_file} vs ${on_file}"
  fi
}

build_multi_config() {
  local out_dir="$1"
  local cfg="${out_dir}/multi_config.csv"
  cat > "${cfg}" <<CFG
[libraries]
fastqs,sample,library_type,feature_types
${GEX_TIER_DIR},A375,Gene Expression,Gene Expression
${CRISPR_TIER_DIR},A375,CRISPR Guide Capture,CRISPR Guide Capture

[feature]
ref,${FEATURE_REF}
CFG
  printf '%s\n' "${cfg}"
}

value_in_sequence_set() {
  local value="$1"
  local sequence="$2"
  local sequence_spaced="${sequence//,/ }"
  local seq_value
  for seq_value in ${sequence_spaced}; do
    if [[ "${value}" == "${seq_value}" ]]; then
      return 0
    fi
  done
  return 1
}

run_mode() {
  local mode_name="$1"
  local dynamic_interface="$2"
  local run_threads="$3"
  local map_permits="$4"
  local variable_threads="$5"
  local retune_every="$6"
  local permit_sequence="$7"
  local read_map_number="$8"
  local pf_controller_mode="${9:-off}"
  local pf_controller_interval_ms="${10:-0}"
  local pf_controller_sequence="${11:-}"
  local pf_controller_max_updates="${12:-0}"
  local mode_dir="${OUT_BASE}/${mode_name}"
  mkdir -p "${mode_dir}"

  local multi_cfg
  multi_cfg="$(build_multi_config "${mode_dir}")"

  local gex_r1 gex_r2
  gex_r1=$(ls "${GEX_TIER_DIR}"/*R1*.fastq.gz | paste -sd, -)
  gex_r2=$(ls "${GEX_TIER_DIR}"/*R2*.fastq.gz | paste -sd, -)

  local extra_args=()
  if [[ -n "${read_map_number}" ]]; then
    extra_args+=(--readMapNumber "${read_map_number}")
  fi
  if [[ "${variable_threads}" == "1" && "${retune_every}" -gt 0 && -n "${permit_sequence}" ]]; then
    local permit_sequence_spaced="${permit_sequence//,/ }"
    # shellcheck disable=SC2206
    local permit_sequence_array=(${permit_sequence_spaced})
    extra_args+=(
      --variableThreadsRetuneEveryAcquires "${retune_every}"
      --variableThreadsPermitSequence
      "${permit_sequence_array[@]}"
    )
  fi
  if [[ "${pf_controller_mode}" != "off" ]]; then
    extra_args+=(
      --dynamicThreadPfControllerMode "${pf_controller_mode}"
      --dynamicThreadPfControllerIntervalMs "${pf_controller_interval_ms}"
    )
    if [[ -n "${pf_controller_sequence}" ]]; then
      local pf_controller_sequence_spaced="${pf_controller_sequence//,/ }"
      # shellcheck disable=SC2206
      local pf_controller_sequence_array=(${pf_controller_sequence_spaced})
      extra_args+=(
        --dynamicThreadPfControllerSequence
        "${pf_controller_sequence_array[@]}"
      )
    fi
    if [[ "${pf_controller_max_updates}" -gt 0 ]]; then
      extra_args+=(--dynamicThreadPfControllerMaxUpdates "${pf_controller_max_updates}")
    fi
  fi

  log "=== Mode: ${mode_name} (dynamicThreadInterface=${dynamic_interface}) ==="
  log "out=${mode_dir}"
  if [[ "${pf_controller_mode}" != "off" ]]; then
    log "pfController mode=${pf_controller_mode} intervalMs=${pf_controller_interval_ms} sequence=${pf_controller_sequence} maxUpdates=${pf_controller_max_updates}"
  fi

  set +e
  /usr/bin/time -p "${STAR_BIN}" \
    --runThreadN "${run_threads}" \
    --genomeDir "${GENOME_DIR}" \
    --readFilesIn "${gex_r2}" "${gex_r1}" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${mode_dir}/" \
    --outSAMtype None \
    --clipAdapterType CellRanger4 \
    --alignEndsType Local \
    --chimSegmentMin 1000000 \
    --soloType CB_UMI_Simple \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 \
    --soloCBwhitelist "${WHITELIST}" \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR \
    --soloMultiMappers Unique \
    --soloCellFilter EmptyDrops_CR \
    --soloCbUbRequireTogether no \
    --soloStrand Unstranded \
    --soloFeatures GeneFull \
    --soloCrGexFeature genefull \
    --pfMultiConfig "${multi_cfg}" \
    --crFeatureRef "${FEATURE_REF}" \
    --crWhitelist "${WHITELIST}" \
    --crMinUmi 3 \
    --dynamicThreadInterface "${dynamic_interface}" \
    --dynamicThreadConstMapPermits "${map_permits}" \
    --dynamicThreadTelemetry 1 \
    --variableThreads "${variable_threads}" \
    "${extra_args[@]}" \
    > "${mode_dir}/STAR.stdout.log" 2> "${mode_dir}/STAR.stderr.log"
  local rc=$?
  set -e

  [[ ${rc} -eq 0 ]] || die "STAR failed for ${mode_name} (exit=${rc})"

  local api_run
  api_run=$(find "${mode_dir}/cr_assign" -type f -name assignBarcodes.api_run.txt | head -n 1 || true)
  [[ -n "${api_run}" ]] || die "missing assignBarcodes.api_run.txt for ${mode_name}"
  LAST_API_RUN="${api_run}"

  log "api_run=${api_run}"

  local hooks_enabled
  hooks_enabled="$(get_kv_value "${api_run}" "enableStarDynamicPermitHooks")"

  if [[ "${dynamic_interface}" == "1" ]]; then
    [[ "${hooks_enabled}" == "1" ]] || die "expected enableStarDynamicPermitHooks=1 in dynamic mode"

    local agg_acq agg_units feat_acq feat_units feat_wait
    local map_acq map_units map_wait
    agg_acq="$(get_kv_value "${api_run}" "dynamicPermitDelta.acquires")"
    agg_units="$(get_kv_value "${api_run}" "dynamicPermitDelta.workUnits")"
    feat_acq="$(get_kv_value "${api_run}" "dynamicPermitDelta.feature.acquires")"
    feat_units="$(get_kv_value "${api_run}" "dynamicPermitDelta.feature.workUnits")"
    feat_wait="$(get_kv_value "${api_run}" "dynamicPermitDelta.feature.waitNs")"
    map_acq="$(get_kv_value "${api_run}" "dynamicPermitDelta.map.acquires")"
    map_units="$(get_kv_value "${api_run}" "dynamicPermitDelta.map.workUnits")"
    map_wait="$(get_kv_value "${api_run}" "dynamicPermitDelta.map.waitNs")"

    num_gt_zero "${agg_acq}" || die "dynamicPermitDelta.acquires not > 0"
    num_gt_zero "${agg_units}" || die "dynamicPermitDelta.workUnits not > 0"
    [[ "${agg_acq}" == "${agg_units}" ]] || die "aggregate acquires/workUnits mismatch: ${agg_acq} vs ${agg_units}"
    num_gt_zero "${feat_acq}" || die "dynamicPermitDelta.feature.acquires not > 0"
    num_gt_zero "${feat_units}" || die "dynamicPermitDelta.feature.workUnits not > 0"
    [[ "${feat_acq}" == "${feat_units}" ]] || die "feature acquires/workUnits mismatch: ${feat_acq} vs ${feat_units}"
    num_gt_zero "${feat_wait}" || die "dynamicPermitDelta.feature.waitNs not > 0"
    num_le "${feat_wait}" "${MAX_FEATURE_WAIT_NS}" || die "dynamicPermitDelta.feature.waitNs exceeds threshold: ${feat_wait} > ${MAX_FEATURE_WAIT_NS}"
    num_eq_zero_or_empty "${map_acq}" || die "expected dynamicPermitDelta.map.acquires=0 in pf stage"
    num_eq_zero_or_empty "${map_units}" || die "expected dynamicPermitDelta.map.workUnits=0 in pf stage"
    num_eq_zero_or_empty "${map_wait}" || die "expected dynamicPermitDelta.map.waitNs=0 in pf stage"
  else
    [[ "${hooks_enabled}" == "0" ]] || die "expected enableStarDynamicPermitHooks=0 in baseline mode"

    local feat_acq map_acq
    feat_acq="$(get_kv_value "${api_run}" "dynamicPermitDelta.feature.acquires")"
    map_acq="$(get_kv_value "${api_run}" "dynamicPermitDelta.map.acquires")"
    num_eq_zero_or_empty "${feat_acq}" || die "baseline should not have feature permit acquires"
    num_eq_zero_or_empty "${map_acq}" || die "baseline should not have map permit acquires"
  fi
}

check_variable_case() {
  local case_name="$1"
  local api_run="$2"
  local permit_sequence="$3"
  local retunes target configured
  retunes="$(get_kv_value "${api_run}" "dynamicPermitDelta.retunes")"
  target="$(get_kv_value "${api_run}" "dynamicPermitAfter.targetPermits")"
  configured="$(get_kv_value "${api_run}" "dynamicPermitAfter.configuredPermits")"

  num_gt_zero "${retunes}" || die "${case_name}: expected dynamicPermitDelta.retunes > 0"
  [[ -n "${target}" ]] || die "${case_name}: missing dynamicPermitAfter.targetPermits"
  [[ -n "${configured}" ]] || die "${case_name}: missing dynamicPermitAfter.configuredPermits"
  value_in_sequence_set "${target}" "${permit_sequence}" || die "${case_name}: targetPermits ${target} not in permit sequence [${permit_sequence}]"
  value_in_sequence_set "${configured}" "${permit_sequence}" || die "${case_name}: configuredPermits ${configured} not in permit sequence [${permit_sequence}]"
}

run_variable_cases() {
  log "=== Variable sequence mode: 3->2->4 ==="
  run_mode \
    "variable_3_2_4" \
    1 \
    "${CASE_3_2_4_THREADS}" \
    "${CASE_3_2_4_INIT_PERMITS}" \
    1 \
    "${VAR_RETUNE_EVERY_ACQUIRES}" \
    "${CASE_3_2_4_SEQUENCE}" \
    "${VAR_READ_MAP_NUMBER}"
  check_variable_case "variable_3_2_4" "${LAST_API_RUN}" "${CASE_3_2_4_SEQUENCE}"

  log "=== Variable sequence mode: 1->2->1 ==="
  run_mode \
    "variable_1_2_1" \
    1 \
    "${CASE_1_2_1_THREADS}" \
    "${CASE_1_2_1_INIT_PERMITS}" \
    1 \
    "${VAR_RETUNE_EVERY_ACQUIRES}" \
    "${CASE_1_2_1_SEQUENCE}" \
    "${VAR_READ_MAP_NUMBER}"
  check_variable_case "variable_1_2_1" "${LAST_API_RUN}" "${CASE_1_2_1_SEQUENCE}"
}

check_pf_controller_case() {
  local case_name="$1"
  local api_run="$2"
  local expected_mode="$3"
  local expected_sequence="$4"
  local retunes target configured
  local controller_mode controller_ticks controller_applied controller_last_target controller_sequence expected_sequence_csv
  retunes="$(get_kv_value "${api_run}" "dynamicPermitDelta.retunes")"
  target="$(get_kv_value "${api_run}" "dynamicPermitAfter.targetPermits")"
  configured="$(get_kv_value "${api_run}" "dynamicPermitAfter.configuredPermits")"
  controller_mode="$(get_kv_value "${api_run}" "pfController.mode")"
  controller_ticks="$(get_kv_value "${api_run}" "pfController.ticks")"
  controller_applied="$(get_kv_value "${api_run}" "pfController.applied")"
  controller_last_target="$(get_kv_value "${api_run}" "pfController.lastTarget")"
  controller_sequence="$(get_kv_value "${api_run}" "pfController.sequence")"
  expected_sequence_csv="${expected_sequence// /,}"
  expected_sequence_csv="${expected_sequence_csv//,,/,}"

  [[ -n "${target}" ]] || die "${case_name}: missing dynamicPermitAfter.targetPermits"
  [[ -n "${configured}" ]] || die "${case_name}: missing dynamicPermitAfter.configuredPermits"
  [[ -n "${controller_mode}" ]] || die "${case_name}: missing pfController.mode"
  [[ "${controller_mode}" == "${expected_mode}" ]] || die "${case_name}: expected pfController.mode=${expected_mode}, got ${controller_mode}"
  num_gt_zero "${controller_ticks}" || die "${case_name}: expected pfController.ticks > 0"
  [[ -n "${controller_sequence}" ]] || die "${case_name}: missing pfController.sequence"
  [[ "${controller_sequence}" == "${expected_sequence_csv}" ]] || die "${case_name}: expected pfController.sequence=${expected_sequence_csv}, got ${controller_sequence}"
  [[ -n "${controller_last_target}" ]] || die "${case_name}: missing pfController.lastTarget"
  value_in_sequence_set "${controller_last_target}" "${expected_sequence}" || die "${case_name}: lastTarget ${controller_last_target} not in controller sequence [${expected_sequence}]"

  if [[ "${expected_mode}" == "active" ]]; then
    num_gt_zero "${controller_applied}" || die "${case_name}: expected pfController.applied > 0 in active mode"
    value_in_sequence_set "${target}" "${expected_sequence}" || die "${case_name}: targetPermits ${target} not in controller sequence [${expected_sequence}]"
    value_in_sequence_set "${configured}" "${expected_sequence}" || die "${case_name}: configuredPermits ${configured} not in controller sequence [${expected_sequence}]"
  else
    num_eq_zero_or_empty "${controller_applied}" || die "${case_name}: expected pfController.applied=0 in shadow mode"
    num_eq_zero_or_empty "${retunes}" || die "${case_name}: expected dynamicPermitDelta.retunes=0 in shadow mode"
  fi
}

run_pf_controller_cases() {
  log "=== PF controller mode: shadow ==="
  run_mode \
    "pf_controller_shadow" \
    1 \
    "${THREADS}" \
    "${MAP_PERMITS}" \
    1 \
    0 \
    "" \
    "${PF_CONTROLLER_READ_MAP_NUMBER}" \
    "shadow" \
    "${PF_CONTROLLER_INTERVAL_MS}" \
    "${PF_CONTROLLER_SEQUENCE}" \
    "${PF_CONTROLLER_MAX_UPDATES}"
  check_pf_controller_case "pf_controller_shadow" "${LAST_API_RUN}" "shadow" "${PF_CONTROLLER_SEQUENCE}"

  log "=== PF controller mode: active ==="
  run_mode \
    "pf_controller_active" \
    1 \
    "${THREADS}" \
    "${MAP_PERMITS}" \
    1 \
    0 \
    "" \
    "${PF_CONTROLLER_READ_MAP_NUMBER}" \
    "active" \
    "${PF_CONTROLLER_INTERVAL_MS}" \
    "${PF_CONTROLLER_SEQUENCE}" \
    "${PF_CONTROLLER_MAX_UPDATES}"
  check_pf_controller_case "pf_controller_active" "${LAST_API_RUN}" "active" "${PF_CONTROLLER_SEQUENCE}"
}

run_forced_exit_probe() {
  local probe_dir="${OUT_BASE}/forced_exit_probe"
  mkdir -p "${probe_dir}"
  local multi_cfg
  multi_cfg="$(build_multi_config "${probe_dir}")"
  local gex_r1 gex_r2
  gex_r1=$(ls "${GEX_TIER_DIR}"/*R1*.fastq.gz | paste -sd, -)
  gex_r2=$(ls "${GEX_TIER_DIR}"/*R2*.fastq.gz | paste -sd, -)

  log "=== Forced-exit probe (timeout ${FORCED_EXIT_TIMEOUT_SEC}s) ==="
  set +e
  timeout --signal=TERM "${FORCED_EXIT_TIMEOUT_SEC}s" \
    "${STAR_BIN}" \
    --runThreadN "${THREADS}" \
    --genomeDir "${GENOME_DIR}" \
    --readFilesIn "${gex_r2}" "${gex_r1}" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${probe_dir}/killed/" \
    --outSAMtype None \
    --clipAdapterType CellRanger4 \
    --alignEndsType Local \
    --chimSegmentMin 1000000 \
    --soloType CB_UMI_Simple \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 \
    --soloCBwhitelist "${WHITELIST}" \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR \
    --soloMultiMappers Unique \
    --soloCellFilter EmptyDrops_CR \
    --soloCbUbRequireTogether no \
    --soloStrand Unstranded \
    --soloFeatures GeneFull \
    --soloCrGexFeature genefull \
    --pfMultiConfig "${multi_cfg}" \
    --crFeatureRef "${FEATURE_REF}" \
    --crWhitelist "${WHITELIST}" \
    --crMinUmi 3 \
    --dynamicThreadInterface 1 \
    --dynamicThreadConstMapPermits "${MAP_PERMITS}" \
    --dynamicThreadTelemetry 1 \
    --variableThreads 1 \
    --variableThreadsRetuneEveryAcquires 1 \
    --variableThreadsPermitSequence 2 4 \
    > "${probe_dir}/forced_exit.stdout.log" 2> "${probe_dir}/forced_exit.stderr.log"
  local rc=$?
  set -e

  if [[ "${rc}" -eq 0 ]]; then
    die "forced-exit probe completed successfully; timeout did not trigger"
  fi
  if [[ "${rc}" -ne 124 && "${rc}" -ne 137 && "${rc}" -ne 143 ]]; then
    die "forced-exit probe returned unexpected code: ${rc}"
  fi

  log "forced-exit probe exit code=${rc} (expected timeout/termination)"

  log "=== Forced-exit recovery probe ==="
  run_mode \
    "forced_exit_recovery" \
    1 \
    "${THREADS}" \
    "${MAP_PERMITS}" \
    0 \
    0 \
    "" \
    "${FORCED_EXIT_RECOVERY_READ_MAP_NUMBER}"
}

main() {
  require_path "${STAR_BIN}" "Missing STAR binary"
  require_path "${GENOME_DIR}" "Missing genomeDir"
  require_path "${FEATURE_REF}" "Missing feature reference"
  require_path "${WHITELIST}" "Missing whitelist"
  require_path "${GEX_TIER_DIR}" "Missing GEX tier directory"
  require_path "${CRISPR_TIER_DIR}" "Missing CRISPR tier directory"

  ls "${GEX_TIER_DIR}"/*R1*.fastq.gz >/dev/null 2>&1 || die "No GEX R1 fastqs in ${GEX_TIER_DIR}"
  ls "${GEX_TIER_DIR}"/*R2*.fastq.gz >/dev/null 2>&1 || die "No GEX R2 fastqs in ${GEX_TIER_DIR}"

  mkdir -p "${OUT_BASE}"

  log "=== PF Dynamic Permit 100K Smoke ==="
  log "STAR_BIN=${STAR_BIN}"
  log "OUT_BASE=${OUT_BASE}"
  log "TIER=${TIER}"
  log "THREADS=${THREADS}"
  log "MAP_PERMITS=${MAP_PERMITS}"
  log "READ_MAP_NUMBER=${READ_MAP_NUMBER:-<all>}"
  log "RUN_VARIABLE_CASES=${RUN_VARIABLE_CASES}"
  log "VAR_RETUNE_EVERY_ACQUIRES=${VAR_RETUNE_EVERY_ACQUIRES}"
  log "VAR_READ_MAP_NUMBER=${VAR_READ_MAP_NUMBER:-<all>}"
  log "RUN_PF_CONTROLLER_CASES=${RUN_PF_CONTROLLER_CASES}"
  log "PF_CONTROLLER_INTERVAL_MS=${PF_CONTROLLER_INTERVAL_MS}"
  log "PF_CONTROLLER_SEQUENCE=${PF_CONTROLLER_SEQUENCE}"
  log "PF_CONTROLLER_MAX_UPDATES=${PF_CONTROLLER_MAX_UPDATES}"
  log "PF_CONTROLLER_READ_MAP_NUMBER=${PF_CONTROLLER_READ_MAP_NUMBER:-<all>}"
  log "MAX_FEATURE_WAIT_NS=${MAX_FEATURE_WAIT_NS}"
  log "RUN_FORCED_EXIT_PROBE=${RUN_FORCED_EXIT_PROBE}"
  log "FORCED_EXIT_TIMEOUT_SEC=${FORCED_EXIT_TIMEOUT_SEC}"
  log "FORCED_EXIT_RECOVERY_READ_MAP_NUMBER=${FORCED_EXIT_RECOVERY_READ_MAP_NUMBER}"
  log "GEX_TIER_DIR=${GEX_TIER_DIR}"
  log "CRISPR_TIER_DIR=${CRISPR_TIER_DIR}"

  if [[ "${RUN_BASELINE}" == "1" ]]; then
    run_mode "baseline_off" 0 "${THREADS}" "${MAP_PERMITS}" 0 0 "" "${READ_MAP_NUMBER}"
  fi
  run_mode "dynamic_on" 1 "${THREADS}" "${MAP_PERMITS}" 0 0 "" "${READ_MAP_NUMBER}"

  if [[ "${RUN_BASELINE}" == "1" ]]; then
    local off_root on_root
    off_root="${OUT_BASE}/baseline_off/outs"
    on_root="${OUT_BASE}/dynamic_on/outs"

    assert_parity_file "${off_root}/raw_feature_bc_matrix/barcodes.tsv.gz" "${on_root}/raw_feature_bc_matrix/barcodes.tsv.gz" "raw barcodes"
    assert_parity_file "${off_root}/raw_feature_bc_matrix/features.tsv.gz" "${on_root}/raw_feature_bc_matrix/features.tsv.gz" "raw features"
    assert_parity_file "${off_root}/raw_feature_bc_matrix/matrix.mtx.gz" "${on_root}/raw_feature_bc_matrix/matrix.mtx.gz" "raw matrix"
    assert_parity_file "${off_root}/filtered_feature_bc_matrix/barcodes.tsv.gz" "${on_root}/filtered_feature_bc_matrix/barcodes.tsv.gz" "filtered barcodes"
    assert_parity_file "${off_root}/filtered_feature_bc_matrix/features.tsv.gz" "${on_root}/filtered_feature_bc_matrix/features.tsv.gz" "filtered features"
    assert_parity_file "${off_root}/filtered_feature_bc_matrix/matrix.mtx.gz" "${on_root}/filtered_feature_bc_matrix/matrix.mtx.gz" "filtered matrix"
  fi

  if [[ "${RUN_VARIABLE_CASES}" == "1" ]]; then
    run_variable_cases
  fi

  if [[ "${RUN_PF_CONTROLLER_CASES}" == "1" ]]; then
    run_pf_controller_cases
  fi

  if [[ "${RUN_FORCED_EXIT_PROBE}" == "1" ]]; then
    run_forced_exit_probe
  fi

  log "PASS: pf dynamic permit 100K smoke succeeded"
  log "Outputs: ${OUT_BASE}"
}

main "$@"
