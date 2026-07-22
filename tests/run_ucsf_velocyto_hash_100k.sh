#!/usr/bin/env bash
# UCSF 100K Stage 2 Velocyto: CB-bucket / integrated-hash deterministic path vs Stage 1 global sort.
# Phase 6 (runbook): Gene/GeneFull parity + Log.final.out timing headers in SUMMARY.
#
# Requires STAR built with Stage 2 (STAR_VELOCYTO_INTEGRATED_HASH + deterministic replay).
#
# DEFAULT: same baseline gate as Stage 1 exact harness (frozen reference vs stream is Stage 1's job;
#   this script focuses on det_sort vs det_hash and hash thread parity).
#   UCSF_VELOCYTO_BASELINE_OUTDIR or UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY=1
#
# Env: UCSF_VELOCYTO_PARITY_THREADS (default 8), UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB (optional),
#   when set, enforced vs PER_LOG_MAX_VM_RSS_KB[det_hash_t1] and [det_hash_tN] only (not global MAX vs det_sort_t1).
#   UCSF_VELOCYTO_HASH_100K_OUTDIR, STAR_BIN

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck disable=SC1090
source "${REPO_ROOT}/tests/external_fixtures_env.sh"

CANONICAL="${REPO_ROOT}/scripts/run_star_velocyto_canonical.sh"
COMPARE=(python3 "${REPO_ROOT}/scripts/compare_velocyto_mex.py")
RSS_REPORT=(python3 "${REPO_ROOT}/scripts/report_velocyto_sorted_replay_rss.py")
N_THREADS="${UCSF_VELOCYTO_PARITY_THREADS:-8}"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="${UCSF_VELOCYTO_HASH_100K_OUTDIR:-${SCRIPT_DIR}/ucsf_velocyto_hash_100k_output_${STAMP}}"

die() {
  echo "ERROR: $*" >&2
  exit 1
}

# When UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB is set, compare only Stage 2 runs (not MAX_VM_RSS_KB across det_sort + det_hash).
enforce_stage2_rss_cap_vs_per_log() {
  local rss_out="$1"
  local cap="${UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB:-}"
  [[ -z "${cap}" ]] && return 0
  local tag v
  for tag in det_hash_t1 det_hash_tN; do
    v="$(echo "${rss_out}" | sed -n "s/^PER_LOG_MAX_VM_RSS_KB\\[${tag}\\]=//p" | head -n1)"
    [[ -n "${v}" && "${v}" != "NA" ]] || die \
      "Stage 2 RSS cap set but missing PER_LOG_MAX_VM_RSS_KB[${tag}] (need VmRSS markers in ${tag}/Log.out)"
    if [[ "${v}" -gt "${cap}" ]]; then
      die "VmRSS ${v} kB for ${tag} exceeds UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB=${cap} (Stage 2 only; det_sort_t1 not capped here)"
    fi
  done
}

log_final_head() {
  local dir="$1"
  local tag="$2"
  if [[ -f "${dir}/Log.final.out" ]]; then
    echo "--- Log.final.out (${tag}) ---"
    head -n 4 "${dir}/Log.final.out" | sed 's/^/  /'
  fi
}

[[ -f "${CANONICAL}" ]] || die "Missing ${CANONICAL}"

if [[ "${UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY:-0}" == "1" ]]; then
  echo "WARNING: UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY=1 — frozen baseline skipped (dev-only)" >&2
else
  [[ -n "${UCSF_VELOCYTO_BASELINE_OUTDIR:-}" ]] || die \
    "Set UCSF_VELOCYTO_BASELINE_OUTDIR or UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY=1"
fi

mkdir -p "${OUT}"
PREP=(--prepare-mex)

echo "=== UCSF 100K Velocyto Stage 2 (integrated hash / CB-bucket) ==="
echo "outdir=${OUT}"

# Optional: confirm baseline still matches stream (Stage 1); cheap 1-thread stream run.
unset STAR_VELOCYTO_INTEGRATED_HASH || true
unset STAR_VELOCYTO_DETERMINISTIC_REPLAY || true
"${CANONICAL}" --profile 100k --threads 1 --out-prefix "${OUT}/stream_t1" "${PREP[@]}"

if [[ -n "${UCSF_VELOCYTO_BASELINE_OUTDIR:-}" ]]; then
  [[ -d "${UCSF_VELOCYTO_BASELINE_OUTDIR}/Solo.out" ]] || die "UCSF_VELOCYTO_BASELINE_OUTDIR must contain Solo.out/: ${UCSF_VELOCYTO_BASELINE_OUTDIR}"
  "${COMPARE[@]}" --mode all "${UCSF_VELOCYTO_BASELINE_OUTDIR}" "${OUT}/stream_t1"
  echo "PASS: stream_t1 matches frozen baseline"
fi

export STAR_VELOCYTO_DETERMINISTIC_REPLAY=1
unset STAR_VELOCYTO_INTEGRATED_HASH || true
"${CANONICAL}" --profile 100k --threads 1 --out-prefix "${OUT}/det_sort_t1" "${PREP[@]}"

"${COMPARE[@]}" --mode all "${OUT}/stream_t1" "${OUT}/det_sort_t1"
echo "PASS: stream vs Stage-1 deterministic (global sort)"

export STAR_VELOCYTO_INTEGRATED_HASH=1
"${CANONICAL}" --profile 100k --threads 1 --out-prefix "${OUT}/det_hash_t1" "${PREP[@]}"

"${COMPARE[@]}" --mode all "${OUT}/det_sort_t1" "${OUT}/det_hash_t1"
echo "PASS: Stage 1 sorted replay vs Stage 2 CB-bucket (1 thread)"

"${CANONICAL}" --profile 100k --threads "${N_THREADS}" --out-prefix "${OUT}/det_hash_tN" "${PREP[@]}"

"${COMPARE[@]}" --mode all "${OUT}/det_hash_t1" "${OUT}/det_hash_tN"
echo "PASS: Stage 2 deterministic 1-thread vs ${N_THREADS}-thread"

# Phase 6: Velocyto variants must not change Solo Gene / GeneFull outputs.
"${COMPARE[@]}" --mode genes "${OUT}/stream_t1" "${OUT}/det_sort_t1" || die "Phase 6 Gene/GeneFull: stream_t1 vs det_sort_t1"
echo "PASS: Phase 6 Gene/GeneFull stream vs Stage-1 det_sort"
"${COMPARE[@]}" --mode genes "${OUT}/stream_t1" "${OUT}/det_hash_t1" || die "Phase 6 Gene/GeneFull: stream_t1 vs det_hash_t1"
echo "PASS: Phase 6 Gene/GeneFull stream vs Stage-2 det_hash"
"${COMPARE[@]}" --mode genes "${OUT}/det_hash_t1" "${OUT}/det_hash_tN" || die "Phase 6 Gene/GeneFull: det_hash_t1 vs det_hash_tN"
echo "PASS: Phase 6 Gene/GeneFull Stage-2 1t vs Nt"

RSS_OUT="$("${RSS_REPORT[@]}" "${OUT}/det_sort_t1/Log.out" "${OUT}/det_hash_t1/Log.out" "${OUT}/det_hash_tN/Log.out")" || die "RSS report failed"
echo "${RSS_OUT}"
echo "VmRSS: optional cap uses PER_LOG_MAX_VM_RSS_KB[det_hash_t1] and [det_hash_tN] only; MAX_VM_RSS_KB is informational across all passed logs."
enforce_stage2_rss_cap_vs_per_log "${RSS_OUT}"

unset STAR_VELOCYTO_INTEGRATED_HASH || true

{
  echo "OK: ${OUT}"
  log_final_head "${OUT}/stream_t1" "stream_t1"
  log_final_head "${OUT}/det_sort_t1" "det_sort_t1"
  log_final_head "${OUT}/det_hash_t1" "det_hash_t1"
  log_final_head "${OUT}/det_hash_tN" "det_hash_tN"
  echo "${RSS_OUT}"
} | tee "${OUT}/SUMMARY.txt"
