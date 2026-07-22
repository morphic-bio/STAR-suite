#!/usr/bin/env bash
# UCSF 2M Stage 2 Velocyto: same gates as hash_100k + UCSF_2M_PFCONFIG + RSS cap (like exact_2m).
# Phase 6: Gene/GeneFull parity + Log.final.out timing in SUMMARY.

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
OUT="${UCSF_VELOCYTO_HASH_2M_OUTDIR:-${SCRIPT_DIR}/ucsf_velocyto_hash_2m_output_${STAMP}}"

die() {
  echo "ERROR: $*" >&2
  exit 1
}

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

[[ -n "${UCSF_2M_PFCONFIG:-}" ]] || die "Set UCSF_2M_PFCONFIG to the full-sample pfMultiConfig CSV"
[[ -f "${CANONICAL}" ]] || die "Missing ${CANONICAL}"

if [[ "${UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY:-0}" == "1" ]]; then
  echo "WARNING: UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY=1 — frozen baseline skipped (dev-only)" >&2
else
  [[ -n "${UCSF_VELOCYTO_BASELINE_OUTDIR:-}" ]] || die \
    "Set UCSF_VELOCYTO_BASELINE_OUTDIR or UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY=1"
fi

if [[ "${UCSF_VELOCYTO_ALLOW_UNCAPPED_2M:-0}" == "1" ]]; then
  echo "WARNING: UCSF_VELOCYTO_ALLOW_UNCAPPED_2M=1 — RSS budget not enforced (dev-only)" >&2
else
  [[ -n "${UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB:-}" ]] || die \
    "Set UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB or UCSF_VELOCYTO_ALLOW_UNCAPPED_2M=1"
fi

mkdir -p "${OUT}"
PREP=(--prepare-mex)

echo "=== UCSF 2M Velocyto Stage 2 (integrated hash / CB-bucket) ==="
echo "outdir=${OUT}"

unset STAR_VELOCYTO_INTEGRATED_HASH || true
unset STAR_VELOCYTO_DETERMINISTIC_REPLAY || true
"${CANONICAL}" --profile 2m --threads 1 --out-prefix "${OUT}/stream_t1" "${PREP[@]}"

if [[ -n "${UCSF_VELOCYTO_BASELINE_OUTDIR:-}" ]]; then
  [[ -d "${UCSF_VELOCYTO_BASELINE_OUTDIR}/Solo.out" ]] || die "UCSF_VELOCYTO_BASELINE_OUTDIR must contain Solo.out/: ${UCSF_VELOCYTO_BASELINE_OUTDIR}"
  "${COMPARE[@]}" --mode all "${UCSF_VELOCYTO_BASELINE_OUTDIR}" "${OUT}/stream_t1"
  echo "PASS: stream_t1 matches frozen baseline"
fi

export STAR_VELOCYTO_DETERMINISTIC_REPLAY=1
unset STAR_VELOCYTO_INTEGRATED_HASH || true
"${CANONICAL}" --profile 2m --threads 1 --out-prefix "${OUT}/det_sort_t1" "${PREP[@]}"

"${COMPARE[@]}" --mode all "${OUT}/stream_t1" "${OUT}/det_sort_t1"
echo "PASS: stream vs Stage-1 deterministic (global sort)"

export STAR_VELOCYTO_INTEGRATED_HASH=1
"${CANONICAL}" --profile 2m --threads 1 --out-prefix "${OUT}/det_hash_t1" "${PREP[@]}"

"${COMPARE[@]}" --mode all "${OUT}/det_sort_t1" "${OUT}/det_hash_t1"
echo "PASS: Stage 1 sorted replay vs Stage 2 CB-bucket (1 thread)"

"${CANONICAL}" --profile 2m --threads "${N_THREADS}" --out-prefix "${OUT}/det_hash_tN" "${PREP[@]}"

"${COMPARE[@]}" --mode all "${OUT}/det_hash_t1" "${OUT}/det_hash_tN"
echo "PASS: Stage 2 deterministic 1-thread vs ${N_THREADS}-thread"

"${COMPARE[@]}" --mode genes "${OUT}/stream_t1" "${OUT}/det_sort_t1" || die "Phase 6 Gene/GeneFull: stream_t1 vs det_sort_t1"
echo "PASS: Phase 6 Gene/GeneFull stream vs Stage-1 det_sort"
"${COMPARE[@]}" --mode genes "${OUT}/stream_t1" "${OUT}/det_hash_t1" || die "Phase 6 Gene/GeneFull: stream_t1 vs det_hash_t1"
echo "PASS: Phase 6 Gene/GeneFull stream vs Stage-2 det_hash"
"${COMPARE[@]}" --mode genes "${OUT}/det_hash_t1" "${OUT}/det_hash_tN" || die "Phase 6 Gene/GeneFull: det_hash_t1 vs det_hash_tN"
echo "PASS: Phase 6 Gene/GeneFull Stage-2 1t vs Nt"

RSS_OUT="$("${RSS_REPORT[@]}" "${OUT}/det_sort_t1/Log.out" "${OUT}/det_hash_t1/Log.out" "${OUT}/det_hash_tN/Log.out")" || die "RSS report failed"
echo "${RSS_OUT}"
echo "VmRSS: 2M cap applies to PER_LOG_MAX_VM_RSS_KB[det_hash_t1] and [det_hash_tN] only (not global MAX_VM_RSS_KB vs det_sort_t1)."

if [[ "${UCSF_VELOCYTO_ALLOW_UNCAPPED_2M:-0}" != "1" ]]; then
  enforce_stage2_rss_cap_vs_per_log "${RSS_OUT}"
fi

unset STAR_VELOCYTO_INTEGRATED_HASH || true

{
  echo "OK: ${OUT}"
  log_final_head "${OUT}/stream_t1" "stream_t1"
  log_final_head "${OUT}/det_sort_t1" "det_sort_t1"
  log_final_head "${OUT}/det_hash_t1" "det_hash_t1"
  log_final_head "${OUT}/det_hash_tN" "det_hash_tN"
  echo "${RSS_OUT}"
} | tee "${OUT}/SUMMARY.txt"
