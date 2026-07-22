#!/usr/bin/env bash
# UCSF 2M: same gates as run_ucsf_velocyto_exact_100k.sh (UCSF_2M_PFCONFIG, RSS cap, baseline default).

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
OUT="${UCSF_VELOCYTO_EXACT_2M_OUTDIR:-${SCRIPT_DIR}/ucsf_velocyto_exact_2m_output_${STAMP}}"

die() {
  echo "ERROR: $*" >&2
  exit 1
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
  echo "WARNING: UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY=1 — frozen baseline skipped (not sufficient for publication Stage 1)" >&2
else
  [[ -n "${UCSF_VELOCYTO_BASELINE_OUTDIR:-}" ]] || die \
    "Set UCSF_VELOCYTO_BASELINE_OUTDIR or UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY=1 for dev-only"
fi

if [[ "${UCSF_VELOCYTO_ALLOW_UNCAPPED_2M:-0}" == "1" ]]; then
  echo "WARNING: UCSF_VELOCYTO_ALLOW_UNCAPPED_2M=1 — RSS budget not enforced (not acceptable for production 2M acceptance)" >&2
else
  [[ -n "${UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB:-}" ]] || die \
    "Set UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB (VmRSS budget kB) or UCSF_VELOCYTO_ALLOW_UNCAPPED_2M=1"
fi

mkdir -p "${OUT}"
PREP=(--prepare-mex)

echo "=== UCSF 2M Velocyto Stage 1 parity ==="
echo "outdir=${OUT}"

unset STAR_VELOCYTO_DETERMINISTIC_REPLAY || true
"${CANONICAL}" --profile 2m --threads 1 --out-prefix "${OUT}/stream_t1" "${PREP[@]}"

if [[ -n "${UCSF_VELOCYTO_BASELINE_OUTDIR:-}" ]]; then
  [[ -d "${UCSF_VELOCYTO_BASELINE_OUTDIR}/Solo.out" ]] || die "UCSF_VELOCYTO_BASELINE_OUTDIR must contain Solo.out/: ${UCSF_VELOCYTO_BASELINE_OUTDIR}"
  "${COMPARE[@]}" --mode all "${UCSF_VELOCYTO_BASELINE_OUTDIR}" "${OUT}/stream_t1"
  echo "PASS: stream_t1 matches frozen baseline outdir"
fi

export STAR_VELOCYTO_DETERMINISTIC_REPLAY=1
"${CANONICAL}" --profile 2m --threads 1 --out-prefix "${OUT}/det_t1" "${PREP[@]}"

"${COMPARE[@]}" --mode all "${OUT}/stream_t1" "${OUT}/det_t1"
echo "PASS: stream vs deterministic (1 thread)"

"${CANONICAL}" --profile 2m --threads "${N_THREADS}" --out-prefix "${OUT}/det_tN" "${PREP[@]}"

"${COMPARE[@]}" --mode all "${OUT}/det_t1" "${OUT}/det_tN"
echo "PASS: deterministic 1-thread vs ${N_THREADS}-thread"

"${COMPARE[@]}" --mode genes "${OUT}/stream_t1" "${OUT}/det_t1" || die "Phase 6 Gene/GeneFull: stream_t1 vs det_t1"
echo "PASS: Phase 6 Gene/GeneFull stream vs deterministic (1t)"
"${COMPARE[@]}" --mode genes "${OUT}/det_t1" "${OUT}/det_tN" || die "Phase 6 Gene/GeneFull: det_t1 vs det_tN"
echo "PASS: Phase 6 Gene/GeneFull deterministic 1t vs Nt"

RSS_OUT="$("${RSS_REPORT[@]}" "${OUT}/det_t1/Log.out" "${OUT}/det_tN/Log.out")" || die "RSS report failed"
echo "${RSS_OUT}"
echo "VmRSS: use PER_LOG_MAX_VM_RSS_KB lines to compare det_t1 vs det_tN when MAX_VM_RSS_KB is ambiguous."
MAX_RSS="$(echo "${RSS_OUT}" | sed -n 's/^MAX_VM_RSS_KB=//p')"

if [[ "${UCSF_VELOCYTO_ALLOW_UNCAPPED_2M:-0}" != "1" ]]; then
  if [[ -z "${MAX_RSS}" ]]; then
    die "Could not parse MAX_VM_RSS_KB from Log.out"
  fi
  if [[ "${MAX_RSS}" -gt "${UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB}" ]]; then
    die "VmRSS ${MAX_RSS} kB exceeds UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB=${UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB}"
  fi
fi

if [[ "${UCSF_VELOCYTO_CHECK_LEGACY_THREAD_PARITY:-0}" == "1" ]]; then
  unset STAR_VELOCYTO_DETERMINISTIC_REPLAY || true
  "${CANONICAL}" --profile 2m --threads "${N_THREADS}" --out-prefix "${OUT}/stream_tN" "${PREP[@]}"
  "${COMPARE[@]}" --mode all "${OUT}/stream_t1" "${OUT}/stream_tN"
  echo "PASS: stream 1t vs Nt"
else
  echo "OPTIONAL SKIP: stream 1t vs Nt (UCSF_VELOCYTO_CHECK_LEGACY_THREAD_PARITY=1)"
fi

{
  echo "OK: ${OUT}"
  log_final_head "${OUT}/stream_t1" "stream_t1"
  log_final_head "${OUT}/det_t1" "det_t1"
  log_final_head "${OUT}/det_tN" "det_tN"
  echo "${RSS_OUT}"
} | tee "${OUT}/SUMMARY.txt"
