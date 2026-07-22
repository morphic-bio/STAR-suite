#!/usr/bin/env bash
# Small-fixture GEX-only Velocyto exactness harness.
# Intended as the first correctness/debugging loop before full-sample perturb runs.
#
# Flow:
#   1. build or reuse a 100K downsampled GEX-only fixture from corrected EBs2_2
#   2. run stream_t1
#   3. run det_t1
#   4. run det_tN
#   5. compare Velocyto raw+packaged outputs and Gene/GeneFull outputs
#
# Baseline compare:
#   - publication-grade: set UCSF_VELOCYTO_BASELINE_OUTDIR to a frozen pre-refactor GEX-only run
#   - dev-only: export UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY=1 to skip baseline

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck disable=SC1090
source "${REPO_ROOT}/tests/external_fixtures_env.sh"

SETUP_FIXTURE="${REPO_ROOT}/tests/setup_ucsf_ebs2_2_gex_100k_fixture.sh"
CANONICAL="${REPO_ROOT}/scripts/run_star_velocyto_gexonly_canonical.sh"
COMPARE=(python3 "${REPO_ROOT}/scripts/compare_velocyto_mex.py")
RSS_REPORT=(python3 "${REPO_ROOT}/scripts/report_velocyto_sorted_replay_rss.py")

SOURCE_GEX_DIR="${UCSF_VELOCYTO_GEX_SOURCE_DIR:-${UCSF_2M_GEX_DIR:-/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2/GEX}}"
FIXTURE_ROOT="${UCSF_VELOCYTO_GEX_100K_FIXTURE_ROOT:-/storage/ucsf-velocyto-validation/fixtures/ebs2_2_gexonly_100k}"
FIXTURE_GEX_DIR="${UCSF_VELOCYTO_GEX_100K_GEX_DIR:-${FIXTURE_ROOT}/GEX}"
TOTAL_READS="${UCSF_VELOCYTO_GEX_100K_TOTAL_READS:-100000}"
FIXTURE_MODE="${UCSF_VELOCYTO_GEX_100K_MODE:-head}"
N_THREADS="${UCSF_VELOCYTO_PARITY_THREADS:-8}"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="${UCSF_VELOCYTO_GEXONLY_EXACT_100K_OUTDIR:-${SCRIPT_DIR}/ucsf_velocyto_gexonly_exact_100k_output_${STAMP}}"

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

[[ -x "${SETUP_FIXTURE}" ]] || die "Missing ${SETUP_FIXTURE}"
[[ -x "${CANONICAL}" ]] || die "Missing ${CANONICAL}"

if [[ ! -d "${FIXTURE_GEX_DIR}" || "${UCSF_VELOCYTO_REBUILD_GEX_100K_FIXTURE:-0}" == "1" ]]; then
  bash "${SETUP_FIXTURE}" \
    --source-gex-dir "${SOURCE_GEX_DIR}" \
    --fixture-root "${FIXTURE_ROOT}" \
    --total-reads "${TOTAL_READS}" \
    --mode "${FIXTURE_MODE}"
fi

[[ -d "${FIXTURE_GEX_DIR}" ]] || die "Fixture GEX dir missing after setup: ${FIXTURE_GEX_DIR}"

export UCSF_VELOCYTO_GEX_GENOME_DIR="${UCSF_VELOCYTO_GEX_GENOME_DIR:-${UCSF_2M_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}}"
export UCSF_VELOCYTO_GEX_CB_WHITELIST="${UCSF_VELOCYTO_GEX_CB_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-february-2018_TRU.txt}"

if [[ "${UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY:-0}" == "1" ]]; then
  echo "WARNING: UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY=1 — frozen baseline skipped (debug only)" >&2
else
  [[ -n "${UCSF_VELOCYTO_BASELINE_OUTDIR:-}" ]] || die \
    "Set UCSF_VELOCYTO_BASELINE_OUTDIR to frozen GEX-only baseline or UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY=1 for debug-only"
fi

mkdir -p "${OUT}"
PREP=(--prepare-mex)

echo "=== UCSF GEX-only 100K Velocyto exact debug ==="
echo "source_gex=${SOURCE_GEX_DIR}"
echo "fixture_gex=${FIXTURE_GEX_DIR}"
echo "outdir=${OUT}"

unset STAR_VELOCYTO_DETERMINISTIC_REPLAY || true
"${CANONICAL}" --gex-dir "${FIXTURE_GEX_DIR}" --threads 1 --out-prefix "${OUT}/stream_t1" "${PREP[@]}"

if [[ -n "${UCSF_VELOCYTO_BASELINE_OUTDIR:-}" ]]; then
  [[ -d "${UCSF_VELOCYTO_BASELINE_OUTDIR}/Solo.out" ]] || die "UCSF_VELOCYTO_BASELINE_OUTDIR must contain Solo.out/: ${UCSF_VELOCYTO_BASELINE_OUTDIR}"
  "${COMPARE[@]}" --mode all "${UCSF_VELOCYTO_BASELINE_OUTDIR}" "${OUT}/stream_t1"
  echo "PASS: stream_t1 matches frozen baseline outdir"
fi

export STAR_VELOCYTO_DETERMINISTIC_REPLAY=1
"${CANONICAL}" --gex-dir "${FIXTURE_GEX_DIR}" --threads 1 --out-prefix "${OUT}/det_t1" "${PREP[@]}"

"${COMPARE[@]}" --mode all "${OUT}/stream_t1" "${OUT}/det_t1"
echo "PASS: stream vs deterministic (1 thread)"

"${CANONICAL}" --gex-dir "${FIXTURE_GEX_DIR}" --threads "${N_THREADS}" --out-prefix "${OUT}/det_tN" "${PREP[@]}"

"${COMPARE[@]}" --mode all "${OUT}/det_t1" "${OUT}/det_tN"
echo "PASS: deterministic 1-thread vs ${N_THREADS}-thread"

"${COMPARE[@]}" --mode genes "${OUT}/stream_t1" "${OUT}/det_t1" || die "Gene/GeneFull: stream_t1 vs det_t1"
echo "PASS: Gene/GeneFull stream vs deterministic (1t)"
"${COMPARE[@]}" --mode genes "${OUT}/det_t1" "${OUT}/det_tN" || die "Gene/GeneFull: det_t1 vs det_tN"
echo "PASS: Gene/GeneFull deterministic 1t vs Nt"

RSS_OUT="$("${RSS_REPORT[@]}" "${OUT}/det_t1/Log.out" "${OUT}/det_tN/Log.out")" || die "RSS report failed"
echo "${RSS_OUT}"
MAX_RSS="$(echo "${RSS_OUT}" | sed -n 's/^MAX_VM_RSS_KB=//p')"
if [[ -n "${UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB:-}" && -n "${MAX_RSS}" ]]; then
  if [[ "${MAX_RSS}" -gt "${UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB}" ]]; then
    die "VmRSS ${MAX_RSS} kB exceeds UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB=${UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB}"
  fi
fi

if [[ "${UCSF_VELOCYTO_CHECK_LEGACY_THREAD_PARITY:-0}" == "1" ]]; then
  unset STAR_VELOCYTO_DETERMINISTIC_REPLAY || true
  "${CANONICAL}" --gex-dir "${FIXTURE_GEX_DIR}" --threads "${N_THREADS}" --out-prefix "${OUT}/stream_tN" "${PREP[@]}"
  "${COMPARE[@]}" --mode all "${OUT}/stream_t1" "${OUT}/stream_tN"
  echo "PASS: stream 1-thread vs ${N_THREADS}-thread"
else
  echo "OPTIONAL SKIP: stream 1t vs Nt (UCSF_VELOCYTO_CHECK_LEGACY_THREAD_PARITY=1)"
fi

{
  echo "OK: ${OUT}"
  echo "FIXTURE_ROOT=${FIXTURE_ROOT}"
  [[ -f "${FIXTURE_ROOT}/FIXTURE_MANIFEST.txt" ]] && cat "${FIXTURE_ROOT}/FIXTURE_MANIFEST.txt"
  log_final_head "${OUT}/stream_t1" "stream_t1"
  log_final_head "${OUT}/det_t1" "det_t1"
  log_final_head "${OUT}/det_tN" "det_tN"
  echo "${RSS_OUT}"
} | tee "${OUT}/SUMMARY.txt"
