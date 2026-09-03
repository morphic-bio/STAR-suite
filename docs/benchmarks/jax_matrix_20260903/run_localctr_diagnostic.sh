#!/usr/bin/env bash
# DIAGNOSTIC: hash-only CBQ arm with the per-thread local counters
# (FlexLocalCounters) replacing the per-read shared atomics. Waits for the
# hash diagnostic to finish so the binary is not rebuilt under a running STAR,
# rebuilds, then runs once with a cold cache on a quiet box. Compare its phase
# stamps with the recorded cbq_noalign arm (mapping 230 s, tail 139 s).
set -uo pipefail
ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
MATRIX="${ROOT_DIR}/docs/benchmarks/jax_matrix_20260903/run_jax_matrix.sh"
DIAG=/home/lhhung/jax_matrix_20260903/diag
OUT="${DIAG}/localctr"
SENT="${OUT}/localctr.done"
say() { echo "[$(date -u +%FT%TZ)] $*"; }
mkdir -p "${OUT}"; rm -f "${SENT}"

for i in $(seq 1 240); do [[ -f "${DIAG}/hash_diagnostic.done" ]] && break; sleep 15; done
[[ -f "${DIAG}/hash_diagnostic.done" ]] || { echo "hash diagnostic never finished" > "${SENT}"; exit 1; }
if pgrep -f "[s]ource/STAR --runThreadN" > /dev/null; then
    echo "refusing: a STAR run is in progress" > "${SENT}"; exit 1
fi

say "rebuilding STAR with per-thread local counters"
make -C "${ROOT_DIR}/core/legacy/source" -j24 STAR > "${OUT}/build.log" 2>&1 \
    || { echo "build failed; see ${OUT}/build.log" > "${SENT}"; exit 1; }
[[ "$(grep -a -c "Fused producers aligned" "${ROOT_DIR}/core/legacy/source/STAR")" -gt 0 ]] \
    || { echo "unexpected binary" > "${SENT}"; exit 1; }
say "build OK; running cbq_noalign with local counters (cold cache)"

ALLOW_DIRTY=1 REPS=1 ARTIFACT_ROOT="${OUT}" MATRIX_LOG_DIR="${OUT}/logs" \
    bash "${MATRIX}" cbq_noalign > "${OUT}/driver.log" 2>&1
rc=$?

stamps() { grep -hE "started mapping|finished mapping|finished Solo counting" "$1" 2>/dev/null | sed 's/^/    /'; }
{
    echo "finished=$(date -u +%FT%TZ)"
    echo "run_exit=${rc}"
    echo "--- recorded cbq_noalign (shared atomics) ---"
    f="$(ls -t "${ROOT_DIR}"/docs/benchmarks/jax_matrix_20260903/cbq_noalign_rep1_*.stdout.txt | head -1)"; stamps "$f"
    grep -h "Elapsed (wall clock)\|Percent of CPU" "${f%.stdout.txt}.time.txt" | sed 's/^/    /'
    echo "--- local counters ---"
    f="$(ls -t "${OUT}"/logs/cbq_noalign_rep1_*.stdout.txt 2>/dev/null | head -1)"; stamps "$f"
    grep -h "Elapsed (wall clock)\|Percent of CPU" "${f%.stdout.txt}.time.txt" 2>/dev/null | sed 's/^/    /'
    echo "--- counters still exact? ---"
    grep -h "Number of input reads\|Hash screen: KEEP \|Hash screen: DENY \|Hash screen: PASS " "${OUT}"/cbq_noalign_rep1_*/Log.final.out 2>/dev/null | sed 's/^/    /'
    grep -h "Flex pipeline complete" "${OUT}"/cbq_noalign_rep1_*/Log.out 2>/dev/null | sed 's/^/    /'
} > "${SENT}"
say "done"; cat "${SENT}"
