#!/usr/bin/env bash
# DIAGNOSTIC: break down the ~26 s serial sumThreads block.
#
# sumThreads walks the 32 per-thread Solo structures one at a time and, for each,
# merges the inline hash and pending-ambiguous sidecars, adds stats, adds counts,
# then destroys that thread's tables. Which of those four costs the 26 s has not
# been isolated, and the redesign (each thread merging itself as its final act,
# versus a parallel reduction over the dense count arrays) depends on the answer.
# This run reports the four sub-totals; no behaviour changes.
set -uo pipefail
ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
MATRIX="${ROOT_DIR}/docs/benchmarks/jax_matrix_20260903/run_jax_matrix.sh"
BENCH="${ROOT_DIR}/docs/benchmarks/jax_matrix_20260903"
OUT=/home/lhhung/jax_matrix_20260903/diag/sumthreads
SENT="${OUT}/sumthreads.done"
say() { echo "[$(date -u +%FT%TZ)] $*"; }
mkdir -p "${OUT}"; rm -f "${SENT}"

for i in $(seq 1 240); do pgrep -f "[s]ource/STAR --runThreadN" > /dev/null || break; sleep 15; done
pgrep -f "[s]ource/STAR --runThreadN" > /dev/null \
    && { echo "refusing: a STAR run is still in progress" > "${SENT}"; exit 1; }

say "rebuilding STAR with sumThreads sub-timers"
make -C "${ROOT_DIR}/core/legacy/source" -j24 STAR > "${OUT}/build.log" 2>&1 \
    || { echo "build failed; see ${OUT}/build.log" > "${SENT}"; exit 1; }
say "build OK; running cbq_noalign (cold cache)"

ALLOW_DIRTY=1 REPS=1 ARTIFACT_ROOT="${OUT}" MATRIX_LOG_DIR="${OUT}/logs" \
    bash "${MATRIX}" cbq_noalign > "${OUT}/driver.log" 2>&1
rc=$?

logout="$(ls -dt "${OUT}"/cbq_noalign_rep1_*/ 2>/dev/null | head -1)Log.out"
stdout="$(ls -t "${OUT}"/logs/cbq_noalign_rep1_*.stdout.txt 2>/dev/null | head -1)"
{
    echo "finished=$(date -u +%FT%TZ)"; echo "run_exit=${rc}"
    echo "--- sumThreads breakdown ---"
    grep -h "Solo timing: sumThreads" "${logout}" 2>/dev/null | sed 's/^/    /'
    echo "--- surrounding phases ---"
    grep -hE "Starting Solo post-map|Starting bucket-parallel|Solo timing:" "${logout}" 2>/dev/null | sed 's/^/    /'
    awk '/started mapping/ {m0=$3} /finished mapping/ {m1=$3} /finished Solo counting/ {s1=$3}
         END { split(m0,a,":"); split(m1,b,":"); split(s1,c,":");
               t0=a[1]*3600+a[2]*60+a[3]; t1=b[1]*3600+b[2]*60+b[3]; t2=c[1]*3600+c[2]*60+c[3];
               printf "    mapping=%ds tail=%ds\n", t1-t0, t2-t1 }' "${stdout}" 2>/dev/null
    grep -h "Elapsed (wall clock)\|Percent of CPU" "${stdout%.stdout.txt}.time.txt" 2>/dev/null | sed 's/^/    /'
    a="$(ls -t "${BENCH}"/cbq_noalign_rep1_*.manifest.tsv | head -1)"
    b="$(ls -t "${OUT}"/logs/cbq_noalign_rep1_*.manifest.tsv 2>/dev/null | head -1)"
    if [[ -s "$a" && -s "$b" ]]; then
        diff -q "$a" "$b" >/dev/null && echo "  per-sample outputs: BYTE-IDENTICAL" \
                                     || { echo "  per-sample outputs: DIFFER"; diff "$a" "$b" | head -6; }
    fi
} > "${SENT}"
say "done"; cat "${SENT}"
