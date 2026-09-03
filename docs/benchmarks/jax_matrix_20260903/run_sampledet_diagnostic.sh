#!/usr/bin/env bash
# DIAGNOSTIC: hash-only CBQ arm with (a) per-thread local counters and (b) the
# sample-detector token registration made lock-free. The DWARF profile of the
# recorded arm attributed ~30% of all samples to a kernel spinlock reached from
# SampleDetector::detectSampleFromPackedTag -> registerSampleToken, which took
# a process-wide mutex on every read to write an idempotent 32-entry table.
# Waits for the local-counter diagnostic so the binary is never rebuilt under a
# running STAR, rebuilds, runs once cold on a quiet box.
set -uo pipefail
ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
MATRIX="${ROOT_DIR}/docs/benchmarks/jax_matrix_20260903/run_jax_matrix.sh"
DIAG=/home/lhhung/jax_matrix_20260903/diag
OUT="${DIAG}/sampledet"
SENT="${OUT}/sampledet.done"
say() { echo "[$(date -u +%FT%TZ)] $*"; }
mkdir -p "${OUT}"; rm -f "${SENT}"

for i in $(seq 1 240); do [[ -f "${DIAG}/localctr/localctr.done" ]] && break; sleep 15; done
[[ -f "${DIAG}/localctr/localctr.done" ]] || { echo "local-counter diagnostic never finished" > "${SENT}"; exit 1; }
if pgrep -f "[s]ource/STAR --runThreadN" > /dev/null; then
    echo "refusing: a STAR run is in progress" > "${SENT}"; exit 1
fi

say "rebuilding STAR with the lock-free sample-token table"
make -C "${ROOT_DIR}/core/legacy/source" -j24 STAR > "${OUT}/build.log" 2>&1 \
    || { echo "build failed; see ${OUT}/build.log" > "${SENT}"; exit 1; }
say "build OK; running cbq_noalign (cold cache)"

ALLOW_DIRTY=1 REPS=1 ARTIFACT_ROOT="${OUT}" MATRIX_LOG_DIR="${OUT}/logs" \
    bash "${MATRIX}" cbq_noalign > "${OUT}/driver.log" 2>&1
rc=$?

stamps() { grep -hE "started mapping|finished mapping|finished Solo counting" "$1" 2>/dev/null | sed 's/^/    /'; }
phase() { awk '/started mapping/ {m0=$3} /finished mapping/ {m1=$3} /finished Solo counting/ {s1=$3}
    END { split(m0,a,":"); split(m1,b,":"); split(s1,c,":");
          t0=a[1]*3600+a[2]*60+a[3]; t1=b[1]*3600+b[2]*60+b[3]; t2=c[1]*3600+c[2]*60+c[3];
          if (t1>0 && t0>0) printf "    mapping=%ds tail=%ds\n", t1-t0, t2-t1 }' "$1"; }
report() {  # report <label> <stdout.txt>
    echo "--- $1 ---"; stamps "$2"; phase "$2"
    grep -h "Elapsed (wall clock)\|Percent of CPU\|Maximum resident" "${2%.stdout.txt}.time.txt" 2>/dev/null | sed 's/^/    /'
}
{
    echo "finished=$(date -u +%FT%TZ)"; echo "run_exit=${rc}"
    report "recorded cbq_noalign (shared atomics, per-read mutex)" \
        "$(ls -t "${ROOT_DIR}"/docs/benchmarks/jax_matrix_20260903/cbq_noalign_rep1_*.stdout.txt | head -1)"
    report "local counters only" "$(ls -t "${DIAG}"/localctr/logs/cbq_noalign_rep1_*.stdout.txt 2>/dev/null | head -1)"
    report "local counters + lock-free sample table" "$(ls -t "${OUT}"/logs/cbq_noalign_rep1_*.stdout.txt 2>/dev/null | head -1)"
    echo "--- functionality unchanged? (must match the recorded arm) ---"
    for d in "${ROOT_DIR}/docs/benchmarks/jax_matrix_20260903" "${OUT}/logs"; do
        f="$(ls -t "$d"/cbq_noalign_rep1_*.Log.final.out 2>/dev/null | head -1)"
        echo "  $(basename "$f")"; grep -h "Number of input reads\|Hash screen: KEEP \|Hash screen: DENY \|Hash screen: PASS " "$f" | sed 's/^/    /'
    done
    a="$(ls -t "${ROOT_DIR}"/docs/benchmarks/jax_matrix_20260903/cbq_noalign_rep1_*.manifest.tsv | head -1)"
    b="$(ls -t "${OUT}"/logs/cbq_noalign_rep1_*.manifest.tsv 2>/dev/null | head -1)"
    if [[ -s "$a" && -s "$b" ]]; then
        diff -q "$a" "$b" >/dev/null && echo "  per-sample outputs: BYTE-IDENTICAL to the recorded arm" \
                                     || { echo "  per-sample outputs: DIFFER"; diff "$a" "$b" | head -6; }
    fi
} > "${SENT}"
say "done"; cat "${SENT}"
