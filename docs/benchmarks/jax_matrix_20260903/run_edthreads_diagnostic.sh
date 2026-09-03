#!/usr/bin/env bash
# DIAGNOSTIC: hash-only CBQ arm with the EmptyDrops Monte Carlo given the whole
# thread budget (FlexFilter split it evenly across 16 tags, leaving the 4 tags
# that actually run simulations on 2 threads each while 24 cores idled).
# flexfilter took 68 s of the 137 s Flex tail. Each simulation seeds from its
# iteration index alone and the per-thread tallies are summed as integers, so
# the p-values are independent of the thread split: this must be a pure speedup
# with byte-identical outputs. That is exactly what this checks.
set -uo pipefail
ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
MATRIX="${ROOT_DIR}/docs/benchmarks/jax_matrix_20260903/run_jax_matrix.sh"
BENCH="${ROOT_DIR}/docs/benchmarks/jax_matrix_20260903"
DIAG=/home/lhhung/jax_matrix_20260903/diag
OUT="${DIAG}/edthreads"
SENT="${OUT}/edthreads.done"
say() { echo "[$(date -u +%FT%TZ)] $*"; }
mkdir -p "${OUT}"; rm -f "${SENT}"

for i in $(seq 1 240); do pgrep -f "[s]ource/STAR --runThreadN" > /dev/null || break; sleep 15; done
if pgrep -f "[s]ource/STAR --runThreadN" > /dev/null; then
    echo "refusing: a STAR run is still in progress" > "${SENT}"; exit 1
fi

say "rebuilding STAR with a threaded EmptyDrops Monte Carlo"
make -C "${ROOT_DIR}/core/legacy/source" -j24 STAR > "${OUT}/build.log" 2>&1 \
    || { echo "build failed; see ${OUT}/build.log" > "${SENT}"; exit 1; }
say "build OK; running cbq_noalign (cold cache)"

ALLOW_DIRTY=1 REPS=1 ARTIFACT_ROOT="${OUT}" MATRIX_LOG_DIR="${OUT}/logs" \
    bash "${MATRIX}" cbq_noalign > "${OUT}/driver.log" 2>&1
rc=$?

flexsec() {  # seconds between the flexfilter start and completion lines
    awk '/Running flexfilter inline/ {split($3,a,":"); t0=a[1]*3600+a[2]*60+a[3]}
         /Flexfilter pipeline complete/ {split($3,b,":"); t1=b[1]*3600+b[2]*60+b[3]}
         END {if (t0 && t1) printf "    flexfilter=%ds\n", t1-t0}' "$1"
}
phase() { awk '/started mapping/ {m0=$3} /finished mapping/ {m1=$3} /finished Solo counting/ {s1=$3}
    END { split(m0,a,":"); split(m1,b,":"); split(s1,c,":");
          t0=a[1]*3600+a[2]*60+a[3]; t1=b[1]*3600+b[2]*60+b[3]; t2=c[1]*3600+c[2]*60+c[3];
          if (t1>0 && t0>0) printf "    mapping=%ds tail=%ds\n", t1-t0, t2-t1 }' "$1"; }
report() {  # report <label> <stdout.txt> <Log.out>
    echo "--- $1 ---"; phase "$2"; flexsec "$3"
    grep -h "Elapsed (wall clock)\|Percent of CPU" "${2%.stdout.txt}.time.txt" 2>/dev/null | sed 's/^/    /'
}
{
    echo "finished=$(date -u +%FT%TZ)"; echo "run_exit=${rc}"
    prev_s="$(ls -t "${DIAG}"/sampledet/logs/cbq_noalign_rep1_*.stdout.txt 2>/dev/null | head -1)"
    prev_l="$(ls -dt "${DIAG}"/sampledet/cbq_noalign_rep1_*/ 2>/dev/null | head -1)Log.out"
    new_s="$(ls -t "${OUT}"/logs/cbq_noalign_rep1_*.stdout.txt 2>/dev/null | head -1)"
    new_l="$(ls -dt "${OUT}"/cbq_noalign_rep1_*/ 2>/dev/null | head -1)Log.out"
    report "lock-free sample table (serial EmptyDrops)" "${prev_s}" "${prev_l}"
    report "plus threaded EmptyDrops" "${new_s}" "${new_l}"
    echo "--- cell calls per tag: must be identical ---"
    paste <(grep -oE "\[BC0[0-9]{2}\] .*Final=[0-9]+" "${prev_l}" 2>/dev/null | sed 's/.*\(\[BC[0-9]*\]\).*Final=\([0-9]*\).*/\1 \2/') \
          <(grep -oE "Final=[0-9]+" "${new_l}" 2>/dev/null | cut -d= -f2) 2>/dev/null | sed 's/^/    /'
    a="$(ls -t "${BENCH}"/cbq_noalign_rep1_*.manifest.tsv | head -1)"
    b="$(ls -t "${OUT}"/logs/cbq_noalign_rep1_*.manifest.tsv 2>/dev/null | head -1)"
    if [[ -s "$a" && -s "$b" ]]; then
        diff -q "$a" "$b" >/dev/null && echo "  per-sample outputs: BYTE-IDENTICAL to the recorded benchmark arm" \
                                     || { echo "  per-sample outputs: DIFFER from the recorded arm"; diff "$a" "$b" | head -8; }
    fi
} > "${SENT}"
say "done"; cat "${SENT}"
