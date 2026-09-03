#!/usr/bin/env bash
# Confirm the new FlexFilter default (tags parallel + full-budget Monte Carlo)
# with NO environment overrides set, so the shipped code path is what is timed.
# Expect flexfilter ~19 s and byte-identical outputs.
set -uo pipefail
ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
MATRIX="${ROOT_DIR}/docs/benchmarks/jax_matrix_20260903/run_jax_matrix.sh"
BENCH="${ROOT_DIR}/docs/benchmarks/jax_matrix_20260903"
OUT=/home/lhhung/jax_matrix_20260903/diag/confirm
SENT="${OUT}/confirm.done"
say() { echo "[$(date -u +%FT%TZ)] $*"; }
mkdir -p "${OUT}"; rm -f "${SENT}"

for i in $(seq 1 240); do pgrep -f "[s]ource/STAR --runThreadN" > /dev/null || break; sleep 15; done
pgrep -f "[s]ource/STAR --runThreadN" > /dev/null \
    && { echo "refusing: a STAR run is still in progress" > "${SENT}"; exit 1; }

say "rebuilding STAR with the hybrid default"
make -C "${ROOT_DIR}/core/legacy/source" -j24 STAR > "${OUT}/build.log" 2>&1 \
    || { echo "build failed; see ${OUT}/build.log" > "${SENT}"; exit 1; }

unset STAR_FLEXFILTER_TAG_THREADS STAR_FLEXFILTER_MC_THREADS
say "running cbq_noalign with no overrides (cold cache)"
ALLOW_DIRTY=1 REPS=1 ARTIFACT_ROOT="${OUT}" MATRIX_LOG_DIR="${OUT}/logs" \
    bash "${MATRIX}" cbq_noalign > "${OUT}/driver.log" 2>&1
rc=$?

logout="$(ls -dt "${OUT}"/cbq_noalign_rep1_*/ 2>/dev/null | head -1)Log.out"
stdout="$(ls -t "${OUT}"/logs/cbq_noalign_rep1_*.stdout.txt 2>/dev/null | head -1)"
{
    echo "finished=$(date -u +%FT%TZ)"; echo "run_exit=${rc}"
    grep -hE "FlexFilter\] Processing|FlexFilter\] EmptyDrops MC" "${stdout%.stdout.txt}.stderr.txt" 2>/dev/null | sed 's/^/    /'
    awk '/started mapping/ {m0=$3} /finished mapping/ {m1=$3} /finished Solo counting/ {s1=$3}
         END { split(m0,a,":"); split(m1,b,":"); split(s1,c,":");
               t0=a[1]*3600+a[2]*60+a[3]; t1=b[1]*3600+b[2]*60+b[3]; t2=c[1]*3600+c[2]*60+c[3];
               printf "    mapping=%ds tail=%ds\n", t1-t0, t2-t1 }' "${stdout}" 2>/dev/null
    awk '/Running flexfilter inline/ {split($3,a,":"); t0=a[1]*3600+a[2]*60+a[3]}
         /Flexfilter pipeline complete/ {split($3,b,":"); t1=b[1]*3600+b[2]*60+b[3]}
         END {if (t0 && t1) printf "    flexfilter=%ds\n", t1-t0}' "${logout}" 2>/dev/null
    grep -h "Elapsed (wall clock)\|Percent of CPU\|Maximum resident" "${stdout%.stdout.txt}.time.txt" 2>/dev/null | sed 's/^/    /'
    a="$(ls -t "${BENCH}"/cbq_noalign_rep1_*.manifest.tsv | head -1)"
    b="$(ls -t "${OUT}"/logs/cbq_noalign_rep1_*.manifest.tsv 2>/dev/null | head -1)"
    if [[ -s "$a" && -s "$b" ]]; then
        diff -q "$a" "$b" >/dev/null && echo "  per-sample outputs: BYTE-IDENTICAL to the recorded benchmark arm" \
                                     || { echo "  per-sample outputs: DIFFER"; diff "$a" "$b" | head -8; }
    fi
} > "${SENT}"
say "confirm complete"; cat "${SENT}"
