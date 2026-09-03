#!/usr/bin/env bash
# DIAGNOSTIC: decompose the flexfilter speedup.
#
# Changing FlexFilter's thread allocation changed TWO things at once: the Monte
# Carlo went from 2 threads to 32, and tag processing went from 16 threads to 1.
# The 68 s -> 27 s gain cannot be attributed to either on its own. This runs the
# 2x2 so the split is measured:
#
#   tags=16 mc=2    the original allocation                (expect ~68 s)
#   tags=1  mc=32   the current default                    (measured 27 s)
#   tags=1  mc=2    isolates the cost of serializing tags
#   tags=16 mc=32   hybrid, oversubscribed; the likely optimum
#
# Only flexfilter wall time matters here; mapping is identical in all four and
# each run costs ~5.5 min. Outputs are checked for byte-identity every time.
set -uo pipefail
ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
MATRIX="${ROOT_DIR}/docs/benchmarks/jax_matrix_20260903/run_jax_matrix.sh"
BENCH="${ROOT_DIR}/docs/benchmarks/jax_matrix_20260903"
OUT=/home/lhhung/jax_matrix_20260903/diag/ff2x2
SENT="${OUT}/ff2x2.done"
say() { echo "[$(date -u +%FT%TZ)] $*"; }
mkdir -p "${OUT}"; rm -f "${SENT}" "${SENT}.rows"

for i in $(seq 1 240); do pgrep -f "[s]ource/STAR --runThreadN" > /dev/null || break; sleep 15; done
pgrep -f "[s]ource/STAR --runThreadN" > /dev/null \
    && { echo "refusing: a STAR run is still in progress" > "${SENT}"; exit 1; }

say "rebuilding STAR with the flexfilter thread-split overrides"
make -C "${ROOT_DIR}/core/legacy/source" -j24 STAR > "${OUT}/build.log" 2>&1 \
    || { echo "build failed; see ${OUT}/build.log" > "${SENT}"; exit 1; }
[[ "$(grep -a -c STAR_FLEXFILTER_MC_THREADS "${ROOT_DIR}/core/legacy/source/STAR")" -gt 0 ]] \
    || { echo "override hook absent from the binary" > "${SENT}"; exit 1; }
say "build OK"

REF_MANIFEST="$(ls -t "${BENCH}"/cbq_noalign_rep1_*.manifest.tsv | head -1)"

run_cell() {  # run_cell <tagThreads> <mcThreads>
    local tags="$1" mc="$2" label="tags${1}_mc${2}"
    say "${label} starting"
    STAR_FLEXFILTER_TAG_THREADS="${tags}" STAR_FLEXFILTER_MC_THREADS="${mc}" \
    ALLOW_DIRTY=1 REPS=1 ARTIFACT_ROOT="${OUT}/${label}" MATRIX_LOG_DIR="${OUT}/${label}/logs" \
        bash "${MATRIX}" cbq_noalign > "${OUT}/${label}.driver.log" 2>&1
    local rc=$? logout stdout ff wall cpu ident
    logout="$(ls -dt "${OUT}/${label}"/cbq_noalign_rep1_*/ 2>/dev/null | head -1)Log.out"
    stdout="$(ls -t "${OUT}/${label}"/logs/cbq_noalign_rep1_*.stdout.txt 2>/dev/null | head -1)"
    ff="$(awk '/Running flexfilter inline/ {split($3,a,":"); t0=a[1]*3600+a[2]*60+a[3]}
               /Flexfilter pipeline complete/ {split($3,b,":"); t1=b[1]*3600+b[2]*60+b[3]}
               END {if (t0 && t1) print t1-t0}' "${logout}" 2>/dev/null)"
    wall="$(awk -F': ' '/Elapsed \(wall clock\)/ {print $NF}' "${stdout%.stdout.txt}.time.txt" 2>/dev/null)"
    cpu="$(awk -F': ' '/Percent of CPU/ {print $NF}' "${stdout%.stdout.txt}.time.txt" 2>/dev/null)"
    local newman; newman="$(ls -t "${OUT}/${label}"/logs/cbq_noalign_rep1_*.manifest.tsv 2>/dev/null | head -1)"
    if [[ -s "${newman}" && -s "${REF_MANIFEST}" ]]; then
        diff -q "${newman}" "${REF_MANIFEST}" >/dev/null && ident=identical || ident=DIFFERS
    else
        ident=no_manifest
    fi
    printf '%-8s %-8s %-14s %-12s %-9s %s\n' "${tags}" "${mc}" "${ff:-NA}s" "${wall:-NA}" "${cpu:-NA}" "${ident}" \
        >> "${SENT}.rows"
    say "${label} done: flexfilter=${ff:-NA}s wall=${wall:-NA} outputs=${ident} (exit ${rc})"
}

run_cell 16 2
run_cell 1 2
run_cell 16 32
run_cell 1 32

{
    echo "finished=$(date -u +%FT%TZ)"
    echo "flexfilter wall time by thread split (mapping is unchanged across all four):"
    printf '%-8s %-8s %-14s %-12s %-9s %s\n' tagThreads mcThreads flexfilter wall cpu outputs
    cat "${SENT}.rows" 2>/dev/null
    echo
    echo "Reading it: comparing 16/2 with 1/2 gives the cost of serializing tag work."
    echo "Comparing 1/2 with 1/32 gives the Monte Carlo threading gain alone."
    echo "16/32 is the hybrid; if it beats 1/32 the default should change."
} > "${SENT}"
say "2x2 complete"; cat "${SENT}"
