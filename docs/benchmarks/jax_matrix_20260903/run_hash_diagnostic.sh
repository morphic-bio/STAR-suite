#!/usr/bin/env bash
# DIAGNOSTIC, not a benchmark: where does STAR-Flex spend its time relative to
# cyto on the same CBQ lanes? Two runs of the hash-only CBQ arm on the SSD stage:
#
#   1. the unmodified path, with system-wide `perf record` sampling in the
#      background, so the breakdown is measured rather than argued. The wall
#      time of this run is perturbed by sampling and is not a benchmark number.
#   2. STAR_FLEX_HASH_H0_ONLY=1, the hypothesis under test: exact-match tier
#      only, no H1/deny lookup on a miss. Cold cache, quiet box; compare its
#      phase stamps (started/finished mapping, finished Solo) with the recorded
#      cbq_noalign arm.
#
# Both runs go through the matrix script so every path, flag and locus check is
# identical to the recorded arm. Their logs and results rows are kept under
# MATRIX_LOG_DIR inside the diagnostic tree, never in the benchmark record.
set -uo pipefail
ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
MATRIX="${ROOT_DIR}/docs/benchmarks/jax_matrix_20260903/run_jax_matrix.sh"
BENCH_LOG_DIR="${ROOT_DIR}/docs/benchmarks/jax_matrix_20260903"
OUT=/home/lhhung/jax_matrix_20260903/diag
SENT="${OUT}/hash_diagnostic.done"
PERF_DATA="${OUT}/perf_cbq_noalign.data"
say() { echo "[$(date -u +%FT%TZ)] $*"; }
mkdir -p "${OUT}"
rm -f "$SENT" "$SENT.partial"

if pgrep -f "[s]ource/STAR --runThreadN" > /dev/null; then
    echo "refusing: a STAR run is in progress" > "$SENT"; exit 1
fi

# Rebuild only now that no STAR process holds the binary.
say "rebuilding STAR with the H0-only diagnostic hook"
make -C "${ROOT_DIR}/core/legacy/source" -j24 STAR > "${OUT}/build_diag.log" 2>&1 \
    || { echo "build failed; see ${OUT}/build_diag.log" > "$SENT"; exit 1; }
# grep the binary directly: `strings | grep -q` trips pipefail when grep closes
# the pipe early, which reported a present hook as absent.
[[ "$(grep -a -c STAR_FLEX_HASH_H0_ONLY "${ROOT_DIR}/core/legacy/source/STAR")" -gt 0 ]] \
    || { echo "hook absent from binary" > "$SENT"; exit 1; }
say "build OK"

stamps() { grep -hE "started mapping|finished mapping|finished Solo counting" "$1" 2>/dev/null | sed 's/^/    /'; }
phase_seconds() {  # phase_seconds <stdout.txt>: mapping and tail durations
    awk '/started mapping/ {m0=$3} /finished mapping/ {m1=$3} /finished Solo counting/ {s1=$3}
         END { split(m0,a,":"); split(m1,b,":"); split(s1,c,":");
               t0=a[1]*3600+a[2]*60+a[3]; t1=b[1]*3600+b[2]*60+b[3]; t2=c[1]*3600+c[2]*60+c[3];
               if (t1>0 && t0>0) printf "    mapping=%ds tail=%ds\n", t1-t0, t2-t1 }' "$1"
}

# ---- 1. profile the unmodified path ----------------------------------------
if sudo -n perf stat -e cycles true > /dev/null 2>&1; then
    say "run 1: unmodified path under system-wide perf sampling"
    sudo -n perf record -F 99 -a -g -o "${PERF_DATA}" > "${OUT}/perf_record.log" 2>&1 &
    PERF_PID=$!
    sleep 2
    ALLOW_DIRTY=1 REPS=1 ARTIFACT_ROOT="${OUT}/perf_run" MATRIX_LOG_DIR="${OUT}/perf_run/logs" \
        bash "${MATRIX}" cbq_noalign > "${OUT}/perf_run.driver.log" 2>&1
    echo "perf_run_exit=$?" >> "$SENT.partial"
    sudo -n kill -INT "${PERF_PID}" 2>/dev/null; wait "${PERF_PID}" 2>/dev/null
    sudo -n chown "$(id -u):$(id -g)" "${PERF_DATA}" 2>/dev/null
    perf report -i "${PERF_DATA}" --no-children --percent-limit 1 --stdio --comm STAR -g none 2>/dev/null \
        | grep -v "^#" | grep -v "^$" | head -50 > "${OUT}/perf_self.txt"
    perf report -i "${PERF_DATA}" --children --percent-limit 2 --stdio --comm STAR -g none 2>/dev/null \
        | grep -v "^#" | grep -v "^$" | head -50 > "${OUT}/perf_inclusive.txt"
    perf report -i "${PERF_DATA}" --no-children --sort dso --stdio --comm STAR -g none 2>/dev/null \
        | grep -v "^#" | grep -v "^$" | head -12 > "${OUT}/perf_dso.txt"
else
    say "perf is not usable under sudo -n; skipping the profile"
    echo "perf_run_exit=skipped" >> "$SENT.partial"
fi

# ---- 2. H0-only hypothesis --------------------------------------------------
say "run 2: H0-only (no H1/deny lookup), cold cache, unprofiled"
STAR_FLEX_HASH_H0_ONLY=1 ALLOW_DIRTY=1 REPS=1 ARTIFACT_ROOT="${OUT}/h0only" MATRIX_LOG_DIR="${OUT}/h0only/logs" \
    bash "${MATRIX}" cbq_noalign > "${OUT}/h0only.driver.log" 2>&1
echo "h0only_run_exit=$?" >> "$SENT.partial"

# ---- summary ----------------------------------------------------------------
{
    echo "finished=$(date -u +%FT%TZ)"
    cat "$SENT.partial" 2>/dev/null
    echo "--- recorded cbq_noalign benchmark (H0+H1, unprofiled) ---"
    f="$(ls -t "${BENCH_LOG_DIR}"/cbq_noalign_rep1_*.stdout.txt 2>/dev/null | head -1)"; stamps "$f"; phase_seconds "$f"
    grep -h "Elapsed (wall clock)" "${f%.stdout.txt}.time.txt" 2>/dev/null | sed 's/^/    /'
    echo "--- run 1: H0+H1 under perf (perturbed) ---"
    f="$(ls -t "${OUT}"/perf_run/logs/cbq_noalign_rep1_*.stdout.txt 2>/dev/null | head -1)"; stamps "$f"; phase_seconds "$f"
    grep -h "Elapsed (wall clock)" "${f%.stdout.txt}.time.txt" 2>/dev/null | sed 's/^/    /'
    echo "--- run 2: H0-only, unprofiled ---"
    f="$(ls -t "${OUT}"/h0only/logs/cbq_noalign_rep1_*.stdout.txt 2>/dev/null | head -1)"; stamps "$f"; phase_seconds "$f"
    grep -h "Elapsed (wall clock)" "${f%.stdout.txt}.time.txt" 2>/dev/null | sed 's/^/    /'
    grep -h "Hash screen: \(KEEP\|DENY\|PASS\) " "${OUT}"/h0only/cbq_noalign_rep1_*/Log.final.out 2>/dev/null | sed 's/^/    /'
    echo "--- perf: self time by DSO ---"; cat "${OUT}/perf_dso.txt" 2>/dev/null
    echo "--- perf: top self-time symbols ---"; head -40 "${OUT}/perf_self.txt" 2>/dev/null
} > "$SENT"
say "diagnostic complete"; cat "$SENT"
