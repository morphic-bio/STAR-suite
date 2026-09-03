#!/usr/bin/env bash
# DIAGNOSTIC: re-profile the CBQ hash-only arm now that the per-read sample-token
# mutex is gone (mapping 230s -> 150s). Finds the next bottleneck in both the
# mapping phase and the Solo tail. Wall time here is perturbed by sampling.
set -uo pipefail
ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
MATRIX="${ROOT_DIR}/docs/benchmarks/jax_matrix_20260903/run_jax_matrix.sh"
OUT=/home/lhhung/jax_matrix_20260903/diag/reprofile
SENT="${OUT}/reprofile.done"
DATA="${OUT}/perf.data"
say() { echo "[$(date -u +%FT%TZ)] $*"; }
mkdir -p "${OUT}"; rm -f "${SENT}"

if pgrep -f "[s]ource/STAR --runThreadN" > /dev/null; then
    echo "refusing: a STAR run is in progress" > "${SENT}"; exit 1
fi

say "profiling the fixed binary on cbq_noalign"
sudo -n perf record -F 99 -a -g -o "${DATA}" > "${OUT}/perf_record.log" 2>&1 &
PERF_PID=$!
sleep 2
ALLOW_DIRTY=1 REPS=1 ARTIFACT_ROOT="${OUT}/run" MATRIX_LOG_DIR="${OUT}/logs" \
    bash "${MATRIX}" cbq_noalign > "${OUT}/driver.log" 2>&1
echo "run_exit=$?" >> "${SENT}.partial"
real_perf=$(pgrep -f "^/usr/lib/linux-tools/.*/perf record" | head -1)
[ -n "${real_perf}" ] && sudo -n kill -INT "${real_perf}"
for i in $(seq 1 60); do pgrep -f "^/usr/lib/linux-tools/.*/perf record" >/dev/null || break; sleep 1; done
wait "${PERF_PID}" 2>/dev/null
sudo -n chown "$(id -u):$(id -g)" "${DATA}" 2>/dev/null

# Phase windows: the mapping phase is roughly the first half, the Solo tail the rest.
for w in "mapping:2%-48%" "tail:50%-99%"; do
    n=${w%%:*}; r=${w##*:}
    sudo -n perf report -i "${DATA}" --force --no-children --stdio -g none --comm STAR \
        --sort sym --time "$r" --percent-limit 1 2>/dev/null \
        | grep -v "^#" | grep -v "^$" \
        | sed 's/std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >/string/g; s/std::allocator<[^>]*>//g; s/(.*$//' \
        | head -22 > "${OUT}/self_${n}.txt"
done
sudo -n perf report -i "${DATA}" --force --no-children --stdio --comm STAR \
    -S native_queued_spin_lock_slowpath -G --percent-limit 2 2>/dev/null \
    | grep -v "^#" | grep -v "^$" | head -40 > "${OUT}/locks.txt"

stamps() { grep -hE "started mapping|finished mapping|finished Solo counting" "$1" 2>/dev/null | sed 's/^/    /'; }
{
    echo "finished=$(date -u +%FT%TZ)"; cat "${SENT}.partial" 2>/dev/null
    echo "--- phase stamps (perturbed by sampling) ---"
    stamps "$(ls -t "${OUT}"/logs/cbq_noalign_rep1_*.stdout.txt 2>/dev/null | head -1)"
    echo "--- MAPPING: top self time ---"; cat "${OUT}/self_mapping.txt" 2>/dev/null
    echo "--- TAIL: top self time ---"; cat "${OUT}/self_tail.txt" 2>/dev/null
    echo "--- remaining lock contention ---"; head -24 "${OUT}/locks.txt" 2>/dev/null
} > "${SENT}"
say "re-profile complete"; cat "${SENT}"
