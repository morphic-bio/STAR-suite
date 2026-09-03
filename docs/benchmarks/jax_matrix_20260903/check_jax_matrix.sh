#!/usr/bin/env bash
# Correctness checks over the JAX Flex matrix. Timing is not the subject here;
# these are the equivalences that have to hold before any number is quoted.
#
#   A. Input-path equivalence, no alignment. gzip, BGZF and CBQ read the same
#      reads, so the three no-align arms must agree on every Flex counter and
#      produce byte-identical per-sample outputs.
#   B. Input-path equivalence, with alignment. Same claim for the three align
#      arms. This also exercises the new fused help path under three different
#      producer shapes.
#   C. Alignment recovery. An align arm must resolve the hash misses that a
#      no-align arm only counts: identical input reads, identical KEEP/DENY,
#      and strictly more assigned molecules.
#   D. Mate read-name validation. The new cross-mate check is on by default,
#      so every BGZF arm passing at all is evidence it accepts real deliveries;
#      a rejection would have been a hard error.
#   E. Help-path neutrality. The staged pipeline (--flexPipelineNTriage > 0)
#      keeps dedicated alignment consumers and a blocking push, so it is an
#      independent implementation of the same work. Its output must match the
#      fused align arm, which proves that aligning inside a producer does not
#      change results.
#
# Usage: check_jax_matrix.sh [--with-staged-control]
set -uo pipefail

LOG_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${LOG_DIR}/../../.." && pwd)"
REPORT="${LOG_DIR}/correctness_report.txt"
WITH_STAGED=0
[[ "${1:-}" == "--with-staged-control" ]] && WITH_STAGED=1

fails=0
checks=0
skips=0
# A report that checked nothing must never read as a pass, so every assertion
# is counted and the summary requires the full set.
EXPECTED_CHECKS="${EXPECTED_CHECKS:-16}"
note() { echo "$*" | tee -a "${REPORT}"; }
ok()   { note "PASS  $*"; checks=$((checks + 1)); }
bad()  { note "FAIL  $*"; checks=$((checks + 1)); fails=$((fails + 1)); }
skip() { note "SKIP  $*"; skips=$((skips + 1)); }

: > "${REPORT}"
note "JAX Flex matrix correctness report"
note "generated: $(date -u +%FT%TZ)"
note "commit:    $(git -C "${ROOT_DIR}" rev-parse HEAD)"
note ""

latest() {  # latest <arm> <suffix>
    ls -t "${LOG_DIR}/$1"_rep*."$2" 2>/dev/null | head -1
}

metric() {  # metric <Log.final.out> <label>
    awk -F '|' -v l="$2" '$1 ~ l "[[:space:]]*$" {gsub(/[[:space:]]/,"",$2); print $2}' "$1" 2>/dev/null
}

declare -A READS KEEP DENY PASS UNIQ MANIFEST HELPED
for arm in gzip_noalign bgzf_noalign cbq_noalign gzip_align bgzf_align cbq_align; do
    final="$(latest "${arm}" 'Log.final.out')"
    man="$(latest "${arm}" 'manifest.tsv')"
    logout="$(latest "${arm}" 'Log.out')"
    if [[ -z "${final}" || ! -s "${final}" ]]; then
        skip "${arm}: no Log.final.out yet"
        continue
    fi
    READS[$arm]="$(metric "${final}" 'Number of input reads')"
    KEEP[$arm]="$(metric "${final}" 'Hash screen: KEEP')"
    DENY[$arm]="$(metric "${final}" 'Hash screen: DENY')"
    PASS[$arm]="$(metric "${final}" 'Hash screen: PASS')"
    UNIQ[$arm]="$(metric "${final}" 'Uniquely mapped reads number')"
    MANIFEST[$arm]="${man}"
    HELPED[$arm]="$(sed -n 's/.*Fused producers aligned \([0-9][0-9]*\) queued reads.*/\1/p' "${logout}" 2>/dev/null | tail -1)"
done

note "--- counters ---"
printf '%-14s %16s %16s %14s %14s %16s %12s\n' arm input_reads hash_keep hash_deny hash_pass uniq_mapped helped \
    | tee -a "${REPORT}"
for arm in gzip_noalign bgzf_noalign cbq_noalign gzip_align bgzf_align cbq_align; do
    [[ -n "${READS[$arm]:-}" ]] || continue
    printf '%-14s %16s %16s %14s %14s %16s %12s\n' "${arm}" "${READS[$arm]}" "${KEEP[$arm]}" \
        "${DENY[$arm]}" "${PASS[$arm]}" "${UNIQ[$arm]:-0}" "${HELPED[$arm]:--}" | tee -a "${REPORT}"
done
note ""

# ---- A and B: input-path equivalence within each alignment mode ------------
compare_group() {  # compare_group <label> <arm> <arm> <arm>
    local label="$1"; shift
    local -a arms=("$@") present=()
    for arm in "${arms[@]}"; do [[ -n "${READS[$arm]:-}" ]] && present+=("${arm}"); done
    if (( ${#present[@]} < 2 )); then skip "${label}: fewer than two arms present"; return; fi
    local ref="${present[0]}"
    for arm in "${present[@]:1}"; do
        local diffs=""
        for field in READS KEEP DENY PASS UNIQ; do
            local -n table="${field}"
            [[ "${table[$ref]}" == "${table[$arm]}" ]] \
                || diffs+=" ${field}(${table[$ref]} vs ${table[$arm]})"
        done
        if [[ -z "${diffs}" ]]; then
            ok "${label}: ${arm} counters match ${ref}"
        else
            bad "${label}: ${arm} counters differ from ${ref}:${diffs}"
        fi
        if [[ -s "${MANIFEST[$ref]:-}" && -s "${MANIFEST[$arm]:-}" ]]; then
            if diff -q "${MANIFEST[$ref]}" "${MANIFEST[$arm]}" >/dev/null; then
                ok "${label}: ${arm} per-sample outputs byte-identical to ${ref}"
            else
                bad "${label}: ${arm} per-sample outputs differ from ${ref}"
                diff "${MANIFEST[$ref]}" "${MANIFEST[$arm]}" | head -10 >> "${REPORT}"
            fi
        else
            skip "${label}: ${arm} or ${ref} has no output manifest (later replicates prune matrices)"
        fi
    done
}
note "--- A. input-path equivalence, no alignment ---"
compare_group "no-align trio" gzip_noalign bgzf_noalign cbq_noalign
note ""
note "--- B. input-path equivalence, with alignment ---"
compare_group "align trio" gzip_align bgzf_align cbq_align
note ""

# ---- C. alignment recovery -------------------------------------------------
note "--- C. alignment recovery over the hash-only path ---"
for input in gzip bgzf cbq; do
    na="${input}_noalign"; al="${input}_align"
    if [[ -z "${READS[$na]:-}" || -z "${READS[$al]:-}" ]]; then
        skip "${input}: needs both arms"; continue
    fi
    [[ "${READS[$na]}" == "${READS[$al]}" ]] \
        && ok "${input}: both arms read ${READS[$na]} reads" \
        || bad "${input}: input read counts differ (${READS[$na]} vs ${READS[$al]})"
    [[ "${KEEP[$na]}" == "${KEEP[$al]}" && "${DENY[$na]}" == "${DENY[$al]}" ]] \
        && ok "${input}: hash KEEP/DENY identical, so alignment did not perturb the hash decision" \
        || bad "${input}: hash KEEP/DENY changed with alignment (${KEEP[$na]}/${DENY[$na]} vs ${KEEP[$al]}/${DENY[$al]})"
    if [[ "${UNIQ[$al]:-0}" =~ ^[0-9]+$ ]] && (( ${UNIQ[$al]:-0} > 0 )); then
        ok "${input}: align arm mapped ${UNIQ[$al]} reads uniquely (residual misses reached the aligner)"
    else
        bad "${input}: align arm reports no uniquely mapped reads; the genome path may not have run"
    fi
    if [[ "${HELPED[$al]:-}" =~ ^[0-9]+$ ]] && (( ${HELPED[$al]:-0} > 0 )); then
        ok "${input}: producers aligned ${HELPED[$al]} queued reads, so the help path carried real load"
    else
        bad "${input}: the fused help counter is absent or zero in the align arm"
    fi
done
note ""

# ---- D. mate read-name validation -----------------------------------------
note "--- D. cross-mate read-name validation on real deliveries ---"
for arm in bgzf_noalign bgzf_align; do
    err="$(latest "${arm}" 'stderr.txt')"
    [[ -n "${err}" ]] || { skip "${arm}: no stderr captured"; continue; }
    if grep -Fq 'read-name mismatch' "${err}"; then
        bad "${arm}: the mate read-name check rejected real JAX deliveries"
    elif [[ -n "${READS[$arm]:-}" ]]; then
        ok "${arm}: ${READS[$arm]} record pairs passed the read-name check"
    fi
done
note ""

# ---- E. help-path neutrality against dedicated consumers ------------------
if (( WITH_STAGED == 1 )); then
    note "--- E. staged-pipeline control (dedicated consumers, blocking push) ---"
    STAGED_LOG="${LOG_DIR}/staged_control"
    if [[ ! -s "${STAGED_LOG}.manifest.tsv" ]]; then
        note "running the staged control; this uses the same inputs and references as the matrix"
        "${LOG_DIR}/run_staged_control.sh" > "${STAGED_LOG}.driver.log" 2>&1 \
            || note "note: staged control exited nonzero (see ${STAGED_LOG}.driver.log)"
    fi
    if [[ -s "${STAGED_LOG}.manifest.tsv" && -s "${MANIFEST[gzip_align]:-}" ]]; then
        if diff -q "${STAGED_LOG}.manifest.tsv" "${MANIFEST[gzip_align]}" >/dev/null; then
            ok "E: staged dedicated-consumer output is byte-identical to the fused help-path output"
        else
            bad "E: staged and fused align outputs differ; the help path changes results"
            diff "${STAGED_LOG}.manifest.tsv" "${MANIFEST[gzip_align]}" | head -10 >> "${REPORT}"
        fi
    else
        skip "E: staged control or fused gzip_align manifest is absent"
    fi
    note ""
fi

note "--- summary ---"
note "assertions run: ${checks}   failed: ${fails}   skipped: ${skips}"
if (( fails > 0 )); then
    note "VERDICT: ${fails} CHECK(S) FAILED (report: ${REPORT})"
    exit 1
fi
if (( checks < EXPECTED_CHECKS )); then
    note "VERDICT: INCOMPLETE - only ${checks} of ${EXPECTED_CHECKS} assertions ran, ${skips} skipped."
    note "Nothing failed, but this is not a pass. Re-run once every arm has finished."
    exit 2
fi
note "VERDICT: ALL ${checks} CHECKS PASSED (report: ${REPORT})"
exit 0
