#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_cbq_e2e_module_regression_${STAMP}}"
RUN_NETWORK="${RUN_NETWORK:-0}"
RUN_CHROMAP_MAPPING_SMOKE="${RUN_CHROMAP_MAPPING_SMOKE:-auto}"

die() {
    echo "ERROR: $*" >&2
    exit 1
}

resolve_bqtools() {
    if [[ -n "${BQTOOLS:-}" ]]; then
        if [[ -x "$BQTOOLS" ]]; then
            printf '%s\n' "$BQTOOLS"
            return 0
        fi
        if command -v "$BQTOOLS" >/dev/null 2>&1; then
            command -v "$BQTOOLS"
            return 0
        fi
        return 1
    fi
    if [[ -x /tmp/star_suite_bqtools/bin/bqtools ]]; then
        printf '%s\n' /tmp/star_suite_bqtools/bin/bqtools
        return 0
    fi
    command -v bqtools >/dev/null 2>&1 || return 1
    command -v bqtools
}

BQTOOLS_BIN="$(resolve_bqtools)" || die "bqtools not found; set BQTOOLS=/path/to/bqtools"

rm -rf "$OUT_ROOT"
mkdir -p "$OUT_ROOT"
SUMMARY="$OUT_ROOT/summary.tsv"
printf 'status\tcase\tlog_dir\n' > "$SUMMARY"

run_case() {
    local case_name="$1"
    local script="$2"
    local case_dir="$OUT_ROOT/$case_name"
    local work_dir="$case_dir/work"

    [[ -f "$ROOT_DIR/$script" ]] || die "missing test script: $script"
    mkdir -p "$case_dir"

    printf 'RUN\t%s\t%s\n' "$case_name" "$case_dir" | tee -a "$SUMMARY"
    set +e
    (
        cd "$ROOT_DIR"
        env \
            BQTOOLS="$BQTOOLS_BIN" \
            OUT_ROOT="$work_dir" \
            RUN_CHROMAP_MAPPING_SMOKE="$RUN_CHROMAP_MAPPING_SMOKE" \
            bash "$ROOT_DIR/$script"
    ) > "$case_dir/stdout.log" 2> "$case_dir/stderr.log"
    local status=$?
    set -e

    printf '%s\n' "$status" > "$case_dir/exit_status"
    if [[ "$status" -ne 0 ]]; then
        printf 'FAIL\t%s\t%s\n' "$case_name" "$case_dir" | tee -a "$SUMMARY"
        return "$status"
    fi

    if grep -q '^SKIP:' "$case_dir/stdout.log"; then
        printf 'FAIL\t%s\t%s\n' "$case_name" "$case_dir" | tee -a "$SUMMARY"
        echo "ERROR: required CBQ regression case skipped: $case_name" >&2
        return 1
    fi

    printf 'PASS\t%s\t%s\n' "$case_name" "$case_dir" | tee -a "$SUMMARY"
}

run_case binseq_probe_contract tests/run_binseq_probe_smoke.sh
run_case cbq_cpp_reader_contract tests/run_cbq_cpp_reader_smoke.sh
run_case cbq_star_mapper_e2e tests/run_cbq_star_input_smoke.sh
run_case cbq_starsolo_e2e tests/run_cbq_solo_e2e_smoke.sh
run_case cbq_process_features_adapter tests/run_cbq_pf_adapter_smoke.sh
run_case cbq_chromap_adapter tests/run_cbq_chromap_adapter_smoke.sh

if [[ "$RUN_NETWORK" == "1" ]]; then
    run_case binseq_upstream_arc_fixture tests/run_binseq_upstream_fixture_smoke.sh
    run_case cbq_flex_tiny_public tests/run_cbq_flex_tiny_public_smoke.sh
else
    printf 'SKIP\tbinseq_upstream_arc_fixture\tRUN_NETWORK=0\n' | tee -a "$SUMMARY"
    printf 'SKIP\tcbq_flex_tiny_public\tRUN_NETWORK=0\n' | tee -a "$SUMMARY"
fi

echo "PASS: CBQ E2E/module regression suite completed at $OUT_ROOT"
