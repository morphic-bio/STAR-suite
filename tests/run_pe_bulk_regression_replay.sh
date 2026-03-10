#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
Y_VALIDATOR="${SCRIPT_DIR}/validate_bulk_yremove_output.py"
VB_VALIDATOR="${SCRIPT_DIR}/transcriptvb/validate_pe_autodetect_output.py"

BAD_Y_OUTDIR="${BAD_Y_OUTDIR:-/tmp/pe_bulk_benchmark_real_v3/downsampled/integrated}"
BAD_VB_OUTDIR="${BAD_VB_OUTDIR:-/tmp/pe_bulk_benchmark_real_v4/downsampled/integrated}"
GOOD_VB_OUTDIR="${GOOD_VB_OUTDIR:-${SCRIPT_DIR}/transcriptvb_auto_fix_output_20260309_235417}"
GOOD_Y_OUTDIR="${GOOD_Y_OUTDIR:-${SCRIPT_DIR}/transcriptvb_auto_fix_output_20260309_235417}"

expect_fail() {
    local label="$1"
    shift
    if "$@"; then
        echo "FAIL: ${label} unexpectedly passed" >&2
        return 1
    fi
    echo "PASS: ${label} failed as expected"
}

expect_pass() {
    local label="$1"
    shift
    "$@"
    echo "PASS: ${label}"
}

echo "=== PE Bulk Regression Replay ==="
echo "Bad Y output: ${BAD_Y_OUTDIR}"
echo "Bad VB output: ${BAD_VB_OUTDIR}"
echo "Good Y output: ${GOOD_Y_OUTDIR}"
echo "Good VB output: ${GOOD_VB_OUTDIR}"
echo

expect_fail \
    "old Y/noY FASTQ emission artifact validation" \
    python3 "${Y_VALIDATOR}" --outdir "${BAD_Y_OUTDIR}"

expect_pass \
    "fixed Y/noY FASTQ emission artifact validation" \
    python3 "${Y_VALIDATOR}" --outdir "${GOOD_Y_OUTDIR}" --require-y-reads

expect_fail \
    "old TranscriptVB auto-detect artifact validation" \
    python3 "${VB_VALIDATOR}" --outdir "${BAD_VB_OUTDIR}" --expected-format IU

expect_pass \
    "fixed TranscriptVB auto-detect artifact validation" \
    python3 "${VB_VALIDATOR}" --outdir "${GOOD_VB_OUTDIR}" --expected-format IU

echo
echo "PASS: archived bad outputs fail and fixed outputs pass"
