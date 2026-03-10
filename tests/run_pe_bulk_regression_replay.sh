#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
Y_VALIDATOR="${SCRIPT_DIR}/validate_bulk_yremove_output.py"
VB_VALIDATOR="${SCRIPT_DIR}/transcriptvb/validate_pe_autodetect_output.py"

BAD_Y_OUTDIR="${BAD_Y_OUTDIR:-/tmp/pe_bulk_benchmark_real_v3/downsampled/integrated}"
BAD_VB_OUTDIR="${BAD_VB_OUTDIR:-/tmp/pe_bulk_benchmark_real_v4/downsampled/integrated}"
GOOD_VB_OUTDIR="${GOOD_VB_OUTDIR:-${SCRIPT_DIR}/transcriptvb_auto_fix_output_20260309_235417}"
GOOD_Y_OUTDIR="${GOOD_Y_OUTDIR:-${SCRIPT_DIR}/transcriptvb_auto_fix_output_20260309_235417}"
BAD_Y_ROUTE_ROOT="${BAD_Y_ROUTE_ROOT:-/tmp/pe_bulk_feature_benchmark_rerun_20260310_010323/downsampled}"
GOOD_Y_ROUTE_ROOT="${GOOD_Y_ROUTE_ROOT:-/tmp/pe_bulk_feature_benchmark_yfix_20260310_015945/downsampled}"
BAD_Y_ROUTE_OUTDIR="${BAD_Y_ROUTE_OUTDIR:-${BAD_Y_ROUTE_ROOT}/integrated}"
GOOD_Y_ROUTE_OUTDIR="${GOOD_Y_ROUTE_OUTDIR:-${GOOD_Y_ROUTE_ROOT}/integrated}"
BAD_Y_ROUTE_REF_DIR="${BAD_Y_ROUTE_REF_DIR:-${BAD_Y_ROUTE_ROOT}/external/y_fastq_split}"
GOOD_Y_ROUTE_REF_DIR="${GOOD_Y_ROUTE_REF_DIR:-${GOOD_Y_ROUTE_ROOT}/external/y_fastq_split}"

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

locate_ref_fastq() {
    local dir="$1"
    local pattern="$2"
    local matches=()
    mapfile -t matches < <(find "${dir}" -maxdepth 1 -type f -name "${pattern}" | sort)
    if [[ "${#matches[@]}" -ne 1 ]]; then
        echo "ERROR: expected exactly one match for ${pattern} under ${dir}, found ${#matches[@]}" >&2
        return 2
    fi
    printf '%s\n' "${matches[0]}"
}

run_reference_y_validation() {
    local outdir="$1"
    local refdir="$2"
    local ref_y_r1 ref_y_r2 ref_noy_r1 ref_noy_r2
    ref_y_r1="$(locate_ref_fastq "${refdir}" '*_R1_Y.fastq*')"
    ref_y_r2="$(locate_ref_fastq "${refdir}" '*_R2_Y.fastq*')"
    ref_noy_r1="$(locate_ref_fastq "${refdir}" '*_R1_noY.fastq*')"
    ref_noy_r2="$(locate_ref_fastq "${refdir}" '*_R2_noY.fastq*')"
    python3 "${Y_VALIDATOR}" \
        --outdir "${outdir}" \
        --require-y-reads \
        --reference-y-r1 "${ref_y_r1}" \
        --reference-y-r2 "${ref_y_r2}" \
        --reference-noy-r1 "${ref_noy_r1}" \
        --reference-noy-r2 "${ref_noy_r2}"
}

echo "=== PE Bulk Regression Replay ==="
echo "Bad Y output: ${BAD_Y_OUTDIR}"
echo "Bad VB output: ${BAD_VB_OUTDIR}"
echo "Good Y output: ${GOOD_Y_OUTDIR}"
echo "Good VB output: ${GOOD_VB_OUTDIR}"
echo "Bad benchmark Y-route output: ${BAD_Y_ROUTE_OUTDIR}"
echo "Bad benchmark Y-route reference: ${BAD_Y_ROUTE_REF_DIR}"
echo "Good benchmark Y-route output: ${GOOD_Y_ROUTE_OUTDIR}"
echo "Good benchmark Y-route reference: ${GOOD_Y_ROUTE_REF_DIR}"
echo

expect_fail \
    "old Y/noY FASTQ emission artifact validation" \
    python3 "${Y_VALIDATOR}" --outdir "${BAD_Y_OUTDIR}"

expect_pass \
    "fixed Y/noY FASTQ emission artifact validation" \
    python3 "${Y_VALIDATOR}" --outdir "${GOOD_Y_OUTDIR}" --require-y-reads

expect_fail \
    "old bulk PE Y-route benchmark validation" \
    run_reference_y_validation "${BAD_Y_ROUTE_OUTDIR}" "${BAD_Y_ROUTE_REF_DIR}"

expect_pass \
    "fixed bulk PE Y-route benchmark validation" \
    run_reference_y_validation "${GOOD_Y_ROUTE_OUTDIR}" "${GOOD_Y_ROUTE_REF_DIR}"

expect_fail \
    "old TranscriptVB auto-detect artifact validation" \
    python3 "${VB_VALIDATOR}" --outdir "${BAD_VB_OUTDIR}" --expected-format IU

expect_pass \
    "fixed TranscriptVB auto-detect artifact validation" \
    python3 "${VB_VALIDATOR}" --outdir "${GOOD_VB_OUTDIR}" --expected-format IU

echo
echo "PASS: archived bad outputs fail and fixed outputs pass"
