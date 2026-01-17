#!/bin/bash
# Run STAR-Slam on the default fixture and compare against GRAND-SLAM reference.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

SLAM_FIXTURE_ROOT="${SLAM_FIXTURE_ROOT:-$ROOT_DIR/test/fixtures/slam}"
SLAM_WORK="${SLAM_WORK:-$ROOT_DIR/test/tmp_slam_fixture}"
STAR_BIN="${STAR_BIN:-$ROOT_DIR/core/legacy/source/STAR}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

FASTQ="${FASTQ:-$SLAM_FIXTURE_ROOT/raw/slam_100000_reads_SRR32576116.fastq.gz}"
STAR_INDEX="${STAR_INDEX:-$SLAM_FIXTURE_ROOT/ref/star_index}"
SNPS_BED="${SNPS_BED:-$SLAM_FIXTURE_ROOT/ref/snps.bed}"
REF_TSV="${REF_TSV:-$SLAM_FIXTURE_ROOT/expected/fixture_ref_human.tsv.gz}"
OUT_PREFIX="${OUT_PREFIX:-$SLAM_WORK/star_slam_}"
SLAM_OUT="${SLAM_OUT:-${OUT_PREFIX}SlamQuant.out}"
GRAND_SLAM_OUT="${GRAND_SLAM_OUT:-${OUT_PREFIX}SlamQuant.grandslam.tsv}"

COMPARE_SCRIPT="$SCRIPT_DIR/slam/compare_fixture.py"

mkdir -p "$SLAM_WORK"

if [[ ! -x "$STAR_BIN" ]]; then
    echo "FAIL: STAR binary not found: $STAR_BIN"
    exit 1
fi
if [[ ! -f "$FASTQ" ]]; then
    echo "FAIL: fixture FASTQ not found: $FASTQ"
    exit 1
fi
if [[ ! -d "$STAR_INDEX" ]]; then
    echo "FAIL: STAR index not found: $STAR_INDEX"
    exit 1
fi
if [[ ! -f "$SNPS_BED" ]]; then
    echo "FAIL: SNP BED not found: $SNPS_BED"
    exit 1
fi
if [[ ! -f "$REF_TSV" ]]; then
    echo "FAIL: reference fixture not found: $REF_TSV"
    exit 1
fi
if [[ ! -f "$COMPARE_SCRIPT" ]]; then
    echo "FAIL: compare script not found: $COMPARE_SCRIPT"
    exit 1
fi

if [[ "${RUN_STAR_SLAM:-0}" -eq 1 ]]; then
    if [[ -z "${STAR_SLAM_ARGS:-}" ]]; then
        echo "FAIL: RUN_STAR_SLAM=1 requires STAR_SLAM_ARGS (slam flags)"
        exit 1
    fi
    echo "=== Running STAR-Slam on fixture ==="
    # Use adapter clipping (same as original fixture generation)
    # NOT EndToEnd alignment - that was incorrect
    "$STAR_BIN" \
        --runThreadN 4 \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$FASTQ" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$OUT_PREFIX" \
        --outSAMtype None \
        --clip3pAdapterSeq AGATCGGAAGAG \
        --clip3pAdapterMMp 0.1 \
        ${STAR_SLAM_ARGS} \
        > "${OUT_PREFIX}slam.log" 2>&1
fi

if [[ ! -f "$SLAM_OUT" ]]; then
    echo "FAIL: STAR-Slam output not found: $SLAM_OUT"
    echo "Set RUN_STAR_SLAM=1 and STAR_SLAM_ARGS to generate it."
    exit 1
fi

if [[ ! -f "$GRAND_SLAM_OUT" ]]; then
    echo "FAIL: GRAND-SLAM-style output not found: $GRAND_SLAM_OUT"
    exit 1
fi

prefix_base="${OUT_PREFIX##*/}"
while [[ -n "$prefix_base" && ( "${prefix_base: -1}" == "_" || "${prefix_base: -1}" == "." ) ]]; do
    prefix_base="${prefix_base%?}"
done
if [[ -z "$prefix_base" ]]; then
    prefix_base="STAR"
fi
expected_header=$(printf "Gene\tSymbol\t%s Readcount\t%s 0.05 quantile\t%s Mean\t%s MAP\t%s 0.95 quantile\t%s alpha\t%s beta\t%s Conversions\t%s Coverage\t%s Double-Hits\t%s Double-Hit Coverage\t%s min2\tLength" \
    "$prefix_base" "$prefix_base" "$prefix_base" "$prefix_base" "$prefix_base" "$prefix_base" "$prefix_base" "$prefix_base" "$prefix_base" "$prefix_base" "$prefix_base" "$prefix_base")
actual_header="$(head -n 1 "$GRAND_SLAM_OUT")"
if [[ "$actual_header" != "$expected_header" ]]; then
    echo "FAIL: GRAND-SLAM header mismatch"
    echo "Expected: $expected_header"
    echo "Actual:   $actual_header"
    exit 1
fi

if [[ "${DIRECT_EM_COMPARE:-0}" -eq 1 ]]; then
    echo "=== Comparing STAR-Slam vs GRAND-SLAM ==="
    "$PYTHON_BIN" "$COMPARE_SCRIPT" \
        --reference "$REF_TSV" \
        --test "$SLAM_OUT"
else
    echo "Skipping STAR-Slam vs GRAND-SLAM comparison (set DIRECT_EM_COMPARE=1 to enable)."
fi
