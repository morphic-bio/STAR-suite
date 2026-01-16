#!/bin/bash
# Quick test to verify SNP mask build segfault fix

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
STAR_BIN="${PROJECT_ROOT}/source/STAR"

# Use production index
STAR_INDEX="${STAR_INDEX:-/storage/autoindex_110_44/bulk_index}"

# Test with 0h FASTQ (same as original run)
WT_FASTQ="/storage/SLAM-Seq-prod-compare-20260109/input/WDHD1-0h-3_S201_R1_001.fastq.gz"

# Test output directory
TEST_DIR="/storage/slam_e2e_test_fix_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$TEST_DIR"/mask

MASK_BED="$TEST_DIR/mask/wt0.mask.bed.gz"
MASK_SUMMARY="$TEST_DIR/mask/wt0.mask.summary.tsv"
FOFN="$TEST_DIR/mask/wt0.fofn"

echo "Testing SNP mask build fix..."
echo "STAR binary: $STAR_BIN"
echo "STAR index: $STAR_INDEX"
echo "Input FASTQ: $WT_FASTQ"
echo "Output dir: $TEST_DIR"
echo

# Create FOFN
echo "$WT_FASTQ" > "$FOFN"

# Run SNP mask build (this should not segfault now)
echo "Running SNP mask build..."
"$STAR_BIN" \
    --runThreadN 4 \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "$WT_FASTQ" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$TEST_DIR/mask/wt0_" \
    --outSAMtype None \
    --slamQuantMode 1 \
    --slamSnpMaskBuildFastqs "$FOFN" \
    --slamSnpMaskBedOut "$MASK_BED" \
    --slamSnpMaskSummaryOut "$MASK_SUMMARY" \
    --slamSnpMaskOnly 1 \
    > "$TEST_DIR/mask_build.log" 2>&1

EXIT_CODE=$?

if [[ $EXIT_CODE -eq 0 ]]; then
    if [[ -f "$MASK_BED" ]]; then
        echo "✅ SUCCESS: SNP mask build completed without segfault!"
        echo "   Output: $MASK_BED"
        echo "   Summary: $MASK_SUMMARY"
        ls -lh "$MASK_BED" "$MASK_SUMMARY" 2>/dev/null || true
        exit 0
    else
        echo "⚠️  WARNING: Process exited successfully but mask file not created"
        echo "   Check log: $TEST_DIR/mask_build.log"
        exit 1
    fi
else
    echo "❌ FAILURE: Process exited with code $EXIT_CODE"
    echo "   Check log: $TEST_DIR/mask_build.log"
    tail -50 "$TEST_DIR/mask_build.log"
    exit 1
fi
