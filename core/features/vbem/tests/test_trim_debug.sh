#!/bin/bash
# Quick test script to run STAR with trimming debug instrumentation
# Tests on nfcore_real dataset to compare with Trim Galore results

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${SCRIPT_DIR}/../source/STAR"
GENOME_DIR="/storage/flex_filtered_reference/star_index"
TEST_DIR="/tmp/star_trim_debug_test"
INPUT_R1="${SCRIPT_DIR}/../test/integration/trim/nfcore_real/input_R1.fastq"
INPUT_R2="${SCRIPT_DIR}/../test/integration/trim/nfcore_real/input_R2.fastq"

# Clean up old test dir
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"

echo "=========================================="
echo "STAR Trimming Debug Test"
echo "=========================================="
echo "Dataset: nfcore_real (same as Trim Galore test)"
echo "Debug logging: First 100 pairs"
echo "Output directory: $TEST_DIR"
echo ""

# Set debug environment variable
export STAR_TRIM_DEBUG_N=100

# Run STAR with trimming
echo "Running STAR with trimming..."
"$STAR_BIN" \
    --runThreadN 4 \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$INPUT_R1" "$INPUT_R2" \
    --trimCutadapt Yes \
    --trimCutadaptQuality 20 \
    --trimCutadaptMinLength 20 \
    --trimCutadaptAdapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix "$TEST_DIR/" \
    --outTmpDir "$TEST_DIR/tmp" \
    2>&1 | tee "$TEST_DIR/star.log"

echo ""
echo "=========================================="
echo "Results Summary"
echo "=========================================="

# Extract key stats
if [ -f "$TEST_DIR/Log.final.out" ]; then
    echo ""
    echo "From Log.final.out:"
    grep -E "Input reads|Pairs processed|Pairs dropped|Pairs kept|Reads dropped.*minimum length|Number of reads unmapped: too short" "$TEST_DIR/Log.final.out" || true
fi

echo ""
echo "=========================================="
echo "Debug Logs (QUAL_CHECK, TRIM_DEBUG, TRIM_DROP, MAP_SHORT)"
echo "=========================================="
echo ""

# Extract debug logs
if [ -f "$TEST_DIR/Log.out" ]; then
    echo "QUAL_CHECK logs (quality encoding before trimming):"
    grep "QUAL_CHECK:" "$TEST_DIR/Log.out" | head -10
    echo ""
    
    echo "TRIM_DEBUG logs (trimming details):"
    grep "TRIM_DEBUG:" "$TEST_DIR/Log.out" | head -10
    echo ""
    
    echo "TRIM_DROP logs (pairs dropped by trimming):"
    grep "TRIM_DROP:" "$TEST_DIR/Log.out" | head -10
    echo ""
    
    echo "MAP_SHORT logs (reads unmapped due to mapping filters):"
    grep "MAP_SHORT:" "$TEST_DIR/Log.out" | head -10
    echo ""
fi

echo ""
echo "Full logs available in:"
echo "  - $TEST_DIR/Log.out (main log with debug output)"
echo "  - $TEST_DIR/Log.final.out (final statistics)"
echo "  - $TEST_DIR/star.log (STDERR output)"

