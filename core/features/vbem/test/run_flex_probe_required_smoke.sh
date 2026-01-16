#!/bin/bash
# Smoke test to verify Flex mode requires probe list

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-$SCRIPT_DIR/../source/STAR}"

# Use test fixtures
TEST_FASTA="${TEST_FASTA:-$SCRIPT_DIR/fixtures/cellranger_format/test_genome.fa}"
TEST_GTF="${TEST_GTF:-$SCRIPT_DIR/fixtures/cellranger_format/test_genes.gtf}"
TEST_PROBE_LIST="${TEST_PROBE_LIST:-$SCRIPT_DIR/fixtures/flex_probe_required/test_probe_list.txt}"
TEST_CB_WHITELIST="${TEST_CB_WHITELIST:-$SCRIPT_DIR/fixtures/flex_probe_required/test_cb_whitelist.txt}"

TEST_DIR=$(mktemp -d)
trap "rm -rf $TEST_DIR" EXIT

if [ ! -x "$STAR_BIN" ]; then
    echo "ERROR: STAR binary not found: $STAR_BIN"
    exit 1
fi

if [ ! -f "$TEST_FASTA" ] || [ ! -f "$TEST_GTF" ] || [ ! -f "$TEST_PROBE_LIST" ] || [ ! -f "$TEST_CB_WHITELIST" ]; then
    echo "ERROR: Test fixtures not found"
    echo "  TEST_FASTA: $TEST_FASTA"
    echo "  TEST_GTF: $TEST_GTF"
    echo "  TEST_PROBE_LIST: $TEST_PROBE_LIST"
    echo "  TEST_CB_WHITELIST: $TEST_CB_WHITELIST"
    exit 1
fi

echo "=============================================="
echo "Flex Probe Requirement Smoke Test"
echo "=============================================="
echo "STAR binary: $STAR_BIN"
echo "Test dir: $TEST_DIR"
echo ""

# First, build an index without probe list (for Test A)
echo ">>> Building test index (without probe list)..."
"$STAR_BIN" --runMode genomeGenerate \
     --genomeFastaFiles "$TEST_FASTA" \
     --sjdbGTFfile "$TEST_GTF" \
     --cellrangerStyleIndex Yes \
     --genomeDir "$TEST_DIR/test_index_no_probe" \
     --runThreadN 1 \
     --genomeSAindexNbases 4 \
     --genomeChrBinNbits 10 \
     --sjdbOverhang 50 2>&1 | grep -v "^working on" || true

# Test A: Flex mode without probe list → should fail
echo ""
echo ">>> Test A: Flex mode without probe list (should fail)..."
mkdir -p "$TEST_DIR/test_a_no_probe"
OUTPUT_A="$(timeout 10 "$STAR_BIN" --runMode alignReads \
     --genomeDir "$TEST_DIR/test_index_no_probe" \
     --readFilesIn /dev/null \
     --soloType CB_UMI_Simple \
     --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
     --soloCBwhitelist "$TEST_CB_WHITELIST" \
     --soloFeatures Gene \
     --flex yes \
     --soloFlexExpectedCellsPerTag 1000 \
     --outFileNamePrefix "$TEST_DIR/test_a_no_probe/" \
     --runThreadN 1 2>&1 || true)"

if echo "$OUTPUT_A" | grep -E "(Flex mode requires probe list|EXITING because of fatal PARAMETERS error.*probe)" > /dev/null; then
    echo "    ✓ PASSED: Correctly failed with probe requirement error"
    echo "$OUTPUT_A" | grep -E "(Flex mode requires|probe)" | head -3
else
    echo "    ✗ FAILED: Should have failed with probe requirement error"
    echo "    Output:"
    echo "$OUTPUT_A" | head -10
    exit 1
fi

# Test B: Flex mode with --soloProbeList → should proceed (will fail later due to no reads, but should pass probe check)
echo ""
echo ">>> Test B: Flex mode with --soloProbeList (should pass probe check)..."
mkdir -p "$TEST_DIR/test_b_with_probe"
OUTPUT_B="$(timeout 10 "$STAR_BIN" --runMode alignReads \
     --genomeDir "$TEST_DIR/test_index_no_probe" \
     --readFilesIn /dev/null \
     --soloType CB_UMI_Simple \
     --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
     --soloCBwhitelist "$TEST_CB_WHITELIST" \
     --soloFeatures Gene \
     --soloProbeList "$TEST_PROBE_LIST" \
     --flex yes \
     --soloFlexExpectedCellsPerTag 1000 \
     --outFileNamePrefix "$TEST_DIR/test_b_with_probe/" \
     --runThreadN 1 2>&1 || true)"

if echo "$OUTPUT_B" | grep -E "(Flex mode requires probe list|EXITING because of fatal PARAMETERS error.*probe)" > /dev/null; then
    echo "    ✗ FAILED: Should have passed probe check"
    echo "    Output:"
    echo "$OUTPUT_B" | head -10
    exit 1
elif echo "$OUTPUT_B" | grep -E "(Flex mode: Using probe list)" > /dev/null; then
    echo "    ✓ PASSED: Probe list accepted, log shows probe list usage"
    echo "$OUTPUT_B" | grep -E "(Flex mode: Using probe list)" | head -2
else
    # Check if it got past the probe check (may fail later due to no reads or other issues)
    if ! echo "$OUTPUT_B" | grep -E "(Flex mode requires probe list|EXITING because of fatal PARAMETERS error.*probe)" > /dev/null; then
        echo "    ✓ PASSED: Probe list accepted (probe check passed, may fail later due to no reads)"
    else
        echo "    ✗ FAILED: Probe check failed unexpectedly"
        echo "$OUTPUT_B" | head -10
        exit 1
    fi
fi

# Test C: Filter-only mode (--soloRunFlexFilter yes with --flex no) without probe list → should NOT fail on probe requirement
echo ""
echo ">>> Test C: Filter-only mode (--soloRunFlexFilter yes, --flex no) without probe list (should pass probe check)..."
mkdir -p "$TEST_DIR/test_c_filter_only"
OUTPUT_C="$(timeout 10 "$STAR_BIN" --runMode alignReads \
     --genomeDir "$TEST_DIR/test_index_no_probe" \
     --readFilesIn /dev/null \
     --soloType CB_UMI_Simple \
     --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
     --soloCBwhitelist "$TEST_CB_WHITELIST" \
     --soloFeatures Gene \
     --flex no \
     --soloRunFlexFilter yes \
     --soloFlexExpectedCellsPerTag 1000 \
     --outFileNamePrefix "$TEST_DIR/test_c_filter_only/" \
     --runThreadN 1 2>&1 || true)"

# Check log file for the message
LOG_C="$TEST_DIR/test_c_filter_only/Log.out"
LOG_MESSAGE=""
if [ -f "$LOG_C" ]; then
    LOG_MESSAGE=$(grep -i "Flex filter-only mode does not require probe list" "$LOG_C" 2>/dev/null || true)
fi

if echo "$OUTPUT_C" | grep -E "(Flex mode requires probe list|EXITING because of fatal PARAMETERS error.*probe)" > /dev/null; then
    echo "    ✗ FAILED: Should NOT require probe list for filter-only mode"
    echo "    Output:"
    echo "$OUTPUT_C" | grep -E "(Flex mode requires|probe|EXITING)" | head -5
    exit 1
elif [ -n "$LOG_MESSAGE" ]; then
    echo "    ✓ PASSED: Filter-only mode correctly does not require probe list"
    echo "    Log message: $LOG_MESSAGE"
elif ! echo "$OUTPUT_C" | grep -E "(Flex mode requires probe list|EXITING because of fatal PARAMETERS error.*probe)" > /dev/null; then
    echo "    ✓ PASSED: Filter-only mode passed probe check (probe requirement not triggered)"
    echo "    (Log message may not appear if run exits early, but probe check passed)"
else
    echo "    ✗ FAILED: Probe check failed unexpectedly for filter-only mode"
    echo "$OUTPUT_C" | head -10
    exit 1
fi

echo ""
echo "=============================================="
echo "✓ FLEX PROBE REQUIREMENT SMOKE TEST PASSED"
echo "=============================================="

