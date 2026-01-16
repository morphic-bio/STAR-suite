#!/bin/bash
# Parity test script: runs trimvalidate on all synthetic fixtures and compares with expected outputs
#
# This script tests the synthetic fixtures under test/fixtures/trim/.
# For integration testing with real-world data, see test/integration/trim/nfcore_smoke/ and nfcore_real/
#
# IMPORTANT: Expected outputs are generated using Trim Galore (cutadapt backend).
# The adapter matching algorithm in source/libtrim/adapter_trim.cpp is synchronized with cutadapt v5.1.
# If cutadapt is upgraded, regenerate expected fixtures and re-run this test suite.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIXTURES_DIR="$SCRIPT_DIR/../../test/fixtures/trim"
INTEGRATION_DIR="$SCRIPT_DIR/../../test/integration/trim"
VALIDATOR="$SCRIPT_DIR/trimvalidate"

if [ ! -f "$VALIDATOR" ]; then
    echo "Error: trimvalidate not found. Build it first with 'make'"
    exit 1
fi

if [ ! -d "$FIXTURES_DIR" ]; then
    echo "Error: Fixtures directory not found: $FIXTURES_DIR"
    exit 1
fi

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

PASSED=0
FAILED=0

echo "Running synthetic fixture tests..."
echo ""

for fixture in "$FIXTURES_DIR"/*/; do
    if [ ! -d "$fixture" ]; then
        continue
    fi
    
    name=$(basename "$fixture")
    echo -n "Testing: $name ... "
    
    if [ ! -f "$fixture/input_R1.fastq" ] || [ ! -f "$fixture/input_R2.fastq" ]; then
        echo "SKIP (missing input files)"
        continue
    fi
    
    if [ ! -f "$fixture/expected_R1.fastq" ] || [ ! -f "$fixture/expected_R2.fastq" ]; then
        echo "SKIP (missing expected files)"
        continue
    fi
    
    # Run trimvalidate
    "$VALIDATOR" -1 "$fixture/input_R1.fastq" \
                 -2 "$fixture/input_R2.fastq" \
                 -o1 "$TMPDIR/out_R1.fastq" \
                 -o2 "$TMPDIR/out_R2.fastq" \
                 --quality 20 --length 20 > /dev/null 2>&1
    
    # Compare outputs
    if diff -q "$TMPDIR/out_R1.fastq" "$fixture/expected_R1.fastq" > /dev/null 2>&1 && \
       diff -q "$TMPDIR/out_R2.fastq" "$fixture/expected_R2.fastq" > /dev/null 2>&1; then
        echo "PASS"
        PASSED=$((PASSED + 1))
    else
        echo "FAIL"
        echo "  Differences in R1:"
        diff "$TMPDIR/out_R1.fastq" "$fixture/expected_R1.fastq" | head -5 || true
        echo "  Differences in R2:"
        diff "$TMPDIR/out_R2.fastq" "$fixture/expected_R2.fastq" | head -5 || true
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "Synthetic fixture summary: $PASSED passed, $FAILED failed"
echo ""

# Run integration tests if available
INTEGRATION_PASSED=0
INTEGRATION_FAILED=0

# Test nfcore_smoke
if [ -d "$INTEGRATION_DIR/nfcore_smoke" ] && \
   [ -f "$INTEGRATION_DIR/nfcore_smoke/input_R1.fastq" ] && \
   [ -f "$INTEGRATION_DIR/nfcore_smoke/input_R2.fastq" ] && \
   [ -f "$INTEGRATION_DIR/nfcore_smoke/expected_R1.fastq" ] && \
   [ -f "$INTEGRATION_DIR/nfcore_smoke/expected_R2.fastq" ]; then
    echo "Running integration test (nfcore_smoke)..."
    echo ""
    
    # Create results directory
    mkdir -p "$INTEGRATION_DIR/nfcore_smoke/results"
    
    # Run trimvalidate
    "$VALIDATOR" -1 "$INTEGRATION_DIR/nfcore_smoke/input_R1.fastq" \
                 -2 "$INTEGRATION_DIR/nfcore_smoke/input_R2.fastq" \
                 -o1 "$TMPDIR/nfcore_smoke_out_R1.fastq" \
                 -o2 "$TMPDIR/nfcore_smoke_out_R2.fastq" \
                 --quality 20 --length 20 > /dev/null 2>&1
    
    # Compare and store diffs
    diff -u "$INTEGRATION_DIR/nfcore_smoke/expected_R1.fastq" "$TMPDIR/nfcore_smoke_out_R1.fastq" > "$INTEGRATION_DIR/nfcore_smoke/results/diff_R1.txt" 2>&1 || true
    diff -u "$INTEGRATION_DIR/nfcore_smoke/expected_R2.fastq" "$TMPDIR/nfcore_smoke_out_R2.fastq" > "$INTEGRATION_DIR/nfcore_smoke/results/diff_R2.txt" 2>&1 || true
    
    # Check if test passed
    if diff -q "$INTEGRATION_DIR/nfcore_smoke/expected_R1.fastq" "$TMPDIR/nfcore_smoke_out_R1.fastq" > /dev/null 2>&1 && \
       diff -q "$INTEGRATION_DIR/nfcore_smoke/expected_R2.fastq" "$TMPDIR/nfcore_smoke_out_R2.fastq" > /dev/null 2>&1; then
        echo "Integration test (nfcore_smoke): PASS"
        echo "PASS" > "$INTEGRATION_DIR/nfcore_smoke/results/status.txt"
        INTEGRATION_PASSED=$((INTEGRATION_PASSED + 1))
    else
        echo "Integration test (nfcore_smoke): FAIL"
        echo "FAIL" > "$INTEGRATION_DIR/nfcore_smoke/results/status.txt"
        INTEGRATION_FAILED=$((INTEGRATION_FAILED + 1))
        echo "  Differences stored in: $INTEGRATION_DIR/nfcore_smoke/results/diff_R*.txt"
        echo "  Differences in R1:"
        head -10 "$INTEGRATION_DIR/nfcore_smoke/results/diff_R1.txt" || true
        echo "  Differences in R2:"
        head -10 "$INTEGRATION_DIR/nfcore_smoke/results/diff_R2.txt" || true
    fi
    echo ""
fi

# Test nfcore_real
if [ -d "$INTEGRATION_DIR/nfcore_real" ] && \
   [ -f "$INTEGRATION_DIR/nfcore_real/input_R1.fastq" ] && \
   [ -f "$INTEGRATION_DIR/nfcore_real/input_R2.fastq" ] && \
   [ -f "$INTEGRATION_DIR/nfcore_real/expected_R1.fastq" ] && \
   [ -f "$INTEGRATION_DIR/nfcore_real/expected_R2.fastq" ]; then
    echo "Running integration test (nfcore_real)..."
    echo ""
    
    # Create results directory
    mkdir -p "$INTEGRATION_DIR/nfcore_real/results"
    
    # Run trimvalidate
    "$VALIDATOR" -1 "$INTEGRATION_DIR/nfcore_real/input_R1.fastq" \
                 -2 "$INTEGRATION_DIR/nfcore_real/input_R2.fastq" \
                 -o1 "$TMPDIR/nfcore_real_out_R1.fastq" \
                 -o2 "$TMPDIR/nfcore_real_out_R2.fastq" \
                 --quality 20 --length 20 > /dev/null 2>&1
    
    # Compare and store diffs
    diff -u "$INTEGRATION_DIR/nfcore_real/expected_R1.fastq" "$TMPDIR/nfcore_real_out_R1.fastq" > "$INTEGRATION_DIR/nfcore_real/results/diff_R1.txt" 2>&1 || true
    diff -u "$INTEGRATION_DIR/nfcore_real/expected_R2.fastq" "$TMPDIR/nfcore_real_out_R2.fastq" > "$INTEGRATION_DIR/nfcore_real/results/diff_R2.txt" 2>&1 || true
    
    # Check if test passed
    if diff -q "$INTEGRATION_DIR/nfcore_real/expected_R1.fastq" "$TMPDIR/nfcore_real_out_R1.fastq" > /dev/null 2>&1 && \
       diff -q "$INTEGRATION_DIR/nfcore_real/expected_R2.fastq" "$TMPDIR/nfcore_real_out_R2.fastq" > /dev/null 2>&1; then
        echo "Integration test (nfcore_real): PASS"
        echo "PASS" > "$INTEGRATION_DIR/nfcore_real/results/status.txt"
        INTEGRATION_PASSED=$((INTEGRATION_PASSED + 1))
    else
        echo "Integration test (nfcore_real): FAIL"
        echo "FAIL" > "$INTEGRATION_DIR/nfcore_real/results/status.txt"
        INTEGRATION_FAILED=$((INTEGRATION_FAILED + 1))
        echo "  Differences stored in: $INTEGRATION_DIR/nfcore_real/results/diff_R*.txt"
        echo "  Differences in R1:"
        head -10 "$INTEGRATION_DIR/nfcore_real/results/diff_R1.txt" || true
        echo "  Differences in R2:"
        head -10 "$INTEGRATION_DIR/nfcore_real/results/diff_R2.txt" || true
    fi
    echo ""
fi

# Overall summary
TOTAL_PASSED=$((PASSED + INTEGRATION_PASSED))
TOTAL_FAILED=$((FAILED + INTEGRATION_FAILED))
echo "Overall summary: $TOTAL_PASSED passed, $TOTAL_FAILED failed"
if [ $INTEGRATION_PASSED -gt 0 ] || [ $INTEGRATION_FAILED -gt 0 ]; then
    echo "  Synthetic fixtures: $PASSED passed, $FAILED failed"
    echo "  Integration tests: $INTEGRATION_PASSED passed, $INTEGRATION_FAILED failed"
fi
echo ""

if [ $TOTAL_FAILED -gt 0 ]; then
    exit 1
fi

exit 0
