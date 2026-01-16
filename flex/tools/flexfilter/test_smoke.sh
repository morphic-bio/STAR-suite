#!/bin/bash
# Smoke test for run_flexfilter_mex standalone tool
# Uses tests/gold_standard/ fixtures if present; skips gracefully if missing

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
FIXTURE_DIR="$REPO_ROOT/tests/gold_standard"
OUTPUT_DIR="$SCRIPT_DIR/test_output"

# Note: Sample names come from the sample whitelist if provided; otherwise they
# are auto-derived from tag sequences (BC001, BC002, ...). Do not assume BC###.

echo "=== FlexFilter Smoke Test ==="
echo "Repo root: $REPO_ROOT"
echo "Fixture dir: $FIXTURE_DIR"

# Skip-if-missing guard
if [ ! -d "$FIXTURE_DIR/raw" ]; then
    echo "WARNING: Gold standard fixtures not found at $FIXTURE_DIR/raw"
    echo "Skipping smoke test (fixtures required for testing)"
    echo "To run this test, ensure tests/gold_standard/raw/ contains:"
    echo "  - matrix.mtx"
    echo "  - barcodes.tsv"
    echo "  - features.tsv"
    exit 0
fi

# Check for required MEX files
for file in matrix.mtx barcodes.tsv features.tsv; do
    if [ ! -f "$FIXTURE_DIR/raw/$file" ]; then
        echo "ERROR: Missing required file: $FIXTURE_DIR/raw/$file"
        exit 1
    fi
done

echo "Found gold standard fixtures"

# Build the tool if needed
TOOL="$SCRIPT_DIR/run_flexfilter_mex"
if [ ! -f "$TOOL" ]; then
    echo "Building run_flexfilter_mex..."
    make -C "$REPO_ROOT/source" flexfilter
fi

if [ ! -f "$TOOL" ]; then
    echo "ERROR: Failed to build run_flexfilter_mex"
    exit 1
fi

echo "Tool built successfully: $TOOL"

# Clean up previous test output
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# Run the tool
echo "Running FlexFilter on gold standard input..."
"$TOOL" \
    --mex-dir "$FIXTURE_DIR/raw" \
    --total-expected 12000 \
    --output-prefix "$OUTPUT_DIR"

# Verify output
echo ""
echo "=== Verifying Output ==="

PASS=true

# Auto-detect produced sample directories (any immediate subdirectory)
PRODUCED_SAMPLES=($(find "$OUTPUT_DIR" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort))

if [ ${#PRODUCED_SAMPLES[@]} -eq 0 ]; then
    echo "FAIL: No sample directories found in output"
    PASS=false
else
    echo "Found ${#PRODUCED_SAMPLES[@]} sample directories: ${PRODUCED_SAMPLES[*]}"
    
    # Check per-sample directories
    for sample in "${PRODUCED_SAMPLES[@]}"; do
        sample_dir="$OUTPUT_DIR/$sample/Gene/filtered"
        if [ ! -d "$sample_dir" ]; then
            echo "FAIL: Missing filtered directory: $sample_dir"
            PASS=false
            continue
        fi
        
        # Check for MEX files
        for file in matrix.mtx barcodes.tsv features.tsv; do
            if [ ! -f "$sample_dir/$file" ]; then
                echo "FAIL: Missing file: $sample_dir/$file"
                PASS=false
            elif [ ! -s "$sample_dir/$file" ]; then
                echo "FAIL: Empty file: $sample_dir/$file"
                PASS=false
            fi
        done
        
        if [ -f "$sample_dir/matrix.mtx" ] && [ -f "$sample_dir/barcodes.tsv" ]; then
            # Count cells and entries
            n_cells=$(wc -l < "$sample_dir/barcodes.tsv")
            n_entries=$(tail -n +3 "$sample_dir/matrix.mtx" | head -1 | awk '{print $3}')
            echo "OK: $sample - $n_cells cells, $n_entries entries"
        fi
    done
fi

# Check summary file exists
if [ -f "$OUTPUT_DIR/flexfilter_summary.tsv" ]; then
    echo "OK: Summary file generated"
    echo ""
    echo "Summary:"
    cat "$OUTPUT_DIR/flexfilter_summary.tsv"
else
    echo "WARN: No summary file generated"
fi

echo ""
if [ "$PASS" = true ]; then
    echo "=== SMOKE TEST PASSED ==="
    exit 0
else
    echo "=== SMOKE TEST FAILED ==="
    exit 1
fi

