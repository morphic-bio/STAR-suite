#!/bin/bash
# Compare MEX outputs between baseline and current runs
# Usage: compare_mex_outputs.sh <baseline_dir> <current_dir>

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <baseline_dir> <current_dir>" >&2
    exit 1
fi

BASELINE_DIR="$1"
CURRENT_DIR="$2"

TOTAL_DIFFS=0
TOTAL_FILES=0

# Function to compare a file and report differences
compare_file() {
    local BASELINE_FILE="$1"
    local CURRENT_FILE="$2"
    local REL_PATH="$3"
    
    TOTAL_FILES=$((TOTAL_FILES + 1))
    
    if [ ! -f "$BASELINE_FILE" ]; then
        echo "  ✗ Missing in baseline: $REL_PATH"
        TOTAL_DIFFS=$((TOTAL_DIFFS + 1))
        return 1
    fi
    
    if [ ! -f "$CURRENT_FILE" ]; then
        echo "  ✗ Missing in current: $REL_PATH"
        TOTAL_DIFFS=$((TOTAL_DIFFS + 1))
        return 1
    fi
    
    if diff -q "$BASELINE_FILE" "$CURRENT_FILE" > /dev/null 2>&1; then
        echo "  ✓ $REL_PATH"
        return 0
    else
        echo "  ✗ $REL_PATH differs"
        TOTAL_DIFFS=$((TOTAL_DIFFS + 1))
        
        # Show first few differences for matrix.mtx (sparse format)
        if [[ "$REL_PATH" == *"matrix.mtx" ]]; then
            # For matrix.mtx, show header info and first few data differences
            BASELINE_HEADER=$(head -2 "$BASELINE_FILE" | tail -1)
            CURRENT_HEADER=$(head -2 "$CURRENT_FILE" | tail -1)
            if [ "$BASELINE_HEADER" != "$CURRENT_HEADER" ]; then
                echo "    Header differs:"
                echo "      Baseline: $BASELINE_HEADER"
                echo "      Current:  $CURRENT_HEADER"
            fi
            echo "    First data differences:"
            diff "$BASELINE_FILE" "$CURRENT_FILE" | grep -E "^[<>]" | head -5 | sed 's/^/      /'
        else
            # For TSV files, show line count and first few differences
            BASELINE_LINES=$(wc -l < "$BASELINE_FILE")
            CURRENT_LINES=$(wc -l < "$CURRENT_FILE")
            if [ "$BASELINE_LINES" -ne "$CURRENT_LINES" ]; then
                echo "    Line count: baseline=$BASELINE_LINES current=$CURRENT_LINES"
            fi
            echo "    First differences:"
            diff "$BASELINE_FILE" "$CURRENT_FILE" | grep -E "^[<>]" | head -5 | sed 's/^/      /'
        fi
        return 1
    fi
}

echo "Comparing MEX outputs..."
echo ""

# Compare pooled raw outputs
echo "=== Pooled Raw Outputs ==="
compare_file \
    "${BASELINE_DIR}/Solo.out/Gene/raw/matrix.mtx" \
    "${CURRENT_DIR}/Solo.out/Gene/raw/matrix.mtx" \
    "Solo.out/Gene/raw/matrix.mtx"

compare_file \
    "${BASELINE_DIR}/Solo.out/Gene/raw/barcodes.tsv" \
    "${CURRENT_DIR}/Solo.out/Gene/raw/barcodes.tsv" \
    "Solo.out/Gene/raw/barcodes.tsv"

compare_file \
    "${BASELINE_DIR}/Solo.out/Gene/raw/features.tsv" \
    "${CURRENT_DIR}/Solo.out/Gene/raw/features.tsv" \
    "Solo.out/Gene/raw/features.tsv"

echo ""

# Compare per-sample outputs
echo "=== Per-Sample Outputs ==="
for SAMPLE in BC004 BC006 BC007 BC008; do
    echo "Sample: $SAMPLE"
    
    compare_file \
        "${BASELINE_DIR}/per_sample/${SAMPLE}/Gene/filtered/matrix.mtx" \
        "${CURRENT_DIR}/per_sample/${SAMPLE}/Gene/filtered/matrix.mtx" \
        "per_sample/${SAMPLE}/Gene/filtered/matrix.mtx"
    
    compare_file \
        "${BASELINE_DIR}/per_sample/${SAMPLE}/Gene/filtered/barcodes.tsv" \
        "${CURRENT_DIR}/per_sample/${SAMPLE}/Gene/filtered/barcodes.tsv" \
        "per_sample/${SAMPLE}/Gene/filtered/barcodes.tsv"
    
    compare_file \
        "${BASELINE_DIR}/per_sample/${SAMPLE}/Gene/filtered/features.tsv" \
        "${CURRENT_DIR}/per_sample/${SAMPLE}/Gene/filtered/features.tsv" \
        "per_sample/${SAMPLE}/Gene/filtered/features.tsv"
    
    echo ""
done

# Summary
echo "=== Summary ==="
echo "Files compared: $TOTAL_FILES"
echo "Differences found: $TOTAL_DIFFS"

if [ $TOTAL_DIFFS -eq 0 ]; then
    echo "✓ All MEX outputs match!"
    exit 0
else
    echo "⚠ $TOTAL_DIFFS file(s) differ"
    exit 1
fi

