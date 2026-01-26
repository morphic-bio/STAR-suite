#!/bin/bash
# test_assignbarcodes_regression.sh
# 
# Regression test for assignBarcodes that compares output against a saved baseline.
# This catches byte-level changes in output files to detect regressions.
#
# Usage: ./test_assignbarcodes_regression.sh [--update-baseline]
#
# Options:
#   --update-baseline   Regenerate baseline instead of comparing

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIXTURE_DIR="${SCRIPT_DIR}/fixtures/assignbarcodes_baseline"
INPUT_DIR="${FIXTURE_DIR}/input"
EXPECTED_DIR="${FIXTURE_DIR}/expected/input"
ACTUAL_DIR="${FIXTURE_DIR}/actual"
PROCESS_FEATURES_DIR="${SCRIPT_DIR}/.."
ASSIGNBARCODES="${PROCESS_FEATURES_DIR}/assignBarcodes"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Files to compare (text files that should be byte-identical)
COMPARE_FILES=(
    "barcodes.txt"
    "features.txt"
    "matrix.mtx"
    "stats.txt"
    "feature_per_cell.csv"
    "deduped_counts_histograms.txt"
)

# Check if updating baseline
UPDATE_BASELINE=false
if [ "${1:-}" = "--update-baseline" ]; then
    UPDATE_BASELINE=true
fi

# Clean up actual dir
rm -rf "${ACTUAL_DIR}"
mkdir -p "${ACTUAL_DIR}"

echo "=============================================="
echo "assignBarcodes Regression Test"
echo "=============================================="
echo ""

# Check that assignBarcodes exists
if [ ! -x "${ASSIGNBARCODES}" ]; then
    echo -e "${RED}ERROR: assignBarcodes not found at ${ASSIGNBARCODES}${NC}"
    echo "Run 'make' in ${PROCESS_FEATURES_DIR} first."
    exit 1
fi

# Run assignBarcodes
echo "Running assignBarcodes on fixture data..."
"${ASSIGNBARCODES}" \
    -w "${INPUT_DIR}/whitelist.txt" \
    -f "${INPUT_DIR}/features.csv" \
    -d "${ACTUAL_DIR}" \
    "${INPUT_DIR}" \
    -b 16 -u 12 2>&1 | grep -v "^$" | head -20

echo ""

# Get actual output subdirectory (matches input directory name)
ACTUAL_OUTPUT="${ACTUAL_DIR}/input"

if [ ! -d "${ACTUAL_OUTPUT}" ]; then
    echo -e "${RED}ERROR: Expected output directory not created: ${ACTUAL_OUTPUT}${NC}"
    exit 1
fi

# Update baseline if requested
if [ "$UPDATE_BASELINE" = true ]; then
    echo -e "${YELLOW}Updating baseline...${NC}"
    rm -rf "${EXPECTED_DIR}"
    mkdir -p "$(dirname "${EXPECTED_DIR}")"
    cp -r "${ACTUAL_OUTPUT}" "${EXPECTED_DIR}"
    echo -e "${GREEN}Baseline updated at: ${EXPECTED_DIR}${NC}"
    exit 0
fi

# Compare files
echo "Comparing outputs against baseline..."
echo ""

TOTAL_DIFFS=0
TOTAL_FILES=0

compare_file() {
    local rel_path="$1"
    local expected_file="${EXPECTED_DIR}/${rel_path}"
    local actual_file="${ACTUAL_OUTPUT}/${rel_path}"
    
    TOTAL_FILES=$((TOTAL_FILES + 1))
    
    if [ ! -f "${expected_file}" ]; then
        echo -e "  ${YELLOW}SKIP${NC} ${rel_path} (no baseline)"
        return 0
    fi
    
    if [ ! -f "${actual_file}" ]; then
        echo -e "  ${RED}FAIL${NC} ${rel_path} (missing in output)"
        TOTAL_DIFFS=$((TOTAL_DIFFS + 1))
        return 1
    fi
    
    if diff -q "${expected_file}" "${actual_file}" > /dev/null 2>&1; then
        echo -e "  ${GREEN}PASS${NC} ${rel_path}"
        return 0
    else
        echo -e "  ${RED}FAIL${NC} ${rel_path} (differs from baseline)"
        TOTAL_DIFFS=$((TOTAL_DIFFS + 1))
        
        # Show first few differences
        echo "       First differences:"
        diff "${expected_file}" "${actual_file}" 2>&1 | head -10 | sed 's/^/         /'
        return 1
    fi
}

for file in "${COMPARE_FILES[@]}"; do
    compare_file "${file}" || true
done

echo ""
echo "=============================================="
echo "Summary"
echo "=============================================="
echo "Files compared: ${TOTAL_FILES}"
echo "Differences:    ${TOTAL_DIFFS}"

if [ "${TOTAL_DIFFS}" -eq 0 ]; then
    echo -e "${GREEN}All files match baseline!${NC}"
    exit 0
else
    echo -e "${RED}${TOTAL_DIFFS} file(s) differ from baseline${NC}"
    echo ""
    echo "To update the baseline (if changes are intentional):"
    echo "  $0 --update-baseline"
    exit 1
fi
