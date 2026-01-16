#!/bin/bash
# Run the SNP threshold auto-estimation unit tests

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="${SCRIPT_DIR}/../source"
TEST_DIR="${SCRIPT_DIR}/slam"

cd "${TEST_DIR}"

echo "Compiling test_snp_threshold.cpp..."
g++ -std=c++17 -O2 \
    -I"${SOURCE_DIR}" \
    -o test_snp_threshold \
    test_snp_threshold.cpp \
    "${SOURCE_DIR}/SlamQuant.cpp" \
    "${SOURCE_DIR}/SlamSolver.cpp"

echo "Running test_snp_threshold..."
./test_snp_threshold

# Clean up
rm -f test_snp_threshold

echo ""
echo "Test completed successfully."

