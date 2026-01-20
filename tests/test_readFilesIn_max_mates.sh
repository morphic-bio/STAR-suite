#!/usr/bin/env bash
set -euo pipefail

# Regression test for MAX_N_MATES error handling
# Verifies that STAR fails cleanly (no segfault) when >3 readFilesIn arguments are provided

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAR_BIN="${ROOT_DIR}/core/legacy/source/STAR"
IDX_DIR="${ROOT_DIR}/tests/solo_smoke/ref/star_index"
R1="${ROOT_DIR}/tests/solo_smoke/fastq/R1.fastq"
R2="${ROOT_DIR}/tests/solo_smoke/fastq/R2.fastq"
WL="${ROOT_DIR}/tests/solo_smoke/whitelist.txt"

if [[ ! -x "${STAR_BIN}" ]]; then
    echo "ERROR: STAR binary not found at ${STAR_BIN}. Build STAR first." >&2
    exit 1
fi

if [[ ! -d "${IDX_DIR}" ]]; then
    echo "ERROR: Genome index not found at ${IDX_DIR}" >&2
    exit 1
fi

OUT_DIR="/tmp/test_max_mates"
rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

echo "=== Testing MAX_N_MATES error handling ==="
echo "Providing 4 readFilesIn arguments (exceeds MAX_N_MATES=3)..."
echo ""

# Run STAR with 4 readFilesIn arguments (should fail with error, not segfault)
set +e
"${STAR_BIN}" \
    --runThreadN 1 \
    --genomeDir "${IDX_DIR}" \
    --readFilesIn "${R2}" "${R1}" "${R2}" "${R1}" \
    --outFileNamePrefix "${OUT_DIR}/" \
    --outSAMtype None \
    --soloType CB_UMI_Simple \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloBarcodeReadLength 0 --soloCBwhitelist "${WL}" \
    --soloFeatures Gene --soloCellFilter None \
    2>&1 | tee "${OUT_DIR}/star.log"
exit_code=$?
set -e

echo ""
echo "=== Checking results ==="

# Check if STAR exited with error code (expected)
if [[ ${exit_code} -eq 0 ]]; then
    echo "ERROR: STAR succeeded unexpectedly (should have failed with >3 mates)" >&2
    exit 1
fi

# Check if the error message contains the expected text
if ! grep -q "exceeds MAX_N_MATES" "${OUT_DIR}/star.log"; then
    echo "ERROR: Expected error message 'exceeds MAX_N_MATES' not found in log" >&2
    echo "Log contents:" >&2
    cat "${OUT_DIR}/star.log" >&2
    exit 1
fi

# Check if the error message mentions the solution
if ! grep -q "comma-separated per mate" "${OUT_DIR}/star.log"; then
    echo "ERROR: Expected solution text 'comma-separated per mate' not found in log" >&2
    exit 1
fi

# Check if the error message mentions index reads
if ! grep -q "Do not include index reads" "${OUT_DIR}/star.log"; then
    echo "ERROR: Expected warning about index reads not found in log" >&2
    exit 1
fi

# Verify no segfault occurred (check for ASAN/SEGV messages)
if grep -qiE "(segmentation|segfault|SIGSEGV|AddressSanitizer.*SEGV)" "${OUT_DIR}/star.log"; then
    echo "ERROR: Segfault detected! This should not happen with the new error handling." >&2
    exit 1
fi

echo "✅ Test PASSED: STAR failed cleanly with expected error message (no segfault)"
echo ""
echo "Error message excerpt:"
grep -A3 "exceeds MAX_N_MATES" "${OUT_DIR}/star.log" | head -5

exit 0
