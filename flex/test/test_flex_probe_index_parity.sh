#!/bin/bash
#
# test_flex_probe_index_parity.sh
# Automated parity test: runs shell script and C++ helper on chr21 fixture,
# then diffs all output files to ensure identical results.
#
# Usage: ./test/test_flex_probe_index_parity.sh
#
# Exit codes:
#   0 - All outputs identical
#   1 - Diff found or test failure
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Test fixtures - adjust paths as needed
PROBE_CSV="${PROBE_CSV:-/mnt/pikachu/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv}"
GTF_PATH="${GTF_PATH:-/home/lhhung/cellranger-9.0.1/external/cellranger_tiny_ref/genes/genes.gtf.gz}"
BASE_FASTA="${BASE_FASTA:-/home/lhhung/cellranger-9.0.1/external/cellranger_tiny_ref/fasta/genome.fa}"

# Tools
SHELL_SCRIPT="$REPO_ROOT/scripts/filter_probes_to_gtf.sh"
CPP_HELPER="$REPO_ROOT/test/test_flex_probe_index"

# Output directories
TEST_DIR=$(mktemp -d)
SCRIPT_OUT="$TEST_DIR/script_output"
HELPER_OUT="$TEST_DIR/helper_output"

trap "rm -rf $TEST_DIR" EXIT

echo "=============================================="
echo "Flex Probe Index Parity Test"
echo "=============================================="
echo "Probe CSV: $PROBE_CSV"
echo "GTF: $GTF_PATH"
echo "Base FASTA: $BASE_FASTA"
echo "Test dir: $TEST_DIR"
echo ""

# Check prerequisites
if [ ! -f "$PROBE_CSV" ]; then
    echo "ERROR: Probe CSV not found: $PROBE_CSV"
    exit 1
fi
if [ ! -f "$GTF_PATH" ]; then
    echo "ERROR: GTF not found: $GTF_PATH"
    exit 1
fi
if [ ! -f "$BASE_FASTA" ]; then
    echo "ERROR: Base FASTA not found: $BASE_FASTA"
    exit 1
fi
if [ ! -x "$SHELL_SCRIPT" ]; then
    echo "ERROR: Shell script not found or not executable: $SHELL_SCRIPT"
    exit 1
fi
if [ ! -x "$CPP_HELPER" ]; then
    echo "ERROR: C++ helper not found. Build with: cd source && make test_flex_probe_index"
    exit 1
fi

# Run shell script
echo ">>> Running shell script..."
mkdir -p "$SCRIPT_OUT"
"$SHELL_SCRIPT" \
    --probe-set "$PROBE_CSV" \
    --gtf "$GTF_PATH" \
    --base-fasta "$BASE_FASTA" \
    --output-dir "$SCRIPT_OUT" \
    --quiet

echo "    Script completed."

# Run C++ helper
echo ">>> Running C++ helper..."
mkdir -p "$HELPER_OUT"
"$CPP_HELPER" \
    "$PROBE_CSV" \
    "$GTF_PATH" \
    "$BASE_FASTA" \
    "$HELPER_OUT" > /dev/null

echo "    Helper completed."

# Compare outputs
echo ""
echo ">>> Comparing outputs..."

FILES_TO_COMPARE=(
    "filtered_probe_set.csv"
    "probes_only.fa"
    "probes_only.gtf"
    "genome.filtered.fa"
    "genes.filtered.gtf"
    "probe_genes_exons.bed"
    "probe_list.txt"
)

PASS_COUNT=0
FAIL_COUNT=0

for f in "${FILES_TO_COMPARE[@]}"; do
    SCRIPT_FILE="$SCRIPT_OUT/$f"
    HELPER_FILE="$HELPER_OUT/$f"
    
    if [ ! -f "$SCRIPT_FILE" ]; then
        echo "    ✗ $f: MISSING in script output"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi
    
    if [ ! -f "$HELPER_FILE" ]; then
        echo "    ✗ $f: MISSING in helper output"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi
    
    if diff -q "$SCRIPT_FILE" "$HELPER_FILE" > /dev/null 2>&1; then
        echo "    ✓ $f: IDENTICAL"
        PASS_COUNT=$((PASS_COUNT + 1))
    else
        echo "    ✗ $f: DIFFERS"
        echo "      First diff:"
        diff "$SCRIPT_FILE" "$HELPER_FILE" | head -10 || true
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
done

# Summary
echo ""
echo "=============================================="
echo "Results: $PASS_COUNT passed, $FAIL_COUNT failed"
echo "=============================================="

if [ "$FAIL_COUNT" -eq 0 ]; then
    echo "✓ PARITY TEST PASSED"
    exit 0
else
    echo "✗ PARITY TEST FAILED"
    exit 1
fi

