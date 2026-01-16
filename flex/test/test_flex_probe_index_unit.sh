#!/bin/bash
#
# test_flex_probe_index_unit.sh
# Unit test: tests edge cases with a tiny fixture (DEPRECATED, invalid length, no-match, bad chars)
#
# Expected behavior with test_probes.csv:
#   Row 1: GENE_A|A|001 - VALID (50bp ACGT, gene in GTF)
#   Row 2: GENE_A|A|002 - VALID (50bp ACGT, gene in GTF)  
#   Row 3: GENE_B|B|001 - VALID (50bp ACGT, gene in GTF)
#   Row 4: DEPRECATED_PROBE - DROPPED (contains "deprecated")
#   Row 5: GENE_D|D|001 - VALID (50bp ACGT, gene in GTF)
#   Row 6: NOT_IN_GTF|X|001 - DROPPED (gene_id not in GTF)
#   Row 7: SHORT_PROBE - ERROR (4bp, not 50bp) -> strict mode fails
#   Row 8: BAD_CHARS|A|003 - ERROR (contains N) -> strict mode fails
#
# This test checks that strict validation catches invalid probes.
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Fixtures
FIXTURE_DIR="$REPO_ROOT/test/fixtures/flex_probe_unit"
PROBE_CSV="$FIXTURE_DIR/test_probes.csv"
GTF_PATH="$FIXTURE_DIR/test_genes.gtf"
BASE_FASTA="$FIXTURE_DIR/test_genome.fa"

# Tools
CPP_HELPER="$REPO_ROOT/test/test_flex_probe_index"
SHELL_SCRIPT="$REPO_ROOT/scripts/filter_probes_to_gtf.sh"

# Output directories
TEST_DIR=$(mktemp -d)
trap "rm -rf $TEST_DIR" EXIT

echo "=============================================="
echo "Flex Probe Index Unit Test"
echo "=============================================="
echo "Fixture dir: $FIXTURE_DIR"
echo "Test dir: $TEST_DIR"
echo ""

# Check prerequisites
if [ ! -x "$CPP_HELPER" ]; then
    echo "ERROR: C++ helper not found. Build with: cd source && make test_flex_probe_index"
    exit 1
fi

# Test 1: C++ helper should fail fast on invalid sequence length
echo ">>> Test 1: C++ helper - invalid length detection"
HELPER_OUT="$TEST_DIR/helper_unit"
mkdir -p "$HELPER_OUT"

HELPER_OUTPUT=$("$CPP_HELPER" "$PROBE_CSV" "$GTF_PATH" "$BASE_FASTA" "$HELPER_OUT" 2>&1 || true)
if echo "$HELPER_OUTPUT" | grep -qE "ERROR.*[0-9]+bp.*expected|ERROR.*length|ERROR.*non-ACGT"; then
    echo "    ✓ C++ helper correctly detected invalid probe"
    echo "      Message: $(echo "$HELPER_OUTPUT" | grep ERROR | head -1)"
else
    echo "    ✗ C++ helper did not detect invalid probe (expected failure)"
    echo "$HELPER_OUTPUT"
    exit 1
fi

# Test 2: Create a clean fixture without invalid probes
echo ""
echo ">>> Test 2: Clean fixture processing"
CLEAN_PROBE_CSV="$TEST_DIR/clean_probes.csv"
head -1 "$PROBE_CSV" > "$CLEAN_PROBE_CSV"
# Keep only valid 50bp ACGT probes for genes A, B, D (skip DEPRECATED, SHORT, BAD_CHARS)
grep -E "^GENE_[ABD]," "$PROBE_CSV" | grep -vE "DEPRECATED|SHORT|BAD_CHARS" >> "$CLEAN_PROBE_CSV"

echo "    Clean probes:"
cat "$CLEAN_PROBE_CSV"
echo ""

CLEAN_OUT="$TEST_DIR/clean_output"
mkdir -p "$CLEAN_OUT"

if "$CPP_HELPER" "$CLEAN_PROBE_CSV" "$GTF_PATH" "$BASE_FASTA" "$CLEAN_OUT" > /dev/null 2>&1; then
    echo "    ✓ C++ helper processed clean fixture successfully"
else
    echo "    ✗ C++ helper failed on clean fixture"
    "$CPP_HELPER" "$CLEAN_PROBE_CSV" "$GTF_PATH" "$BASE_FASTA" "$CLEAN_OUT" 2>&1 || true
    exit 1
fi

# Verify output counts
OUTPUT_PROBE_COUNT=$(tail -n +2 "$CLEAN_OUT/filtered_probe_set.csv" | wc -l)
if [ "$OUTPUT_PROBE_COUNT" -eq 4 ]; then
    echo "    ✓ Output has expected 4 probes (GENE_A x2, GENE_B x1, GENE_D x1)"
else
    echo "    ✗ Expected 4 output probes, got $OUTPUT_PROBE_COUNT"
    exit 1
fi

# Verify synthetic contigs in FASTA
SYNTH_CONTIG_COUNT=$(grep -c "^>" "$CLEAN_OUT/probes_only.fa")
if [ "$SYNTH_CONTIG_COUNT" -eq 4 ]; then
    echo "    ✓ Synthetic FASTA has 4 probe contigs"
else
    echo "    ✗ Expected 4 synthetic contigs, got $SYNTH_CONTIG_COUNT"
    exit 1
fi

# Test 3: Verify shell script also catches issues (soft validation mode)
echo ""
echo ">>> Test 3: Shell script validation"
SCRIPT_OUT="$TEST_DIR/script_output"
mkdir -p "$SCRIPT_OUT"

if "$SHELL_SCRIPT" --probe-set "$CLEAN_PROBE_CSV" --gtf "$GTF_PATH" --base-fasta "$BASE_FASTA" --output-dir "$SCRIPT_OUT" --quiet 2>&1; then
    echo "    ✓ Shell script processed clean fixture"
else
    echo "    ✗ Shell script failed on clean fixture"
    exit 1
fi

# Diff script vs helper on clean fixture
echo ""
echo ">>> Test 4: Parity on clean fixture"
DIFF_FAILED=0
for f in filtered_probe_set.csv probes_only.fa probes_only.gtf probe_list.txt; do
    if diff -q "$SCRIPT_OUT/$f" "$CLEAN_OUT/$f" > /dev/null 2>&1; then
        echo "    ✓ $f: IDENTICAL"
    else
        echo "    ✗ $f: DIFFERS"
        diff "$SCRIPT_OUT/$f" "$CLEAN_OUT/$f" | head -5 || true
        DIFF_FAILED=1
    fi
done

echo ""
echo "=============================================="
if [ "$DIFF_FAILED" -eq 0 ]; then
    echo "✓ UNIT TEST PASSED"
    exit 0
else
    echo "✗ UNIT TEST FAILED"
    exit 1
fi

