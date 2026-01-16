#!/bin/bash
set -euo pipefail

# Test: cellrangerStyle + FlexProbeIndex integration
# This test verifies that both features work together

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-$SCRIPT_DIR/../source/STAR}"

if [ ! -x "$STAR_BIN" ]; then
    echo "ERROR: STAR binary not found: $STAR_BIN"
    echo "Build STAR first: cd source && make"
    exit 1
fi

# Use test fixtures if available, otherwise require paths
TEST_FASTA="${TEST_FASTA:-$SCRIPT_DIR/fixtures/cellranger_format/test_genome.fa}"
TEST_GTF="${TEST_GTF:-$SCRIPT_DIR/fixtures/cellranger_format/test_genes.gtf}"
TEST_PROBES="${TEST_PROBES:-}"

if [ ! -f "$TEST_FASTA" ]; then
    echo "ERROR: Test FASTA not found: $TEST_FASTA"
    echo "Set TEST_FASTA environment variable or ensure fixtures exist"
    exit 1
fi

if [ ! -f "$TEST_GTF" ]; then
    echo "ERROR: Test GTF not found: $TEST_GTF"
    echo "Set TEST_GTF environment variable or ensure fixtures exist"
    exit 1
fi

TEST_DIR=$(mktemp -d)
trap "rm -rf $TEST_DIR" EXIT

echo "=============================================="
echo "CellRanger + Flex Integration Smoke Test"
echo "=============================================="
echo "STAR binary: $STAR_BIN"
echo "Test FASTA: $TEST_FASTA"
echo "Test GTF: $TEST_GTF"
echo "Test dir: $TEST_DIR"
echo ""

# Test 1: cellrangerStyle alone
echo ">>> Test 1: cellrangerStyle formatting..."
"$STAR_BIN" --runMode genomeGenerate \
     --cellrangerStyleIndex Yes \
     --genomeFastaFiles "$TEST_FASTA" \
     --sjdbGTFfile "$TEST_GTF" \
     --genomeDir "$TEST_DIR/test_cr_only" \
     --runThreadN 1 \
     --genomeSAindexNbases 4 \
     --genomeChrBinNbits 10 \
     --sjdbOverhang 50 2>&1 | grep -v "^working on" || true

if [ ! -f "$TEST_DIR/test_cr_only/cellranger_ref/genome.fa" ]; then
    echo "    ✗ FAILED: cellranger_ref/genome.fa not created"
    exit 1
fi
if [ ! -f "$TEST_DIR/test_cr_only/cellranger_ref/genes.gtf" ]; then
    echo "    ✗ FAILED: cellranger_ref/genes.gtf not created"
    exit 1
fi
    echo "    ✓ PASSED: cellrangerStyle formatting works"

# Test 2: cellrangerStyle + FlexProbeIndex (if probes available)
if [ -n "$TEST_PROBES" ] && [ -f "$TEST_PROBES" ]; then
    echo ""
    echo ">>> Test 2: cellrangerStyle + FlexProbeIndex..."
    "$STAR_BIN" --runMode genomeGenerate \
         --cellrangerStyleIndex Yes \
         --flexGeneProbeSet "$TEST_PROBES" \
         --genomeFastaFiles "$TEST_FASTA" \
         --sjdbGTFfile "$TEST_GTF" \
         --genomeDir "$TEST_DIR/test_cr_flex" \
         --runThreadN 1 \
         --genomeSAindexNbases 4 \
         --genomeChrBinNbits 10 \
         --sjdbOverhang 50 2>&1 | grep -v "^working on" || true
    
    if [ ! -f "$TEST_DIR/test_cr_flex/cellranger_ref/genome.fa" ]; then
        echo "    ✗ FAILED: cellranger_ref/genome.fa not created"
        exit 1
    fi
    if [ ! -d "$TEST_DIR/test_cr_flex/flex_probe_artifacts" ]; then
        echo "    ✗ FAILED: flex_probe_artifacts directory not created"
        exit 1
    fi
    echo "    ✓ PASSED: cellrangerStyle + FlexProbeIndex integration works"
else
    echo ""
    echo ">>> Test 2: SKIPPED (no probe CSV provided)"
    echo "    Set TEST_PROBES environment variable to enable FlexProbeIndex test"
fi

echo ""
echo "=============================================="
echo "✓ FLEX SMOKE TEST PASSED"
echo "=============================================="

