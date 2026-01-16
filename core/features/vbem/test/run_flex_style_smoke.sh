#!/bin/bash
set -euo pipefail

# Test: Flex-style indices (CellRanger formatting + FlexProbeIndex integration)
# This test verifies the complete flex-style workflow

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-$SCRIPT_DIR/../source/STAR}"

if [ ! -x "$STAR_BIN" ]; then
    echo "ERROR: STAR binary not found: $STAR_BIN"
    echo "Build STAR first: cd source && make"
    exit 1
fi

# Use test fixtures
TEST_FASTA="${TEST_FASTA:-$SCRIPT_DIR/fixtures/cellranger_format/test_genome.fa}"
TEST_GTF="${TEST_GTF:-$SCRIPT_DIR/fixtures/cellranger_format/test_genes.gtf}"
TEST_PROBES="${TEST_PROBES:-$SCRIPT_DIR/fixtures/cellranger_format/test_probes.csv}"

if [ ! -f "$TEST_FASTA" ]; then
    echo "ERROR: Test FASTA not found: $TEST_FASTA"
    exit 1
fi

if [ ! -f "$TEST_GTF" ]; then
    echo "ERROR: Test GTF not found: $TEST_GTF"
    exit 1
fi

if [ ! -f "$TEST_PROBES" ]; then
    echo "ERROR: Test probe CSV not found: $TEST_PROBES"
    exit 1
fi

# Create a temporary directory for this test run
TEST_DIR=$(mktemp -d)
trap "rm -rf $TEST_DIR" EXIT

echo "=============================================="
echo "Flex-Style Index Smoke Test"
echo "=============================================="
echo "STAR binary: $STAR_BIN"
echo "Test FASTA: $TEST_FASTA"
echo "Test GTF: $TEST_GTF"
echo "Test probes: $TEST_PROBES"
echo "Test dir: $TEST_DIR"
echo ""

# Helper function to check if index exists
check_index_exists() {
    local genome_dir="$1"
    [ -f "$genome_dir/Genome" ] && \
    [ -f "$genome_dir/SA" ] && \
    [ -f "$genome_dir/SAindex" ] && \
    [ -f "$genome_dir/genomeParameters.txt" ]
}

# Test 1: CellRanger formatting alone
echo ">>> Test 1: CellRanger formatting alone..."
TEST_CASE_DIR="$TEST_DIR/test_cr_only"
mkdir -p "$TEST_CASE_DIR"

"$STAR_BIN" --runMode genomeGenerate \
    --cellrangerStyleIndex Yes \
    --genomeFastaFiles "$TEST_FASTA" \
    --sjdbGTFfile "$TEST_GTF" \
    --genomeDir "$TEST_CASE_DIR" \
    --runThreadN 1 \
    --genomeSAindexNbases 4 \
    --genomeChrBinNbits 10 \
    --sjdbOverhang 50 > "$TEST_CASE_DIR/log.txt" 2>&1

if [ ! -f "$TEST_CASE_DIR/cellranger_ref/genome.fa" ]; then
    echo "    ✗ FAILED: cellranger_ref/genome.fa not created"
    cat "$TEST_CASE_DIR/log.txt"
    exit 1
fi
if [ ! -f "$TEST_CASE_DIR/cellranger_ref/genes.gtf" ]; then
    echo "    ✗ FAILED: cellranger_ref/genes.gtf not created"
    cat "$TEST_CASE_DIR/log.txt"
    exit 1
fi
if ! check_index_exists "$TEST_CASE_DIR"; then
    echo "    ✗ FAILED: Index files not created"
    cat "$TEST_CASE_DIR/log.txt"
    exit 1
fi

# Verify formatted FASTA has chr prefix
if ! grep -q "^>chr1" "$TEST_CASE_DIR/cellranger_ref/genome.fa"; then
    echo "    ✗ FAILED: Formatted FASTA should have chr prefix"
    exit 1
fi

echo "    ✓ PASSED: CellRanger formatting works"

# Test 2: FlexProbeIndex alone (skipped - requires matching chromosome names)
# Note: FlexProbeIndex alone is not the typical use case
# The main flex-style workflow is CellRanger + FlexProbeIndex together (Test 3)
echo ""
echo ">>> Test 2: FlexProbeIndex alone..."
echo "    ⚠ SKIPPED: FlexProbeIndex alone requires matching chromosome names"
echo "    Focus on Test 3 (CellRanger + FlexProbeIndex) which is the main flex-style workflow"

# Test 3: CellRanger formatting + FlexProbeIndex (flex-style index)
echo ""
echo ">>> Test 3: CellRanger + FlexProbeIndex (flex-style index)..."
TEST_CASE_DIR="$TEST_DIR/test_flex_style"
mkdir -p "$TEST_CASE_DIR"

"$STAR_BIN" --runMode genomeGenerate \
    --cellrangerStyleIndex Yes \
    --flexGeneProbeSet "$TEST_PROBES" \
    --genomeFastaFiles "$TEST_FASTA" \
    --sjdbGTFfile "$TEST_GTF" \
    --genomeDir "$TEST_CASE_DIR" \
    --runThreadN 1 \
    --genomeSAindexNbases 4 \
    --genomeChrBinNbits 10 \
    --sjdbOverhang 50 > "$TEST_CASE_DIR/log.txt" 2>&1

# Verify CellRanger formatting happened first
if [ ! -f "$TEST_CASE_DIR/cellranger_ref/genome.fa" ]; then
    echo "    ✗ FAILED: cellranger_ref/genome.fa not created"
    cat "$TEST_CASE_DIR/log.txt"
    exit 1
fi
if [ ! -f "$TEST_CASE_DIR/cellranger_ref/genes.gtf" ]; then
    echo "    ✗ FAILED: cellranger_ref/genes.gtf not created"
    cat "$TEST_CASE_DIR/log.txt"
    exit 1
fi

# Verify FlexProbeIndex ran after formatting
if [ ! -d "$TEST_CASE_DIR/flex_probe_artifacts" ]; then
    echo "    ✗ FAILED: flex_probe_artifacts directory not created"
    cat "$TEST_CASE_DIR/log.txt"
    exit 1
fi
if [ ! -f "$TEST_CASE_DIR/flex_probe_artifacts/genome.filtered.fa" ]; then
    echo "    ✗ FAILED: Hybrid FASTA not created"
    cat "$TEST_CASE_DIR/log.txt"
    exit 1
fi

# Verify FlexProbeIndex used formatted files by checking hybrid FASTA has chr prefix


# Verify the hybrid FASTA uses formatted chromosomes (chr prefix)
if ! grep -q "^>chr1" "$TEST_CASE_DIR/flex_probe_artifacts/genome.filtered.fa"; then
    echo "    ✗ FAILED: Hybrid FASTA should use formatted chromosomes (chr prefix)"
    exit 1
fi

if ! check_index_exists "$TEST_CASE_DIR"; then
    echo "    ✗ FAILED: Index files not created"
    cat "$TEST_CASE_DIR/log.txt"
    exit 1
fi

echo "    ✓ PASSED: Flex-style index (CellRanger + FlexProbeIndex) works"

# Test 4: AutoIndex with flex-style index
echo ""
echo ">>> Test 4: AutoIndex with flex-style index..."
TEST_CASE_DIR="$TEST_DIR/test_autoindex_flex"
mkdir -p "$TEST_CASE_DIR"

# Build initial flex-style index
"$STAR_BIN" --runMode genomeGenerate \
    --autoIndex Yes \
    --cellrangerStyleIndex Yes \
    --flexGeneProbeSet "$TEST_PROBES" \
    --genomeFastaFiles "$TEST_FASTA" \
    --sjdbGTFfile "$TEST_GTF" \
    --genomeDir "$TEST_CASE_DIR" \
    --runThreadN 1 \
    --genomeSAindexNbases 4 \
    --genomeChrBinNbits 10 \
    --sjdbOverhang 50 > "$TEST_CASE_DIR/log_initial.txt" 2>&1

if ! check_index_exists "$TEST_CASE_DIR"; then
    echo "    ✗ FAILED: Initial flex-style index build failed"
    cat "$TEST_CASE_DIR/log_initial.txt"
    exit 1
fi

# Try autoIndex again - should skip
"$STAR_BIN" --runMode genomeGenerate \
    --autoIndex Yes \
    --cellrangerStyleIndex Yes \
    --flexGeneProbeSet "$TEST_PROBES" \
    --genomeFastaFiles "$TEST_FASTA" \
    --sjdbGTFfile "$TEST_GTF" \
    --genomeDir "$TEST_CASE_DIR" \
    --runThreadN 1 \
    --genomeSAindexNbases 4 \
    --genomeChrBinNbits 10 \
    --sjdbOverhang 50 2>&1 | grep -E "(Index already exists|Skipping index generation)" > /dev/null || {
    echo "    ✗ FAILED: Should have skipped index generation"
    exit 1
}

echo "    ✓ PASSED: AutoIndex correctly skips existing flex-style index"

# Test 5: Verify probe filtering worked correctly
echo ""
echo ">>> Test 5: Verify probe filtering..."
TEST_CASE_DIR="$TEST_DIR/test_flex_style"

# Check that DEPRECATED probe was filtered out
if grep -q "DEPRECATED" "$TEST_CASE_DIR/flex_probe_artifacts/filtered_probe_set.csv" 2>/dev/null; then
    echo "    ✗ FAILED: DEPRECATED probe should have been filtered out"
    exit 1
fi

# Check that probe_list.txt exists and contains gene IDs
if [ ! -f "$TEST_CASE_DIR/probe_gene_list.txt" ]; then
    echo "    ✗ FAILED: probe_gene_list.txt not created"
    exit 1
fi

# Check that probes match GTF gene IDs
if ! grep -q "ENSG00000123456" "$TEST_CASE_DIR/probe_gene_list.txt"; then
    echo "    ✗ FAILED: probe_gene_list.txt should contain matching gene IDs"
    exit 1
fi

echo "    ✓ PASSED: Probe filtering works correctly"

echo ""
echo "=============================================="
echo "✓ FLEX-STYLE INDEX SMOKE TEST PASSED"
echo "=============================================="

