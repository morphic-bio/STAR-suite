#!/bin/bash
# Smoke test to verify FlexProbeIndex excludes probes with included=FALSE

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-$SCRIPT_DIR/../source/STAR}"

# Use test fixtures
TEST_FASTA="${TEST_FASTA:-$SCRIPT_DIR/fixtures/cellranger_format/test_genome.fa}"
TEST_GTF="${TEST_GTF:-$SCRIPT_DIR/fixtures/cellranger_format/test_genes.gtf}"
TEST_PROBE_CSV="${TEST_PROBE_CSV:-$SCRIPT_DIR/fixtures/cellranger_format/test_probes_included_false.csv}"

TEST_DIR=$(mktemp -d)
trap "rm -rf $TEST_DIR" EXIT

if [ ! -x "$STAR_BIN" ]; then
    echo "ERROR: STAR binary not found: $STAR_BIN"
    exit 1
fi

if [ ! -f "$TEST_FASTA" ] || [ ! -f "$TEST_GTF" ] || [ ! -f "$TEST_PROBE_CSV" ]; then
    echo "ERROR: Test fixtures not found"
    echo "  TEST_FASTA: $TEST_FASTA"
    echo "  TEST_GTF: $TEST_GTF"
    echo "  TEST_PROBE_CSV: $TEST_PROBE_CSV"
    exit 1
fi

echo "=============================================="
echo "Included=FALSE Filter Smoke Test"
echo "=============================================="
echo "STAR binary: $STAR_BIN"
echo "Test dir: $TEST_DIR"
echo ""

# Build Flex index with probe CSV containing included=FALSE probe
echo ">>> Building Flex index with included=FALSE probe..."
GENOME_DIR="$TEST_DIR/genome"
mkdir -p "$GENOME_DIR"

"$STAR_BIN" --runMode genomeGenerate \
     --genomeFastaFiles "$TEST_FASTA" \
     --sjdbGTFfile "$TEST_GTF" \
     --genomeDir "$GENOME_DIR" \
     --runThreadN 1 \
     --genomeSAindexNbases 4 \
     --genomeChrBinNbits 10 \
     --sjdbOverhang 50 \
     --cellrangerStyleIndex Yes \
     --flexGeneProbeSet "$TEST_PROBE_CSV" \
     2>&1 | grep -v "^working on" || true

# Check probe_list.txt does NOT include gene ID from excluded probe
PROBE_LIST="$GENOME_DIR/flex_probe_artifacts/probe_list.txt"
if [ ! -f "$PROBE_LIST" ]; then
    echo "✗ FAILED: probe_list.txt not found at $PROBE_LIST"
    exit 1
fi

echo ""
echo ">>> Test 1: probe_list.txt should NOT contain ENSG00000123457 (included=FALSE)..."
if grep -q "^ENSG00000123457$" "$PROBE_LIST"; then
    echo "✗ FAILED: probe_list.txt contains ENSG00000123457 (should be excluded)"
    echo "  probe_list.txt contents:"
    cat "$PROBE_LIST"
    exit 1
else
    echo "✓ PASSED: probe_list.txt does not contain ENSG00000123457"
    echo "  probe_list.txt contents:"
    cat "$PROBE_LIST"
fi

# Check probe_list.txt DOES include gene IDs from included=TRUE probes
echo ""
echo ">>> Test 2: probe_list.txt should contain ENSG00000123456 (included=TRUE)..."
if grep -q "^ENSG00000123456$" "$PROBE_LIST"; then
    echo "✓ PASSED: probe_list.txt contains ENSG00000123456"
else
    echo "✗ FAILED: probe_list.txt does not contain ENSG00000123456 (should be included)"
    exit 1
fi

# Check genome.filtered.fa does NOT include probe contig for excluded probe
HYBRID_FASTA="$GENOME_DIR/flex_probe_artifacts/genome.filtered.fa"
if [ ! -f "$HYBRID_FASTA" ]; then
    echo "✗ FAILED: genome.filtered.fa not found at $HYBRID_FASTA"
    exit 1
fi

echo ""
echo ">>> Test 3: genome.filtered.fa should NOT contain PROBE_002 (included=FALSE)..."
if grep -q "^>PROBE_002$" "$HYBRID_FASTA"; then
    echo "✗ FAILED: genome.filtered.fa contains PROBE_002 (should be excluded)"
    echo "  Checking for PROBE_002:"
    grep -A 2 "^>PROBE_002$" "$HYBRID_FASTA" || true
    exit 1
else
    echo "✓ PASSED: genome.filtered.fa does not contain PROBE_002"
fi

# Check genome.filtered.fa DOES include probe contigs for included=TRUE probes
echo ""
echo ">>> Test 4: genome.filtered.fa should contain PROBE_001 and PROBE_003 (included=TRUE)..."
PROBE_COUNT=$(grep -c "^>PROBE_" "$HYBRID_FASTA" || echo "0")
if [ "$PROBE_COUNT" -ge "2" ]; then
    echo "✓ PASSED: genome.filtered.fa contains at least 2 probe contigs"
    echo "  Probe contigs found:"
    grep "^>PROBE_" "$HYBRID_FASTA" | head -5
else
    echo "✗ FAILED: Expected at least 2 probe contigs, found $PROBE_COUNT"
    exit 1
fi

# Check probes_only.gtf does NOT include excluded probe
PROBES_GTF="$GENOME_DIR/flex_probe_artifacts/probes_only.gtf"
if [ -f "$PROBES_GTF" ]; then
    echo ""
    echo ">>> Test 5: probes_only.gtf should NOT contain PROBE_002 (included=FALSE)..."
    if grep -q "probe_id \"PROBE_002\"" "$PROBES_GTF"; then
        echo "✗ FAILED: probes_only.gtf contains PROBE_002 (should be excluded)"
        exit 1
    else
        echo "✓ PASSED: probes_only.gtf does not contain PROBE_002"
    fi
fi

# Check log for dropped count
BUILD_LOG="$GENOME_DIR/Log.out"
if [ -f "$BUILD_LOG" ]; then
    echo ""
    echo ">>> Test 6: Build log should report dropped_included_false count..."
    if grep -q "Dropped (included=FALSE)" "$BUILD_LOG"; then
        DROPPED_COUNT=$(grep "Dropped (included=FALSE)" "$BUILD_LOG" | grep -o "[0-9]*" | head -1 || echo "0")
        echo "✓ PASSED: Build log reports dropped_included_false: $DROPPED_COUNT"
        if [ "$DROPPED_COUNT" -ge "1" ]; then
            echo "  ✓ Expected at least 1 dropped probe (found $DROPPED_COUNT)"
        else
            echo "  ⚠ Warning: Expected at least 1 dropped probe, but found $DROPPED_COUNT"
        fi
    else
        echo "⚠ Warning: Build log does not contain 'Dropped (included=FALSE)' message"
    fi
fi

echo ""
echo "=============================================="
echo "✓ INCLUDED=FALSE FILTER SMOKE TEST PASSED"
echo "=============================================="

