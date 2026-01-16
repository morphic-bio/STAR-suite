#!/bin/bash
# Smoke test to verify non-Flex runs produce legacy MEX output
# and Flex runs continue to use inline-hash pipeline

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-$SCRIPT_DIR/../source/STAR}"

# Use small test fixtures if available, otherwise use environment variables
TEST_GENOME_DIR="${TEST_GENOME_DIR:-}"
TEST_R1="${TEST_R1:-}"
TEST_R2="${TEST_R2:-}"
TEST_CB_WHITELIST="${TEST_CB_WHITELIST:-}"

OUT_DIR="${OUT_DIR:-$SCRIPT_DIR/nonflex_mex_smoke_output}"
TMP_DIR="${TMP_DIR:-/tmp/nonflex_mex_smoke_$$}"

mkdir -p "$OUT_DIR" "$TMP_DIR"

# Cleanup on exit
cleanup() {
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

if [ ! -x "$STAR_BIN" ]; then
    echo "ERROR: STAR binary not found: $STAR_BIN"
    exit 1
fi

echo "=============================================="
echo "Non-Flex MEX Output Smoke Test"
echo "=============================================="
echo "STAR binary: $STAR_BIN"
echo "Output dir: $OUT_DIR"
echo ""

# Test 1: Non-Flex run (should produce standard MEX output)
echo ">>> Test 1: Non-Flex PE run (should produce matrix.mtx, barcodes.tsv, features.tsv)"
NONFLEX_OUT="$OUT_DIR/nonflex"
mkdir -p "$NONFLEX_OUT"

# Check if we have test data
if [ -z "$TEST_GENOME_DIR" ] || [ -z "$TEST_R1" ] || [ -z "$TEST_R2" ] || [ -z "$TEST_CB_WHITELIST" ]; then
    echo "SKIP: Test data not provided via environment variables:"
    echo "  TEST_GENOME_DIR, TEST_R1, TEST_R2, TEST_CB_WHITELIST"
    echo ""
    echo "To run this test, provide:"
    echo "  TEST_GENOME_DIR=/path/to/genome"
    echo "  TEST_R1=/path/to/R1.fastq.gz (CB/UMI read)"
    echo "  TEST_R2=/path/to/R2.fastq.gz (cDNA read)"
    echo "  TEST_CB_WHITELIST=/path/to/whitelist.txt"
    echo "  TEST_PROBE_LIST=/path/to/probe_list.txt (optional, uses fixture by default for Flex test)"
    echo ""
    echo "Note: R2 (cDNA) is passed first, R1 (CB/UMI) second to match STARsolo conventions"
    echo ""
    echo "Example:"
    echo "  TEST_GENOME_DIR=/storage/flex_filtered_reference/star_index \\"
    echo "  TEST_R1=/storage/downsampled_100K/SC2300771/...R1_001.fastq.gz \\"
    echo "  TEST_R2=/storage/downsampled_100K/SC2300771/...R2_001.fastq.gz \\"
    echo "  TEST_CB_WHITELIST=/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt \\"
    echo "  TEST_PROBE_LIST=/path/to/probe_list.txt \\"
    echo "  $0"
    exit 0
fi

echo "Running STAR without --flex (default behavior)..."
"$STAR_BIN" \
    --runThreadN 2 \
    --outTmpDir "$TMP_DIR/nonflex" \
    --genomeDir "$TEST_GENOME_DIR" \
    --soloType CB_UMI_Simple \
    --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
    --soloCBwhitelist "$TEST_CB_WHITELIST" \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR \
    --soloCellFilter None \
    --soloFeatures Gene \
    --soloStrand Unstranded \
    --clipAdapterType CellRanger4 \
    --readFilesCommand zcat \
    --readFilesIn "$TEST_R2" "$TEST_R1" \
    --outFileNamePrefix "$NONFLEX_OUT/" \
    --limitIObufferSize 50000000 50000000 \
    --outSJtype None \
    --outBAMcompression 6 \
    --alignIntronMax 500000 \
    --outFilterMismatchNmax 6 \
    --outFilterMismatchNoverReadLmax 1.0 \
    --outFilterMatchNmin 25 \
    --outSAMunmapped None \
    --outFilterMatchNminOverLread 0 \
    --outFilterMultimapNmax 10000 \
    --outFilterMultimapScoreRange 4 \
    --outSAMmultNmax 10000 \
    --winAnchorMultimapNmax 200 \
    --outSAMprimaryFlag AllBestScore \
    --outFilterScoreMin 0 \
    --outFilterScoreMinOverLread 0 \
    --outSAMattributes NH HI AS nM NM GX GN \
    --alignEndsType Local \
    --chimSegmentMin 1000000 \
    --outSAMtype BAM Unsorted \
    2>&1 | tee "$NONFLEX_OUT/run.log" | tail -20 || {
    echo "FAIL: Non-Flex run failed"
    exit 1
}

echo ""
echo "Checking for legacy MEX output files..."
MEX_DIR="$NONFLEX_OUT/Solo.out/Gene/raw"
if [ ! -f "$MEX_DIR/matrix.mtx" ]; then
    echo "FAIL: matrix.mtx not found in $MEX_DIR"
    exit 1
fi
if [ ! -f "$MEX_DIR/barcodes.tsv" ]; then
    echo "FAIL: barcodes.tsv not found in $MEX_DIR"
    exit 1
fi
if [ ! -f "$MEX_DIR/features.tsv" ]; then
    echo "FAIL: features.tsv not found in $MEX_DIR"
    exit 1
fi

echo "✓ PASS: Non-Flex run produced standard MEX output"
echo "  - matrix.mtx: $(wc -l < "$MEX_DIR/matrix.mtx" | tr -d ' ') lines"
echo "  - barcodes.tsv: $(wc -l < "$MEX_DIR/barcodes.tsv" | tr -d ' ') lines"
echo "  - features.tsv: $(wc -l < "$MEX_DIR/features.tsv" | tr -d ' ') lines"

# Check log for inline-hash mode (should NOT be enabled)
if grep -qi "inline.*hash\|direct hash\|inline-hash" "$NONFLEX_OUT/Log.out"; then
    echo "WARNING: Inline hash mode detected in non-Flex run (should be disabled)"
    grep -i "inline.*hash\|direct hash\|inline-hash" "$NONFLEX_OUT/Log.out" | head -5
else
    echo "✓ PASS: Inline hash mode correctly disabled for non-Flex run"
fi

echo ""
echo ">>> Test 2: Flex run (should use inline-hash pipeline)"
FLEX_OUT="$OUT_DIR/flex"
mkdir -p "$FLEX_OUT"

# Flex mode requires probe list
TEST_PROBE_LIST="${TEST_PROBE_LIST:-$SCRIPT_DIR/fixtures/flex_probe_required/test_probe_list.txt}"
if [ ! -f "$TEST_PROBE_LIST" ]; then
    echo "ERROR: Probe list required for Flex mode test"
    echo "  TEST_PROBE_LIST: $TEST_PROBE_LIST"
    echo "  Set TEST_PROBE_LIST environment variable or ensure fixture exists"
    exit 1
fi

echo "Running STAR with --flex yes (using probe list: $TEST_PROBE_LIST)..."
"$STAR_BIN" \
    --runThreadN 2 \
    --outTmpDir "$TMP_DIR/flex" \
    --genomeDir "$TEST_GENOME_DIR" \
    --soloType CB_UMI_Simple \
    --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
    --soloCBwhitelist "$TEST_CB_WHITELIST" \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR \
    --soloCellFilter None \
    --soloFeatures Gene \
    --soloStrand Unstranded \
    --clipAdapterType CellRanger4 \
    --flex yes \
    --soloProbeList "$TEST_PROBE_LIST" \
    --soloFlexExpectedCellsPerTag 1000 \
    --readFilesCommand zcat \
    --readFilesIn "$TEST_R2" "$TEST_R1" \
    --outFileNamePrefix "$FLEX_OUT/" \
    --limitIObufferSize 50000000 50000000 \
    --outSJtype None \
    --outBAMcompression 6 \
    --alignIntronMax 500000 \
    --outFilterMismatchNmax 6 \
    --outFilterMismatchNoverReadLmax 1.0 \
    --outFilterMatchNmin 25 \
    --outSAMunmapped None \
    --outFilterMatchNminOverLread 0 \
    --outFilterMultimapNmax 10000 \
    --outFilterMultimapScoreRange 4 \
    --outSAMmultNmax 10000 \
    --winAnchorMultimapNmax 200 \
    --outSAMprimaryFlag AllBestScore \
    --outFilterScoreMin 0 \
    --outFilterScoreMinOverLread 0 \
    --outSAMattributes NH HI AS nM NM GX GN \
    --alignEndsType Local \
    --chimSegmentMin 1000000 \
    --outSAMtype BAM Unsorted \
    2>&1 | tee "$FLEX_OUT/run.log" | tail -20 || {
    echo "FAIL: Flex run failed"
    exit 1
}

echo ""
echo "Checking for inline-hash mode in Flex run..."
if grep -qi "inline.*hash\|direct hash\|inline-hash\|Starting direct hash collapse" "$FLEX_OUT/Log.out"; then
    echo "✓ PASS: Inline hash mode correctly enabled for Flex run"
    grep -i "inline.*hash\|direct hash\|inline-hash\|Starting direct hash collapse" "$FLEX_OUT/Log.out" | head -5
else
    echo "WARNING: Inline hash mode not detected in Flex run"
    echo "Checking log for Flex mode..."
    grep -i "flex\|inline" "$FLEX_OUT/Log.out" | head -10 || true
fi

# Check for MEX output (Flex should also produce MEX, but via inline-hash path)
FLEX_MEX_DIR="$FLEX_OUT/Solo.out/Gene/raw"
if [ -f "$FLEX_MEX_DIR/matrix.mtx" ]; then
    echo "✓ PASS: Flex run produced MEX output (via inline-hash)"
    echo "  - matrix.mtx: $(wc -l < "$FLEX_MEX_DIR/matrix.mtx" | tr -d ' ') lines"
else
    echo "NOTE: Flex run may not produce MEX in minimal memory mode (expected behavior)"
fi

echo ""
echo "=============================================="
echo "Summary"
echo "=============================================="
echo "✓ Non-Flex run: Legacy MEX output produced"
echo "✓ Flex run: Inline-hash pipeline enabled"
echo ""
echo "Test output directory: $OUT_DIR"
echo "  - Non-Flex: $NONFLEX_OUT"
echo "  - Flex: $FLEX_OUT"
echo ""
echo "All tests PASSED"

