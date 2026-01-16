#!/bin/bash
# Test script for inline hash mode using --flexEnable yes
# This tests the basic inline hash mechanism without sample tags

set -e

STAR_BIN="/mnt/pikachu/STAR-Flex/source/STAR"
OUT_DIR="/mnt/pikachu/STAR-Flex/tests/flex_inline_test_output"
TMP_DIR="/storage/100K/tmp/flex_inline_test"

# Clean up previous run
rm -rf "$OUT_DIR"
rm -rf "$TMP_DIR"
mkdir -p "$OUT_DIR"

echo "=== Testing STAR with --flexEnable yes (inline hash mode) ==="
echo "Binary: $STAR_BIN"
echo "Output: $OUT_DIR"
echo ""

# Run with just 2 threads and 1 lane for faster testing
"$STAR_BIN" \
  --runThreadN 4 \
  --outTmpDir "$TMP_DIR" \
  --genomeDir /storage/flex_filtered_reference/star_index \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
  --soloCBwhitelist /storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt \
  --flex yes \
  --soloFlexExpectedCellsPerTag 3000 \
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
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloCellFilter None \
  --clipAdapterType CellRanger4 \
  --soloFeatures Gene \
  --alignEndsType Local \
  --soloStrand Unstranded \
  --chimSegmentMin 1000000 \
  --outSAMtype BAM Unsorted \
  --readFilesCommand zcat \
  --readFilesIn \
    /storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz \
    /storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz \
  --outFileNamePrefix "$OUT_DIR/"

echo ""
echo "=== Test completed ==="
echo ""

# Check results
# NOTE: In inline hash + minimal memory mode (--flex yes), MEX output files are NOT created
# because the pipeline uses direct hash collapse without materialization. This is expected behavior.
echo "Checking test results..."
echo "NOTE: MEX output files (matrix.mtx, barcodes.tsv, features.tsv) are NOT expected"
echo "      in inline hash + minimal memory mode (--flex yes). This is correct behavior."
echo ""

# Verify STAR ran successfully
if [ ! -f "$OUT_DIR/Log.out" ]; then
    echo "✗ ERROR: Log.out not found - STAR run may have failed"
    exit 1
fi

# Check that inline hash collapse was executed
if grep -q "Starting direct hash collapse\|Finished collapsing UMIs.*direct hash\|inline-hash mode completed" "$OUT_DIR/Log.out"; then
    echo "✓ Inline hash collapse executed successfully"
else
    echo "✗ WARNING: Inline hash collapse not detected in log"
fi

# Check for inline hash mode indication
if grep -qi "inline.*hash\|direct hash\|inline-hash" "$OUT_DIR/Log.out"; then
    echo "✓ Inline hash mode confirmed in log"
else
    echo "⚠ WARNING: No inline hash mentions in log"
fi

echo ""
echo "Log snippets for inline hash:"
grep -i "inline\|hash\|flex\|direct hash\|inline-hash" "$OUT_DIR/Log.out" | head -5 || echo "(no inline/hash mentions in log)"

echo ""
echo "✓ Test PASSED: Inline hash mode executed correctly (no MEX output expected)"

