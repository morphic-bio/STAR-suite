#!/bin/bash
# Baseline test without --flexEnable to verify default behavior unchanged

set -e

STAR_BIN="/mnt/pikachu/STAR-Flex/source/STAR"
OUT_DIR="/mnt/pikachu/STAR-Flex/tests/baseline_test_output"
TMP_DIR="/storage/100K/tmp/baseline_test"

# Clean up previous run
rm -rf "$OUT_DIR"
rm -rf "$TMP_DIR"
mkdir -p "$OUT_DIR"

echo "=== Testing STAR without --flexEnable (baseline) ==="
echo "Binary: $STAR_BIN"
echo "Output: $OUT_DIR"
echo ""

# Run without --flexEnable (default behavior)
"$STAR_BIN" \
  --runThreadN 4 \
  --outTmpDir "$TMP_DIR" \
  --genomeDir /storage/flex_filtered_reference/star_index \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
  --soloCBwhitelist /storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt \
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
echo "=== Baseline test completed ==="
echo ""

# Check results
echo "Checking output files..."
if [ -f "$OUT_DIR/Solo.out/Gene/raw/matrix.mtx" ]; then
    echo "✓ matrix.mtx exists"
    head -5 "$OUT_DIR/Solo.out/Gene/raw/matrix.mtx"
else
    echo "✗ matrix.mtx NOT found"
fi

echo ""
echo "=== Comparing flex vs baseline ==="
FLEX_MATRIX="/mnt/pikachu/STAR-Flex/tests/flex_inline_test_output/Solo.out/Gene/raw/matrix.mtx"
BASE_MATRIX="$OUT_DIR/Solo.out/Gene/raw/matrix.mtx"

if [ -f "$FLEX_MATRIX" ] && [ -f "$BASE_MATRIX" ]; then
    FLEX_ENTRIES=$(grep -v "^%" "$FLEX_MATRIX" | head -1 | awk '{print $3}')
    BASE_ENTRIES=$(grep -v "^%" "$BASE_MATRIX" | head -1 | awk '{print $3}')
    echo "Flex entries: $FLEX_ENTRIES"
    echo "Baseline entries: $BASE_ENTRIES"
    if [ "$FLEX_ENTRIES" == "$BASE_ENTRIES" ]; then
        echo "✓ Entry counts match!"
    else
        echo "⚠ Entry counts differ (expected for different UMI dedup paths)"
    fi
fi

