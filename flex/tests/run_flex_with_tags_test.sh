#!/bin/bash
# Test flex pipeline with sample tag detection
set -e

STAR_BIN="/mnt/pikachu/STAR-Flex/source/STAR"
OUT_DIR="/mnt/pikachu/STAR-Flex/tests/flex_tags_test_output"
GENOME_DIR="/storage/flex_filtered_reference/star_index"
WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
READS_R2="/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz"
READS_R1="/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz"
TMP_DIR="/storage/100K/tmp/flex_tags_test"

# Sample tag whitelist (derived from the data)
SAMPLE_WHITELIST="/storage/scRNAseq_output/sample_whitelists/sample_whitelist.tsv"
SAMPLE_PROBES="/storage/scRNAseq_output/sample_whitelists/sample_probes.tsv"
ALLOWED_TAGS="/storage/scRNAseq_output/sample_whitelists/allowed_tags.tsv"

echo "=== Testing STAR flex pipeline with sample tag detection ==="
echo "Binary: $STAR_BIN"
echo "Output: $OUT_DIR"

# Check if sample files exist
if [ ! -f "$SAMPLE_WHITELIST" ]; then
    echo "WARNING: Sample whitelist not found at $SAMPLE_WHITELIST"
    echo "Running without sample tag detection (basic flex mode)"
    SAMPLE_FLAGS=""
else
    echo "Sample whitelist: $SAMPLE_WHITELIST"
    SAMPLE_FLAGS="--soloSampleWhitelist $SAMPLE_WHITELIST"
    if [ -f "$SAMPLE_PROBES" ]; then
        SAMPLE_FLAGS="$SAMPLE_FLAGS --soloSampleProbes $SAMPLE_PROBES"
    fi
    if [ -f "$ALLOWED_TAGS" ]; then
        SAMPLE_FLAGS="$SAMPLE_FLAGS --flexAllowedTags $ALLOWED_TAGS"
    fi
fi

rm -rf "$OUT_DIR" "$TMP_DIR"
mkdir -p "$OUT_DIR" "$TMP_DIR"

LOG_FILE="$OUT_DIR/Log.out"

"$STAR_BIN" \
  --runThreadN 4 \
  --outTmpDir "$TMP_DIR" \
  --genomeDir "$GENOME_DIR" \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
  --soloCBwhitelist "$WHITELIST" \
  --flexEnable yes \
  $SAMPLE_FLAGS \
  --flexExpectedCellsTotal 5000 \
  --flexEdNiters 5000 \
  --flexEdFdr 0.01 \
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
  --readFilesIn "$READS_R2" "$READS_R1" \
  --outFileNamePrefix "$OUT_DIR/" \
  > "$LOG_FILE" 2>&1

echo "=== Test completed ==="
echo ""
echo "Checking output files..."

# Check standard Solo output
if [ -f "$OUT_DIR/Solo.out/Gene/raw/matrix.mtx" ]; then
    echo "✓ Standard matrix.mtx exists"
    head -5 "$OUT_DIR/Solo.out/Gene/raw/matrix.mtx"
else
    echo "✗ Standard matrix.mtx missing"
fi

# Check FlexFilter output (if sample tags were used)
if [ -d "$OUT_DIR/Solo.out/FlexFilter" ]; then
    echo ""
    echo "✓ FlexFilter output directory exists"
    ls -la "$OUT_DIR/Solo.out/FlexFilter/"
    
    if [ -f "$OUT_DIR/Solo.out/FlexFilter/flexfilter_summary.tsv" ]; then
        echo ""
        echo "FlexFilter summary:"
        cat "$OUT_DIR/Solo.out/FlexFilter/flexfilter_summary.tsv"
    fi
else
    echo ""
    echo "Note: FlexFilter output not generated (sample tags may not be configured)"
fi

echo ""
echo "Log snippets for flex/sample detection:"
grep -E "Flex pipeline enabled|SampleDetector|inlineHashMode|detected sample|Building flex matrix|FlexFilter" "$LOG_FILE" 2>/dev/null | head -20 || true

echo ""
echo "Test script completed."

