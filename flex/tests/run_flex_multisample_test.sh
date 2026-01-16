#!/bin/bash
# Test script for 100K downsampled dataset using --flex omnibus flag
# Based on reference/runSTAR.sh
# Compares against gold standard bundled in tests/gold_standard/

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${SCRIPT_DIR}/../source/STAR"
OUT_DIR="${SCRIPT_DIR}/flex_multisample_output"
TMP_DIR="/storage/100K/tmp/flex_multisample_test"
GOLD_DIR="${SCRIPT_DIR}/gold_standard"
GOLD_RAW="${GOLD_DIR}/raw"
GOLD_PER_SAMPLE="${GOLD_DIR}/per_sample"

# Clean up previous run
rm -rf "$OUT_DIR"
rm -rf "$TMP_DIR"
mkdir -p "$OUT_DIR"

echo "=== Testing STAR-Flex with multi-sample output ==="
echo "Binary: $STAR_BIN"
echo "Output: $OUT_DIR"
echo "Gold standard (raw): $GOLD_RAW"
echo "Gold standard (per-sample): $GOLD_PER_SAMPLE"
echo ""

# Run with flex parameters matching reference/runSTAR.sh
"$STAR_BIN" \
  --runThreadN 8 \
  --outTmpDir "$TMP_DIR" \
  --genomeDir /storage/flex_filtered_reference/star_index \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
  --soloCBwhitelist /storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt \
  --flex yes \
  --soloFlexExpectedCellsPerTag 3000 \
  --soloSampleWhitelist /storage/SC2300771_filtered_2M/sample_whitelist.tsv \
  --soloProbeList /storage/flex_filtered_reference/filtered_reference/probe_list.txt \
  --soloSampleProbes /mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt \
  --soloSampleProbeOffset 68 \
  --soloFlexAllowedTags /storage/SC2300771_filtered_2M/sample_whitelist.tsv \
  --soloFlexOutputPrefix "${OUT_DIR}/per_sample" \
  --limitIObufferSize 50000000 50000000 \
  --outSJtype None \
  --outBAMcompression 6 \
  --soloMultiMappers Rescue \
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
  --outSAMattributes NH HI AS nM NM GX GN ZG \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloCellFilter None \
  --clipAdapterType CellRanger4 \
  --soloFeatures Gene \
  --alignEndsType Local \
  --soloAddTagsToUnsorted no \
  --soloStrand Unstranded \
  --chimSegmentMin 1000000 \
  --soloKeysCompat cr \
  --outSAMtype BAM Unsorted \
  --soloSampleSearchNearby no \
  --readFilesCommand zcat \
  --readFilesIn \
    /storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R2_001.fastq.gz \
    /storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R1_001.fastq.gz \
  --outFileNamePrefix "$OUT_DIR/"

echo ""
echo "=== Test completed ==="
echo ""

# Check for output
echo "Checking output files..."

# Check standard Solo output
if [ -f "$OUT_DIR/Solo.out/Gene/raw/matrix.mtx" ]; then
    echo "✓ Solo.out/Gene/raw/matrix.mtx exists"
    ENTRIES=$(tail -n +3 "$OUT_DIR/Solo.out/Gene/raw/matrix.mtx" | head -1 | awk '{print $3}')
    echo "  Entries: $ENTRIES"
else
    echo "✗ Solo.out/Gene/raw/matrix.mtx NOT found"
fi

# Check per-sample output
echo ""
echo "Checking per-sample output..."
if [ -d "$OUT_DIR/per_sample" ]; then
    echo "✓ per_sample directory exists"
    for SAMPLE_LABEL in BC004 BC006 BC007 BC008; do
        OUR_MATRIX="$OUT_DIR/per_sample/${SAMPLE_LABEL}/Gene/filtered/matrix.mtx"
        GOLD_MATRIX="$GOLD_PER_SAMPLE/${SAMPLE_LABEL}/Gene/filtered/matrix.mtx"
        if [ -f "$OUR_MATRIX" ]; then
            OUR_ENTRIES=$(tail -n +3 "$OUR_MATRIX" | head -1 | awk '{print $3}')
            if [ -f "$GOLD_MATRIX" ]; then
                GOLD_ENTRIES=$(tail -n +3 "$GOLD_MATRIX" | head -1 | awk '{print $3}')
                echo "  ✓ ${SAMPLE_LABEL}: matrix.mtx exists (${OUR_ENTRIES} entries)"
                if [ "$OUR_ENTRIES" -ne "$GOLD_ENTRIES" ]; then
                    echo "    ⚠ Entry count differs (flex: ${OUR_ENTRIES}, gold: ${GOLD_ENTRIES})"
                else
                    echo "    ✓ Entry count matches gold standard"
                fi
            else
                echo "  ✓ ${SAMPLE_LABEL}: matrix.mtx exists (${OUR_ENTRIES} entries)"
                echo "    (gold standard not found for comparison)"
            fi
        else
            echo "  ✗ ${SAMPLE_LABEL}: matrix.mtx NOT found"
        fi
    done
else
    echo "✗ per_sample directory NOT found"
fi

echo ""
echo "=== Comparing raw MEX with gold standard ==="
GOLD_RAW_MATRIX="$GOLD_RAW/matrix.mtx"
OUR_RAW_MATRIX="$OUT_DIR/Solo.out/Gene/raw/matrix.mtx"
if [ -f "$GOLD_RAW_MATRIX" ] && [ -f "$OUR_RAW_MATRIX" ]; then
    GOLD_ENTRIES=$(tail -n +3 "$GOLD_RAW_MATRIX" | head -1 | awk '{print $3}')
    OUR_ENTRIES=$(tail -n +3 "$OUR_RAW_MATRIX" | head -1 | awk '{print $3}')
    echo "Gold raw entries: $GOLD_ENTRIES"
    echo "Our raw entries: $OUR_ENTRIES"
    if [ "$OUR_ENTRIES" -eq "$GOLD_ENTRIES" ]; then
        echo "✓ Raw matrix entry counts match!"
    else
        echo "⚠ Raw matrix entry counts differ"
    fi
    
    # Full diff comparison
    echo ""
    echo "=== Full MEX comparison ==="
    if diff -q "$GOLD_RAW/barcodes.tsv" "$OUT_DIR/Solo.out/Gene/raw/barcodes.tsv" > /dev/null 2>&1; then
        echo "✓ barcodes.tsv matches gold standard"
    else
        echo "⚠ barcodes.tsv differs from gold standard"
    fi
    if diff -q "$GOLD_RAW/features.tsv" "$OUT_DIR/Solo.out/Gene/raw/features.tsv" > /dev/null 2>&1; then
        echo "✓ features.tsv matches gold standard"
    else
        echo "⚠ features.tsv differs from gold standard"
    fi
    if diff -q "$GOLD_RAW_MATRIX" "$OUR_RAW_MATRIX" > /dev/null 2>&1; then
        echo "✓ matrix.mtx matches gold standard"
    else
        echo "⚠ matrix.mtx differs from gold standard"
    fi
else
    echo "⚠ Could not compare: gold=$GOLD_RAW_MATRIX exists=$([ -f \"$GOLD_RAW_MATRIX\" ] && echo yes || echo no)"
fi

echo ""
echo "Log snippets:"
grep -E "flex|Flex|probe|hash|Hash" "$OUT_DIR/Log.out" 2>/dev/null | tail -20 || echo "(no flex mentions in log)"
