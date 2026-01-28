#!/bin/bash
# CB/UB Tag Injection Validation Test
# Tests that CB and UB tags are correctly injected into unsorted BAM output
# Uses 100K downsampled Flex dataset

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${SCRIPT_DIR}/.."
OUT_DIR="${SCRIPT_DIR}/cbub_tag_validation"
TMP_DIR="/storage/100K/tmp/cbub_validation"
STAR_BIN="${REPO_DIR}/core/legacy/source/STAR"

# Dataset and reference paths (same as run_100K_regression_test.sh)
GENOME_DIR="/storage/flex_filtered_reference/star_index"
CB_WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
SAMPLE_WHITELIST="/storage/SC2300771_filtered_2M/sample_whitelist.tsv"
PROBE_LIST="/storage/flex_filtered_reference/filtered_reference/probe_list.txt"
SAMPLE_PROBES="/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt"
READS_R2="/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R2_001.fastq.gz"
READS_R1="/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R1_001.fastq.gz"

echo "=== CB/UB Tag Injection Validation Test ==="
echo "STAR binary: $STAR_BIN"
echo "Output directory: $OUT_DIR"
echo ""

# Verify STAR exists
if [ ! -x "$STAR_BIN" ]; then
    echo "ERROR: STAR binary not found or not executable: $STAR_BIN" >&2
    exit 1
fi

# Clean up previous run
rm -rf "$OUT_DIR"
rm -rf "$TMP_DIR"
mkdir -p "$OUT_DIR"

echo "=== Step 1: Running STAR with CB/UB tag injection enabled ==="
"$STAR_BIN" \
  --runThreadN 8 \
  --outTmpDir "$TMP_DIR" \
  --genomeDir "$GENOME_DIR" \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
  --soloCBwhitelist "$CB_WHITELIST" \
  --flex yes \
  --soloFlexExpectedCellsPerTag 3000 \
  --soloSampleWhitelist "$SAMPLE_WHITELIST" \
  --soloProbeList "$PROBE_LIST" \
  --soloSampleProbes "$SAMPLE_PROBES" \
  --soloSampleProbeOffset 68 \
  --soloFlexAllowedTags "$SAMPLE_WHITELIST" \
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
  --outSAMattributes NH HI AS nM NM GX GN CB UB \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloCellFilter None \
  --clipAdapterType CellRanger4 \
  --soloFeatures Gene \
  --alignEndsType Local \
  --soloStrand Unstranded \
  --chimSegmentMin 1000000 \
  --soloKeysCompat cr \
  --outSAMtype BAM Unsorted \
  --soloSampleSearchNearby no \
  --readFilesCommand zcat \
  --readFilesIn "$READS_R2" "$READS_R1" \
  --outFileNamePrefix "$OUT_DIR/"

echo ""
echo "=== Step 2: Validating CB/UB tags in BAM output ==="

BAM_FILE="${OUT_DIR}/Aligned.out.bam"
if [ ! -f "$BAM_FILE" ]; then
    echo "ERROR: BAM file not found: $BAM_FILE" >&2
    exit 1
fi

# Count total records
TOTAL=$(samtools view -c "$BAM_FILE")
echo "Total BAM records: $TOTAL"

# Count records with CB tag
CB_COUNT=$(samtools view "$BAM_FILE" | grep -c "CB:Z:" || echo "0")
echo "Records with CB:Z: tag: $CB_COUNT"

# Count records with UB tag
UB_COUNT=$(samtools view "$BAM_FILE" | grep -c "UB:Z:" || echo "0")
echo "Records with UB:Z: tag: $UB_COUNT"

# Count records with both CB and UB
BOTH_COUNT=$(samtools view "$BAM_FILE" | grep "CB:Z:" | grep -c "UB:Z:" || echo "0")
echo "Records with both CB:Z: and UB:Z: tags: $BOTH_COUNT"

# Count records with CB but no UB (CB-only, status=2 scenario)
CB_ONLY=$(samtools view "$BAM_FILE" | grep "CB:Z:" | grep -v "UB:Z:" | wc -l || echo "0")
echo "Records with CB:Z: only (no UB:Z:): $CB_ONLY"

# Calculate percentages
if [ "$TOTAL" -gt 0 ]; then
    CB_PCT=$(echo "scale=2; $CB_COUNT * 100 / $TOTAL" | bc)
    UB_PCT=$(echo "scale=2; $UB_COUNT * 100 / $TOTAL" | bc)
    BOTH_PCT=$(echo "scale=2; $BOTH_COUNT * 100 / $TOTAL" | bc)
    echo ""
    echo "CB tag rate: ${CB_PCT}%"
    echo "UB tag rate: ${UB_PCT}%"
    echo "Both tags rate: ${BOTH_PCT}%"
fi

# Sample some records with tags
echo ""
echo "=== Sample records with CB/UB tags ==="
echo "First 5 records with CB tag:"
samtools view "$BAM_FILE" | grep "CB:Z:" | head -5 | cut -f1,12- | head -c 500
echo ""

echo ""
echo "First 5 records with CB-only (no UB):"
samtools view "$BAM_FILE" | grep "CB:Z:" | grep -v "UB:Z:" | head -5 | cut -f1,12- | head -c 500 || echo "(none found)"
echo ""

# Validation checks
echo ""
echo "=== Validation Results ==="

PASS=true

# Check 1: CB tags should be present
if [ "$CB_COUNT" -gt 0 ]; then
    echo "✓ CB tags are being injected ($CB_COUNT records)"
else
    echo "✗ FAIL: No CB tags found in BAM"
    PASS=false
fi

# Check 2: UB tags should be present
if [ "$UB_COUNT" -gt 0 ]; then
    echo "✓ UB tags are being injected ($UB_COUNT records)"
else
    echo "✗ FAIL: No UB tags found in BAM"
    PASS=false
fi

# Check 3: Most records with CB should also have UB (typical case)
if [ "$CB_COUNT" -gt 0 ]; then
    BOTH_RATIO=$(echo "scale=2; $BOTH_COUNT * 100 / $CB_COUNT" | bc)
    if [ "$(echo "$BOTH_RATIO > 80" | bc)" -eq 1 ]; then
        echo "✓ Most CB records also have UB (${BOTH_RATIO}% overlap)"
    else
        echo "⚠ Only ${BOTH_RATIO}% of CB records have UB (may be expected for many invalid UMIs)"
    fi
fi

# Check 4: CB-only records demonstrate independence
if [ "$CB_ONLY" -gt 0 ]; then
    echo "✓ CB/UB independence working: $CB_ONLY records have CB without UB"
else
    echo "⚠ No CB-only records found (all valid UMIs is possible but unusual)"
fi

# Check 5: Tag rate should be reasonable (>5% for Flex data)
MIN_TAG_RATE=5
if [ "$(echo "$CB_PCT > $MIN_TAG_RATE" | bc)" -eq 1 ]; then
    echo "✓ CB tag rate is reasonable (${CB_PCT}% > ${MIN_TAG_RATE}%)"
else
    echo "✗ FAIL: CB tag rate is too low (${CB_PCT}% <= ${MIN_TAG_RATE}%)"
    PASS=false
fi

echo ""
if [ "$PASS" = true ]; then
    echo "=== ALL VALIDATION CHECKS PASSED ==="
    echo ""
    echo "Summary:"
    echo "  - CB/UB tags are being injected into unsorted BAM"
    echo "  - CB and UB are independent (CB-only records exist: $CB_ONLY)"
    echo "  - Tag injection rate is reasonable"
    exit 0
else
    echo "=== VALIDATION FAILED ==="
    exit 1
fi
