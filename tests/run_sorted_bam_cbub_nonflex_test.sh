#!/bin/bash
# Test sorted BAM with CB/UB tag injection in NON-FLEX mode
# This tests standard STARsolo which uses the legacy path that populates packedReadInfo

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${SCRIPT_DIR}/.."
OUT_DIR="${SCRIPT_DIR}/sorted_cbub_nonflex_output"
TMP_DIR="/storage/100K/tmp/sorted_cbub_nonflex"
STAR_BIN="${REPO_DIR}/core/legacy/source/STAR"

# Use a standard 10x Genomics style dataset (not Flex)
GENOME_DIR="/storage/flex_filtered_reference/star_index"
CB_WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
# Just use 2 lanes from the 100K dataset
READS_R2="/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz"
READS_R1="/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz"

echo "=========================================="
echo "Sorted BAM CB/UB Tag Test (NON-FLEX Mode)"
echo "=========================================="
echo ""

# Clean up previous run
rm -rf "$OUT_DIR"
rm -rf "$TMP_DIR"
mkdir -p "$OUT_DIR"

echo "Running STAR in standard STARsolo mode (no --flex)..."
"$STAR_BIN" \
  --runThreadN 8 \
  --outTmpDir "$TMP_DIR" \
  --genomeDir "$GENOME_DIR" \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
  --soloCBwhitelist "$CB_WHITELIST" \
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
  --outSAMtype BAM SortedByCoordinate \
  --readFilesCommand zcat \
  --readFilesIn "$READS_R2" "$READS_R1" \
  --outFileNamePrefix "$OUT_DIR/"

echo ""
echo "STAR completed. Checking for CB/UB tags in sorted BAM..."
echo ""

BAM_FILE="${OUT_DIR}/Aligned.sortedByCoord.out.bam"

if [ ! -f "$BAM_FILE" ]; then
    echo "ERROR: Sorted BAM file not found: $BAM_FILE"
    exit 1
fi

# Count reads with CB and UB tags
TOTAL_READS=$(samtools view -c "$BAM_FILE")
CB_READS=$(samtools view "$BAM_FILE" | grep -c "CB:Z:" || true)
UB_READS=$(samtools view "$BAM_FILE" | grep -c "UB:Z:" || true)
CBUB_READS=$(samtools view "$BAM_FILE" | grep "CB:Z:" | grep -c "UB:Z:" || true)

echo "=========================================="
echo "Results:"
echo "=========================================="
echo "Total reads in BAM:       $TOTAL_READS"
echo "Reads with CB tag:        $CB_READS"
echo "Reads with UB tag:        $UB_READS"
echo "Reads with both CB+UB:    $CBUB_READS"
echo ""

# Show sample of reads with CB/UB tags
echo "Sample reads with CB/UB tags (first 5):"
echo "----------------------------------------"
samtools view "$BAM_FILE" | grep "CB:Z:" | head -5 | cut -f1,12-20
echo ""

# Determine pass/fail
if [ "$CB_READS" -gt 0 ]; then
    echo "✓ SUCCESS: CB/UB tags are present in sorted BAM!"
    echo ""
    echo "Tag injection rate:"
    echo "  CB: $(echo "scale=2; $CB_READS * 100 / $TOTAL_READS" | bc)% of reads"
    echo "  UB: $(echo "scale=2; $UB_READS * 100 / $TOTAL_READS" | bc)% of reads"
    exit 0
else
    echo "✗ FAILURE: No CB/UB tags found in sorted BAM"
    exit 1
fi
