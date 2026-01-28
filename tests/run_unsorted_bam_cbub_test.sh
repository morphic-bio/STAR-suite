#!/bin/bash
# Test unsorted BAM with CB/UB tag injection via buffered mode
# This tests the new SamtoolsSorter noSort mode for unsorted BAM tag injection

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${SCRIPT_DIR}/.."
OUT_DIR="${SCRIPT_DIR}/unsorted_cbub_test_output"
TMP_DIR="/storage/100K/tmp/unsorted_cbub_test"
STAR_BIN="${REPO_DIR}/core/legacy/source/STAR"

# Dataset and reference paths (same as sorted test)
GENOME_DIR="/storage/flex_filtered_reference/star_index"
CB_WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
SAMPLE_WHITELIST="/storage/SC2300771_filtered_2M/sample_whitelist.tsv"
PROBE_LIST="/storage/flex_filtered_reference/filtered_reference/probe_list.txt"
SAMPLE_PROBES="/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt"
READS_R2="/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R2_001.fastq.gz"
READS_R1="/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R1_001.fastq.gz"

echo "=========================================="
echo "Unsorted BAM CB/UB Tag Injection Test"
echo "=========================================="
echo ""

# Clean up previous run
rm -rf "$OUT_DIR"
rm -rf "$TMP_DIR"
mkdir -p "$OUT_DIR"

echo "Running STAR with unsorted BAM + CB/UB tags..."
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
echo "STAR completed. Checking for CB/UB tags in unsorted BAM..."
echo ""

BAM_FILE="${OUT_DIR}/Aligned.out.bam"

if [ ! -f "$BAM_FILE" ]; then
    echo "ERROR: Unsorted BAM file not found: $BAM_FILE"
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
samtools view "$BAM_FILE" | grep "CB:Z:" | head -5 | cut -f1,12-20 || true
echo ""

# Determine pass/fail
if [ "$CB_READS" -gt 0 ]; then
    echo "✓ SUCCESS: CB/UB tags are present in unsorted BAM!"
    echo ""
    echo "Tag injection rate:"
    echo "  CB: $(echo "scale=2; $CB_READS * 100 / $TOTAL_READS" | bc)% of reads"
    echo "  UB: $(echo "scale=2; $UB_READS * 100 / $TOTAL_READS" | bc)% of reads"
    exit 0
else
    echo "✗ FAILURE: No CB/UB tags found in unsorted BAM"
    exit 1
fi
