#!/bin/bash
# Integration test for CB/UB tag validation in sorted BAM output
# Tests readId-based validation infrastructure with Flex dataset

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${SCRIPT_DIR}/../core/legacy/source/STAR"
OUT_DIR="${SCRIPT_DIR}/flex_cbub_validation_output"
TMP_DIR="/tmp/flex_cbub_validation_$$"

# Clean up previous run
rm -rf "$OUT_DIR"
rm -rf "$TMP_DIR"
mkdir -p "$OUT_DIR"
# Don't create TMP_DIR - STAR will create it

# Check prerequisites
if [ ! -f "$STAR_BIN" ]; then
    echo "ERROR: STAR binary not found: $STAR_BIN" >&2
    exit 1
fi

if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found in PATH" >&2
    exit 1
fi

if [ ! -d "/storage/flex_filtered_reference/star_index" ]; then
    echo "ERROR: STAR index not found: /storage/flex_filtered_reference/star_index" >&2
    exit 1
fi

echo "=== Testing CB/UB tag validation with sorted BAM ==="
echo "Binary: $STAR_BIN"
echo "Output: $OUT_DIR"
echo ""

# Run STAR with sorted BAM output via samtools sorter
# Enable ZI tag and CB/UB table export via environment variables
# 
# IMPORTANT: We use --flex yes to enable Flex mode, but override two flags that would
#            empty packedReadInfo (required for CB/UB table export):
#            - --soloFlexMinimalMemory no (keeps packedReadInfo populated)
#            - --soloInlineCBCorrection no (keeps packedReadInfo populated)
# 
# This allows us to test the Flex path (sample detection, probe matching, etc.)
# while still having packedReadInfo available for validation.
# 
# NOTE: In production Flex runs with --flex yes defaults, packedReadInfo will be
#       empty and CB/UB table export will hard-fail. This is intentional - table
#       export is a validation/testing feature, not for production Flex runs.
STAR_EMIT_READID_TAG=1 STAR_EMIT_CBUB_TABLE=1 "$STAR_BIN" \
  --runThreadN 4 \
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
  --soloFlexMinimalMemory no \
  --soloInlineCBCorrection no \
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
  --soloAddTagsToUnsorted no \
  --soloStrand Unstranded \
  --chimSegmentMin 1000000 \
  --soloKeysCompat cr \
  --outSAMtype BAM SortedByCoordinate \
  --outBAMsortMethod samtools \
  --soloSampleSearchNearby no \
  --readFilesCommand zcat \
  --readFilesIn \
    /storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz \
    /storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz \
  --outFileNamePrefix "$OUT_DIR/"

# Verify output files exist
BAM_FILE="${OUT_DIR}/Aligned.sortedByCoord.out.bam"
TSV_FILE="${OUT_DIR}/Aligned.out.cb_ub.readId.tsv"

if [ ! -f "$BAM_FILE" ]; then
    echo "ERROR: Sorted BAM file not found: $BAM_FILE" >&2
    exit 1
fi

if [ ! -f "$TSV_FILE" ]; then
    echo "ERROR: CB/UB TSV table not found: $TSV_FILE" >&2
    exit 1
fi

echo ""
echo "=== Running CB/UB tag comparison ==="
echo "BAM: $BAM_FILE"
echo "TSV: $TSV_FILE"
echo ""

# Run comparison script
"${SCRIPT_DIR}/compare_cb_ub_tags.sh" "$BAM_FILE" "$TSV_FILE"

echo ""
echo "=== Test PASSED ==="

# Cleanup
rm -rf "$TMP_DIR"

