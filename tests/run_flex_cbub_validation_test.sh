#!/bin/bash
# Integration test for CB/UB tag validation in sorted BAM output
# Tests readId-based validation infrastructure with Flex dataset

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-${SCRIPT_DIR}/../core/legacy/source/STAR}"
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}/flex_cbub_validation_output}"
TMP_DIR="${TMP_DIR:-/tmp/flex_cbub_validation_$$}"
FLEX_INDEX="${FLEX_INDEX:-/storage/flex_filtered_reference/star_index}"
FLEX_WHITELIST="${FLEX_WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}"
FLEX_SAMPLE_WHITELIST="${FLEX_SAMPLE_WHITELIST:-/storage/SC2300771_filtered_2M/sample_whitelist.tsv}"
FLEX_PROBE_LIST="${FLEX_PROBE_LIST:-/storage/flex_filtered_reference/filtered_reference/probe_list.txt}"
FLEX_SAMPLE_PROBES="${FLEX_SAMPLE_PROBES:-/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}"
FLEX_ALLOWED_TAGS="${FLEX_ALLOWED_TAGS:-/storage/SC2300771_filtered_2M/sample_whitelist.tsv}"
FLEX_FASTQ_R2="${FLEX_FASTQ_R2:-/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz}"
FLEX_FASTQ_R1="${FLEX_FASTQ_R1:-/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz}"

skip() {
    echo "SKIP: $*"
    exit 0
}

# Clean up previous run
rm -rf "$OUT_DIR"
rm -rf "$TMP_DIR"
mkdir -p "$OUT_DIR"
# Don't create TMP_DIR - STAR will create it

# Check prerequisites
if [ ! -f "$STAR_BIN" ]; then
    skip "STAR binary not found: $STAR_BIN"
fi

if ! command -v samtools &> /dev/null; then
    skip "samtools not found in PATH"
fi

if [ ! -d "$FLEX_INDEX" ]; then
    skip "STAR index not found: $FLEX_INDEX"
fi
if [ ! -f "$FLEX_WHITELIST" ]; then
    skip "Whitelist not found: $FLEX_WHITELIST"
fi
if [ ! -f "$FLEX_SAMPLE_WHITELIST" ]; then
    skip "Sample whitelist not found: $FLEX_SAMPLE_WHITELIST"
fi
if [ ! -f "$FLEX_PROBE_LIST" ]; then
    skip "Probe list not found: $FLEX_PROBE_LIST"
fi
if [ ! -f "$FLEX_SAMPLE_PROBES" ]; then
    skip "Sample probes not found: $FLEX_SAMPLE_PROBES"
fi
if [ ! -f "$FLEX_ALLOWED_TAGS" ]; then
    skip "Allowed tags file not found: $FLEX_ALLOWED_TAGS"
fi
if [ ! -f "$FLEX_FASTQ_R1" ]; then
    skip "FASTQ not found: $FLEX_FASTQ_R1"
fi
if [ ! -f "$FLEX_FASTQ_R2" ]; then
    skip "FASTQ not found: $FLEX_FASTQ_R2"
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
  --genomeDir "$FLEX_INDEX" \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
  --soloCBwhitelist "$FLEX_WHITELIST" \
  --flex yes \
  --soloFlexExpectedCellsPerTag 3000 \
  --soloSampleWhitelist "$FLEX_SAMPLE_WHITELIST" \
  --soloProbeList "$FLEX_PROBE_LIST" \
  --soloSampleProbes "$FLEX_SAMPLE_PROBES" \
  --soloSampleProbeOffset 68 \
  --soloFlexAllowedTags "$FLEX_ALLOWED_TAGS" \
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
    "$FLEX_FASTQ_R2" \
    "$FLEX_FASTQ_R1" \
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
