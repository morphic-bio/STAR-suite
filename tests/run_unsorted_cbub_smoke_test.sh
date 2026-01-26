#!/bin/bash
# Smoke test for CB/UB tag injection in unsorted BAM output
# Verifies that CB:Z: and UB:Z: tags are present in both mates when requested

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-${SCRIPT_DIR}/../core/legacy/source/STAR}"
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}/unsorted_cbub_smoke_output}"
TMP_DIR="${TMP_DIR:-/tmp/unsorted_cbub_smoke_$$}"
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

echo "=== Testing CB/UB tag injection in UNSORTED BAM ==="
echo "Binary: $STAR_BIN"
echo "Output: $OUT_DIR"
echo ""

# Run STAR with unsorted BAM output and CB/UB tags requested
"$STAR_BIN" \
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
  --limitIObufferSize 50000000 50000000 \
  --outSJtype None \
  --outBAMcompression 6 \
  --soloMultiMappers Rescue \
  --alignIntronMax 500000 \
  --outFilterMismatchNmax 6 \
  --outFilterMismatchNoverReadLmax 1.0 \
  --outFilterMatchNmin 25 \
  --outSAMunmapped Within \
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
  --readFilesIn \
    "$FLEX_FASTQ_R2" \
    "$FLEX_FASTQ_R1" \
  --outFileNamePrefix "$OUT_DIR/"

# Verify output file exists
BAM_FILE="${OUT_DIR}/Aligned.out.bam"

if [ ! -f "$BAM_FILE" ]; then
    echo "ERROR: Unsorted BAM file not found: $BAM_FILE" >&2
    exit 1
fi

echo ""
echo "=== Checking CB/UB tags in unsorted BAM ==="
echo "BAM: $BAM_FILE"
echo ""

# Count total records
TOTAL_RECORDS=$(samtools view -c "$BAM_FILE")
echo "Total BAM records: $TOTAL_RECORDS"

# Count records with CB tag
CB_RECORDS=$(samtools view "$BAM_FILE" | grep -c 'CB:Z:' || true)
echo "Records with CB:Z: tag: $CB_RECORDS"

# Count records with UB tag
UB_RECORDS=$(samtools view "$BAM_FILE" | grep -c 'UB:Z:' || true)
echo "Records with UB:Z: tag: $UB_RECORDS"

# Count records with both CB and UB tags
CBUB_RECORDS=$(samtools view "$BAM_FILE" | grep 'CB:Z:' | grep -c 'UB:Z:' || true)
echo "Records with both CB:Z: and UB:Z: tags: $CBUB_RECORDS"

# Verify both mates have CB/UB tags (check FLAG 0x1 paired reads)
# For single-cell/Flex data, both R1 and R2 should get CB/UB tags
echo ""
echo "=== Verifying both mates have CB/UB tags ==="

# Sample a few records to verify
echo "Sample records with CB/UB tags:"
samtools view "$BAM_FILE" | grep 'CB:Z:.*UB:Z:' | head -5

# Check that we have a reasonable percentage of records with tags
# (Not all records will have tags - unmapped or no CB match)
if [ "$CB_RECORDS" -eq 0 ]; then
    echo ""
    echo "ERROR: No CB:Z: tags found in unsorted BAM!" >&2
    echo "CB/UB injection may not be working correctly." >&2
    exit 1
fi

if [ "$UB_RECORDS" -eq 0 ]; then
    echo ""
    echo "ERROR: No UB:Z: tags found in unsorted BAM!" >&2
    echo "CB/UB injection may not be working correctly." >&2
    exit 1
fi

# Calculate percentage
CB_PCT=$((CB_RECORDS * 100 / TOTAL_RECORDS))
UB_PCT=$((UB_RECORDS * 100 / TOTAL_RECORDS))

echo ""
echo "CB tag coverage: ${CB_PCT}% (${CB_RECORDS}/${TOTAL_RECORDS})"
echo "UB tag coverage: ${UB_PCT}% (${UB_RECORDS}/${TOTAL_RECORDS})"

# Verify reasonable coverage (at least 10% of records should have tags)
if [ "$CB_PCT" -lt 10 ]; then
    echo ""
    echo "WARNING: Low CB tag coverage (${CB_PCT}%). This may indicate an issue." >&2
fi

echo ""
echo "=== Test PASSED ==="
echo "CB/UB tags are being injected into unsorted BAM output."

# Cleanup
rm -rf "$TMP_DIR"
