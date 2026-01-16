#!/bin/bash
# Generate BAM file for alignment-mode parity testing
# Requires: minimap2, samtools

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FIXTURE_DIR="${SCRIPT_DIR}/../../test/fixtures/salmon_eq"
TRANSCRIPTOME="${FIXTURE_DIR}/chr22_trans.fa"
R1="${FIXTURE_DIR}/synthetic_reads1.fq"
R2="${FIXTURE_DIR}/synthetic_reads2.fq"
OUTPUT_BAM="${FIXTURE_DIR}/alignments.bam"

# Check if inputs exist
if [ ! -f "${TRANSCRIPTOME}" ]; then
    echo "Error: Transcriptome not found: ${TRANSCRIPTOME}"
    exit 1
fi

if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
    echo "Error: Read files not found: ${R1} or ${R2}"
    exit 1
fi

# Check for required tools
if ! command -v minimap2 &> /dev/null; then
    echo "Error: minimap2 not found. Install with: conda install -c bioconda minimap2"
    exit 1
fi

if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found. Install with: conda install -c bioconda samtools"
    exit 1
fi

echo "=== Generating BAM file for parity testing ==="
echo "Transcriptome: ${TRANSCRIPTOME}"
echo "Reads: ${R1}, ${R2}"
echo ""

# Align reads to transcriptome with minimap2 (splice-aware mode disabled for transcriptome)
# -ax sr: short-read paired-end alignment
# -t 4: use 4 threads
echo "Aligning reads..."
minimap2 -ax sr -t 4 "${TRANSCRIPTOME}" "${R1}" "${R2}" | \
    samtools view -bS - | \
    samtools sort -o "${OUTPUT_BAM}" -

# Index the BAM
echo "Indexing BAM..."
samtools index "${OUTPUT_BAM}"

echo ""
echo "Created: ${OUTPUT_BAM}"
echo "  $(samtools view -c ${OUTPUT_BAM}) alignments"
echo "  Index: ${OUTPUT_BAM}.bai"
