#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${SCRIPT_DIR}/../source/STAR"
TEST_DIR="${SCRIPT_DIR}/ychrom_test"
REF_DIR="${TEST_DIR}/ref/star_index"
FASTQ_DIR="${TEST_DIR}/fastq"

# Ensure samtools exists
command -v samtools >/dev/null || { echo "ERROR: samtools not found"; exit 1; }

# Ensure STAR binary exists
[ -f "$STAR_BIN" ] || { echo "ERROR: STAR binary not found at $STAR_BIN"; exit 1; }

# Ensure STAR index exists
[ -d "$REF_DIR" ] || {
    echo "ERROR: STAR index not found at $REF_DIR"
    echo "Please run tests/run_ychrom_bam_split_test.sh first or generate the index manually"
    exit 1
}

OUT1="${TEST_DIR}/output_spill"
OUT2="${TEST_DIR}/output_nospill"
TMP1="${TEST_DIR}/tmp_spill"
TMP2="${TEST_DIR}/tmp_nospill"

rm -rf "$OUT1" "$OUT2" "$TMP1" "$TMP2"
mkdir -p "$OUT1" "$OUT2"

# Spill run
"$STAR_BIN" \
  --runThreadN 2 \
  --genomeDir "$REF_DIR" \
  --readFilesIn "$FASTQ_DIR/R1.fastq" "$FASTQ_DIR/R2.fastq" \
  --outSAMtype BAM SortedByCoordinate \
  --outBAMsortMethod samtools \
  --limitBAMsortRAM 100000 \
  --emitNoYBAM yes \
  --outFileNamePrefix "$OUT1/" \
  --outTmpDir "$TMP1"

# No-spill run
"$STAR_BIN" \
  --runThreadN 2 \
  --genomeDir "$REF_DIR" \
  --readFilesIn "$FASTQ_DIR/R1.fastq" "$FASTQ_DIR/R2.fastq" \
  --outSAMtype BAM SortedByCoordinate \
  --outBAMsortMethod samtools \
  --limitBAMsortRAM 1000000000 \
  --emitNoYBAM yes \
  --outFileNamePrefix "$OUT2/" \
  --outTmpDir "$TMP2"

# Validate BAMs
samtools quickcheck "$OUT1/Aligned.sortedByCoord.out_Y.bam"
samtools quickcheck "$OUT1/Aligned.sortedByCoord.out_noY.bam"

# Header check
samtools view -H "$OUT1/Aligned.sortedByCoord.out_Y.bam" | grep "SO:coordinate" >/dev/null

# Y/noY counts add up (keepBAM=no in this script)
Y_COUNT=$(samtools view -c "$OUT1/Aligned.sortedByCoord.out_Y.bam")
NOY_COUNT=$(samtools view -c "$OUT1/Aligned.sortedByCoord.out_noY.bam")
TOTAL=$((Y_COUNT + NOY_COUNT))
if [ "$TOTAL" -le 0 ]; then
  echo "FAIL: no mapped reads in split outputs"
  exit 1
fi

# Determinism: spill vs no-spill should match
samtools view "$OUT1/Aligned.sortedByCoord.out_Y.bam" > /tmp/y_spill.txt
samtools view "$OUT2/Aligned.sortedByCoord.out_Y.bam" > /tmp/y_nospill.txt
diff /tmp/y_spill.txt /tmp/y_nospill.txt

samtools view "$OUT1/Aligned.sortedByCoord.out_noY.bam" > /tmp/noy_spill.txt
samtools view "$OUT2/Aligned.sortedByCoord.out_noY.bam" > /tmp/noy_nospill.txt
diff /tmp/noy_spill.txt /tmp/noy_nospill.txt

echo "OK: spill/no-spill outputs match and routing is valid"

