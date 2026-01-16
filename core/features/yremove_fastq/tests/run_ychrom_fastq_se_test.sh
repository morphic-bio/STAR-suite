#!/bin/bash
# Validate Y/noY FASTQ emission for single-end reads and naming.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-$SCRIPT_DIR/../source/STAR}"
TEST_DIR="${TEST_DIR:-/tmp/yfastq_se_test_$$}"
GENOME_DIR="$TEST_DIR/genome"
OUT_DIR="$TEST_DIR/out"

cleanup() {
    rm -rf "$TEST_DIR"
}
trap cleanup EXIT

mkdir -p "$GENOME_DIR" "$OUT_DIR"

echo "=== Creating synthetic reference ==="
cat > "$TEST_DIR/ref.fa" << 'EOF'
>chr1
ACGTTGCAAGCTTAGCTAGCTTGACCTTGACGATCGATCGTAGCTAGCTAG
>chrY
TTGACCGTACGATCGTTCGATGCTAGCTTACGATCGGATCGAATTCGATCGT
EOF

echo "=== Building STAR index ==="
"$STAR_BIN" --runMode genomeGenerate \
    --genomeDir "$GENOME_DIR" \
    --genomeFastaFiles "$TEST_DIR/ref.fa" \
    --genomeSAindexNbases 2 \
    --runThreadN 1 \
    --outFileNamePrefix "$TEST_DIR/" \
    > "$TEST_DIR/genome.log" 2>&1

echo "=== Creating synthetic single-end FASTQ ==="
cat > "$TEST_DIR/sample_R1_001.fastq" << 'EOF'
@readY
TTGACCGTACGATCGTTCGA
+
IIIIIIIIIIIIIIIIIIII
@readC
ACGTTGCAAGCTTAGCTAGC
+
IIIIIIIIIIIIIIIIIIII
EOF

echo "=== Running STAR (single-end) ==="
"$STAR_BIN" \
    --runThreadN 2 \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$TEST_DIR/sample_R1_001.fastq" \
    --alignIntronMax 1 \
    --alignEndsType EndToEnd \
    --outFilterMismatchNmax 0 \
    --outFilterMatchNmin 1 \
    --outFilterMatchNminOverLread 0.0 \
    --emitYNoYFastq yes \
    --emitYNoYFastqCompression gz \
    --outSAMtype None \
    --outFileNamePrefix "$OUT_DIR/" \
    > "$OUT_DIR/star.log" 2>&1

Y_R1="$OUT_DIR/sample_Y_R1_001.fastq.gz"
NOY_R1="$OUT_DIR/sample_noY_R1_001.fastq.gz"

echo "=== Validating outputs ==="
if [[ ! -f "$Y_R1" || ! -f "$NOY_R1" ]]; then
    echo "FAIL: Missing output FASTQ (expected $Y_R1 and $NOY_R1)"
    exit 1
fi

if [[ -f "$OUT_DIR/sample_Y_R2_001.fastq.gz" || -f "$OUT_DIR/sample_noY_R2_001.fastq.gz" ]]; then
    echo "FAIL: R2 outputs should not be created for single-end reads"
    exit 1
fi

count_reads() {
    gzip -dc "$1" | grep -c "^@" || true
}

Y_COUNT=$(count_reads "$Y_R1")
NOY_COUNT=$(count_reads "$NOY_R1")

if [[ "$Y_COUNT" -ne 1 || "$NOY_COUNT" -ne 1 ]]; then
    echo "FAIL: Unexpected read counts (Y=$Y_COUNT, noY=$NOY_COUNT)"
    exit 1
fi

if ! gzip -dc "$Y_R1" | grep -q "^@readY"; then
    echo "FAIL: readY not found in Y FASTQ"
    exit 1
fi
if gzip -dc "$Y_R1" | grep -q "^@readC"; then
    echo "FAIL: readC incorrectly found in Y FASTQ"
    exit 1
fi
if ! gzip -dc "$NOY_R1" | grep -q "^@readC"; then
    echo "FAIL: readC not found in noY FASTQ"
    exit 1
fi
if gzip -dc "$NOY_R1" | grep -q "^@readY"; then
    echo "FAIL: readY incorrectly found in noY FASTQ"
    exit 1
fi

echo "=== PASS: single-end Y/noY FASTQ emission ==="
