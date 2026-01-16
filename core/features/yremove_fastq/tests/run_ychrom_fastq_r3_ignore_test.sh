#!/bin/bash
# Validate Y/noY FASTQ emission ignores barcode read (R3) in Solo mode.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-$SCRIPT_DIR/../source/STAR}"
TEST_DIR="${TEST_DIR:-/tmp/yfastq_r3_ignore_test_$$}"
GENOME_DIR="$TEST_DIR/genome"
OUT_DIR="$TEST_DIR/out"
KEEP_TEST_DIR="${KEEP_TEST_DIR:-0}"

cleanup() {
    if [[ "$KEEP_TEST_DIR" -ne 1 ]]; then
        rm -rf "$TEST_DIR"
    else
        echo "Keeping test directory: $TEST_DIR"
    fi
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

cat > "$TEST_DIR/genes.gtf" << 'EOF'
chr1	sim	exon	1	50	.	+	.	gene_id "gene1"; transcript_id "tx1";
chrY	sim	exon	1	50	.	+	.	gene_id "geneY"; transcript_id "ty1";
EOF

echo "=== Building STAR index ==="
"$STAR_BIN" --runMode genomeGenerate \
    --genomeDir "$GENOME_DIR" \
    --genomeFastaFiles "$TEST_DIR/ref.fa" \
    --sjdbGTFfile "$TEST_DIR/genes.gtf" \
    --sjdbOverhang 19 \
    --genomeSAindexNbases 2 \
    --runThreadN 1 \
    --outFileNamePrefix "$TEST_DIR/" \
    > "$TEST_DIR/genome.log" 2>&1

echo "=== Creating synthetic FASTQ (R1/R2 + barcode R3) ==="
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

cat > "$TEST_DIR/sample_R2_001.fastq" << 'EOF'
@readY
GATCCGATCGTAAGCTAGCA
+
IIIIIIIIIIIIIIIIIIII
@readC
CGATCGATCGTCAAGGTCAA
+
IIIIIIIIIIIIIIIIIIII
EOF

cat > "$TEST_DIR/sample_R3_001.fastq" << 'EOF'
@readY
TTTTTTTTTTTTTTTT
+
IIIIIIIIIIIIIIII
@readC
TTTTTTTTTTTTTTTT
+
IIIIIIIIIIIIIIII
EOF

echo "=== Running STAR with Solo + emitYNoYFastq ==="
"$STAR_BIN" \
    --runThreadN 2 \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$TEST_DIR/sample_R1_001.fastq" "$TEST_DIR/sample_R2_001.fastq" "$TEST_DIR/sample_R3_001.fastq" \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist None \
    --soloCBstart 1 --soloCBlen 8 \
    --soloUMIstart 9 --soloUMIlen 8 \
    --alignIntronMax 1 \
    --alignMatesGapMax 100 \
    --alignEndsType EndToEnd \
    --outFilterMismatchNmax 0 \
    --outFilterMismatchNoverReadLmax 0.0 \
    --emitYNoYFastq yes \
    --emitYNoYFastqCompression gz \
    --outSAMtype None \
    --outFileNamePrefix "$OUT_DIR/" \
    > "$OUT_DIR/star.log" 2>&1

Y_R1="$OUT_DIR/sample_Y_R1_001.fastq.gz"
NOY_R1="$OUT_DIR/sample_noY_R1_001.fastq.gz"
Y_R2="$OUT_DIR/sample_Y_R2_001.fastq.gz"
NOY_R2="$OUT_DIR/sample_noY_R2_001.fastq.gz"

echo "=== Validating outputs ==="
for f in "$Y_R1" "$NOY_R1" "$Y_R2" "$NOY_R2"; do
    if [[ ! -f "$f" ]]; then
        echo "FAIL: Missing output FASTQ: $f"
        exit 1
    fi
done

if [[ -f "$OUT_DIR/sample_Y_R3_001.fastq.gz" || -f "$OUT_DIR/sample_noY_R3_001.fastq.gz" ]]; then
    echo "FAIL: R3 FASTQ outputs should not be created"
    exit 1
fi

count_reads() {
    gzip -dc "$1" | grep -c "^@" || true
}

Y_R1_COUNT=$(count_reads "$Y_R1")
NOY_R1_COUNT=$(count_reads "$NOY_R1")
Y_R2_COUNT=$(count_reads "$Y_R2")
NOY_R2_COUNT=$(count_reads "$NOY_R2")

if [[ "$Y_R1_COUNT" -ne 1 || "$NOY_R1_COUNT" -ne 1 || "$Y_R2_COUNT" -ne 1 || "$NOY_R2_COUNT" -ne 1 ]]; then
    echo "FAIL: Unexpected read counts (Y_R1=$Y_R1_COUNT, NOY_R1=$NOY_R1_COUNT, Y_R2=$Y_R2_COUNT, NOY_R2=$NOY_R2_COUNT)"
    exit 1
fi

if gzip -dc "$Y_R1" "$NOY_R1" "$Y_R2" "$NOY_R2" | grep -q "TTTTTTTTTTTTTTTT"; then
    echo "FAIL: Barcode sequence leaked into Y/noY FASTQ output"
    exit 1
fi

echo "=== PASS: R3 barcode read ignored in Y/noY FASTQ emission ==="
