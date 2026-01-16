#!/bin/bash
# Edge case validation for Y/noY FASTQ emission.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-$SCRIPT_DIR/../source/STAR}"
TEST_DIR="${TEST_DIR:-/tmp/yfastq_edge_cases_$$}"
KEEP_TEST_DIR="${KEEP_TEST_DIR:-0}"

cleanup() {
    if [[ "$KEEP_TEST_DIR" -ne 1 ]]; then
        rm -rf "$TEST_DIR"
    else
        echo "Keeping test directory: $TEST_DIR"
    fi
}
trap cleanup EXIT

mkdir -p "$TEST_DIR"

fail() {
    echo "FAIL: $*" >&2
    exit 1
}

count_fastq_reads() {
    local f="$1"
    gzip -dc "$f" | grep -c "^@" || true
}

count_fastq_reads_plain() {
    local f="$1"
    grep -c "^@" "$f" || true
}

count_fasta_records() {
    local f="$1"
    gzip -dc "$f" | grep -c "^>" || true
}

build_index() {
    local ref_fa="$1"
    local genome_dir="$2"
    mkdir -p "$genome_dir"
    "$STAR_BIN" --runMode genomeGenerate \
        --genomeDir "$genome_dir" \
        --genomeFastaFiles "$ref_fa" \
        --genomeSAindexNbases 2 \
        --runThreadN 1 \
        --outFileNamePrefix "$genome_dir/" \
        > "$genome_dir/genome.log" 2>&1
}

echo "=== Preparing references ==="
REF_Y="$TEST_DIR/ref_with_y.fa"
cat > "$REF_Y" << 'EOF'
>chr1
ACGTTGCAAGCTTAGCTAGCTTGACCTTGACGATCGATCGTAGCTAGCTAG
>chrY
TTGACCGTACGATCGTTCGATGCTAGCTTACGATCGGATCGAATTCGATCGT
EOF

REF_NOY="$TEST_DIR/ref_no_y.fa"
cat > "$REF_NOY" << 'EOF'
>chr1
ACGTTGCAAGCTTAGCTAGCTTGACCTTGACGATCGATCGTAGCTAGCTAG
EOF

GENOME_Y="$TEST_DIR/genome_y"
GENOME_NOY="$TEST_DIR/genome_noy"
build_index "$REF_Y" "$GENOME_Y"
build_index "$REF_NOY" "$GENOME_NOY"

echo "=== Case 1: No Y contigs (Y outputs empty) ==="
CASE1="$TEST_DIR/case_no_y_contigs"
mkdir -p "$CASE1/out"
cat > "$CASE1/sample_R1_001.fastq" << 'EOF'
@readC
ACGTTGCAAGCTTAGCTAGC
+
IIIIIIIIIIIIIIIIIIII
EOF

"$STAR_BIN" \
    --runThreadN 2 \
    --genomeDir "$GENOME_NOY" \
    --readFilesIn "$CASE1/sample_R1_001.fastq" \
    --alignIntronMax 1 \
    --alignEndsType EndToEnd \
    --outFilterMismatchNmax 0 \
    --outFilterMatchNmin 1 \
    --outFilterMatchNminOverLread 0.0 \
    --emitYNoYFastq yes \
    --emitYNoYFastqCompression gz \
    --outSAMtype None \
    --outFileNamePrefix "$CASE1/out/" \
    > "$CASE1/out/star.log" 2>&1

if ! grep -q "No Y-chromosome contigs found" "$CASE1/out/Log.out"; then
    fail "Missing no-Y contig warning"
fi

Y1="$CASE1/out/sample_Y_R1_001.fastq.gz"
NOY1="$CASE1/out/sample_noY_R1_001.fastq.gz"
[[ -f "$Y1" && -f "$NOY1" ]] || fail "Missing outputs for no-Y contigs"
[[ "$(count_fastq_reads "$Y1")" -eq 0 ]] || fail "Expected empty Y output for no-Y contigs"
[[ "$(count_fastq_reads "$NOY1")" -eq 1 ]] || fail "Expected noY output for no-Y contigs"

echo "=== Case 2: Multiple input files per mate (name from first file) ==="
CASE2="$TEST_DIR/case_multi_inputs"
mkdir -p "$CASE2/out"
cat > "$CASE2/multi_R1_001.fastq" << 'EOF'
@readY
TTGACCGTACGATCGTTCGA
+
IIIIIIIIIIIIIIIIIIII
EOF
cat > "$CASE2/multi2_R1_001.fastq" << 'EOF'
@readC
ACGTTGCAAGCTTAGCTAGC
+
IIIIIIIIIIIIIIIIIIII
EOF

"$STAR_BIN" \
    --runThreadN 2 \
    --genomeDir "$GENOME_Y" \
    --readFilesIn "$CASE2/multi_R1_001.fastq,$CASE2/multi2_R1_001.fastq" \
    --alignIntronMax 1 \
    --alignEndsType EndToEnd \
    --outFilterMismatchNmax 0 \
    --outFilterMatchNmin 1 \
    --outFilterMatchNminOverLread 0.0 \
    --emitYNoYFastq yes \
    --emitYNoYFastqCompression gz \
    --outSAMtype None \
    --outFileNamePrefix "$CASE2/out/" \
    > "$CASE2/out/star.log" 2>&1

Y2="$CASE2/out/multi_Y_R1_001.fastq.gz"
NOY2="$CASE2/out/multi_noY_R1_001.fastq.gz"
[[ -f "$Y2" && -f "$NOY2" ]] || fail "Missing outputs for multi-input case"
[[ ! -f "$CASE2/out/multi2_Y_R1_001.fastq.gz" ]] || fail "Unexpected output name from second file"
[[ "$(count_fastq_reads "$Y2")" -eq 1 ]] || fail "Expected 1 Y read for multi-input case"
[[ "$(count_fastq_reads "$NOY2")" -eq 1 ]] || fail "Expected 1 noY read for multi-input case"

echo "=== Case 3: No _R1/_R2 token (fallback naming) ==="
CASE3="$TEST_DIR/case_fallback_names"
mkdir -p "$CASE3/out"
cat > "$CASE3/sample.fastq" << 'EOF'
@readY
TTGACCGTACGATCGTTCGA
+
IIIIIIIIIIIIIIIIIIII
@readC
ACGTTGCAAGCTTAGCTAGC
+
IIIIIIIIIIIIIIIIIIII
EOF

"$STAR_BIN" \
    --runThreadN 2 \
    --genomeDir "$GENOME_Y" \
    --readFilesIn "$CASE3/sample.fastq" \
    --alignIntronMax 1 \
    --alignEndsType EndToEnd \
    --outFilterMismatchNmax 0 \
    --outFilterMatchNmin 1 \
    --outFilterMatchNminOverLread 0.0 \
    --emitYNoYFastq yes \
    --emitYNoYFastqCompression gz \
    --outSAMtype None \
    --outFileNamePrefix "$CASE3/out/" \
    > "$CASE3/out/star.log" 2>&1

Y3="$CASE3/out/Y_reads.mate1.fastq.gz"
NOY3="$CASE3/out/noY_reads.mate1.fastq.gz"
[[ -f "$Y3" && -f "$NOY3" ]] || fail "Missing fallback outputs"
[[ "$(count_fastq_reads "$Y3")" -eq 1 ]] || fail "Expected 1 Y read for fallback naming case"
[[ "$(count_fastq_reads "$NOY3")" -eq 1 ]] || fail "Expected 1 noY read for fallback naming case"

echo "=== Case 4: FASTA input ==="
CASE4="$TEST_DIR/case_fasta_input"
mkdir -p "$CASE4/out"
cat > "$CASE4/sample_R1_001.fa" << 'EOF'
>readY
TTGACCGTACGATCGTTCGA
>readC
ACGTTGCAAGCTTAGCTAGC
EOF

"$STAR_BIN" \
    --runThreadN 2 \
    --genomeDir "$GENOME_Y" \
    --readFilesIn "$CASE4/sample_R1_001.fa" \
    --alignIntronMax 1 \
    --alignEndsType EndToEnd \
    --outFilterMismatchNmax 0 \
    --outFilterMatchNmin 1 \
    --outFilterMatchNminOverLread 0.0 \
    --emitYNoYFastq yes \
    --emitYNoYFastqCompression gz \
    --outSAMtype None \
    --outFileNamePrefix "$CASE4/out/" \
    > "$CASE4/out/star.log" 2>&1

Y4="$CASE4/out/sample_Y_R1_001.fa.gz"
NOY4="$CASE4/out/sample_noY_R1_001.fa.gz"
[[ -f "$Y4" && -f "$NOY4" ]] || fail "Missing FASTA outputs"
[[ "$(count_fasta_records "$Y4")" -eq 1 ]] || fail "Expected 1 FASTA record in Y output"
[[ "$(count_fasta_records "$NOY4")" -eq 1 ]] || fail "Expected 1 FASTA record in noY output"
if gzip -dc "$Y4" | grep -q "^+"; then
    fail "FASTA output contains FASTQ '+' separator"
fi
if ! gzip -dc "$Y4" | head -n 1 | grep -q "^>"; then
    fail "FASTA output header missing '>' prefix"
fi

echo "=== Case 5: Uncompressed output ==="
CASE5="$TEST_DIR/case_uncompressed"
mkdir -p "$CASE5/out"
cat > "$CASE5/sample_R1_001.fastq" << 'EOF'
@readY
TTGACCGTACGATCGTTCGA
+
IIIIIIIIIIIIIIIIIIII
@readC
ACGTTGCAAGCTTAGCTAGC
+
IIIIIIIIIIIIIIIIIIII
EOF

"$STAR_BIN" \
    --runThreadN 2 \
    --genomeDir "$GENOME_Y" \
    --readFilesIn "$CASE5/sample_R1_001.fastq" \
    --alignIntronMax 1 \
    --alignEndsType EndToEnd \
    --outFilterMismatchNmax 0 \
    --outFilterMatchNmin 1 \
    --outFilterMatchNminOverLread 0.0 \
    --emitYNoYFastq yes \
    --emitYNoYFastqCompression none \
    --outSAMtype None \
    --outFileNamePrefix "$CASE5/out/" \
    > "$CASE5/out/star.log" 2>&1

Y5="$CASE5/out/sample_Y_R1_001.fastq"
NOY5="$CASE5/out/sample_noY_R1_001.fastq"
[[ -f "$Y5" && -f "$NOY5" ]] || fail "Missing uncompressed outputs"
[[ "$(count_fastq_reads_plain "$Y5")" -eq 1 ]] || fail "Expected 1 Y read for uncompressed output"
[[ "$(count_fastq_reads_plain "$NOY5")" -eq 1 ]] || fail "Expected 1 noY read for uncompressed output"

echo "=== Case 6: Unmapped reads route to noY ==="
CASE6="$TEST_DIR/case_unmapped"
mkdir -p "$CASE6/out"
cat > "$CASE6/sample_R1_001.fastq" << 'EOF'
@readX
GGGGGGGGGGGGGGGGGGGG
+
IIIIIIIIIIIIIIIIIIII
@readZ
CCCCCCCCCCCCCCCCCCCC
+
IIIIIIIIIIIIIIIIIIII
EOF

"$STAR_BIN" \
    --runThreadN 2 \
    --genomeDir "$GENOME_Y" \
    --readFilesIn "$CASE6/sample_R1_001.fastq" \
    --alignIntronMax 1 \
    --alignEndsType EndToEnd \
    --outFilterMismatchNmax 0 \
    --outFilterMatchNmin 0 \
    --outFilterMatchNminOverLread 0.0 \
    --emitYNoYFastq yes \
    --emitYNoYFastqCompression gz \
    --outSAMtype None \
    --outFileNamePrefix "$CASE6/out/" \
    > "$CASE6/out/star.log" 2>&1

Y6="$CASE6/out/sample_Y_R1_001.fastq.gz"
NOY6="$CASE6/out/sample_noY_R1_001.fastq.gz"
[[ -f "$Y6" && -f "$NOY6" ]] || fail "Missing outputs for unmapped case"
[[ "$(count_fastq_reads "$Y6")" -eq 0 ]] || fail "Expected 0 Y reads for unmapped case"
[[ "$(count_fastq_reads "$NOY6")" -eq 2 ]] || fail "Expected 2 noY reads for unmapped case"

echo "=== PASS: all Y/noY FASTQ edge cases ==="
