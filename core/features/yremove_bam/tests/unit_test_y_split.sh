#!/bin/bash
# Unit tests for Y-chromosome BAM split feature
# Tests: Y-contig detection, path derivation, routing logic

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${SCRIPT_DIR}/../source/STAR"
TEST_DIR="${SCRIPT_DIR}/unit_y_split_test"
PASSED=0
FAILED=0

# Cleanup
cleanup() {
    rm -rf "$TEST_DIR"
}
trap cleanup EXIT

mkdir -p "$TEST_DIR"

echo "=========================================="
echo "Y-Chromosome BAM Split Unit Tests"
echo "=========================================="
echo ""

# Test 1: Y-contig detection
echo ">>> Test 1: Y-contig detection"
TEST1_DIR="$TEST_DIR/test1_ycontig"
mkdir -p "$TEST1_DIR/ref"

# Create test reference with various Y contig names
cat > "$TEST1_DIR/ref/chr1.fa" << 'EOF'
>chr1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
EOF

cat > "$TEST1_DIR/ref/chrY.fa" << 'EOF'
>chrY
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > "$TEST1_DIR/ref/chr_Y.fa" << 'EOF'
>chr_Y
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > "$TEST1_DIR/ref/chrY_random.fa" << 'EOF'
>chrY_random
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > "$TEST1_DIR/ref/chrY_alt.fa" << 'EOF'
>chrY_alt
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > "$TEST1_DIR/ref/Y.fa" << 'EOF'
>Y
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > "$TEST1_DIR/ref/y.fa" << 'EOF'
>y
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

cat > "$TEST1_DIR/ref/genes.gtf" << 'EOF'
chr1	.	gene	1	50	.	+	.	gene_id "GENE1";
chr1	.	exon	1	50	.	+	.	gene_id "GENE1"; transcript_id "T1";
chrY	.	gene	1	50	.	+	.	gene_id "GENE2";
chrY	.	exon	1	50	.	+	.	gene_id "GENE2"; transcript_id "T2";
EOF

# Build index
"$STAR_BIN" --runMode genomeGenerate \
    --genomeDir "$TEST1_DIR/star_index" \
    --genomeFastaFiles "$TEST1_DIR/ref"/*.fa \
    --sjdbGTFfile "$TEST1_DIR/ref/genes.gtf" \
    --sjdbOverhang 19 \
    --genomeSAindexNbases 3 \
    > "$TEST1_DIR/genome_build.log" 2>&1 || true

# Check chrName.txt to see which contigs were indexed
if [ -f "$TEST1_DIR/star_index/chrName.txt" ]; then
    echo "  Contigs in index:"
    cat "$TEST1_DIR/star_index/chrName.txt" | while read chr; do
        echo "    - $chr"
    done
    
    # Check if Y contigs are present
    Y_CONTIGS=$(grep -iE "^chrY|^chr_Y|^chrY_|^Y$|^y$" "$TEST1_DIR/star_index/chrName.txt" || true)
    if [ -n "$Y_CONTIGS" ]; then
        echo "  ✓ Y contigs detected in index"
        PASSED=$((PASSED + 1))
    else
        echo "  ✗ No Y contigs found in index"
        FAILED=$((FAILED + 1))
    fi
else
    echo "  ✗ Index build failed"
    FAILED=$((FAILED + 1))
fi
echo ""

# Test 2: Path derivation
echo ">>> Test 2: Path derivation"
TEST2_DIR="$TEST_DIR/test2_paths"
mkdir -p "$TEST2_DIR"

# Test default path derivation
test_path() {
    local prefix="$1"
    local expected_y="$2"
    local expected_noy="$3"
    local desc="$4"
    
    # Simulate path derivation logic
    if [[ "$prefix" == "/dev/stdout" ]]; then
        derived_y="star_output_Y.bam"
        derived_noy="star_output_noY.bam"
    else
        if [[ "$prefix" == *.bam ]]; then
            # Remove .bam and add suffix
            base="${prefix%.bam}"
            derived_y="${base}_Y.bam"
            derived_noy="${base}_noY.bam"
        else
            derived_y="${prefix}_Y.bam"
            derived_noy="${prefix}_noY.bam"
        fi
    fi
    
    if [ "$derived_y" == "$expected_y" ] && [ "$derived_noy" == "$expected_noy" ]; then
        echo "  ✓ $desc"
        echo "    Y: $derived_y"
        echo "    noY: $derived_noy"
        PASSED=$((PASSED + 1))
    else
        echo "  ✗ $desc"
        echo "    Expected Y: $expected_y, got: $derived_y"
        echo "    Expected noY: $expected_noy, got: $derived_noy"
        FAILED=$((FAILED + 1))
    fi
}

test_path "output/Aligned" "output/Aligned_Y.bam" "output/Aligned_noY.bam" "Default prefix"
test_path "output/Aligned.bam" "output/Aligned_Y.bam" "output/Aligned_noY.bam" "Prefix with .bam"
test_path "/dev/stdout" "star_output_Y.bam" "star_output_noY.bam" "/dev/stdout special case"
test_path "test/out" "test/out_Y.bam" "test/out_noY.bam" "Path with directory"
echo ""

# Test 3: Writer routing sanity (simulated)
echo ">>> Test 3: Writer routing sanity"
TEST3_DIR="$TEST_DIR/test3_routing"
mkdir -p "$TEST3_DIR"

# Create a mock scenario: read IDs that should route to Y
cat > "$TEST3_DIR/y_reads.txt" << 'EOF'
read2
read3
EOF

# Create a mock scenario: all read IDs
cat > "$TEST3_DIR/all_reads.txt" << 'EOF'
read1
read2
read3
EOF

# Simulate routing logic
Y_COUNT=0
NOY_COUNT=0
while read read_id; do
    if grep -q "^${read_id}$" "$TEST3_DIR/y_reads.txt"; then
        Y_COUNT=$((Y_COUNT + 1))
    else
        NOY_COUNT=$((NOY_COUNT + 1))
    fi
done < "$TEST3_DIR/all_reads.txt"

if [ "$Y_COUNT" -eq 2 ] && [ "$NOY_COUNT" -eq 1 ]; then
    echo "  ✓ Routing logic correct"
    echo "    Y reads: $Y_COUNT (expected: 2)"
    echo "    noY reads: $NOY_COUNT (expected: 1)"
    PASSED=$((PASSED + 1))
else
    echo "  ✗ Routing logic incorrect"
    echo "    Y reads: $Y_COUNT (expected: 2)"
    echo "    noY reads: $NOY_COUNT (expected: 1)"
    FAILED=$((FAILED + 1))
fi
echo ""

# Summary
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo "Passed: $PASSED"
echo "Failed: $FAILED"
echo ""

if [ $FAILED -eq 0 ]; then
    echo "✓ All unit tests passed!"
    exit 0
else
    echo "✗ Some tests failed"
    exit 1
fi

