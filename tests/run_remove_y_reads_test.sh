#!/bin/bash
# tests/run_remove_y_reads_test.sh
# Self-contained test with synthetic data

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOL="$SCRIPT_DIR/../tools/remove_y_reads/remove_y_reads"
TEST_DIR="/tmp/remove_y_reads_test_$$"

mkdir -p "$TEST_DIR"
trap "rm -rf $TEST_DIR" EXIT

echo "=== Creating synthetic test data ==="

# Create minimal reference and index
cat > "$TEST_DIR/ref.fa" << 'EOF'
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>chrY
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
EOF

# Create synthetic BAM with Y reads (read1, read3 on Y)
cat > "$TEST_DIR/y_reads.sam" << 'EOF'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:52
@SQ	SN:chrY	LN:52
read1	0	chrY	1	60	20M	*	0	0	GGCCGGCCGGCCGGCCGGCC	IIIIIIIIIIIIIIIIIIII
read3	0	chrY	10	60	20M	*	0	0	GGCCGGCCGGCCGGCCGGCC	IIIIIIIIIIIIIIIIIIII
EOF

if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found. Please install samtools to run this test."
    exit 1
fi

samtools view -bS "$TEST_DIR/y_reads.sam" > "$TEST_DIR/y_reads.bam"

# Create test FASTQ with 5 reads (read1, read3 are Y; read2, read4, read5 are noY)
cat > "$TEST_DIR/test.fastq" << 'EOF'
@read1 comment here
GGCCGGCCGGCCGGCCGGCC
+
IIIIIIIIIIIIIIIIIIII
@read2
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
@read3/1
GGCCGGCCGGCCGGCCGGCC
+
IIIIIIIIIIIIIIIIIIII
@read4/2
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
@read5
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
EOF

gzip -c "$TEST_DIR/test.fastq" > "$TEST_DIR/test.fastq.gz"

echo "=== Running remove_y_reads ==="
"$TOOL" -y "$TEST_DIR/y_reads.bam" -o "$TEST_DIR" "$TEST_DIR/test.fastq.gz"

echo "=== Validating outputs ==="

# Decompress outputs
gunzip -c "$TEST_DIR/test_Y.fastq.gz" > "$TEST_DIR/test_Y.fastq"
gunzip -c "$TEST_DIR/test_noY.fastq.gz" > "$TEST_DIR/test_noY.fastq"

# Count reads
Y_COUNT=$(grep -c "^@read" "$TEST_DIR/test_Y.fastq" || echo 0)
NOY_COUNT=$(grep -c "^@read" "$TEST_DIR/test_noY.fastq" || echo 0)
TOTAL=$((Y_COUNT + NOY_COUNT))

echo "Y reads: $Y_COUNT (expected: 2)"
echo "noY reads: $NOY_COUNT (expected: 3)"
echo "Total: $TOTAL (expected: 5)"

# Assertion 1: Y + noY = original
if [ "$TOTAL" -ne 5 ]; then
    echo "FAIL: Total count mismatch"
    exit 1
fi

# Assertion 2: Y count correct
if [ "$Y_COUNT" -ne 2 ]; then
    echo "FAIL: Y count mismatch"
    exit 1
fi

# Assertion 3: noY count correct
if [ "$NOY_COUNT" -ne 3 ]; then
    echo "FAIL: noY count mismatch"
    exit 1
fi

# Assertion 4: Y file contains only Y reads (read1, read3)
if ! grep -q "^@read1" "$TEST_DIR/test_Y.fastq"; then
    echo "FAIL: read1 not in Y output"
    exit 1
fi
if ! grep -q "^@read3" "$TEST_DIR/test_Y.fastq"; then
    echo "FAIL: read3 not in Y output"
    exit 1
fi

# Assertion 5: noY file contains no Y reads
if grep -q "^@read1" "$TEST_DIR/test_noY.fastq"; then
    echo "FAIL: read1 found in noY output"
    exit 1
fi
if grep -q "^@read3" "$TEST_DIR/test_noY.fastq"; then
    echo "FAIL: read3 found in noY output"
    exit 1
fi

# Assertion 6: Order preserved (read2 before read4 before read5 in noY)
READ2_LINE=$(grep -n "^@read2" "$TEST_DIR/test_noY.fastq" | cut -d: -f1)
READ4_LINE=$(grep -n "^@read4" "$TEST_DIR/test_noY.fastq" | cut -d: -f1)
READ5_LINE=$(grep -n "^@read5" "$TEST_DIR/test_noY.fastq" | cut -d: -f1)

if [ "$READ2_LINE" -gt "$READ4_LINE" ] || [ "$READ4_LINE" -gt "$READ5_LINE" ]; then
    echo "FAIL: Order not preserved in noY output"
    exit 1
fi

echo ""
echo "=== All tests passed ==="

