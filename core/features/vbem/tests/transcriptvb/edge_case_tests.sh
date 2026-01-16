#!/bin/bash
# edge_case_tests.sh - Test TranscriptVB edge cases
#
# Usage: ./edge_case_tests.sh [STAR_BIN] [GENOME_DIR]

set -e

STAR_BIN=${1:-/mnt/pikachu/STAR-Flex/source/STAR}
GENOME_DIR=${2:-/tmp/star_vb_test/star_new_index}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

OUTDIR=/tmp/edge_case_tests_$$
READS_DIR=$OUTDIR/test_reads
PASSED=0
FAILED=0
SKIPPED=0

mkdir -p $OUTDIR $READS_DIR
cd $OUTDIR

echo "========================================"
echo "TranscriptVB Edge Case Tests"
echo "========================================"
echo "STAR: $STAR_BIN"
echo "Genome: $GENOME_DIR"
echo "Output: $OUTDIR"
echo "========================================"
echo ""

# Helper function to run test
run_test() {
    local name=$1
    local reads1=$2
    local reads2=$3
    local extra_args=$4
    local expected=$5
    
    echo "--- Test: $name ---"
    
    local prefix="${name// /_}_"
    
    # Build command
    local cmd="$STAR_BIN --runMode alignReads --genomeDir $GENOME_DIR"
    
    if [ -n "$reads2" ]; then
        cmd="$cmd --readFilesIn $reads1 $reads2"
    else
        cmd="$cmd --readFilesIn $reads1"
    fi
    
    # Auto-detect compression
    if [[ "$reads1" == *.gz ]]; then
        cmd="$cmd --readFilesCommand zcat"
    fi
    
    cmd="$cmd --quantMode TranscriptVB --outFileNamePrefix $prefix $extra_args"
    
    # Run and capture output
    set +e
    output=$($cmd 2>&1)
    exit_code=$?
    set -e
    
    # Check result
    case "$expected" in
        success)
            if [ $exit_code -eq 0 ] && [ -f "${prefix}quant.sf" ]; then
                echo "  ✓ PASS: Completed successfully"
                PASSED=$((PASSED + 1))
            else
                echo "  ✗ FAIL: Expected success, got exit code $exit_code"
                echo "$output" | tail -5
                FAILED=$((FAILED + 1))
            fi
            ;;
        empty_output)
            if [ -f "${prefix}quant.sf" ]; then
                zeros=$(tail -n +2 "${prefix}quant.sf" | awk '$5 > 0' | wc -l)
                if [ "$zeros" -eq 0 ]; then
                    echo "  ✓ PASS: All transcripts have zero counts (expected)"
                    PASSED=$((PASSED + 1))
                else
                    echo "  ✗ FAIL: Expected zero counts, got $zeros expressed"
                    FAILED=$((FAILED + 1))
                fi
            else
                echo "  ✗ FAIL: quant.sf not generated"
                FAILED=$((FAILED + 1))
            fi
            ;;
        no_crash)
            if [ $exit_code -eq 0 ]; then
                echo "  ✓ PASS: No crash"
                PASSED=$((PASSED + 1))
            else
                echo "  ✗ FAIL: Crashed with exit code $exit_code"
                FAILED=$((FAILED + 1))
            fi
            ;;
        *)
            echo "  ? SKIP: Unknown expected result"
            SKIPPED=$((SKIPPED + 1))
            ;;
    esac
    echo ""
}

# Generate synthetic test reads
generate_test_reads() {
    echo "Generating synthetic test reads..."
    
    # Empty file
    touch $READS_DIR/empty.fq
    gzip -f $READS_DIR/empty.fq
    
    # Tiny file (1 read)
    cat > $READS_DIR/tiny_1.fq << 'EOF'
@read1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
    cat > $READS_DIR/tiny_2.fq << 'EOF'
@read1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
    gzip -f $READS_DIR/tiny_1.fq $READS_DIR/tiny_2.fq
    
    # Random sequences (won't map)
    python3 << PYEOF
import random
random.seed(42)
bases = 'ACGT'

with open('${READS_DIR}/random_1.fq', 'w') as f1, open('${READS_DIR}/random_2.fq', 'w') as f2:
    for i in range(1000):
        seq = ''.join(random.choice(bases) for _ in range(75))
        qual = 'I' * 75
        f1.write(f"@random_{i}\n{seq}\n+\n{qual}\n")
        seq2 = ''.join(random.choice(bases) for _ in range(75))
        f2.write(f"@random_{i}\n{seq2}\n+\n{qual}\n")
PYEOF
    gzip -f $READS_DIR/random_1.fq $READS_DIR/random_2.fq
    
    # Single-end reads (copy from test data if available)
    if [ -f /mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz ]; then
        # Use first 10000 reads for SE test
        zcat /mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz | head -40000 | gzip > $READS_DIR/se_reads.fq.gz
    else
        echo "Warning: Test data not found, some tests may be skipped"
        touch $READS_DIR/se_reads.fq.gz
    fi
    
    echo "✓ Test reads generated"
    echo ""
}

# Generate test reads
generate_test_reads

# === Edge Case Tests ===

# Test: Empty input file
echo "=== Test Group: Empty/Minimal Input ==="
# Note: Empty input will likely cause STAR to fail before TranscriptVB
# This tests graceful error handling
run_test "empty_input" "$READS_DIR/empty.fq.gz" "" "--runThreadN 1" "no_crash"

# Test: Very few reads
run_test "tiny_input" "$READS_DIR/tiny_1.fq.gz" "$READS_DIR/tiny_2.fq.gz" "--runThreadN 1" "no_crash"

# Test: Unmappable reads
run_test "unmappable_reads" "$READS_DIR/random_1.fq.gz" "$READS_DIR/random_2.fq.gz" "--runThreadN 1" "empty_output"

# Test: Single-end reads
echo "=== Test Group: Single-End Processing ==="
if [ -s "$READS_DIR/se_reads.fq.gz" ]; then
    run_test "single_end" "$READS_DIR/se_reads.fq.gz" "" "--runThreadN 1" "success"
else
    echo "  SKIP: SE test data not available"
    SKIPPED=$((SKIPPED + 1))
fi

# Test: Different parameters
echo "=== Test Group: Parameter Variations ==="

# Low prior
run_test "low_prior" \
    "/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz" \
    "/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz" \
    "--runThreadN 1 --quantVBprior 0.0001" \
    "success"

# High prior
run_test "high_prior" \
    "/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz" \
    "/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz" \
    "--runThreadN 1 --quantVBprior 10.0" \
    "success"

# Test: Thread variations
echo "=== Test Group: Multi-Threading ==="
for threads in 1 2 4; do
    run_test "threads_$threads" \
        "/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz" \
        "/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz" \
        "--runThreadN $threads" \
        "success"
done

# Summary
echo "========================================"
echo "Edge Case Test Summary"
echo "========================================"
echo "Passed:  $PASSED"
echo "Failed:  $FAILED"
echo "Skipped: $SKIPPED"
echo ""

if [ $FAILED -eq 0 ]; then
    echo "✓ All edge case tests passed!"
    rm -rf $OUTDIR
    exit 0
else
    echo "✗ Some tests failed"
    echo "Output preserved in: $OUTDIR"
    exit 1
fi

