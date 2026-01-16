#!/bin/bash
# quick_test.sh - Run basic TranscriptVB validation tests
#
# Usage: ./quick_test.sh [STAR_BIN] [GENOME_DIR] [READS_1] [READS_2]
#
# If arguments not provided, uses defaults for STAR-Flex test data

set -e

# Default paths (adjust for your environment)
STAR_BIN=${1:-/mnt/pikachu/STAR-Flex/source/STAR}
GENOME_DIR=${2:-/tmp/star_vb_test/star_new_index}
READS_1=${3:-/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz}
READS_2=${4:-/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz}

OUTDIR=/tmp/transcriptvb_quick_test_$$
PASSED=0
FAILED=0

mkdir -p $OUTDIR
cd $OUTDIR

echo "========================================"
echo "TranscriptVB Quick Validation Tests"
echo "========================================"
echo "STAR: $STAR_BIN"
echo "Genome: $GENOME_DIR"
echo "Reads: $READS_1, $READS_2"
echo "Output: $OUTDIR"
echo "========================================"
echo ""

# Test 1: Basic TranscriptVB
echo "=== Test 1: Basic TranscriptVB ==="
$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --runThreadN 4 \
    --outFileNamePrefix basic_ \
    2>&1 | tail -5

if [ -f basic_quant.sf ]; then
    echo "✓ quant.sf created"
    
    # Check TPM sum
    TPM_SUM=$(tail -n +2 basic_quant.sf | awk '{sum+=$4} END {printf "%.0f", sum}')
    echo "  TPM sum: $TPM_SUM (expected ~1000000)"
    
    # Check for expressed transcripts
    EXPRESSED=$(tail -n +2 basic_quant.sf | awk '$5 > 0 {count++} END {print count}')
    echo "  Expressed transcripts: $EXPRESSED"
    
    # Check convergence
    if grep -q "converged: yes" basic_Log.out 2>/dev/null; then
        echo "  ✓ Quantification converged"
        PASSED=$((PASSED + 1))
    else
        echo "  ✗ Convergence status unknown"
        FAILED=$((FAILED + 1))
    fi
    
    # Test 1b: Gene-level output (quant.genes.sf)
    if [ -f basic_quant.genes.sf ]; then
        echo "  ✓ quant.genes.sf created"
        
        # Check header
        HEADER=$(head -1 basic_quant.genes.sf)
        if [ "$HEADER" = "Name	Length	EffectiveLength	TPM	NumReads" ]; then
            echo "    ✓ Header format correct"
        else
            echo "    ✗ Header format incorrect: $HEADER"
            FAILED=$((FAILED + 1))
        fi
        
        # Check line count (should be genes + header)
        GENE_COUNT=$(tail -n +2 basic_quant.genes.sf | wc -l)
        echo "    Gene count: $GENE_COUNT"
        
        # Check gene TPM sum (should match transcript TPM sum approximately)
        GENE_TPM_SUM=$(tail -n +2 basic_quant.genes.sf | awk '{sum+=$4} END {printf "%.0f", sum}')
        echo "    Gene TPM sum: $GENE_TPM_SUM (should match transcript TPM sum)"
        
        # Check for no scientific notation
        if grep -qE '[eE][+-]' basic_quant.genes.sf; then
            echo "    ✗ Scientific notation found in output"
            FAILED=$((FAILED + 1))
        else
            echo "    ✓ No scientific notation"
        fi
        
        # Check for NaN/inf
        if grep -qiE '(nan|inf)' basic_quant.genes.sf; then
            echo "    ✗ NaN or Inf values found"
            FAILED=$((FAILED + 1))
        else
            echo "    ✓ No NaN/Inf values"
        fi
        
        # Check decimal precision (spot check first few lines)
        # Length should have max 3 decimals, EffectiveLength max 4, TPM max 6, NumReads max 3
        PRECISION_OK=true
        while IFS=$'\t' read -r name len efflen tpm counts; do
            # Check Length (column 2) - max 3 decimals
            if echo "$len" | grep -qE '\.'; then
                DECIMALS=$(echo "$len" | sed 's/.*\.//' | wc -c)
                if [ $DECIMALS -gt 4 ]; then  # 4 = 3 decimals + newline
                    echo "    ✗ Length precision too high: $len"
                    PRECISION_OK=false
                fi
            fi
            # Check EffectiveLength (column 3) - max 4 decimals
            if echo "$efflen" | grep -qE '\.'; then
                DECIMALS=$(echo "$efflen" | sed 's/.*\.//' | wc -c)
                if [ $DECIMALS -gt 5 ]; then  # 5 = 4 decimals + newline
                    echo "    ✗ EffectiveLength precision too high: $efflen"
                    PRECISION_OK=false
                fi
            fi
            # Check NumReads (column 5) - max 3 decimals
            if echo "$counts" | grep -qE '\.'; then
                DECIMALS=$(echo "$counts" | sed 's/.*\.//' | wc -c)
                if [ $DECIMALS -gt 4 ]; then  # 4 = 3 decimals + newline
                    echo "    ✗ NumReads precision too high: $counts"
                    PRECISION_OK=false
                fi
            fi
        done < <(tail -n +2 basic_quant.genes.sf | head -5)
        
        if [ "$PRECISION_OK" = true ]; then
            echo "    ✓ Decimal precision within limits"
        else
            echo "    ✗ Decimal precision check failed"
            FAILED=$((FAILED + 1))
        fi
        
        PASSED=$((PASSED + 1))
    else
        echo "  ✗ quant.genes.sf NOT created (expected by default)"
        FAILED=$((FAILED + 1))
    fi
else
    echo "✗ quant.sf NOT created"
    FAILED=$((FAILED + 1))
fi
echo ""

# Test 2: GC Bias Collection
echo "=== Test 2: GC Bias Collection ==="
$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --quantVBgcBias 1 \
    --runThreadN 4 \
    --outFileNamePrefix gc_ \
    2>&1 | tail -5

if grep -q "GC bias: collected" gc_Log.out 2>/dev/null; then
    GC_OBS=$(grep "GC bias: collected" gc_Log.out | sed -n 's/.*collected \([0-9]*\) fragment.*/\1/p')
    echo "✓ GC observations collected: $GC_OBS"
    PASSED=$((PASSED + 1))
else
    echo "✗ No GC observations found in log"
    FAILED=$((FAILED + 1))
fi

if grep -q "FLD-adjusted" gc_Log.out 2>/dev/null; then
    echo "✓ FLD-adjusted effective lengths"
else
    echo "  (FLD adjustment not logged)"
fi
echo ""

# Test 3: EM Mode
echo "=== Test 3: EM Mode ==="
$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --quantVBem 1 \
    --runThreadN 4 \
    --outFileNamePrefix em_ \
    2>&1 | tail -5

if [ -f em_quant.sf ]; then
    echo "✓ EM mode completed"
    PASSED=$((PASSED + 1))
else
    echo "✗ EM mode failed"
    FAILED=$((FAILED + 1))
fi
echo ""

# Test 4: Single Thread Determinism
echo "=== Test 4: Single Thread Determinism ==="
$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --runThreadN 1 \
    --outFileNamePrefix st1_ \
    2>&1 | tail -3

$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --runThreadN 1 \
    --outFileNamePrefix st2_ \
    2>&1 | tail -3

if diff -q st1_quant.sf st2_quant.sf > /dev/null 2>&1; then
    echo "✓ Deterministic results (single thread)"
    PASSED=$((PASSED + 1))
else
    echo "✗ Non-deterministic results"
    FAILED=$((FAILED + 1))
fi
echo ""

# Test 5: Gene quant disabled
echo "=== Test 5: Gene Quant Disabled (--quantVBgenes 0) ==="
$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --quantVBgenes 0 \
    --runThreadN 4 \
    --outFileNamePrefix nogene_ \
    2>&1 | tail -3

if [ -f nogene_quant.sf ] && [ ! -f nogene_quant.genes.sf ]; then
    echo "✓ quant.genes.sf correctly NOT created when disabled"
    PASSED=$((PASSED + 1))
else
    if [ -f nogene_quant.genes.sf ]; then
        echo "✗ quant.genes.sf created when it should be disabled"
    else
        echo "✗ quant.sf not created"
    fi
    FAILED=$((FAILED + 1))
fi
echo ""

# Summary
echo "========================================"
echo "Test Summary"
echo "========================================"
echo "Passed: $PASSED"
echo "Failed: $FAILED"
echo ""

if [ $FAILED -eq 0 ]; then
    echo "✓ All tests passed!"
    echo ""
    echo "Output files in: $OUTDIR"
    exit 0
else
    echo "✗ Some tests failed"
    echo ""
    echo "Check output in: $OUTDIR"
    exit 1
fi

