#!/bin/bash
# regression_test.sh - TranscriptVB regression testing
#
# This script compares current TranscriptVB output against golden references.
#
# Usage:
#   ./regression_test.sh generate  - Generate new golden references
#   ./regression_test.sh test      - Test against golden references
#   ./regression_test.sh update    - Update golden references after review

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GOLDEN_DIR="$SCRIPT_DIR/golden"
STAR_BIN=${STAR_BIN:-/mnt/pikachu/STAR-Flex/source/STAR}

# Default test data
GENOME_DIR=${GENOME_DIR:-/tmp/star_vb_test/star_new_index}
READS_1=${READS_1:-/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz}
READS_2=${READS_2:-/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz}

TMPDIR=/tmp/regression_test_$$

# Correlation thresholds for regression test
SPEARMAN_THRESHOLD=0.9999
PEARSON_THRESHOLD=0.9999
MAX_TPM_DIFF=1.0

generate_golden() {
    echo "=== Generating Golden References ==="
    mkdir -p "$GOLDEN_DIR"
    mkdir -p $TMPDIR
    cd $TMPDIR
    
    # Generate VB output
    echo "Generating VB golden output..."
    $STAR_BIN --runMode alignReads \
        --genomeDir $GENOME_DIR \
        --readFilesIn $READS_1 $READS_2 \
        --readFilesCommand zcat \
        --quantMode TranscriptVB \
        --runThreadN 1 \
        --outFileNamePrefix vb_
    
    cp vb_quant.sf "$GOLDEN_DIR/vb_quant.sf"
    echo "✓ Saved: $GOLDEN_DIR/vb_quant.sf"
    
    # Generate EM output
    echo "Generating EM golden output..."
    $STAR_BIN --runMode alignReads \
        --genomeDir $GENOME_DIR \
        --readFilesIn $READS_1 $READS_2 \
        --readFilesCommand zcat \
        --quantMode TranscriptVB \
        --quantVBem 1 \
        --runThreadN 1 \
        --outFileNamePrefix em_
    
    cp em_quant.sf "$GOLDEN_DIR/em_quant.sf"
    echo "✓ Saved: $GOLDEN_DIR/em_quant.sf"
    
    # Generate GC bias output
    echo "Generating GC bias golden output..."
    $STAR_BIN --runMode alignReads \
        --genomeDir $GENOME_DIR \
        --readFilesIn $READS_1 $READS_2 \
        --readFilesCommand zcat \
        --quantMode TranscriptVB \
        --quantVBgcBias 1 \
        --runThreadN 1 \
        --outFileNamePrefix gc_
    
    cp gc_quant.sf "$GOLDEN_DIR/gc_quant.sf"
    echo "✓ Saved: $GOLDEN_DIR/gc_quant.sf"
    
    # Save metadata
    cat > "$GOLDEN_DIR/metadata.txt" << EOF
Generated: $(date -Iseconds)
STAR version: $($STAR_BIN --version 2>&1 | head -1)
STAR path: $STAR_BIN
Genome: $GENOME_DIR
Reads: $READS_1 $READS_2
EOF
    
    echo ""
    echo "✓ Golden references generated in: $GOLDEN_DIR"
    rm -rf $TMPDIR
}

run_regression_test() {
    echo "=== Running Regression Tests ==="
    
    if [ ! -d "$GOLDEN_DIR" ]; then
        echo "ERROR: Golden directory not found: $GOLDEN_DIR"
        echo "Run '$0 generate' first."
        exit 1
    fi
    
    mkdir -p $TMPDIR
    cd $TMPDIR
    
    PASSED=0
    FAILED=0
    
    # Test VB
    echo ""
    echo "--- Test: VB Mode ---"
    $STAR_BIN --runMode alignReads \
        --genomeDir $GENOME_DIR \
        --readFilesIn $READS_1 $READS_2 \
        --readFilesCommand zcat \
        --quantMode TranscriptVB \
        --runThreadN 1 \
        --outFileNamePrefix vb_ 2>&1 | tail -3
    
    if compare_quant "$GOLDEN_DIR/vb_quant.sf" "vb_quant.sf" "VB"; then
        PASSED=$((PASSED + 1))
    else
        FAILED=$((FAILED + 1))
    fi
    
    # Test EM
    echo ""
    echo "--- Test: EM Mode ---"
    $STAR_BIN --runMode alignReads \
        --genomeDir $GENOME_DIR \
        --readFilesIn $READS_1 $READS_2 \
        --readFilesCommand zcat \
        --quantMode TranscriptVB \
        --quantVBem 1 \
        --runThreadN 1 \
        --outFileNamePrefix em_ 2>&1 | tail -3
    
    if compare_quant "$GOLDEN_DIR/em_quant.sf" "em_quant.sf" "EM"; then
        PASSED=$((PASSED + 1))
    else
        FAILED=$((FAILED + 1))
    fi
    
    # Test GC bias
    echo ""
    echo "--- Test: GC Bias Mode ---"
    $STAR_BIN --runMode alignReads \
        --genomeDir $GENOME_DIR \
        --readFilesIn $READS_1 $READS_2 \
        --readFilesCommand zcat \
        --quantMode TranscriptVB \
        --quantVBgcBias 1 \
        --runThreadN 1 \
        --outFileNamePrefix gc_ 2>&1 | tail -3
    
    if compare_quant "$GOLDEN_DIR/gc_quant.sf" "gc_quant.sf" "GC"; then
        PASSED=$((PASSED + 1))
    else
        FAILED=$((FAILED + 1))
    fi
    
    echo ""
    echo "========================================"
    echo "Regression Test Summary"
    echo "========================================"
    echo "Passed: $PASSED"
    echo "Failed: $FAILED"
    echo ""
    
    rm -rf $TMPDIR
    
    if [ $FAILED -eq 0 ]; then
        echo "✓ All regression tests PASSED"
        return 0
    else
        echo "✗ Some regression tests FAILED"
        return 1
    fi
}

compare_quant() {
    local golden=$1
    local test=$2
    local name=$3
    
    python3 << EOF
import pandas as pd
from scipy.stats import spearmanr, pearsonr
import sys

golden = pd.read_csv('$golden', sep='\t')
test = pd.read_csv('$test', sep='\t')

merged = golden.merge(test, on='Name', suffixes=('_golden', '_test'))

# Check correlations
r_spear, _ = spearmanr(merged['NumReads_golden'], merged['NumReads_test'])
r_pear, _ = pearsonr(merged['NumReads_golden'], merged['NumReads_test'])

# Check max TPM difference
merged['tpm_diff'] = abs(merged['TPM_golden'] - merged['TPM_test'])
max_tpm_diff = merged['tpm_diff'].max()

print(f"  Spearman: {r_spear:.6f} (threshold: $SPEARMAN_THRESHOLD)")
print(f"  Pearson:  {r_pear:.6f} (threshold: $PEARSON_THRESHOLD)")
print(f"  Max TPM diff: {max_tpm_diff:.4f} (threshold: $MAX_TPM_DIFF)")

if r_spear >= $SPEARMAN_THRESHOLD and r_pear >= $PEARSON_THRESHOLD and max_tpm_diff <= $MAX_TPM_DIFF:
    print("  ✓ $name PASS")
    sys.exit(0)
else:
    print("  ✗ $name FAIL")
    sys.exit(1)
EOF
}

update_golden() {
    echo "=== Updating Golden References ==="
    
    if [ ! -d "$GOLDEN_DIR" ]; then
        echo "No existing golden directory. Use 'generate' instead."
        exit 1
    fi
    
    # Back up old golden
    BACKUP="$GOLDEN_DIR.backup.$(date +%Y%m%d_%H%M%S)"
    cp -r "$GOLDEN_DIR" "$BACKUP"
    echo "Backed up old golden to: $BACKUP"
    
    # Generate new
    generate_golden
    
    echo ""
    echo "✓ Golden references updated"
    echo "  Old golden backed up to: $BACKUP"
}

# Main
case "${1:-test}" in
    generate)
        generate_golden
        ;;
    test)
        run_regression_test
        ;;
    update)
        update_golden
        ;;
    *)
        echo "Usage: $0 {generate|test|update}"
        echo ""
        echo "Commands:"
        echo "  generate - Create new golden references"
        echo "  test     - Test against golden references"
        echo "  update   - Update golden after review (backs up old)"
        exit 1
        ;;
esac

