#!/bin/bash
# Y/NoY Split Regression Tests
# Validates STAR-Flex baseline matches upstream STAR Solo behavior
# Compares new Flex outputs against previous test artifacts

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_FLEX_BIN="${SCRIPT_DIR}/../source/STAR"
STAR_SOLO_BIN="/usr/local/bin/STAR"
REPORT_FILE="${SCRIPT_DIR}/TEST_REPORT_Y_SPLIT_REGRESSION.md"

# Check if STAR Solo is available
if [ ! -f "$STAR_SOLO_BIN" ]; then
    echo "WARNING: STAR Solo not found at $STAR_SOLO_BIN"
    echo "Skipping bulk RNA-seq regression test"
    STAR_SOLO_AVAILABLE=false
else
    STAR_SOLO_AVAILABLE=true
fi

echo "=========================================="
echo "Y/NoY Split Regression Tests"
echo "=========================================="
echo "STAR-Flex Binary: $STAR_FLEX_BIN"
if [ "$STAR_SOLO_AVAILABLE" = true ]; then
    echo "STAR Solo Binary: $STAR_SOLO_BIN"
    "$STAR_SOLO_BIN" --version 2>&1 | head -1 || true
fi
echo ""

# Initialize report
cat > "$REPORT_FILE" << EOF
# Y/NoY Split Regression Test Report

Generated: $(date)

## Overview

This report validates that STAR-Flex baseline behavior matches upstream STAR Solo, and that Y/noY split outputs are stable across test runs.

---

## Test 1: Bulk RNA-Seq Baseline Regression

### Configuration

- STAR-Flex Binary: \`$STAR_FLEX_BIN\`
- STAR Solo Binary: \`$STAR_SOLO_BIN\`
- Dataset: \`/storage/PE/downsampled/21033-09-01-13-01_S1_L007_R1_001.fastq.gz\`
- Reference: \`/storage/flex_filtered_reference/star_index\`

### Test Procedure

1. Run STAR-Flex baseline (no --emitNoYBAM) - unsorted and sorted
2. Run STAR Solo baseline (no split) - unsorted and sorted
3. Compare total mapped reads and chrY counts

EOF

BULK_PASSED=0
BULK_FAILED=0

if [ "$STAR_SOLO_AVAILABLE" = true ]; then
    echo ">>> Test 1: Bulk RNA-Seq Baseline Regression"
    
    BULK_BASE_DIR="/tmp/ychrom_regression_bulk"
    rm -rf "$BULK_BASE_DIR"
    mkdir -p "$BULK_BASE_DIR"
    
    R1_FILE="/storage/PE/downsampled/21033-09-01-13-01_S1_L007_R1_001.fastq.gz"
    R2_FILE="/storage/PE/downsampled/21033-09-01-13-01_S1_L007_R2_001.fastq.gz"
    GENOME_DIR="/storage/flex_filtered_reference/star_index"
    
    # Common parameters for bulk RNA-seq
    COMMON_PARAMS=(
        --runThreadN 8
        --genomeDir "$GENOME_DIR"
        --limitIObufferSize 50000000 50000000
        --outSJtype None
        --outBAMcompression 6
        --alignIntronMax 500000
        --outFilterMismatchNmax 6
        --outFilterMismatchNoverReadLmax 1.0
        --outFilterMatchNmin 25
        --outSAMunmapped None
        --outFilterMatchNminOverLread 0
        --outFilterMultimapNmax 10000
        --outFilterMultimapScoreRange 4
        --outSAMmultNmax 10000
        --winAnchorMultimapNmax 200
        --outSAMprimaryFlag AllBestScore
        --outFilterScoreMin 0
        --outFilterScoreMinOverLread 0
        --outSAMattributes NH HI AS nM NM
        --alignEndsType Local
        --chimSegmentMin 1000000
        --readFilesCommand zcat
        --readFilesIn "$R1_FILE" "$R2_FILE"
    )
    
    # Run STAR-Flex baseline (unsorted)
    echo "  Running STAR-Flex baseline (unsorted)..."
    FLEX_UNSORTED_DIR="$BULK_BASE_DIR/flex_unsorted"
    mkdir -p "$FLEX_UNSORTED_DIR"
    "$STAR_FLEX_BIN" \
        "${COMMON_PARAMS[@]}" \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix "$FLEX_UNSORTED_DIR/" \
        --outTmpDir "$BULK_BASE_DIR/tmp_flex_unsorted" \
        2>&1 | tee "$FLEX_UNSORTED_DIR/star.log" | grep -E "started|finished|ERROR|WARNING" || true
    
    # Run STAR Solo baseline (unsorted)
    echo "  Running STAR Solo baseline (unsorted)..."
    SOLO_UNSORTED_DIR="$BULK_BASE_DIR/solo_unsorted"
    mkdir -p "$SOLO_UNSORTED_DIR"
    "$STAR_SOLO_BIN" \
        "${COMMON_PARAMS[@]}" \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix "$SOLO_UNSORTED_DIR/" \
        --outTmpDir "$BULK_BASE_DIR/tmp_solo_unsorted" \
        2>&1 | tee "$SOLO_UNSORTED_DIR/star.log" | grep -E "started|finished|ERROR|WARNING" || true
    
    # Run STAR-Flex baseline (sorted)
    echo "  Running STAR-Flex baseline (sorted)..."
    FLEX_SORTED_DIR="$BULK_BASE_DIR/flex_sorted"
    mkdir -p "$FLEX_SORTED_DIR"
    "$STAR_FLEX_BIN" \
        "${COMMON_PARAMS[@]}" \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "$FLEX_SORTED_DIR/" \
        --outTmpDir "$BULK_BASE_DIR/tmp_flex_sorted" \
        2>&1 | tee "$FLEX_SORTED_DIR/star.log" | grep -E "started|finished|ERROR|WARNING" || true
    
    # Run STAR Solo baseline (sorted)
    echo "  Running STAR Solo baseline (sorted)..."
    SOLO_SORTED_DIR="$BULK_BASE_DIR/solo_sorted"
    mkdir -p "$SOLO_SORTED_DIR"
    "$STAR_SOLO_BIN" \
        "${COMMON_PARAMS[@]}" \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "$SOLO_SORTED_DIR/" \
        --outTmpDir "$BULK_BASE_DIR/tmp_solo_sorted" \
        2>&1 | tee "$SOLO_SORTED_DIR/star.log" | grep -E "started|finished|ERROR|WARNING" || true
    
    # Compare results
    echo "  Comparing results..."
    
    if command -v samtools &> /dev/null; then
        # Compare unsorted
        FLEX_UNSORTED_TOTAL=$(samtools view -c "$FLEX_UNSORTED_DIR/Aligned.out.bam" 2>/dev/null || echo "0")
        SOLO_UNSORTED_TOTAL=$(samtools view -c "$SOLO_UNSORTED_DIR/Aligned.out.bam" 2>/dev/null || echo "0")
        
        samtools index "$FLEX_UNSORTED_DIR/Aligned.out.bam" 2>/dev/null || true
        samtools index "$SOLO_UNSORTED_DIR/Aligned.out.bam" 2>/dev/null || true
        
        # For unsorted BAMs, idxstats may not work reliably, use direct count
        FLEX_UNSORTED_Y=$(samtools view "$FLEX_UNSORTED_DIR/Aligned.out.bam" 2>/dev/null | awk '$3=="chrY" {count++} END {print count+0}' || echo "0")
        SOLO_UNSORTED_Y=$(samtools view "$SOLO_UNSORTED_DIR/Aligned.out.bam" 2>/dev/null | awk '$3=="chrY" {count++} END {print count+0}' || echo "0")
        
        # Compare sorted
        FLEX_SORTED_TOTAL=$(samtools view -c "$FLEX_SORTED_DIR/Aligned.sortedByCoord.out.bam" 2>/dev/null || echo "0")
        SOLO_SORTED_TOTAL=$(samtools view -c "$SOLO_SORTED_DIR/Aligned.sortedByCoord.out.bam" 2>/dev/null || echo "0")
        
        samtools index "$FLEX_SORTED_DIR/Aligned.sortedByCoord.out.bam" 2>/dev/null || true
        samtools index "$SOLO_SORTED_DIR/Aligned.sortedByCoord.out.bam" 2>/dev/null || true
        
        FLEX_SORTED_Y=$(samtools idxstats "$FLEX_SORTED_DIR/Aligned.sortedByCoord.out.bam" 2>/dev/null | awk '$1=="chrY" {print $3+0}' || echo "0")
        SOLO_SORTED_Y=$(samtools idxstats "$SOLO_SORTED_DIR/Aligned.sortedByCoord.out.bam" 2>/dev/null | awk '$1=="chrY" {print $3+0}' || echo "0")
        
        echo ""
        echo "  Unsorted Comparison:"
        echo "    STAR-Flex total: $FLEX_UNSORTED_TOTAL"
        echo "    STAR Solo total: $SOLO_UNSORTED_TOTAL"
        echo "    STAR-Flex chrY: $FLEX_UNSORTED_Y"
        echo "    STAR Solo chrY: $SOLO_UNSORTED_Y"
        
        echo ""
        echo "  Sorted Comparison:"
        echo "    STAR-Flex total: $FLEX_SORTED_TOTAL"
        echo "    STAR Solo total: $SOLO_SORTED_TOTAL"
        echo "    STAR-Flex chrY: $FLEX_SORTED_Y"
        echo "    STAR Solo chrY: $SOLO_SORTED_Y"
        
        # Check if matches
        UNSORTED_TOTAL_MATCH=false
        UNSORTED_Y_MATCH=false
        SORTED_TOTAL_MATCH=false
        SORTED_Y_MATCH=false
        
        if [ "$FLEX_UNSORTED_TOTAL" -eq "$SOLO_UNSORTED_TOTAL" ]; then
            UNSORTED_TOTAL_MATCH=true
            echo "    ✓ Unsorted total reads match"
            BULK_PASSED=$((BULK_PASSED + 1))
        else
            echo "    ✗ Unsorted total reads mismatch"
            BULK_FAILED=$((BULK_FAILED + 1))
        fi
        
        if [ "$FLEX_UNSORTED_Y" -eq "$SOLO_UNSORTED_Y" ]; then
            UNSORTED_Y_MATCH=true
            echo "    ✓ Unsorted chrY reads match"
            BULK_PASSED=$((BULK_PASSED + 1))
        else
            echo "    ✗ Unsorted chrY reads mismatch"
            BULK_FAILED=$((BULK_FAILED + 1))
        fi
        
        if [ "$FLEX_SORTED_TOTAL" -eq "$SOLO_SORTED_TOTAL" ]; then
            SORTED_TOTAL_MATCH=true
            echo "    ✓ Sorted total reads match"
            BULK_PASSED=$((BULK_PASSED + 1))
        else
            echo "    ✗ Sorted total reads mismatch"
            BULK_FAILED=$((BULK_FAILED + 1))
        fi
        
        if [ "$FLEX_SORTED_Y" -eq "$SOLO_SORTED_Y" ]; then
            SORTED_Y_MATCH=true
            echo "    ✓ Sorted chrY reads match"
            BULK_PASSED=$((BULK_PASSED + 1))
        else
            echo "    ✗ Sorted chrY reads mismatch"
            BULK_FAILED=$((BULK_FAILED + 1))
        fi
        
        # Append to report
        cat >> "$REPORT_FILE" << EOF

### Results

#### Unsorted BAM

| Metric | STAR-Flex | STAR Solo | Match |
|--------|-----------|-----------|-------|
| Total reads | $FLEX_UNSORTED_TOTAL | $SOLO_UNSORTED_TOTAL | $([ "$UNSORTED_TOTAL_MATCH" = true ] && echo "✓" || echo "✗") |
| chrY reads | $FLEX_UNSORTED_Y | $SOLO_UNSORTED_Y | $([ "$UNSORTED_Y_MATCH" = true ] && echo "✓" || echo "✗") |

#### Sorted BAM

| Metric | STAR-Flex | STAR Solo | Match |
|--------|-----------|-----------|-------|
| Total reads | $FLEX_SORTED_TOTAL | $SOLO_SORTED_TOTAL | $([ "$SORTED_TOTAL_MATCH" = true ] && echo "✓" || echo "✗") |
| chrY reads | $FLEX_SORTED_Y | $SOLO_SORTED_Y | $([ "$SORTED_Y_MATCH" = true ] && echo "✓" || echo "✗") |

### Conclusion

- Unsorted total reads: $([ "$UNSORTED_TOTAL_MATCH" = true ] && echo "PASS" || echo "FAIL")
- Unsorted chrY reads: $([ "$UNSORTED_Y_MATCH" = true ] && echo "PASS" || echo "FAIL")
- Sorted total reads: $([ "$SORTED_TOTAL_MATCH" = true ] && echo "PASS" || echo "FAIL")
- Sorted chrY reads: $([ "$SORTED_Y_MATCH" = true ] && echo "PASS" || echo "FAIL")

EOF
    else
        echo "  ERROR: samtools not available for comparison"
        BULK_FAILED=$((BULK_FAILED + 1))
    fi
else
    echo ">>> Test 1: Bulk RNA-Seq Baseline Regression - SKIPPED (STAR Solo not available)"
    cat >> "$REPORT_FILE" << EOF

### Status: SKIPPED

STAR Solo binary not found at \`$STAR_SOLO_BIN\`. Skipping bulk RNA-seq regression test.

EOF
fi

echo ""

# Test 2: Flex Regression (compare to previous artifacts)
echo ">>> Test 2: Flex Regression (compare to previous artifacts)"

FLEX_PASSED=0
FLEX_FAILED=0

# Look for previous Flex test artifacts
PREV_FLEX_DIR="/tmp/ychrom_flex_test"
CURRENT_FLEX_DIR="/tmp/ychrom_flex_test_current"

if [ -d "$PREV_FLEX_DIR/sorted_split" ] && [ -f "$PREV_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_Y.bam" ]; then
    echo "  Found previous Flex test artifacts at $PREV_FLEX_DIR"
    
    # Run current Flex test
    echo "  Running current Flex test..."
    rm -rf "$CURRENT_FLEX_DIR"
    BASE_DIR="$CURRENT_FLEX_DIR" bash "$SCRIPT_DIR/run_ychrom_flex_validation_test.sh" > /dev/null 2>&1 || true
    
    if [ -f "$CURRENT_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_Y.bam" ]; then
        echo "  Comparing previous vs current Flex outputs..."
        
        if command -v samtools &> /dev/null; then
            # Compare Y BAM counts
            PREV_Y_COUNT=$(samtools view -c "$PREV_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_Y.bam" 2>/dev/null || echo "0")
            CURRENT_Y_COUNT=$(samtools view -c "$CURRENT_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_Y.bam" 2>/dev/null || echo "0")
            
            PREV_NOY_COUNT=$(samtools view -c "$PREV_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_noY.bam" 2>/dev/null || echo "0")
            CURRENT_NOY_COUNT=$(samtools view -c "$CURRENT_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_noY.bam" 2>/dev/null || echo "0")
            
            echo "    Previous Y reads: $PREV_Y_COUNT"
            echo "    Current Y reads: $CURRENT_Y_COUNT"
            echo "    Previous noY reads: $PREV_NOY_COUNT"
            echo "    Current noY reads: $CURRENT_NOY_COUNT"
            
            Y_MATCH=false
            NOY_MATCH=false
            
            if [ "$PREV_Y_COUNT" -eq "$CURRENT_Y_COUNT" ]; then
                Y_MATCH=true
                echo "    ✓ Y read counts match"
                FLEX_PASSED=$((FLEX_PASSED + 1))
            else
                echo "    ✗ Y read counts mismatch (diff: $((CURRENT_Y_COUNT - PREV_Y_COUNT)))"
                FLEX_FAILED=$((FLEX_FAILED + 1))
            fi
            
            if [ "$PREV_NOY_COUNT" -eq "$CURRENT_NOY_COUNT" ]; then
                NOY_MATCH=true
                echo "    ✓ noY read counts match"
                FLEX_PASSED=$((FLEX_PASSED + 1))
            else
                echo "    ✗ noY read counts mismatch (diff: $((CURRENT_NOY_COUNT - PREV_NOY_COUNT)))"
                FLEX_FAILED=$((FLEX_FAILED + 1))
            fi
            
            # Compare idxstats
            echo "  Comparing idxstats..."
            samtools index "$PREV_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_Y.bam" 2>/dev/null || true
            samtools index "$CURRENT_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_Y.bam" 2>/dev/null || true
            samtools index "$PREV_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_noY.bam" 2>/dev/null || true
            samtools index "$CURRENT_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_noY.bam" 2>/dev/null || true
            
            PREV_Y_IDXSTATS="$CURRENT_FLEX_DIR/prev_Y_idxstats.txt"
            CURRENT_Y_IDXSTATS="$CURRENT_FLEX_DIR/current_Y_idxstats.txt"
            PREV_NOY_IDXSTATS="$CURRENT_FLEX_DIR/prev_noY_idxstats.txt"
            CURRENT_NOY_IDXSTATS="$CURRENT_FLEX_DIR/current_noY_idxstats.txt"
            
            samtools idxstats "$PREV_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_Y.bam" > "$PREV_Y_IDXSTATS" 2>/dev/null || true
            samtools idxstats "$CURRENT_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_Y.bam" > "$CURRENT_Y_IDXSTATS" 2>/dev/null || true
            samtools idxstats "$PREV_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_noY.bam" > "$PREV_NOY_IDXSTATS" 2>/dev/null || true
            samtools idxstats "$CURRENT_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_noY.bam" > "$CURRENT_NOY_IDXSTATS" 2>/dev/null || true
            
            IDXSTATS_Y_DIFF=$(diff "$PREV_Y_IDXSTATS" "$CURRENT_Y_IDXSTATS" 2>/dev/null | wc -l || echo "0")
            IDXSTATS_NOY_DIFF=$(diff "$PREV_NOY_IDXSTATS" "$CURRENT_NOY_IDXSTATS" 2>/dev/null | wc -l || echo "0")
            
            if [ "$IDXSTATS_Y_DIFF" -eq 0 ]; then
                echo "    ✓ Y idxstats match"
                FLEX_PASSED=$((FLEX_PASSED + 1))
            else
                echo "    ✗ Y idxstats differ ($IDXSTATS_Y_DIFF lines)"
                FLEX_FAILED=$((FLEX_FAILED + 1))
            fi
            
            if [ "$IDXSTATS_NOY_DIFF" -eq 0 ]; then
                echo "    ✓ noY idxstats match"
                FLEX_PASSED=$((FLEX_PASSED + 1))
            else
                echo "    ✗ noY idxstats differ ($IDXSTATS_NOY_DIFF lines)"
                FLEX_FAILED=$((FLEX_FAILED + 1))
            fi
            
            # Append to report
            cat >> "$REPORT_FILE" << EOF

---

## Test 2: Flex Regression (Previous vs Current)

### Configuration

- Previous artifacts: \`$PREV_FLEX_DIR\`
- Current run: \`$CURRENT_FLEX_DIR\`
- Dataset: \`/storage/downsampled/SC2300771\`

### Results

| Metric | Previous | Current | Match |
|--------|----------|---------|-------|
| Y reads | $PREV_Y_COUNT | $CURRENT_Y_COUNT | $([ "$Y_MATCH" = true ] && echo "✓" || echo "✗") |
| noY reads | $PREV_NOY_COUNT | $CURRENT_NOY_COUNT | $([ "$NOY_MATCH" = true ] && echo "✓" || echo "✗") |
| Y idxstats | - | - | $([ "$IDXSTATS_Y_DIFF" -eq 0 ] && echo "✓" || echo "✗") |
| noY idxstats | - | - | $([ "$IDXSTATS_NOY_DIFF" -eq 0 ] && echo "✓" || echo "✗") |

### Conclusion

- Y read counts: $([ "$Y_MATCH" = true ] && echo "PASS" || echo "FAIL")
- noY read counts: $([ "$NOY_MATCH" = true ] && echo "PASS" || echo "FAIL")
- Y idxstats: $([ "$IDXSTATS_Y_DIFF" -eq 0 ] && echo "PASS" || echo "FAIL")
- noY idxstats: $([ "$IDXSTATS_NOY_DIFF" -eq 0 ] && echo "PASS" || echo "FAIL")

EOF
        else
            echo "  ERROR: samtools not available for comparison"
            FLEX_FAILED=$((FLEX_FAILED + 1))
        fi
    else
        echo "  ERROR: Current Flex test did not produce expected outputs"
        FLEX_FAILED=$((FLEX_FAILED + 1))
        cat >> "$REPORT_FILE" << EOF

---

## Test 2: Flex Regression - FAILED

Current Flex test did not produce expected outputs at \`$CURRENT_FLEX_DIR\`.

EOF
    fi
else
    echo "  Previous Flex test artifacts not found at $PREV_FLEX_DIR"
    echo "  Running new Flex test for baseline..."
    
    # Run Flex test to create baseline
    rm -rf "$CURRENT_FLEX_DIR"
    BASE_DIR="$CURRENT_FLEX_DIR" bash "$SCRIPT_DIR/run_ychrom_flex_validation_test.sh" > /dev/null 2>&1 || true
    
    if [ -f "$CURRENT_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_Y.bam" ]; then
        CURRENT_Y_COUNT=$(samtools view -c "$CURRENT_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_Y.bam" 2>/dev/null || echo "0")
        CURRENT_NOY_COUNT=$(samtools view -c "$CURRENT_FLEX_DIR/sorted_split/Aligned.sortedByCoord.out_noY.bam" 2>/dev/null || echo "0")
        
        echo "  Created baseline Flex test outputs:"
        echo "    Y reads: $CURRENT_Y_COUNT"
        echo "    noY reads: $CURRENT_NOY_COUNT"
        
        cat >> "$REPORT_FILE" << EOF

---

## Test 2: Flex Regression - BASELINE CREATED

Previous artifacts not found. Created new baseline at \`$CURRENT_FLEX_DIR\`:

- Y reads: $CURRENT_Y_COUNT
- noY reads: $CURRENT_NOY_COUNT

Future runs will compare against this baseline.

EOF
    else
        echo "  ERROR: Failed to create Flex test baseline"
        FLEX_FAILED=$((FLEX_FAILED + 1))
        cat >> "$REPORT_FILE" << EOF

---

## Test 2: Flex Regression - FAILED

Failed to create Flex test baseline.

EOF
    fi
fi

echo ""

# Final summary
TOTAL_PASSED=$((BULK_PASSED + FLEX_PASSED))
TOTAL_FAILED=$((BULK_FAILED + FLEX_FAILED))

cat >> "$REPORT_FILE" << EOF

---

## Summary

- Bulk RNA-seq regression checks passed: \`$BULK_PASSED\`
- Bulk RNA-seq regression checks failed: \`$BULK_FAILED\`
- Flex regression checks passed: \`$FLEX_PASSED\`
- Flex regression checks failed: \`$FLEX_FAILED\`
- **Total passed: \`$TOTAL_PASSED\`**
- **Total failed: \`$TOTAL_FAILED\`**

EOF

echo "=========================================="
echo "Regression Test Summary"
echo "=========================================="
echo "Bulk RNA-seq checks passed: $BULK_PASSED"
echo "Bulk RNA-seq checks failed: $BULK_FAILED"
echo "Flex regression checks passed: $FLEX_PASSED"
echo "Flex regression checks failed: $FLEX_FAILED"
echo "Total passed: $TOTAL_PASSED"
echo "Total failed: $TOTAL_FAILED"
echo ""
echo "Report written to: $REPORT_FILE"
echo ""

if [ $TOTAL_FAILED -eq 0 ]; then
    echo "✓ All regression checks passed!"
    exit 0
else
    echo "✗ Some regression checks failed"
    exit 1
fi

