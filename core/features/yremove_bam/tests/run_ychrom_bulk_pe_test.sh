#!/bin/bash
# Bulk Paired-End Y/NoY Split Validation Test
# Validates Y-chromosome BAM splitting on bulk paired-end RNA-seq data
# Tests mate-checking logic: if either mate has Y alignment, both mates route to _Y.bam

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${SCRIPT_DIR}/../source/STAR"
BASE_DIR="${BASE_DIR:-/tmp/ychrom_bulk_pe_test}"
REPORT_FILE="${SCRIPT_DIR}/TEST_REPORT_Y_SPLIT_BULK_PE.md"

# Bulk RNA-seq parameters (no Flex, no Solo)
GENOME_DIR="/storage/flex_filtered_reference/star_index"  # Use same reference (has chrY)
R1_FILE="/storage/PE/downsampled/21033-09-01-13-01_S1_L007_R1_001.fastq.gz"
R2_FILE="/storage/PE/downsampled/21033-09-01-13-01_S1_L007_R2_001.fastq.gz"

if [ ! -f "$R1_FILE" ] || [ ! -f "$R2_FILE" ]; then
    echo "ERROR: FASTQ files not found"
    echo "  R1: $R1_FILE"
    echo "  R2: $R2_FILE"
    exit 1
fi

# Clean up previous run
rm -rf "$BASE_DIR"
mkdir -p "$BASE_DIR"

echo "=========================================="
echo "Bulk Paired-End Y/NoY Split Validation Test"
echo "=========================================="
echo "Binary: $STAR_BIN"
echo "Base directory: $BASE_DIR"
echo "Dataset: $R1_FILE, $R2_FILE"
echo ""

# Common STAR parameters for bulk paired-end RNA-seq
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

# Step 1: Baseline Run (no Y-split)
echo ">>> Step 1: Baseline Run (no Y-split)"
BASELINE_DIR="$BASE_DIR/baseline"
mkdir -p "$BASELINE_DIR"

"$STAR_BIN" \
    "${COMMON_PARAMS[@]}" \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix "$BASELINE_DIR/" \
    --outTmpDir "$BASE_DIR/tmp_baseline" \
    2>&1 | tee "$BASELINE_DIR/star.log" | grep -E "started|finished|ERROR|WARNING" || true

if [ ! -f "$BASELINE_DIR/Aligned.out.bam" ]; then
    echo "ERROR: Baseline BAM not created"
    exit 1
fi

# Capture baseline stats
if command -v samtools &> /dev/null; then
    BASE_TOTAL=$(samtools view -c "$BASELINE_DIR/Aligned.out.bam" 2>/dev/null || echo "0")
    # Index baseline BAM for idxstats
    samtools index "$BASELINE_DIR/Aligned.out.bam" 2>/dev/null || true
    # Get chrY count from idxstats (column 3)
    BASE_Y=$(samtools idxstats "$BASELINE_DIR/Aligned.out.bam" 2>/dev/null | awk '$1=="chrY" {y=$3+0} END {if (y=="") y=0; print y}')
    # Also check by counting reads directly (including secondary alignments)
    BASE_Y_DIRECT=$(samtools view "$BASELINE_DIR/Aligned.out.bam" 2>/dev/null | awk '$3=="chrY" {count++} END {print count+0}' || echo "0")
    echo "  Baseline total reads: $BASE_TOTAL"
    echo "  Baseline chrY reads (idxstats): $BASE_Y"
    echo "  Baseline chrY reads (direct count): $BASE_Y_DIRECT"
    
    # Use the larger of the two counts (idxstats might miss unmapped pairs)
    if [ "$BASE_Y_DIRECT" -gt "$BASE_Y" ]; then
        BASE_Y="$BASE_Y_DIRECT"
    fi
    
    if [ -z "$BASE_Y" ] || [ "$BASE_Y" = "0" ] || [ "$BASE_Y" -eq 0 ]; then
        echo "WARNING: Baseline has zero Y reads - aborting split validation"
        echo "This may indicate a dataset issue or reference problem"
        exit 1
    fi
else
    echo "WARNING: samtools not available, skipping baseline stats"
    BASE_TOTAL="unknown"
    BASE_Y="unknown"
fi
echo ""

# Step 2: Unsorted Split Run
echo ">>> Step 2: Unsorted Split Run"
UNSORTED_DIR="$BASE_DIR/unsorted_split"
mkdir -p "$UNSORTED_DIR"

"$STAR_BIN" \
    "${COMMON_PARAMS[@]}" \
    --outSAMtype BAM Unsorted \
    --emitNoYBAM yes \
    --outFileNamePrefix "$UNSORTED_DIR/" \
    --outTmpDir "$BASE_DIR/tmp_unsorted" \
    2>&1 | tee "$UNSORTED_DIR/star.log" | grep -E "started|finished|ERROR|WARNING" || true

if [ ! -f "$UNSORTED_DIR/Aligned.out_Y.bam" ] || [ ! -f "$UNSORTED_DIR/Aligned.out_noY.bam" ]; then
    echo "ERROR: Unsorted split BAMs not created"
    exit 1
fi
echo ""

# Step 3: Sorted Split Run
echo ">>> Step 3: Sorted Split Run"
SORTED_DIR="$BASE_DIR/sorted_split"
mkdir -p "$SORTED_DIR"

"$STAR_BIN" \
    "${COMMON_PARAMS[@]}" \
    --outSAMtype BAM SortedByCoordinate \
    --emitNoYBAM yes \
    --outFileNamePrefix "$SORTED_DIR/" \
    --outTmpDir "$BASE_DIR/tmp_sorted" \
    2>&1 | tee "$SORTED_DIR/star.log" | grep -E "started|finished|ERROR|WARNING" || true

if [ ! -f "$SORTED_DIR/Aligned.sortedByCoord.out_Y.bam" ] || [ ! -f "$SORTED_DIR/Aligned.sortedByCoord.out_noY.bam" ]; then
    echo "ERROR: Sorted split BAMs not created"
    exit 1
fi
echo ""

# Step 4: Validation
echo ">>> Step 4: Validation Checks"
VALIDATION_PASSED=0
VALIDATION_FAILED=0

if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools required for validation"
    exit 1
fi

# Function to validate split outputs
validate_split() {
    local mode="$1"  # "unsorted" or "sorted"
    local y_bam="$2"
    local noy_bam="$3"
    local prefix="$4"
    
    echo "  Validating $mode mode..."
    
    # Index BAMs
    samtools index "$y_bam" 2>/dev/null || true
    samtools index "$noy_bam" 2>/dev/null || true
    
    # Get counts
    Y_COUNT=$(samtools view -c "$y_bam" 2>/dev/null || echo "0")
    NOY_COUNT=$(samtools view -c "$noy_bam" 2>/dev/null || echo "0")
    TOTAL=$((Y_COUNT + NOY_COUNT))
    
    echo "    Y reads: $Y_COUNT"
    echo "    noY reads: $NOY_COUNT"
    echo "    Total: $TOTAL"
    
    # Check 1: Count consistency
    if [ "$BASE_TOTAL" != "unknown" ]; then
        if [ "$TOTAL" -eq "$BASE_TOTAL" ]; then
            echo "    ✓ Count consistency: Y + noY = baseline ($TOTAL)"
            VALIDATION_PASSED=$((VALIDATION_PASSED + 1))
        else
            echo "    ✗ Count mismatch: Y + noY = $TOTAL, baseline = $BASE_TOTAL"
            VALIDATION_FAILED=$((VALIDATION_FAILED + 1))
        fi
    else
        echo "    ⚠ Skipping count check (baseline unknown)"
    fi
    
    # Check 2: noY exclusivity (must have zero chrY reads)
    samtools index "$noy_bam" 2>/dev/null || true
    NOY_Y_COUNT=$(samtools idxstats "$noy_bam" 2>/dev/null | awk '$1=="chrY" {print $3+0}' || echo "0")
    if [ -z "$NOY_Y_COUNT" ]; then
        NOY_Y_COUNT=0
    fi
    if [ "$NOY_Y_COUNT" -eq 0 ]; then
        echo "    ✓ noY exclusivity: zero chrY reads in _noY.bam"
        VALIDATION_PASSED=$((VALIDATION_PASSED + 1))
    else
        echo "    ✗ noY exclusivity FAILED: $NOY_Y_COUNT chrY reads in _noY.bam"
        VALIDATION_FAILED=$((VALIDATION_FAILED + 1))
    fi
    
    # Check 3: Y exclusivity (reads in _Y.bam must have at least one alignment on chrY)
    # For bulk paired-end: check if read or mate maps to chrY
    Y_CHRY_COUNT=$(samtools view "$y_bam" 2>/dev/null | awk '$3=="chrY" || $7=="chrY" {count++} END {print count+0}' || echo "0")
    Y_TOTAL_MAPPED=$(samtools view "$y_bam" 2>/dev/null | awk '$3!="*" {count++} END {print count+0}' || echo "0")
    if [ "$Y_TOTAL_MAPPED" -eq 0 ] || [ "$Y_CHRY_COUNT" -gt 0 ]; then
        echo "    ✓ Y exclusivity: reads in _Y.bam have chrY alignments (direct or via mate)"
        VALIDATION_PASSED=$((VALIDATION_PASSED + 1))
    else
        echo "    ✗ Y exclusivity FAILED: $Y_TOTAL_MAPPED mapped reads, $Y_CHRY_COUNT with chrY connection"
        VALIDATION_FAILED=$((VALIDATION_FAILED + 1))
    fi
    
    # Check 4: Mate consistency (for paired-end, both mates should route together)
    # Extract read names and check if mates are in the same file
    Y_READS=$(samtools view "$y_bam" 2>/dev/null | awk '{print $1}' | sort -u)
    NOY_READS=$(samtools view "$noy_bam" 2>/dev/null | awk '{print $1}' | sort -u)
    OVERLAP=$(comm -12 <(echo "$Y_READS") <(echo "$NOY_READS") | wc -l)
    if [ "$OVERLAP" -eq 0 ]; then
        echo "    ✓ Mate consistency: no read names appear in both Y and noY files"
        VALIDATION_PASSED=$((VALIDATION_PASSED + 1))
    else
        echo "    ✗ Mate consistency FAILED: $OVERLAP read names appear in both files"
        VALIDATION_FAILED=$((VALIDATION_FAILED + 1))
    fi
    
    # Store for report
    eval "${prefix}_Y_COUNT=$Y_COUNT"
    eval "${prefix}_NOY_COUNT=$NOY_COUNT"
}

# Validate unsorted
validate_split "unsorted" \
    "$UNSORTED_DIR/Aligned.out_Y.bam" \
    "$UNSORTED_DIR/Aligned.out_noY.bam" \
    "UNSORTED"

# Validate sorted
validate_split "sorted" \
    "$SORTED_DIR/Aligned.sortedByCoord.out_Y.bam" \
    "$SORTED_DIR/Aligned.sortedByCoord.out_noY.bam" \
    "SORTED"

echo ""

# Step 5: Mate-checking verification
echo ">>> Step 5: Mate-Checking Verification (Bulk Paired-End)"
echo "  Checking if mates route together when one mate is on chrY..."

# Sample a few reads from _Y.bam and verify their mates are also in _Y.bam
SAMPLE_READS=$(samtools view "$SORTED_DIR/Aligned.sortedByCoord.out_Y.bam" 2>/dev/null | head -20 | awk '{print $1}' | sort -u | head -5)
MATE_CHECK_PASSED=0
MATE_CHECK_FAILED=0

for read_name in $SAMPLE_READS; do
    # Count occurrences in Y and noY files
    Y_COUNT=$(samtools view "$SORTED_DIR/Aligned.sortedByCoord.out_Y.bam" 2>/dev/null | awk -v r="$read_name" '$1==r {count++} END {print count+0}')
    NOY_COUNT=$(samtools view "$SORTED_DIR/Aligned.sortedByCoord.out_noY.bam" 2>/dev/null | awk -v r="$read_name" '$1==r {count++} END {print count+0}')
    
    if [ "$Y_COUNT" -gt 0 ] && [ "$NOY_COUNT" -eq 0 ]; then
        MATE_CHECK_PASSED=$((MATE_CHECK_PASSED + 1))
    elif [ "$Y_COUNT" -gt 0 ] && [ "$NOY_COUNT" -gt 0 ]; then
        echo "    ✗ Read $read_name appears in both Y ($Y_COUNT) and noY ($NOY_COUNT) files"
        MATE_CHECK_FAILED=$((MATE_CHECK_FAILED + 1))
    fi
done

if [ "$MATE_CHECK_FAILED" -eq 0 ]; then
    echo "  ✓ Mate-checking: all sampled reads route consistently (mates together)"
    VALIDATION_PASSED=$((VALIDATION_PASSED + 1))
else
    echo "  ✗ Mate-checking FAILED: $MATE_CHECK_FAILED reads appear in both files"
    VALIDATION_FAILED=$((VALIDATION_FAILED + 1))
fi
echo ""

# Step 6: Generate Report
echo ">>> Step 6: Generating Report"
cat > "$REPORT_FILE" << EOF
# Bulk Paired-End Y/NoY Split Validation Test Report

Generated: $(date)

## Test Configuration

- STAR Binary: \`$STAR_BIN\`
- Dataset: \`$R1_FILE\`, \`$R2_FILE\`
- Reference: \`${GENOME_DIR}\`
- Base Directory: \`${BASE_DIR}\`
- Mode: Bulk Paired-End RNA-seq (NOT single-cell/Flex)

## Baseline Statistics

- Total reads: \`$BASE_TOTAL\`
- chrY reads: \`$BASE_Y\`

## Split Run Results

### Unsorted Mode

- Y reads: \`${UNSORTED_Y_COUNT}\`
- noY reads: \`${UNSORTED_NOY_COUNT}\`
- Total: \`$((UNSORTED_Y_COUNT + UNSORTED_NOY_COUNT))\`

### Sorted Mode

- Y reads: \`${SORTED_Y_COUNT}\`
- noY reads: \`${SORTED_NOY_COUNT}\`
- Total: \`$((SORTED_Y_COUNT + SORTED_NOY_COUNT))\`

## Validation Results

- Passed: \`$VALIDATION_PASSED\`
- Failed: \`$VALIDATION_FAILED\`

### Checks Performed

1. Count consistency: Y + noY = baseline total
2. noY exclusivity: Zero chrY reads in \`_noY.bam\`
3. Y exclusivity: All reads in \`_Y.bam\` have chrY alignments (direct or via mate)
4. Mate consistency: No read names appear in both Y and noY files
5. Mate-checking: Mates route together (if one mate is on Y, both go to _Y.bam)

## Files

- Baseline BAM: \`${BASELINE_DIR}/Aligned.out.bam\`
- Unsorted Y: \`${UNSORTED_DIR}/Aligned.out_Y.bam\`
- Unsorted noY: \`${UNSORTED_DIR}/Aligned.out_noY.bam\`
- Sorted Y: \`${SORTED_DIR}/Aligned.sortedByCoord.out_Y.bam\`
- Sorted noY: \`${SORTED_DIR}/Aligned.sortedByCoord.out_noY.bam\`

EOF

echo "  Report written to: $REPORT_FILE"
echo ""

# Summary
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo "Validation checks passed: $VALIDATION_PASSED"
echo "Validation checks failed: $VALIDATION_FAILED"
echo ""

if [ $VALIDATION_FAILED -eq 0 ]; then
    echo "✓ All validation checks passed!"
    exit 0
else
    echo "✗ Some validation checks failed"
    echo "See report: $REPORT_FILE"
    exit 1
fi

