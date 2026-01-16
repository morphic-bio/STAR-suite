#!/bin/bash
# Comprehensive test for remove_y_reads FASTQ Splitter
# Tests single-threaded and multithreaded modes with validation

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOL="$SCRIPT_DIR/../tools/remove_y_reads/remove_y_reads"
REPORT_FILE="$SCRIPT_DIR/TEST_REPORT_REMOVE_Y_FASTQ.md"
TEST_DIR="/tmp/y_removal_comprehensive_test_$$"

mkdir -p "$TEST_DIR"
trap "rm -rf $TEST_DIR" EXIT

echo "=========================================="
echo "Comprehensive remove_y_reads Test"
echo "=========================================="
echo "Tool: $TOOL"
echo "Test directory: $TEST_DIR"
echo ""

# Check if tool exists
if [ ! -f "$TOOL" ]; then
    echo "ERROR: Tool not found at $TOOL"
    echo "Please build it first: cd tools/remove_y_reads && make"
    exit 1
fi

# Check if samtools is available
if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found. Required for test data creation."
    exit 1
fi

# Initialize report
cat > "$REPORT_FILE" << EOF
# remove_y_reads FASTQ Splitter Test Report

Generated: $(date)

## Test Configuration

- Tool: \`$TOOL\`
- Test Directory: \`$TEST_DIR\`

EOF

PASSED=0
FAILED=0

# Helper function to normalize read name (matches tool's normalization exactly)
normalize_qname() {
    local name="$1"
    local raw_len=${#name}
    local start=0
    local len=$raw_len
    
    # Skip leading '@' if present (FASTQ name lines)
    if [ $len -gt 0 ] && [ "${name:0:1}" = "@" ]; then
        start=1
        len=$((len - 1))
    fi
    
    # Find end: stop at whitespace, tab, or newline
    local end=0
    while [ $end -lt $len ]; do
        local char="${name:$((start + end)):1}"
        if [ "$char" = " " ] || [ "$char" = $'\t' ] || [ "$char" = $'\n' ] || [ "$char" = $'\r' ]; then
            break
        fi
        end=$((end + 1))
    done
    len=$end
    
    # Strip trailing /1 or /2 (mate suffix)
    if [ $len -ge 2 ]; then
        local last_two="${name:$((start + len - 2)):2}"
        if [ "$last_two" = "/1" ] || [ "$last_two" = "/2" ]; then
            len=$((len - 2))
        fi
    fi
    
    # Extract normalized name
    echo "${name:$start:$len}"
}

# Helper function to extract read names from FASTQ
extract_read_names() {
    local fastq="$1"
    if [[ "$fastq" == *.gz ]]; then
        gunzip -c "$fastq" | grep "^@" | while read line; do
            normalize_qname "$line"
        done
    else
        grep "^@" "$fastq" | while read line; do
            normalize_qname "$line"
        done
    fi
}

# Helper function to extract read names from BAM
extract_bam_qnames() {
    local bam="$1"
    samtools view "$bam" | awk '{print $1}' | sort -u | while read name; do
        normalize_qname "$name"
    done
}

echo "=== Step 1: Creating Test Data ==="

# Create synthetic BAM with Y reads
cat > "$TEST_DIR/y_reads.sam" << 'EOF'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:100
@SQ	SN:chrY	LN:100
read1	0	chrY	1	60	20M	*	0	0	ATCGATCGATCGATCGATCG	IIIIIIIIIIIIIIIIIIII
read3	0	chrY	10	60	20M	*	0	0	ATCGATCGATCGATCGATCG	IIIIIIIIIIIIIIIIIIII
read7	0	chrY	20	60	20M	*	0	0	ATCGATCGATCGATCGATCG	IIIIIIIIIIIIIIIIIIII
EOF

samtools view -bS "$TEST_DIR/y_reads.sam" > "$TEST_DIR/y_reads.bam"

# Create test FASTQ files (multiple files for multithreaded test)
cat > "$TEST_DIR/test1.fastq" << 'EOF'
@read1 comment here
ATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIII
@read2
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
@read3/1
ATCGATCGATCGATCGATCG
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

cat > "$TEST_DIR/test2.fastq" << 'EOF'
@read6
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
@read7
ATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIII
@read8
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
EOF

# Create gzipped versions
gzip -c "$TEST_DIR/test1.fastq" > "$TEST_DIR/test1.fastq.gz"
gzip -c "$TEST_DIR/test2.fastq" > "$TEST_DIR/test2.fastq.gz"

# Extract Y read names from BAM for validation
Y_READS=$(extract_bam_qnames "$TEST_DIR/y_reads.bam" | sort)
Y_COUNT=$(echo "$Y_READS" | grep -c . || echo "0")
echo "Y BAM contains $Y_COUNT unique read names: $(echo $Y_READS | tr '\n' ' ')"

cat >> "$REPORT_FILE" << EOF

## Test Data

- Y BAM: \`$TEST_DIR/y_reads.bam\` ($Y_COUNT unique read names)
- Test FASTQ 1: \`$TEST_DIR/test1.fastq.gz\` (5 reads)
- Test FASTQ 2: \`$TEST_DIR/test2.fastq.gz\` (3 reads)
- Expected Y reads: read1, read3, read7
- Expected noY reads: read2, read4, read5, read6, read8

---

## Test Results

EOF

# Function to run tests and validate
run_test() {
    local mode="$1"
    local threads="$2"
    local outdir="$TEST_DIR/output_${mode}"
    local fastq_args="$3"
    
    echo ""
    echo "=== Testing $mode mode (--threads $threads) ==="
    
    rm -rf "$outdir"
    mkdir -p "$outdir"
    
    # Run tool
    "$TOOL" -y "$TEST_DIR/y_reads.bam" --threads "$threads" -o "$outdir" $fastq_args 2>&1 | tee "$outdir/tool.log"
    
    # Count output files
    local y_files=$(ls "$outdir"/*_Y.fastq.gz 2>/dev/null | wc -l)
    local noy_files=$(ls "$outdir"/*_noY.fastq.gz 2>/dev/null | wc -l)
    
    echo "  Output files: $y_files Y files, $noy_files noY files"
    
    # Validate outputs exist
    local check_passed=true
    if [ "$y_files" -eq 0 ] || [ "$noy_files" -eq 0 ]; then
        echo "  ✗ FAIL: Output files missing"
        check_passed=false
        FAILED=$((FAILED + 1))
    else
        echo "  ✓ Output files exist"
        PASSED=$((PASSED + 1))
    fi
    
    # Validate each FASTQ file
    local total_y=0
    local total_noy=0
    local total_original=0
    
    for fastq in $fastq_args; do
        local basename=$(basename "$fastq" .gz)
        basename="${basename%.fastq}"
        basename="${basename%.fq}"
        
        local y_file="$outdir/${basename}_Y.fastq.gz"
        local noy_file="$outdir/${basename}_noY.fastq.gz"
        
        if [ ! -f "$y_file" ] || [ ! -f "$noy_file" ]; then
            echo "  ✗ FAIL: Missing outputs for $basename"
            check_passed=false
            FAILED=$((FAILED + 1))
            continue
        fi
        
        # Count reads
        local y_count=$(gunzip -c "$y_file" 2>/dev/null | grep -c "^@" || echo "0")
        local noy_count=$(gunzip -c "$noy_file" 2>/dev/null | grep -c "^@" || echo "0")
        local original_count=$(gunzip -c "$fastq" 2>/dev/null | grep -c "^@" || echo "0")
        
        total_y=$((total_y + y_count))
        total_noy=$((total_noy + noy_count))
        total_original=$((total_original + original_count))
        
        echo "  $basename: Y=$y_count, noY=$noy_count, original=$original_count"
        
        # Check count consistency
        if [ $((y_count + noy_count)) -ne "$original_count" ]; then
            echo "  ✗ FAIL: Count mismatch for $basename (Y+noY=$((y_count + noy_count)) != original=$original_count)"
            check_passed=false
            FAILED=$((FAILED + 1))
        else
            echo "  ✓ Count consistency for $basename"
            PASSED=$((PASSED + 1))
        fi
        
        # Extract read names and validate (preserve order for order checks)
        local y_names=$(extract_read_names "$y_file")
        local noy_names=$(extract_read_names "$noy_file")
        local y_names_sorted=$(echo "$y_names" | sort)
        local noy_names_sorted=$(echo "$noy_names" | sort)
        
        # Enhanced order preservation check
        if [ "$basename" = "test1" ]; then
            # Build ordered list of read names from original FASTQ with line numbers
            local original_order_file="$outdir/${basename}_original_order.txt"
            gunzip -c "$fastq" | grep -n "^@" | while IFS=: read line_num line; do
                local norm_name=$(normalize_qname "$line")
                echo "$line_num $norm_name"
            done > "$original_order_file"
            
            # Build ordered list from Y output
            local y_order_file="$outdir/${basename}_y_order.txt"
            gunzip -c "$y_file" | grep -n "^@" | while IFS=: read line_num line; do
                local norm_name=$(normalize_qname "$line")
                echo "$line_num $norm_name"
            done > "$y_order_file"
            
            # Build ordered list from noY output
            local noy_order_file="$outdir/${basename}_noy_order.txt"
            gunzip -c "$noy_file" | grep -n "^@" | while IFS=: read line_num line; do
                local norm_name=$(normalize_qname "$line")
                echo "$line_num $norm_name"
            done > "$noy_order_file"
            
            # Extract line numbers for known noY reads in original (read2, read4, read5)
            local read2_orig_line=$(grep " read2$" "$original_order_file" | awk '{print $1}')
            local read4_orig_line=$(grep " read4$" "$original_order_file" | awk '{print $1}')
            local read5_orig_line=$(grep " read5$" "$original_order_file" | awk '{print $1}')
            
            # Extract line numbers for same reads in noY output
            local read2_noy_line=$(grep " read2$" "$noy_order_file" | awk '{print $1}')
            local read4_noy_line=$(grep " read4$" "$noy_order_file" | awk '{print $1}')
            local read5_noy_line=$(grep " read5$" "$noy_order_file" | awk '{print $1}')
            
            # Verify relative order is preserved: read2 < read4 < read5 in both original and noY
            local order_preserved=true
            if [ -n "$read2_orig_line" ] && [ -n "$read4_orig_line" ] && [ -n "$read5_orig_line" ]; then
                if [ "$read2_orig_line" -gt "$read4_orig_line" ] || [ "$read4_orig_line" -gt "$read5_orig_line" ]; then
                    echo "  ✗ FAIL: Original order invalid (read2=$read2_orig_line, read4=$read4_orig_line, read5=$read5_orig_line)"
                    order_preserved=false
                fi
            fi
            
            if [ -n "$read2_noy_line" ] && [ -n "$read4_noy_line" ] && [ -n "$read5_noy_line" ]; then
                if [ "$read2_noy_line" -gt "$read4_noy_line" ] || [ "$read4_noy_line" -gt "$read5_noy_line" ]; then
                    echo "  ✗ FAIL: Order not preserved in noY output (read2=$read2_noy_line, read4=$read4_noy_line, read5=$read5_noy_line)"
                    order_preserved=false
                fi
            fi
            
            # Check Y reads order (read1 should come before read3)
            local read1_y_line=$(grep " read1$" "$y_order_file" | awk '{print $1}')
            local read3_y_line=$(grep " read3$" "$y_order_file" | awk '{print $1}')
            local read1_orig_line=$(grep " read1$" "$original_order_file" | awk '{print $1}')
            local read3_orig_line=$(grep " read3$" "$original_order_file" | awk '{print $1}')
            
            if [ -n "$read1_y_line" ] && [ -n "$read3_y_line" ]; then
                if [ "$read1_y_line" -gt "$read3_y_line" ]; then
                    echo "  ✗ FAIL: Y reads order not preserved (read1=$read1_y_line, read3=$read3_y_line)"
                    order_preserved=false
                fi
            fi
            
            if [ "$order_preserved" = true ]; then
                echo "  ✓ Order preservation: Relative order maintained (read2<read4<read5 in noY, read1<read3 in Y)"
                PASSED=$((PASSED + 1))
            else
                check_passed=false
                FAILED=$((FAILED + 1))
            fi
        fi
        
        # Enhanced name membership check using hash-based validation
        # Build Y hash set from BAM (matching tool's logic)
        local y_hash_file="$outdir/y_hash_set.txt"
        extract_bam_qnames "$TEST_DIR/y_reads.bam" | sort -u > "$y_hash_file"
        
        # Verify every Y output read name is in hash set
        local y_membership_fail=false
        local y_missing_count=0
        while IFS= read -r name; do
            if [ -z "$name" ]; then continue; fi
            if ! grep -q "^${name}$" "$y_hash_file"; then
                echo "  ✗ FAIL: $name in Y output but not in Y hash set"
                y_membership_fail=true
                y_missing_count=$((y_missing_count + 1))
            fi
        done <<< "$y_names"
        
        if [ "$y_membership_fail" = true ]; then
            echo "  ✗ FAIL: $y_missing_count read(s) in Y output not in Y hash set"
            check_passed=false
            FAILED=$((FAILED + 1))
        else
            echo "  ✓ Y hash membership: All $y_count reads in Y output verified in hash set"
            PASSED=$((PASSED + 1))
        fi
        
        # Verify no noY output read name is in hash set
        local noy_membership_fail=false
        local noy_contamination_count=0
        while IFS= read -r name; do
            if [ -z "$name" ]; then continue; fi
            if grep -q "^${name}$" "$y_hash_file"; then
                echo "  ✗ FAIL: $name in noY output but present in Y hash set"
                noy_membership_fail=true
                noy_contamination_count=$((noy_contamination_count + 1))
            fi
        done <<< "$noy_names"
        
        if [ "$noy_membership_fail" = true ]; then
            echo "  ✗ FAIL: $noy_contamination_count read(s) in noY output contaminated with Y reads"
            check_passed=false
            FAILED=$((FAILED + 1))
        else
            echo "  ✓ noY hash membership: All $noy_count reads verified absent from Y hash set"
            PASSED=$((PASSED + 1))
        fi
    done
    
    # Overall count check
    echo ""
    echo "  Overall counts: Y=$total_y, noY=$total_noy, original=$total_original"
    if [ $((total_y + total_noy)) -eq "$total_original" ]; then
        echo "  ✓ Overall count consistency"
        PASSED=$((PASSED + 1))
    else
        echo "  ✗ FAIL: Overall count mismatch"
        check_passed=false
        FAILED=$((FAILED + 1))
    fi
    
    # Y count should match Y BAM count (approximately, accounting for multiple alignments)
    if [ "$total_y" -ge "$Y_COUNT" ]; then
        echo "  ✓ Y count >= Y BAM count (expected due to multiple alignments)"
        PASSED=$((PASSED + 1))
    else
        echo "  ⚠ Y count ($total_y) < Y BAM count ($Y_COUNT) - may be expected"
    fi
    
    # Write to report
    cat >> "$REPORT_FILE" << EOF

### $mode Mode (--threads $threads)

| Check | Result |
|-------|--------|
| Output files exist | $([ "$check_passed" = true ] && echo "✓ PASS" || echo "✗ FAIL") |
| Count consistency | $([ $((total_y + total_noy)) -eq "$total_original" ] && echo "✓ PASS" || echo "✗ FAIL") |
| Y hash membership | $([ "$y_membership_fail" != true ] && echo "✓ PASS" || echo "✗ FAIL") |
| noY hash membership | $([ "$noy_membership_fail" != true ] && echo "✓ PASS" || echo "✗ FAIL") |
| Order preservation | $([ "$order_preserved" = true ] && echo "✓ PASS" || echo "✗ FAIL") |
| Overall counts | Y=$total_y, noY=$total_noy, original=$total_original |

EOF
}

# Test 1: Single-threaded
run_test "single_threaded" "1" "$TEST_DIR/test1.fastq.gz"

# Test 2: Multithreaded with single file
run_test "multithreaded_single" "4" "$TEST_DIR/test1.fastq.gz"

# Test 3: Multithreaded with multiple files
run_test "multithreaded_multi" "2" "$TEST_DIR/test1.fastq.gz $TEST_DIR/test2.fastq.gz"

# Final summary
cat >> "$REPORT_FILE" << EOF

---

## Summary

- **Total checks passed**: $PASSED
- **Total checks failed**: $FAILED
- **Status**: $([ $FAILED -eq 0 ] && echo "✓ ALL TESTS PASSED" || echo "✗ SOME TESTS FAILED")

EOF

echo ""
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo "Checks passed: $PASSED"
echo "Checks failed: $FAILED"
echo ""
echo "Report written to: $REPORT_FILE"
echo ""

if [ $FAILED -eq 0 ]; then
    echo "✓ All tests passed!"
    exit 0
else
    echo "✗ Some tests failed"
    exit 1
fi

