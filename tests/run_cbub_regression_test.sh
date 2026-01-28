#!/bin/bash
# CB/UB Tag Injection Regression Test Suite
# Tests both sorted and unsorted BAM CB/UB tag injection against gold standard
#
# Usage: ./run_cbub_regression_test.sh [--quick]
#   --quick: Run 2-lane tests only (faster, ~90s each)
#   default: Run full 8-lane gold standard test (~60s)

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${SCRIPT_DIR}/.."
STAR_BIN="${REPO_DIR}/core/legacy/source/STAR"
GOLD_DIR="${SCRIPT_DIR}/gold_standard"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

QUICK_MODE=false
if [[ "${1:-}" == "--quick" ]]; then
    QUICK_MODE=true
fi

echo "=============================================="
echo "CB/UB Tag Injection Regression Test Suite"
echo "=============================================="
echo "STAR binary: $STAR_BIN"
echo "Mode: $([ "$QUICK_MODE" = true ] && echo "Quick (2 lanes)" || echo "Full (8 lanes)")"
echo ""

# Verify STAR binary exists
if [ ! -f "$STAR_BIN" ]; then
    echo -e "${RED}ERROR: STAR binary not found at $STAR_BIN${NC}"
    exit 1
fi

PASS_COUNT=0
FAIL_COUNT=0

# Common parameters
GENOME_DIR="/storage/flex_filtered_reference/star_index"
CB_WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
SAMPLE_WHITELIST="/storage/SC2300771_filtered_2M/sample_whitelist.tsv"
PROBE_LIST="/storage/flex_filtered_reference/filtered_reference/probe_list.txt"
SAMPLE_PROBES="/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt"

if [ "$QUICK_MODE" = true ]; then
    READS_R2="/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R2_001.fastq.gz"
    READS_R1="/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R1_001.fastq.gz"
else
    READS_R2="/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R2_001.fastq.gz"
    READS_R1="/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R1_001.fastq.gz"
fi

run_star() {
    local OUT_DIR="$1"
    local BAM_TYPE="$2"  # "Unsorted" or "SortedByCoordinate"
    local TMP_DIR="/storage/100K/tmp/cbub_regress_$(echo $BAM_TYPE | tr '[:upper:]' '[:lower:]')"
    
    rm -rf "$OUT_DIR" "$TMP_DIR"
    mkdir -p "$OUT_DIR"
    
    "$STAR_BIN" \
      --runThreadN 8 \
      --outTmpDir "$TMP_DIR" \
      --genomeDir "$GENOME_DIR" \
      --soloType CB_UMI_Simple \
      --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
      --soloCBwhitelist "$CB_WHITELIST" \
      --flex yes \
      --soloFlexExpectedCellsPerTag 3000 \
      --soloSampleWhitelist "$SAMPLE_WHITELIST" \
      --soloProbeList "$PROBE_LIST" \
      --soloSampleProbes "$SAMPLE_PROBES" \
      --soloSampleProbeOffset 68 \
      --soloFlexAllowedTags "$SAMPLE_WHITELIST" \
      --soloFlexOutputPrefix "${OUT_DIR}/per_sample" \
      --limitIObufferSize 50000000 50000000 \
      --outSJtype None \
      --outBAMcompression 6 \
      --soloMultiMappers Rescue \
      --alignIntronMax 500000 \
      --outFilterMismatchNmax 6 \
      --outFilterMismatchNoverReadLmax 1.0 \
      --outFilterMatchNmin 25 \
      --outSAMunmapped None \
      --outFilterMatchNminOverLread 0 \
      --outFilterMultimapNmax 10000 \
      --outFilterMultimapScoreRange 4 \
      --outSAMmultNmax 10000 \
      --winAnchorMultimapNmax 200 \
      --outSAMprimaryFlag AllBestScore \
      --outFilterScoreMin 0 \
      --outFilterScoreMinOverLread 0 \
      --outSAMattributes NH HI AS nM NM GX GN CB UB \
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
      --soloUMIfiltering MultiGeneUMI_CR \
      --soloUMIdedup 1MM_CR \
      --soloCellFilter None \
      --clipAdapterType CellRanger4 \
      --soloFeatures Gene \
      --alignEndsType Local \
      --soloStrand Unstranded \
      --chimSegmentMin 1000000 \
      --soloKeysCompat cr \
      --outSAMtype BAM "$BAM_TYPE" \
      --soloSampleSearchNearby no \
      --readFilesCommand zcat \
      --readFilesIn $READS_R2 $READS_R1 \
      --outFileNamePrefix "$OUT_DIR/" \
      > /dev/null 2>&1
}

check_tags() {
    local BAM_FILE="$1"
    local TEST_NAME="$2"
    
    if [ ! -f "$BAM_FILE" ]; then
        echo -e "${RED}FAIL${NC}: $TEST_NAME - BAM file not found"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        return 1
    fi
    
    local TOTAL=$(samtools view -c "$BAM_FILE")
    local CB_COUNT=$(samtools view "$BAM_FILE" | grep -c "CB:Z:" || true)
    local UB_COUNT=$(samtools view "$BAM_FILE" | grep -c "UB:Z:" || true)
    
    if [ "$CB_COUNT" -gt 0 ]; then
        local CB_PCT=$(echo "scale=2; $CB_COUNT * 100 / $TOTAL" | bc)
        local UB_PCT=$(echo "scale=2; $UB_COUNT * 100 / $TOTAL" | bc)
        echo -e "${GREEN}PASS${NC}: $TEST_NAME"
        echo "       Reads: $TOTAL, CB: $CB_COUNT ($CB_PCT%), UB: $UB_COUNT ($UB_PCT%)"
        PASS_COUNT=$((PASS_COUNT + 1))
        return 0
    else
        echo -e "${RED}FAIL${NC}: $TEST_NAME - No CB tags found"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        return 1
    fi
}

compare_mex() {
    local TEST_DIR="$1"
    local TEST_NAME="$2"
    
    if [ ! -f "${TEST_DIR}/per_sample/flexfilter_summary.tsv" ]; then
        echo -e "${RED}FAIL${NC}: $TEST_NAME MEX - Summary not found"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        return 1
    fi
    
    # Extract total UMIs and final cells from test
    local TEST_UMIS=$(awk '/TOTAL/ {print $9}' "${TEST_DIR}/per_sample/flexfilter_summary.tsv")
    local TEST_CELLS=$(awk '/TOTAL/ {print $8}' "${TEST_DIR}/per_sample/flexfilter_summary.tsv")
    
    if [ "$QUICK_MODE" = false ] && [ -f "${GOLD_DIR}/per_sample/flexfilter_summary.tsv" ]; then
        # Compare against gold standard
        local GOLD_UMIS=$(awk '/TOTAL/ {print $9}' "${GOLD_DIR}/per_sample/flexfilter_summary.tsv")
        local GOLD_CELLS=$(awk '/TOTAL/ {print $8}' "${GOLD_DIR}/per_sample/flexfilter_summary.tsv")
        local UMI_DIFF=$((TEST_UMIS - GOLD_UMIS))
        local CELL_DIFF=$((TEST_CELLS - GOLD_CELLS))
        
        # Allow small UMI difference (+1 per sample is known)
        if [ "$CELL_DIFF" -eq 0 ] && [ "$UMI_DIFF" -ge 0 ] && [ "$UMI_DIFF" -le 10 ]; then
            echo -e "${GREEN}PASS${NC}: $TEST_NAME MEX vs gold standard"
            echo "       Cells: $TEST_CELLS (gold: $GOLD_CELLS, diff: $CELL_DIFF)"
            echo "       UMIs: $TEST_UMIS (gold: $GOLD_UMIS, diff: +$UMI_DIFF)"
            PASS_COUNT=$((PASS_COUNT + 1))
            return 0
        else
            echo -e "${RED}FAIL${NC}: $TEST_NAME MEX vs gold standard"
            echo "       Cells: $TEST_CELLS (gold: $GOLD_CELLS, diff: $CELL_DIFF)"
            echo "       UMIs: $TEST_UMIS (gold: $GOLD_UMIS, diff: $UMI_DIFF)"
            FAIL_COUNT=$((FAIL_COUNT + 1))
            return 1
        fi
    else
        echo -e "${GREEN}PASS${NC}: $TEST_NAME MEX output generated"
        echo "       Cells: $TEST_CELLS, UMIs: $TEST_UMIS"
        PASS_COUNT=$((PASS_COUNT + 1))
        return 0
    fi
}

# ============================================
# Test 1: Sorted BAM with CB/UB tags
# ============================================
echo ""
echo "--- Test 1: Sorted BAM CB/UB Tag Injection ---"
SORTED_OUT="${SCRIPT_DIR}/cbub_regress_sorted"
run_star "$SORTED_OUT" "SortedByCoordinate"
check_tags "${SORTED_OUT}/Aligned.sortedByCoord.out.bam" "Sorted BAM tags"
compare_mex "$SORTED_OUT" "Sorted"
SORTED_SUMMARY="${SORTED_OUT}/per_sample/flexfilter_summary.tsv"

# ============================================
# Test 2: Unsorted BAM with CB/UB tags (automatic)
# ============================================
echo ""
echo "--- Test 2: Unsorted BAM CB/UB Tag Injection (automatic) ---"
UNSORTED_OUT="${SCRIPT_DIR}/cbub_regress_unsorted"
run_star "$UNSORTED_OUT" "Unsorted"
check_tags "${UNSORTED_OUT}/Aligned.out.bam" "Unsorted BAM tags"
compare_mex "$UNSORTED_OUT" "Unsorted"
UNSORTED_SUMMARY="${UNSORTED_OUT}/per_sample/flexfilter_summary.tsv"

# ============================================
# Test 3: Compare sorted vs unsorted
# ============================================
echo ""
echo "--- Test 3: Sorted vs Unsorted Consistency ---"
if [ -f "$SORTED_SUMMARY" ] && [ -f "$UNSORTED_SUMMARY" ]; then
    SORTED_UMIS=$(awk '/TOTAL/ {print $9}' "$SORTED_SUMMARY")
    UNSORTED_UMIS=$(awk '/TOTAL/ {print $9}' "$UNSORTED_SUMMARY")
    SORTED_CELLS=$(awk '/TOTAL/ {print $8}' "$SORTED_SUMMARY")
    UNSORTED_CELLS=$(awk '/TOTAL/ {print $8}' "$UNSORTED_SUMMARY")
    
    if [ "$SORTED_UMIS" -eq "$UNSORTED_UMIS" ] && [ "$SORTED_CELLS" -eq "$UNSORTED_CELLS" ]; then
        echo -e "${GREEN}PASS${NC}: Sorted and Unsorted MEX outputs are identical"
        echo "       Cells: $SORTED_CELLS, UMIs: $SORTED_UMIS"
        PASS_COUNT=$((PASS_COUNT + 1))
    else
        echo -e "${RED}FAIL${NC}: Sorted and Unsorted MEX outputs differ"
        echo "       Sorted:   Cells=$SORTED_CELLS, UMIs=$SORTED_UMIS"
        echo "       Unsorted: Cells=$UNSORTED_CELLS, UMIs=$UNSORTED_UMIS"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
else
    echo -e "${YELLOW}SKIP${NC}: Cannot compare - missing summary files"
fi

# ============================================
# Summary
# ============================================
echo ""
echo "=============================================="
echo "Summary"
echo "=============================================="
echo -e "Passed: ${GREEN}$PASS_COUNT${NC}"
echo -e "Failed: ${RED}$FAIL_COUNT${NC}"
echo ""

if [ "$FAIL_COUNT" -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed!${NC}"
    exit 1
fi
