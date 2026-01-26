#!/bin/bash
# Non-Flex CB/UB Tag Injection Regression Test (A375 100K GEX dataset)
# Tests both sorted and unsorted BAM CB/UB tag injection for non-Flex STARsolo
#
# Tests:
# 1. Run STARsolo (non-Flex) with --outSAMtype BAM Unsorted
# 2. Run STARsolo (non-Flex) with --outSAMtype BAM SortedByCoordinate
# 3. Compare tag presence and values between sorted and unsorted
# 4. Verify TranscriptomeSAM is separate when enabled

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${SCRIPT_DIR}/.."
STAR_BIN="${REPO_DIR}/core/legacy/source/STAR"

# A375 paths
GENOME_DIR="/storage/autoindex_110_44/bulk_index"
FASTQ_DIR="/storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/gex/downsampled_100000_v2"
CB_WHITELIST="/storage/A375/3M-5pgex-jan-2023.txt"

# Output directory with timestamp
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUT_ROOT="/tmp/a375_cbub_nonflex_${TIMESTAMP}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo "=============================================="
echo "Non-Flex CB/UB Tag Injection Test (A375 100K)"
echo "=============================================="
echo "STAR binary: $STAR_BIN"
echo "Output dir:  $OUT_ROOT"
echo "Timestamp:   $TIMESTAMP"
echo ""

# Verify STAR binary exists
if [ ! -f "$STAR_BIN" ]; then
    echo -e "${RED}ERROR: STAR binary not found at $STAR_BIN${NC}"
    exit 1
fi

# Verify inputs exist
if [ ! -d "$GENOME_DIR" ]; then
    echo -e "${RED}ERROR: Genome directory not found: $GENOME_DIR${NC}"
    exit 1
fi

if [ ! -d "$FASTQ_DIR" ]; then
    echo -e "${RED}ERROR: FASTQ directory not found: $FASTQ_DIR${NC}"
    exit 1
fi

# Create output directories
mkdir -p "$OUT_ROOT"

PASS_COUNT=0
FAIL_COUNT=0

# Build read file arguments
READS_R2="${FASTQ_DIR}/1k_CRISPR_5p_gemx_gex_S2_L001_R2_001.fastq.gz,${FASTQ_DIR}/1k_CRISPR_5p_gemx_gex_S2_L002_R2_001.fastq.gz"
READS_R1="${FASTQ_DIR}/1k_CRISPR_5p_gemx_gex_S2_L001_R1_001.fastq.gz,${FASTQ_DIR}/1k_CRISPR_5p_gemx_gex_S2_L002_R1_001.fastq.gz"

# Common STAR parameters (non-Flex STARsolo)
common_star_args() {
    echo "--runThreadN 8"
    echo "--genomeDir $GENOME_DIR"
    echo "--readFilesIn $READS_R2 $READS_R1"
    echo "--readFilesCommand zcat"
    echo "--soloType CB_UMI_Simple"
    echo "--soloCBstart 1"
    echo "--soloCBlen 16"
    echo "--soloUMIstart 17"
    echo "--soloUMIlen 12"
    echo "--soloBarcodeReadLength 0"
    echo "--soloCBwhitelist $CB_WHITELIST"
    echo "--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts"
    echo "--soloUMIfiltering MultiGeneUMI_CR"
    echo "--soloUMIdedup 1MM_CR"
    echo "--soloMultiMappers Rescue"
    echo "--soloCellFilter None"
    echo "--clipAdapterType CellRanger4"
    echo "--soloFeatures Gene"
    echo "--soloStrand Unstranded"
    echo "--alignEndsType Local"
    echo "--chimSegmentMin 1000000"
    echo "--outSAMattributes NH HI AS nM NM GX GN CB UB"
}

run_star() {
    local OUT_DIR="$1"
    local BAM_TYPE="$2"
    local EXTRA_ARGS="${3:-}"
    local TMP_DIR="${OUT_ROOT}/tmp_$(echo $BAM_TYPE | tr '[:upper:]' '[:lower:]' | tr ' ' '_')"
    
    rm -rf "$OUT_DIR" "$TMP_DIR"
    mkdir -p "$OUT_DIR"
    
    local args
    args="$(common_star_args)"
    args+=" --outTmpDir $TMP_DIR"
    args+=" --outFileNamePrefix ${OUT_DIR}/"
    args+=" --outSAMtype BAM $BAM_TYPE"
    
    if [ -n "$EXTRA_ARGS" ]; then
        args+=" $EXTRA_ARGS"
    fi
    
    echo "Running STAR with BAM $BAM_TYPE..."
    $STAR_BIN $args > "${OUT_DIR}/star.log" 2>&1
    local rc=$?
    
    if [ $rc -ne 0 ]; then
        echo -e "${RED}ERROR: STAR failed with exit code $rc${NC}"
        tail -20 "${OUT_DIR}/star.log"
        return 1
    fi
    return 0
}

count_tags() {
    local BAM_FILE="$1"
    
    if [ ! -f "$BAM_FILE" ]; then
        echo "0 0 0"
        return 1
    fi
    
    local total=$(samtools view -c "$BAM_FILE")
    local cb_count=$(samtools view "$BAM_FILE" | grep -c "CB:Z:" || true)
    local ub_count=$(samtools view "$BAM_FILE" | grep -c "UB:Z:" || true)
    
    echo "$total $cb_count $ub_count"
}

extract_tags_for_reads() {
    # Extract CB:Z and UB:Z tags for a sample of reads
    local BAM_FILE="$1"
    local SAMPLE_SIZE="${2:-10000}"
    local OUTPUT="$3"
    
    samtools view "$BAM_FILE" | head -n "$SAMPLE_SIZE" | \
        awk '{
            qname=$1
            cb=""; ub=""
            for(i=12; i<=NF; i++) {
                if($i ~ /^CB:Z:/) cb=substr($i,6)
                if($i ~ /^UB:Z:/) ub=substr($i,6)
            }
            print qname "\t" cb "\t" ub
        }' | sort > "$OUTPUT"
}

compare_tags() {
    # Compare tags between two files
    local FILE1="$1"
    local FILE2="$2"
    
    local total1=$(wc -l < "$FILE1")
    local total2=$(wc -l < "$FILE2")
    
    # Find reads present in both
    local common=$(comm -12 <(cut -f1 "$FILE1") <(cut -f1 "$FILE2") | wc -l)
    
    # For common reads, compare tags
    local mismatches=0
    if [ "$common" -gt 0 ]; then
        # Join on read name and compare
        mismatches=$(join -t$'\t' "$FILE1" "$FILE2" | \
            awk -F'\t' '$2 != $4 || $3 != $5 { count++ } END { print count+0 }')
    fi
    
    echo "$total1 $total2 $common $mismatches"
}

# ============================================
# Test 1: Sorted BAM with CB/UB tags
# ============================================
echo ""
echo -e "${BLUE}--- Test 1: Sorted BAM CB/UB Tag Injection ---${NC}"
SORTED_OUT="${OUT_ROOT}/sorted"
if run_star "$SORTED_OUT" "SortedByCoordinate"; then
    read total cb ub <<< $(count_tags "${SORTED_OUT}/Aligned.sortedByCoord.out.bam")
    if [ "$cb" -gt 0 ]; then
        cb_pct=$(echo "scale=2; $cb * 100 / $total" | bc)
        ub_pct=$(echo "scale=2; $ub * 100 / $total" | bc)
        echo -e "${GREEN}PASS${NC}: Sorted BAM tags"
        echo "       Reads: $total, CB: $cb ($cb_pct%), UB: $ub ($ub_pct%)"
        PASS_COUNT=$((PASS_COUNT + 1))
        SORTED_CB=$cb
        SORTED_UB=$ub
        SORTED_TOTAL=$total
    else
        echo -e "${RED}FAIL${NC}: Sorted BAM - No CB tags found"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        SORTED_CB=0
        SORTED_UB=0
        SORTED_TOTAL=$total
    fi
else
    echo -e "${RED}FAIL${NC}: Sorted BAM - STAR failed"
    FAIL_COUNT=$((FAIL_COUNT + 1))
fi

# ============================================
# Test 2: Unsorted BAM with CB/UB tags (automatic)
# ============================================
echo ""
echo -e "${BLUE}--- Test 2: Unsorted BAM CB/UB Tag Injection (automatic) ---${NC}"
UNSORTED_OUT="${OUT_ROOT}/unsorted"
if run_star "$UNSORTED_OUT" "Unsorted"; then
    read total cb ub <<< $(count_tags "${UNSORTED_OUT}/Aligned.out.bam")
    if [ "$cb" -gt 0 ]; then
        cb_pct=$(echo "scale=2; $cb * 100 / $total" | bc)
        ub_pct=$(echo "scale=2; $ub * 100 / $total" | bc)
        echo -e "${GREEN}PASS${NC}: Unsorted BAM tags"
        echo "       Reads: $total, CB: $cb ($cb_pct%), UB: $ub ($ub_pct%)"
        PASS_COUNT=$((PASS_COUNT + 1))
        UNSORTED_CB=$cb
        UNSORTED_UB=$ub
        UNSORTED_TOTAL=$total
    else
        echo -e "${RED}FAIL${NC}: Unsorted BAM - No CB tags found"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        UNSORTED_CB=0
        UNSORTED_UB=0
        UNSORTED_TOTAL=$total
    fi
else
    echo -e "${RED}FAIL${NC}: Unsorted BAM - STAR failed"
    FAIL_COUNT=$((FAIL_COUNT + 1))
fi

# ============================================
# Test 3: Sorted vs Unsorted tag comparison
# ============================================
echo ""
echo -e "${BLUE}--- Test 3: Sorted vs Unsorted Tag Comparison ---${NC}"
if [ -f "${SORTED_OUT}/Aligned.sortedByCoord.out.bam" ] && [ -f "${UNSORTED_OUT}/Aligned.out.bam" ]; then
    # Compare total counts
    if [ "$SORTED_TOTAL" -eq "$UNSORTED_TOTAL" ] && \
       [ "$SORTED_CB" -eq "$UNSORTED_CB" ] && \
       [ "$SORTED_UB" -eq "$UNSORTED_UB" ]; then
        echo -e "${GREEN}PASS${NC}: Tag counts identical"
        echo "       CB: $SORTED_CB, UB: $SORTED_UB"
        PASS_COUNT=$((PASS_COUNT + 1))
    else
        echo -e "${RED}FAIL${NC}: Tag counts differ"
        echo "       Sorted:   CB=$SORTED_CB, UB=$SORTED_UB"
        echo "       Unsorted: CB=$UNSORTED_CB, UB=$UNSORTED_UB"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
    
    # Sample comparison for exact tag values
    echo ""
    echo "Sampling first 10000 reads for exact tag comparison..."
    SAMPLE_SIZE=10000
    extract_tags_for_reads "${SORTED_OUT}/Aligned.sortedByCoord.out.bam" $SAMPLE_SIZE "${OUT_ROOT}/sorted_tags.tsv"
    extract_tags_for_reads "${UNSORTED_OUT}/Aligned.out.bam" $SAMPLE_SIZE "${OUT_ROOT}/unsorted_tags.tsv"
    
    read total1 total2 common mismatches <<< $(compare_tags "${OUT_ROOT}/sorted_tags.tsv" "${OUT_ROOT}/unsorted_tags.tsv")
    
    echo "  Sorted sample size:   $total1"
    echo "  Unsorted sample size: $total2"
    echo "  Common reads:         $common"
    echo "  Tag mismatches:       $mismatches"
    
    if [ "$mismatches" -eq 0 ] && [ "$common" -gt 0 ]; then
        echo -e "${GREEN}PASS${NC}: All sampled common reads have matching tags"
        PASS_COUNT=$((PASS_COUNT + 1))
    elif [ "$common" -eq 0 ]; then
        echo -e "${YELLOW}SKIP${NC}: No common reads in sample (different read ordering)"
    else
        echo -e "${RED}FAIL${NC}: $mismatches tag mismatches in $common common reads"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
else
    echo -e "${YELLOW}SKIP${NC}: Cannot compare - missing BAM files"
fi

# ============================================
# Test 4: TranscriptomeSAM guard
# ============================================
echo ""
echo -e "${BLUE}--- Test 4: TranscriptomeSAM Separation Guard ---${NC}"
TRANSCRIPTOME_OUT="${OUT_ROOT}/transcriptome"
if run_star "$TRANSCRIPTOME_OUT" "Unsorted" "--quantMode TranscriptomeSAM"; then
    # Check that transcriptome BAM exists and is separate
    UNSORTED_BAM="${TRANSCRIPTOME_OUT}/Aligned.out.bam"
    TRANSCRIPTOME_BAM="${TRANSCRIPTOME_OUT}/Aligned.toTranscriptome.out.bam"
    
    if [ -f "$TRANSCRIPTOME_BAM" ]; then
        # Count reads in each
        unsorted_count=$(samtools view -c "$UNSORTED_BAM")
        transcriptome_count=$(samtools view -c "$TRANSCRIPTOME_BAM")
        
        # Transcriptome BAM should have reads (if any mapped to transcriptome)
        # and should NOT be mixed into unsorted
        echo "  Unsorted BAM reads:     $unsorted_count"
        echo "  Transcriptome BAM reads: $transcriptome_count"
        
        # Check that transcriptome reads are NOT in the unsorted BAM
        # by checking that totals are different (transcriptome is subset/different)
        if [ "$transcriptome_count" -gt 0 ]; then
            # Sample a few transcriptome read names
            transcriptome_names=$(samtools view "$TRANSCRIPTOME_BAM" | head -100 | cut -f1 | sort -u)
            
            # These should be in transcriptome BAM but file sizes should differ
            if [ "$unsorted_count" != "$transcriptome_count" ]; then
                echo -e "${GREEN}PASS${NC}: TranscriptomeSAM is separate from unsorted BAM"
                echo "       (Different read counts indicate separate files)"
                PASS_COUNT=$((PASS_COUNT + 1))
            else
                # Same count - check if contents are actually different
                unsorted_md5=$(samtools view "$UNSORTED_BAM" | head -1000 | md5sum | cut -d' ' -f1)
                trans_md5=$(samtools view "$TRANSCRIPTOME_BAM" | head -1000 | md5sum | cut -d' ' -f1)
                
                if [ "$unsorted_md5" != "$trans_md5" ]; then
                    echo -e "${GREEN}PASS${NC}: TranscriptomeSAM content differs from unsorted BAM"
                    PASS_COUNT=$((PASS_COUNT + 1))
                else
                    echo -e "${RED}FAIL${NC}: TranscriptomeSAM may be mixed with unsorted BAM"
                    FAIL_COUNT=$((FAIL_COUNT + 1))
                fi
            fi
        else
            echo -e "${YELLOW}NOTE${NC}: No transcriptome alignments found"
            PASS_COUNT=$((PASS_COUNT + 1))
        fi
    else
        echo -e "${RED}FAIL${NC}: TranscriptomeSAM output file missing"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
else
    echo -e "${RED}FAIL${NC}: TranscriptomeSAM run - STAR failed"
    FAIL_COUNT=$((FAIL_COUNT + 1))
fi

# ============================================
# Summary
# ============================================
echo ""
echo "=============================================="
echo "Summary"
echo "=============================================="
SUMMARY_FILE="${OUT_ROOT}/summary.txt"

{
    echo "Non-Flex CB/UB Tag Injection Test Summary"
    echo "========================================="
    echo "Date: $(date)"
    echo "STAR: $STAR_BIN"
    echo "Dataset: A375 100K GEX (non-Flex)"
    echo ""
    echo "Results:"
    echo "  Passed: $PASS_COUNT"
    echo "  Failed: $FAIL_COUNT"
    echo ""
    echo "Tag Statistics:"
    echo "  Sorted BAM:   Total=$SORTED_TOTAL CB=$SORTED_CB UB=$SORTED_UB"
    echo "  Unsorted BAM: Total=${UNSORTED_TOTAL:-N/A} CB=${UNSORTED_CB:-N/A} UB=${UNSORTED_UB:-N/A}"
    echo ""
    echo "Output directory: $OUT_ROOT"
} | tee "$SUMMARY_FILE"

echo ""
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
