#!/bin/bash
# Test binomial p-value SNP mask implementation
# This script builds a mask using the new binomial model and compares with EM model

set -euo pipefail

# Configuration
WORK_DIR="/storage/test_binomial_mask_$(date +%Y%m%d)"
STAR="/mnt/pikachu/STAR-Flex/source/STAR"
STAR_INDEX="/storage/autoindex_110_44/bulk_index"
FASTQ="/mnt/pikachu/NW-5-21/SLAM-Seq/ARID1A-no4su_S50_R1_001.fastq.gz"

# Create directories
mkdir -p "$WORK_DIR"/{binom,em,compare,report}

echo "=== Testing Binomial SNP Mask Implementation ==="
echo "Work directory: $WORK_DIR"
echo "FASTQ: $FASTQ"
echo ""

# Check if FASTQ exists
if [ ! -f "$FASTQ" ]; then
    echo "ERROR: FASTQ not found: $FASTQ"
    exit 1
fi

# Create FOFN
echo "$FASTQ" > "$WORK_DIR/wt.fofn"

# === Test 1: Binomial model (default) ===
echo "=== Test 1: Building mask with binomial model (default) ==="
echo "Running STAR SNP mask build with binomial model..."

"$STAR" \
    --runThreadN 8 \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "$FASTQ" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$WORK_DIR/binom/binom_" \
    --outSAMtype None \
    --slamQuantMode 1 \
    --slamSnpMaskBuildFastqs "$WORK_DIR/wt.fofn" \
    --slamSnpMaskBedOut "$WORK_DIR/binom/binom_mask.bed.gz" \
    --slamSnpMaskSummaryOut "$WORK_DIR/binom/binom_summary.tsv" \
    --slamSnpMaskOnly 1 \
    --autoTrim variance \
    --trimScope first \
    --trimSource "$FASTQ" \
    --autoTrimDetectionReads 100000 \
    > "$WORK_DIR/report/binom_build.log" 2>&1

if [ ! -f "$WORK_DIR/binom/binom_mask.bed.gz" ]; then
    echo "ERROR: Binomial mask build failed. Check log: $WORK_DIR/report/binom_build.log"
    exit 1
fi

echo "Binomial mask built successfully"
echo ""

# Display summary
if [ -f "$WORK_DIR/binom/binom_summary.tsv" ]; then
    echo "Binomial Mask Summary:"
    cat "$WORK_DIR/binom/binom_summary.tsv"
    echo ""
fi

# Extract key stats from log
echo "Key stats from log:"
grep -E "building SNP mask|binomial|model|masked|p_err|pval|finished" "$WORK_DIR/report/binom_build.log" | head -20 || true
echo ""

# === Test 2: EM model (for comparison) ===
echo "=== Test 2: Building mask with EM model (for comparison) ==="
echo "Running STAR SNP mask build with EM model..."

"$STAR" \
    --runThreadN 8 \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "$FASTQ" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$WORK_DIR/em/em_" \
    --outSAMtype None \
    --slamQuantMode 1 \
    --slamSnpMaskBuildFastqs "$WORK_DIR/wt.fofn" \
    --slamSnpMaskBedOut "$WORK_DIR/em/em_mask.bed.gz" \
    --slamSnpMaskSummaryOut "$WORK_DIR/em/em_summary.tsv" \
    --slamSnpMaskOnly 1 \
    --slamSnpMaskModel em \
    > "$WORK_DIR/report/em_build.log" 2>&1

if [ ! -f "$WORK_DIR/em/em_mask.bed.gz" ]; then
    echo "ERROR: EM mask build failed. Check log: $WORK_DIR/report/em_build.log"
    exit 1
fi

echo "EM mask built successfully"
echo ""

# Display summary
if [ -f "$WORK_DIR/em/em_summary.tsv" ]; then
    echo "EM Mask Summary:"
    cat "$WORK_DIR/em/em_summary.tsv"
    echo ""
fi

# === Test 3: Compare masks ===
echo "=== Test 3: Comparing binomial vs EM masks ==="

# Extract T→C sites (cov≥20) from both masks
zcat "$WORK_DIR/binom/binom_mask.bed.gz" | tail -n +2 | \
    awk -F'\t' '$4=="T" && $5=="C" && $6>=20 {print $1"\t"$2"\t"$3}' | \
    sort -k1,1 -k2,2n | uniq > "$WORK_DIR/compare/binom_tc_cov20.bed"

zcat "$WORK_DIR/em/em_mask.bed.gz" | tail -n +2 | \
    awk -F'\t' '$4=="T" && $5=="C" && $6>=20 {print $1"\t"$2"\t"$3}' | \
    sort -k1,1 -k2,2n | uniq > "$WORK_DIR/compare/em_tc_cov20.bed"

BINOM_COUNT=$(wc -l < "$WORK_DIR/compare/binom_tc_cov20.bed")
EM_COUNT=$(wc -l < "$WORK_DIR/compare/em_tc_cov20.bed")

echo "Mask counts (T→C, cov≥20):"
echo "  Binomial: $BINOM_COUNT sites"
echo "  EM: $EM_COUNT sites"
echo ""

# Compute overlap
OVERLAP_BED="$WORK_DIR/compare/overlap.bed"
bedtools intersect -a "$WORK_DIR/compare/binom_tc_cov20.bed" \
    -b "$WORK_DIR/compare/em_tc_cov20.bed" -u > "$OVERLAP_BED"
OVERLAP_COUNT=$(wc -l < "$OVERLAP_BED")

# Compute Jaccard
UNION_BED="$WORK_DIR/compare/union.bed"
bedtools merge -i <(cat "$WORK_DIR/compare/binom_tc_cov20.bed" "$WORK_DIR/compare/em_tc_cov20.bed" | sort -k1,1 -k2,2n) > "$UNION_BED"
UNION_COUNT=$(wc -l < "$UNION_BED")

if [ "$UNION_COUNT" -gt 0 ]; then
    JACCARD=$(echo "scale=6; $OVERLAP_COUNT / $UNION_COUNT" | bc)
else
    JACCARD="0.000000"
fi

echo "Overlap metrics:"
echo "  Overlap: $OVERLAP_COUNT sites"
echo "  Binomial in EM: $(echo "scale=2; $OVERLAP_COUNT * 100 / $BINOM_COUNT" | bc)%"
echo "  EM in Binomial: $(echo "scale=2; $OVERLAP_COUNT * 100 / $EM_COUNT" | bc)%"
echo "  Jaccard: $JACCARD"
echo ""

# === Test 4: Verify p_err usage ===
echo "=== Test 4: Verifying p_err usage in binomial model ==="
echo "Checking if p_err was computed and used..."

if grep -q "snp_err_used" "$WORK_DIR/report/binom_build.log"; then
    echo "✓ p_err computation found in log"
    grep "snp_err_used\|p_err\|binomial model parameters" "$WORK_DIR/report/binom_build.log" | head -10
else
    echo "WARNING: p_err not found in log (may not have run detection pass)"
fi
echo ""

# === Summary ===
echo "=== Test Summary ==="
{
    echo "Binomial vs EM SNP Mask Comparison"
    echo "===================================="
    echo ""
    echo "Date: $(date)"
    echo "Work directory: $WORK_DIR"
    echo ""
    echo "Inputs:"
    echo "  FASTQ: $FASTQ"
    echo "  STAR index: $STAR_INDEX"
    echo ""
    echo "Results:"
    echo "  Binomial mask: $WORK_DIR/binom/binom_mask.bed.gz"
    echo "  EM mask: $WORK_DIR/em/em_mask.bed.gz"
    echo ""
    echo "Mask Counts (T→C, cov≥20):"
    echo "  Binomial: $BINOM_COUNT"
    echo "  EM: $EM_COUNT"
    echo ""
    echo "Overlap:"
    echo "  Overlap: $OVERLAP_COUNT sites"
    echo "  Jaccard: $JACCARD"
    echo ""
} > "$WORK_DIR/report/test_summary.txt"

cat "$WORK_DIR/report/test_summary.txt"
echo ""
echo "Full summary: $WORK_DIR/report/test_summary.txt"
echo ""
echo "=== Test completed successfully ==="
