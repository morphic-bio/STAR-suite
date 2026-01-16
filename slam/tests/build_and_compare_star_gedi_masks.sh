#!/bin/bash
# Build new STAR SNP mask and compare with existing GEDI mask
# This script:
# 1. Uses existing GEDI mask (from previous runs)
# 2. Builds a NEW STAR mask using the new methodology
# 3. Compares the SNP lists

set -euo pipefail

# Configuration
WORK_DIR="/storage/slam_e2e_arid1a_$(date +%Y%m%d)"
STAR="/mnt/pikachu/STAR-Flex/source/STAR"
STAR_INDEX="/storage/autoindex_110_44/bulk_index"
FASTQ="/mnt/pikachu/NW-5-21/SLAM-Seq/ARID1A-no4su_S50_R1_001.fastq.gz"

# Create directories
mkdir -p "$WORK_DIR"/{mask,star,compare,report}

echo "=== STAR vs GEDI SNP Mask Comparison ==="
echo "Work directory: $WORK_DIR"
echo "FASTQ: $FASTQ"
echo ""

# === TASK 1: Find existing GEDI mask ===
echo "=== TASK 1: Locating existing GEDI mask ==="
GEDI_MASK=""
GEDI_SNPDATA=""

# Check common locations
for loc in \
    "/mnt/pikachu/STAR-Flex/test/tmp_wt_snp_comparison/gedi_snp.snpdata" \
    "/storage/SLAM-Seq-prod-compare-20260109/gedi/*snpdata" \
    "/storage/slam_e2e_*/gedi/*snpdata" \
    "/storage/test/tmp_wt_snp_comparison/gedi_snp.snpdata"; do
    if ls $loc 2>/dev/null | head -1 | grep -q .; then
        GEDI_SNPDATA=$(ls $loc 2>/dev/null | head -1)
        echo "Found GEDI snpdata: $GEDI_SNPDATA"
        break
    fi
done

# Also check for converted BED
for loc in \
    "/mnt/pikachu/STAR-Flex/test/tmp_wt_snp_comparison/gedi_cov20_matched.bed" \
    "/storage/test/tmp_wt_snp_comparison/gedi_cov20_matched.bed"; do
    if [ -f "$loc" ]; then
        GEDI_MASK="$loc"
        echo "Found GEDI BED: $GEDI_MASK"
        break
    fi
done

# If we have snpdata but no BED, convert it
if [ -n "$GEDI_SNPDATA" ] && [ -z "$GEDI_MASK" ]; then
    echo "Converting GEDI snpdata to BED..."
    GEDI_MASK="$WORK_DIR/compare/gedi_mask.bed"
    tail -n +2 "$GEDI_SNPDATA" | awk -F'\t' '$2 >= 20 {
        split($1, loc, ":");
        chrom = loc[1];
        pos = loc[2];
        if (chrom == "MT") chrom = "chrMT";
        else if (chrom !~ /^chr/) chrom = "chr" chrom;
        print chrom "\t" pos "\t" (pos+1)
    }' | sort -k1,1 -k2,2n | uniq > "$GEDI_MASK"
    echo "Created: $GEDI_MASK"
fi

if [ -z "$GEDI_MASK" ] && [ -z "$GEDI_SNPDATA" ]; then
    echo "WARNING: GEDI mask not found. Will build STAR mask only."
    echo "Searched locations:"
    echo "  - /mnt/pikachu/STAR-Flex/test/tmp_wt_snp_comparison/"
    echo "  - /storage/SLAM-Seq-prod-compare-20260109/gedi/"
    echo "  - /storage/slam_e2e_*/gedi/"
fi

# === TASK 2: Build NEW STAR mask ===
echo ""
echo "=== TASK 2: Building NEW STAR SNP mask ==="
STAR_MASK_BED="$WORK_DIR/mask/star_new.mask.bed.gz"
STAR_MASK_SUMMARY="$WORK_DIR/mask/star_new.mask.summary.tsv"

echo "$FASTQ" > "$WORK_DIR/mask/wt.fofn"

echo "Running STAR SNP mask build (new methodology)..."
"$STAR" \
    --runThreadN 8 \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "$FASTQ" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$WORK_DIR/mask/star_new_" \
    --outSAMtype None \
    --slamQuantMode 1 \
    --slamSnpMaskBuildFastqs "$WORK_DIR/mask/wt.fofn" \
    --slamSnpMaskBedOut "$STAR_MASK_BED" \
    --slamSnpMaskSummaryOut "$STAR_MASK_SUMMARY" \
    --slamSnpMaskOnly 1 \
    > "$WORK_DIR/report/build_star_mask.log" 2>&1

if [ ! -f "$STAR_MASK_BED" ]; then
    echo "ERROR: STAR mask build failed. Check log: $WORK_DIR/report/build_star_mask.log"
    exit 1
fi

echo "STAR mask built successfully: $STAR_MASK_BED"
echo "Summary: $STAR_MASK_SUMMARY"

# Display STAR mask summary
if [ -f "$STAR_MASK_SUMMARY" ]; then
    echo ""
    echo "STAR Mask Summary:"
    cat "$STAR_MASK_SUMMARY"
fi

# === TASK 3: Compare masks ===
if [ -n "$GEDI_MASK" ] || [ -n "$GEDI_SNPDATA" ]; then
    echo ""
    echo "=== TASK 3: Comparing STAR vs GEDI masks ==="
    
    # Extract STAR T→C sites (cov≥20) for comparison
    STAR_TC_BED="$WORK_DIR/compare/star_tc_cov20.bed"
    echo "Extracting STAR T→C sites (cov≥20)..."
    zcat "$STAR_MASK_BED" | tail -n +2 | \
        awk -F'\t' '$4=="T" && $5=="C" && $6>=20 {print $1"\t"$2"\t"$3}' | \
        sort -k1,1 -k2,2n | uniq > "$STAR_TC_BED"
    
    STAR_TC_COUNT=$(wc -l < "$STAR_TC_BED")
    echo "STAR T→C cov≥20 sites: $STAR_TC_COUNT"
    
    if [ -n "$GEDI_MASK" ]; then
        GEDI_COUNT=$(wc -l < "$GEDI_MASK")
        echo "GEDI T→C cov≥20 sites: $GEDI_COUNT"
        
        # Compute overlap
        OVERLAP_BED="$WORK_DIR/compare/overlap.bed"
        bedtools intersect -a "$STAR_TC_BED" -b "$GEDI_MASK" -u > "$OVERLAP_BED"
        OVERLAP_COUNT=$(wc -l < "$OVERLAP_BED")
        
        # Compute Jaccard
        UNION_BED="$WORK_DIR/compare/union.bed"
        bedtools merge -i <(cat "$STAR_TC_BED" "$GEDI_MASK" | sort -k1,1 -k2,2n) > "$UNION_BED"
        UNION_COUNT=$(wc -l < "$UNION_BED")
        
        if [ "$UNION_COUNT" -gt 0 ]; then
            JACCARD=$(echo "scale=6; $OVERLAP_COUNT / $UNION_COUNT" | bc)
        else
            JACCARD="0.000000"
        fi
        
        # STAR-only and GEDI-only
        STAR_ONLY_BED="$WORK_DIR/compare/star_only.bed"
        GEDI_ONLY_BED="$WORK_DIR/compare/gedi_only.bed"
        bedtools subtract -a "$STAR_TC_BED" -b "$GEDI_MASK" > "$STAR_ONLY_BED"
        bedtools subtract -a "$GEDI_MASK" -b "$STAR_TC_BED" > "$GEDI_ONLY_BED"
        STAR_ONLY_COUNT=$(wc -l < "$STAR_ONLY_BED")
        GEDI_ONLY_COUNT=$(wc -l < "$GEDI_ONLY_BED")
        
        # Generate comparison report
        REPORT="$WORK_DIR/report/mask_comparison.txt"
        {
            echo "STAR vs GEDI SNP Mask Comparison"
            echo "=================================="
            echo ""
            echo "Date: $(date)"
            echo ""
            echo "STAR Mask (NEW):"
            echo "  BED: $STAR_MASK_BED"
            echo "  Summary: $STAR_MASK_SUMMARY"
            echo "  T→C cov≥20 sites: $STAR_TC_COUNT"
            echo ""
            echo "GEDI Mask (EXISTING):"
            echo "  BED: $GEDI_MASK"
            if [ -n "$GEDI_SNPDATA" ]; then
                echo "  Source snpdata: $GEDI_SNPDATA"
            fi
            echo "  T→C cov≥20 sites: $GEDI_COUNT"
            echo ""
            echo "Overlap Metrics:"
            echo "  Overlap (STAR ∩ GEDI): $OVERLAP_COUNT"
            echo "  STAR in GEDI: $(echo "scale=2; $OVERLAP_COUNT * 100 / $STAR_TC_COUNT" | bc)%"
            echo "  GEDI in STAR: $(echo "scale=2; $OVERLAP_COUNT * 100 / $GEDI_COUNT" | bc)%"
            echo "  Jaccard similarity: $JACCARD"
            echo ""
            echo "Unique Sites:"
            echo "  STAR-only: $STAR_ONLY_COUNT"
            echo "  GEDI-only: $GEDI_ONLY_COUNT"
            echo ""
            echo "Files:"
            echo "  STAR T→C BED: $STAR_TC_BED"
            echo "  GEDI BED: $GEDI_MASK"
            echo "  Overlap BED: $OVERLAP_BED"
            echo "  STAR-only BED: $STAR_ONLY_BED"
            echo "  GEDI-only BED: $GEDI_ONLY_BED"
        } > "$REPORT"
        
        cat "$REPORT"
        echo ""
        echo "Full report: $REPORT"
    else
        echo "GEDI mask BED not available for comparison"
    fi
else
    echo ""
    echo "=== TASK 3: Skipped (GEDI mask not found) ==="
    echo "STAR mask built: $STAR_MASK_BED"
fi

echo ""
echo "=== Workflow completed ==="
