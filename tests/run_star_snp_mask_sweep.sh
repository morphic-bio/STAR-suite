#!/bin/bash
# STAR SNP-Mask Parameter Sweep
# Compares STAR SNP detection vs GEDI reference across parameter combinations

set -e

OUTDIR="test/tmp_wt_snp_comparison"
WT_FASTQ="/storage/SLAM-Seq-prod-compare-20260109/input/ARID1A-no4su_S50_R1_001.fastq.gz"
STAR_BIN="./core/legacy/source/STAR"
GENOME_DIR="/storage/autoindex_110_44/bulk_index"

# Ensure GEDI reference exists
if [[ ! -f "$OUTDIR/gedi_cov20_matched.bed" ]]; then
    echo "ERROR: GEDI reference BED not found. Run normalization first."
    exit 1
fi
GEDI_COUNT=$(wc -l < "$OUTDIR/gedi_cov20_matched.bed")
echo "GEDI T→C cov≥20 reference: $GEDI_COUNT sites"
echo ""

# Create FOFN
echo "$WT_FASTQ" > "$OUTDIR/wt.fofn"

# Parameter arrays
MIN_COVS=(20 10 6)
POSTERIORS=(0.99 0.95 0.90)
MIN_ALTS=(3 2)

RESULTS_FILE="$OUTDIR/sweep_results.tsv"
echo -e "MinCov\tPosterior\tMinAlt\tSTAR_Total\tSTAR_TC_cov20\tOverlap\tSTAR_in_GEDI_pct\tGEDI_in_STAR_pct\tJaccard" > "$RESULTS_FILE"

COMBINATIONS=0
TOTAL=$(( ${#MIN_COVS[@]} * ${#POSTERIORS[@]} * ${#MIN_ALTS[@]} ))

for COV in "${MIN_COVS[@]}"; do
    for POST in "${POSTERIORS[@]}"; do
        for ALT in "${MIN_ALTS[@]}"; do
            TAG="cov${COV}_post${POST}_alt${ALT}"
            COMBINATIONS=$((COMBINATIONS + 1))
            
            echo "[$COMBINATIONS/$TOTAL] Testing: MinCov=$COV, Posterior=$POST, MinAlt=$ALT"
            
            # Run STAR SNP mask build
            $STAR_BIN \
                --runThreadN 8 \
                --genomeDir "$GENOME_DIR" \
                --readFilesIn "$WT_FASTQ" \
                --readFilesCommand zcat \
                --outFileNamePrefix "$OUTDIR/mask_${TAG}_" \
                --outSAMtype None \
                --slamQuantMode 1 \
                --slamSnpMaskBuildFastqs "$OUTDIR/wt.fofn" \
                --slamSnpMaskBedOut "$OUTDIR/star_snps_${TAG}.bed.gz" \
                --slamSnpMaskSummaryOut "$OUTDIR/star_snps_${TAG}_summary.tsv" \
                --slamSnpMaskOnly 1 \
                --slamSnpMaskMinCov $COV \
                --slamSnpMaskPosterior $POST \
                --slamSnpMaskMinAlt $ALT \
                > "$OUTDIR/mask_${TAG}.log" 2>&1
            
            # Extract T→C SNPs with cov≥20 (using raw coordinates, matching report)
            zcat "$OUTDIR/star_snps_${TAG}.bed.gz" | tail -n +2 | \
                awk -F'\t' '$4=="T" && $5=="C" && $6>=20 {print $1"\t"$2"\t"$3}' | \
                sort -k1,1 -k2,2n > "$OUTDIR/star_tc_${TAG}.bed"
            
            STAR_TC=$(wc -l < "$OUTDIR/star_tc_${TAG}.bed")
            STAR_TOTAL=$(zcat "$OUTDIR/star_snps_${TAG}.bed.gz" | tail -n +2 | wc -l)
            
            # Compute overlap with GEDI
            OVERLAP=$(bedtools intersect -a "$OUTDIR/star_tc_${TAG}.bed" \
                -b "$OUTDIR/gedi_cov20_matched.bed" -u 2>/dev/null | wc -l)
            
            # GEDI in STAR
            GEDI_IN_STAR=$(bedtools intersect -a "$OUTDIR/gedi_cov20_matched.bed" \
                -b "$OUTDIR/star_tc_${TAG}.bed" -u 2>/dev/null | wc -l)
            
            # Calculate percentages and Jaccard
            if [[ $STAR_TC -gt 0 ]]; then
                STAR_IN_GEDI_PCT=$(echo "scale=2; $OVERLAP * 100 / $STAR_TC" | bc)
            else
                STAR_IN_GEDI_PCT=0
            fi
            
            if [[ $GEDI_COUNT -gt 0 ]]; then
                GEDI_IN_STAR_PCT=$(echo "scale=2; $GEDI_IN_STAR * 100 / $GEDI_COUNT" | bc)
            else
                GEDI_IN_STAR_PCT=0
            fi
            
            UNION=$((STAR_TC + GEDI_COUNT - OVERLAP))
            if [[ $UNION -gt 0 ]]; then
                JACCARD=$(echo "scale=4; $OVERLAP / $UNION" | bc)
            else
                JACCARD=0
            fi
            
            # Append to results
            echo -e "$COV\t$POST\t$ALT\t$STAR_TOTAL\t$STAR_TC\t$OVERLAP\t$STAR_IN_GEDI_PCT\t$GEDI_IN_STAR_PCT\t$JACCARD" >> "$RESULTS_FILE"
            
            echo "  → STAR T→C cov≥20: $STAR_TC, Overlap: $OVERLAP, Jaccard: $JACCARD"
            echo ""
        done
    done
done

echo "=== Sweep Complete ==="
echo ""
echo "Results table: $RESULTS_FILE"
echo ""
echo "Summary:"
column -t -s $'\t' "$RESULTS_FILE" | head -20
