#!/bin/bash
# Parity test harness for pileup_snp vs GEDI and STAR

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOL="$SCRIPT_DIR/pileup_snp"

# Check if tool exists
if [ ! -f "$TOOL" ]; then
    echo "ERROR: pileup_snp not found. Run 'make' first."
    exit 1
fi

# Parse arguments
if [ $# -lt 5 ]; then
    echo "Usage: $0 <bam> <ref_fa> <star_bed> <gedi_bed> <candidate_bed> [output_dir]"
    echo ""
    echo "Arguments:"
    echo "  bam: STAR-generated WT BAM file"
    echo "  ref_fa: Reference FASTA (indexed)"
    echo "  star_bed: STAR SNP mask BED (for comparison)"
    echo "  gedi_bed: GEDI SNP mask BED (for comparison)"
    echo "  candidate_bed: Candidate loci BED (subset for testing)"
    echo "  output_dir: Output directory [default: ./parity_test_$(date +%Y%m%d_%H%M%S)]"
    exit 1
fi

BAM="$1"
REF_FA="$2"
STAR_BED="$3"
GEDI_BED="$4"
CANDIDATE_BED="$5"
OUT_DIR="${6:-./parity_test_$(date +%Y%m%d_%H%M%S)}"

mkdir -p "$OUT_DIR"

echo "=== Parity Test: pileup_snp vs STAR/GEDI ==="
echo "BAM: $BAM"
echo "Reference: $REF_FA"
echo "STAR BED: $STAR_BED"
echo "GEDI BED: $GEDI_BED"
echo "Candidates: $CANDIDATE_BED"
echo "Output: $OUT_DIR"
echo ""

# Test 1: STAR-compat mode (conv, primary-only)
echo "[1/4] Running pileup_snp in STAR-compat mode..."
"$TOOL" \
    --bam "$BAM" \
    --bed "$CANDIDATE_BED" \
    --ref "$REF_FA" \
    --output "$OUT_DIR/pileup_star_compat.bed.gz" \
    --debug-tsv "$OUT_DIR/pileup_star_compat_debug.tsv" \
    --kMode conv \
    --includeSecondary 0 \
    --minMapQ 0 \
    --minBaseQ 0 \
    --pErr 0.001 \
    --pval 0.001 \
    --minTcRatio 0.3 \
    --minCov 6 \
    --minAlt 1

# Test 2: GEDI-compat mode (any mismatch, include secondary)
echo "[2/4] Running pileup_snp in GEDI-compat mode..."
"$TOOL" \
    --bam "$BAM" \
    --bed "$CANDIDATE_BED" \
    --ref "$REF_FA" \
    --output "$OUT_DIR/pileup_gedi_compat.bed.gz" \
    --debug-tsv "$OUT_DIR/pileup_gedi_compat_debug.tsv" \
    --kMode any \
    --includeSecondary 1 \
    --minMapQ 0 \
    --minBaseQ 0 \
    --pErr 0.001 \
    --pval 0.001 \
    --minTcRatio 0.3 \
    --minCov 6 \
    --minAlt 1

# Test 3: Compare overlaps
echo "[3/4] Computing overlaps..."

# Intersect with STAR BED
if [ -f "$STAR_BED" ]; then
    if [[ "$STAR_BED" == *.gz ]]; then
        STAR_INTERSECT="$OUT_DIR/star_intersect.bed"
        zcat "$STAR_BED" | bedtools intersect -a "$CANDIDATE_BED" -b stdin -wa -u > "$STAR_INTERSECT" || true
    else
        STAR_INTERSECT="$OUT_DIR/star_intersect.bed"
        bedtools intersect -a "$CANDIDATE_BED" -b "$STAR_BED" -wa -u > "$STAR_INTERSECT" || true
    fi
    
    PILEUP_STAR_INTERSECT="$OUT_DIR/pileup_star_intersect.bed"
    zcat "$OUT_DIR/pileup_star_compat.bed.gz" | bedtools intersect -a "$CANDIDATE_BED" -b stdin -wa -u > "$PILEUP_STAR_INTERSECT" || true
    
    STAR_COUNT=$(wc -l < "$STAR_INTERSECT" || echo 0)
    PILEUP_STAR_COUNT=$(wc -l < "$PILEUP_STAR_INTERSECT" || echo 0)
    
    echo "STAR overlap: $STAR_COUNT sites in candidate set"
    echo "pileup_snp (STAR-compat) overlap: $PILEUP_STAR_COUNT sites in candidate set"
fi

# Intersect with GEDI BED
if [ -f "$GEDI_BED" ]; then
    if [[ "$GEDI_BED" == *.gz ]]; then
        GEDI_INTERSECT="$OUT_DIR/gedi_intersect.bed"
        zcat "$GEDI_BED" | bedtools intersect -a "$CANDIDATE_BED" -b stdin -wa -u > "$GEDI_INTERSECT" || true
    else
        GEDI_INTERSECT="$OUT_DIR/gedi_intersect.bed"
        bedtools intersect -a "$CANDIDATE_BED" -b "$GEDI_BED" -wa -u > "$GEDI_INTERSECT" || true
    fi
    
    PILEUP_GEDI_INTERSECT="$OUT_DIR/pileup_gedi_intersect.bed"
    zcat "$OUT_DIR/pileup_gedi_compat.bed.gz" | bedtools intersect -a "$CANDIDATE_BED" -b stdin -wa -u > "$PILEUP_GEDI_INTERSECT" || true
    
    GEDI_COUNT=$(wc -l < "$GEDI_INTERSECT" || echo 0)
    PILEUP_GEDI_COUNT=$(wc -l < "$PILEUP_GEDI_INTERSECT" || echo 0)
    
    echo "GEDI overlap: $GEDI_COUNT sites in candidate set"
    echo "pileup_snp (GEDI-compat) overlap: $PILEUP_GEDI_COUNT sites in candidate set"
fi

# Test 4: Generate summary report
echo "[4/4] Generating summary report..."
SUMMARY="$OUT_DIR/summary.txt"
{
    echo "=== Parity Test Summary ==="
    echo "Date: $(date)"
    echo ""
    echo "Inputs:"
    echo "  BAM: $BAM"
    echo "  Reference: $REF_FA"
    echo "  Candidate loci: $(wc -l < "$CANDIDATE_BED")"
    echo ""
    echo "Outputs:"
    echo "  pileup_snp (STAR-compat): $(zcat "$OUT_DIR/pileup_star_compat.bed.gz" | grep -v '^#' | wc -l) sites"
    echo "  pileup_snp (GEDI-compat): $(zcat "$OUT_DIR/pileup_gedi_compat.bed.gz" | grep -v '^#' | wc -l) sites"
    echo ""
    if [ -f "$STAR_BED" ]; then
        echo "STAR comparison:"
        echo "  STAR sites in candidates: $STAR_COUNT"
        echo "  pileup_snp (STAR-compat) sites in candidates: $PILEUP_STAR_COUNT"
    fi
    if [ -f "$GEDI_BED" ]; then
        echo "GEDI comparison:"
        echo "  GEDI sites in candidates: $GEDI_COUNT"
        echo "  pileup_snp (GEDI-compat) sites in candidates: $PILEUP_GEDI_COUNT"
    fi
} > "$SUMMARY"

cat "$SUMMARY"
echo ""
echo "Results written to: $OUT_DIR"
