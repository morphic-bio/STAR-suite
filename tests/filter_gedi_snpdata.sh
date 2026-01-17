#!/bin/bash
# filter_gedi_snpdata.sh
# Post-filter GEDI snpdata to match STAR binomial mask criteria
#
# Usage: ./filter_gedi_snpdata.sh <gedi.snpdata> <p_err> <output.bed> [star_mask.bed.gz]
#
# STAR binomial filters mirrored:
#   - minCov = 20
#   - minAlt = 3
#   - minTcRatio = 0.3
#   - pval < 0.001 (using p_err)
#   - T→C only (GEDI SLAM mode already T→C focused)
#
# Note: GEDI doesn't provide MAPQ/BaseQ/junction/indel flank info,
#       so these cannot be matched (residual mismatch source)

set -e

GEDI_SNPDATA="${1:?Usage: $0 <gedi.snpdata> <p_err> <output.bed> [star_mask.bed.gz]}"
P_ERR="${2:?Missing p_err (e.g., 0.001)}"
OUTPUT_BED="${3:?Missing output BED path}"
STAR_BED="${4:-}"  # Optional: STAR mask for overlap comparison

# Filter thresholds (matching STAR binomial defaults)
MIN_COV=20
MIN_ALT=3
MIN_RATIO=0.3
MAX_PVAL=0.001

echo "=== GEDI SNP Post-Filter (STAR-matching) ==="
echo "Input: $GEDI_SNPDATA"
echo "p_err: $P_ERR"
echo "Filters: cov>=$MIN_COV, alt>=$MIN_ALT, ratio>=$MIN_RATIO, pval<$MAX_PVAL"
echo ""

# Step 1: Basic awk filter (cov, alt, ratio)
# Then Python for binomial p-value recalculation
python3 << PYFILTER
import sys
from scipy import stats
import math

# Parameters
p_err = float("$P_ERR")
min_cov = $MIN_COV
min_alt = $MIN_ALT
min_ratio = $MIN_RATIO
max_pval = $MAX_PVAL

input_file = "$GEDI_SNPDATA"
output_file = "$OUTPUT_BED"

def binom_pval(n, k, p):
    """Upper tail binomial p-value: P[X >= k | n, p]"""
    if n == 0 or k == 0:
        return 1.0
    if k > n:
        return 0.0
    # Use survival function (1 - CDF(k-1))
    return stats.binom.sf(k - 1, n, p)

passed = 0
total = 0

with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
    # Write BED header
    fout.write("#chrom\tstart\tend\tcov\talt\tratio\tgedi_pval\tstar_pval\n")
    
    for line in fin:
        if line.startswith("Location"):
            continue  # Skip header
        
        parts = line.strip().split('\t')
        if len(parts) < 4:
            continue
        
        total += 1
        
        # Parse GEDI columns
        # IMPORTANT: GEDI `*.snpdata` Location is 0-based (confirmed by instrumentation of SlamDetectSnps)
        loc = parts[0]  # chr:pos (0-based)
        cov = float(parts[1])
        alt = float(parts[2])
        gedi_pval = float(parts[3])
        
        # Apply filters
        if cov < min_cov:
            continue
        if alt < min_alt:
            continue
        
        ratio = alt / cov if cov > 0 else 0
        if ratio < min_ratio:
            continue
        
        # Compute STAR-style binomial p-value
        n = int(cov)
        k = int(alt)
        star_pval = binom_pval(n, k, p_err)
        
        if star_pval >= max_pval:
            continue
        
        # Convert to BED format (0-based)
        chrom_pos = loc.split(':')
        chrom = "chr" + chrom_pos[0] if not chrom_pos[0].startswith("chr") else chrom_pos[0]
        pos0 = int(chrom_pos[1])  # 0-based
        
        passed += 1
        fout.write(f"{chrom}\t{pos0}\t{pos0+1}\t{n}\t{k}\t{ratio:.4f}\t{gedi_pval:.2e}\t{star_pval:.2e}\n")

print(f"Total GEDI sites: {total}")
print(f"Passed all filters: {passed}")
print(f"Output: {output_file}")
PYFILTER

echo ""

# Step 2: Sort the output BED
sort -k1,1 -k2,2n "$OUTPUT_BED" -o "$OUTPUT_BED"
echo "Sorted output BED"

# Step 3: Compare with STAR mask if provided
if [ -n "$STAR_BED" ] && [ -f "$STAR_BED" ]; then
    echo ""
    echo "=== Overlap with STAR binomial mask ==="
    
    # Extract STAR positions (handle gzipped or plain, skip header)
    STAR_TMP=$(mktemp)
    if [[ "$STAR_BED" == *.gz ]]; then
        zcat "$STAR_BED" | grep -v "^#" | cut -f1-3 | sort -k1,1 -k2,2n > "$STAR_TMP"
    else
        grep -v "^#" "$STAR_BED" | cut -f1-3 | sort -k1,1 -k2,2n > "$STAR_TMP"
    fi
    
    # Extract GEDI positions
    GEDI_TMP=$(mktemp)
    tail -n +2 "$OUTPUT_BED" | cut -f1-3 | sort -k1,1 -k2,2n > "$GEDI_TMP"
    
    STAR_COUNT=$(wc -l < "$STAR_TMP")
    GEDI_COUNT=$(wc -l < "$GEDI_TMP")
    
    echo "STAR sites: $STAR_COUNT"
    echo "GEDI filtered sites: $GEDI_COUNT"
    
    if [ "$STAR_COUNT" -gt 0 ] && [ "$GEDI_COUNT" -gt 0 ]; then
        OVERLAP=$(bedtools intersect -a "$STAR_TMP" -b "$GEDI_TMP" -u | wc -l)
        STAR_IN_GEDI=$OVERLAP
        GEDI_IN_STAR=$(bedtools intersect -a "$GEDI_TMP" -b "$STAR_TMP" -u | wc -l)
        
        UNION=$((STAR_COUNT + GEDI_COUNT - OVERLAP))
        JACCARD=$(echo "scale=4; $OVERLAP / $UNION" | bc)
        
        echo ""
        echo "Overlap: $OVERLAP"
        echo "STAR in GEDI: $STAR_IN_GEDI / $STAR_COUNT ($(echo "scale=1; 100*$STAR_IN_GEDI/$STAR_COUNT" | bc)%)"
        echo "GEDI in STAR: $GEDI_IN_STAR / $GEDI_COUNT ($(echo "scale=1; 100*$GEDI_IN_STAR/$GEDI_COUNT" | bc)%)"
        echo "Jaccard: $JACCARD"
        
        echo ""
        echo "=== STAR-only sites (sample) ==="
        bedtools intersect -a "$STAR_TMP" -b "$GEDI_TMP" -v | head -5
        
        echo ""
        echo "=== GEDI-only sites (sample) ==="
        bedtools intersect -a "$GEDI_TMP" -b "$STAR_TMP" -v | head -5
    fi
    
    rm -f "$STAR_TMP" "$GEDI_TMP"
fi

echo ""
echo "=== Done ==="
