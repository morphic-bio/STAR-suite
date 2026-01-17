#!/bin/bash
# Analyze STAR SNP-Mask Parameter Sweep Results

RESULTS_FILE="test/tmp_wt_snp_comparison/sweep_results.tsv"

if [[ ! -f "$RESULTS_FILE" ]]; then
    echo "ERROR: Results file not found: $RESULTS_FILE"
    echo "Run tests/run_star_snp_mask_sweep.sh first"
    exit 1
fi

echo "=== STAR SNP-Mask Parameter Sweep Analysis ==="
echo ""

# Show full table
echo "Full Results Table:"
echo "==================="
column -t -s $'\t' "$RESULTS_FILE"
echo ""

# Find best Jaccard
BEST_JACCARD=$(tail -n +2 "$RESULTS_FILE" | awk -F'\t' '{print $9}' | sort -rn | head -1)
BEST_LINE=$(tail -n +2 "$RESULTS_FILE" | awk -F'\t' -v best="$BEST_JACCARD" '$9==best')

echo "=== Best Jaccard: $BEST_JACCARD ==="
echo "Parameters:"
echo "$BEST_LINE" | awk -F'\t' '{print "  MinCov: " $1 ", Posterior: " $2 ", MinAlt: " $3}'
echo "$BEST_LINE" | awk -F'\t' '{print "  STAR T→C cov≥20: " $5 ", Overlap: " $6 ", Jaccard: " $9}'
echo ""

# Find best GEDI-in-STAR percentage
BEST_GEDI=$(tail -n +2 "$RESULTS_FILE" | awk -F'\t' '{print $8}' | sort -rn | head -1)
BEST_GEDI_LINE=$(tail -n +2 "$RESULTS_FILE" | awk -F'\t' -v best="$BEST_GEDI" '$8==best')

echo "=== Best GEDI Capture: ${BEST_GEDI}% ==="
echo "Parameters:"
echo "$BEST_GEDI_LINE" | awk -F'\t' '{print "  MinCov: " $1 ", Posterior: " $2 ", MinAlt: " $3}'
echo "$BEST_GEDI_LINE" | awk -F'\t' '{print "  STAR T→C cov≥20: " $5 ", GEDI in STAR: " $8 "%"}'
echo ""

# Summary by parameter
echo "=== Summary by Parameter ==="
echo ""
echo "Effect of MinCov:"
tail -n +2 "$RESULTS_FILE" | awk -F'\t' '{
    cov[$1] += $9
    count[$1]++
}
END {
    for (c in cov) {
        printf "  MinCov %s: avg Jaccard = %.4f (n=%d)\n", c, cov[c]/count[c], count[c]
    }
}' | sort

echo ""
echo "Effect of Posterior:"
tail -n +2 "$RESULTS_FILE" | awk -F'\t' '{
    post[$2] += $9
    count[$2]++
}
END {
    for (p in post) {
        printf "  Posterior %s: avg Jaccard = %.4f (n=%d)\n", p, post[p]/count[p], count[p]
    }
}' | sort -rn

echo ""
echo "Effect of MinAlt:"
tail -n +2 "$RESULTS_FILE" | awk -F'\t' '{
    alt[$3] += $9
    count[$3]++
}
END {
    for (a in alt) {
        printf "  MinAlt %s: avg Jaccard = %.4f (n=%d)\n", a, alt[a]/count[a], count[a]
    }
}' | sort -rn
