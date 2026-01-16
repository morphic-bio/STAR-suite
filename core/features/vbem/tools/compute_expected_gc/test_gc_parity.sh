#!/bin/bash
# Test GC bias expected distribution parity with Salmon

set -e

# Test data paths
GENOME_FASTA="/tmp/nfcore_ec_parity_test/transcriptome.fasta"  # Using transcriptome as "genome" for simplicity
GTF_FILE=""  # We'll need a GTF file or use transcriptome directly
TRANSCRIPTOME="/tmp/nfcore_ec_parity_test/transcriptome.fasta"
SALMON_OUT="/tmp/salmon_gc_test"
OUR_OUT="/tmp/our_gc_test"

echo "=== GC Bias Parity Test ==="
echo ""

# Check if test data exists
if [ ! -f "$TRANSCRIPTOME" ]; then
    echo "Error: Test transcriptome not found at $TRANSCRIPTOME"
    echo "Please run the EC filter parity test first to generate test data"
    exit 1
fi

# Build our tool if needed
cd "$(dirname "$0")"
if [ ! -f "./compute_expected_gc" ]; then
    echo "Building compute_expected_gc..."
    make
fi

echo "Step 1: Running Salmon with GC bias correction..."
echo "----------------------------------------"
# Run Salmon with GC bias to generate expected GC
# Create output directories
mkdir -p "$SALMON_OUT"
mkdir -p "$OUR_OUT"

SALMON_DETERMINISTIC_EM=1 \
/mnt/pikachu/salmon/build/src/salmon quant \
    -t "$TRANSCRIPTOME" \
    -l A \
    -a /tmp/nfcore_ec_parity_test/star_align/Aligned.toTranscriptome.out.bam \
    --gcBias \
    --dumpEqWeights \
    --noEffectiveLengthCorrection \
    --noFragLengthDist \
    -p 1 \
    -o "$SALMON_OUT" 2>&1 | tail -10

# Check if Salmon created expected GC file
if [ -f "$SALMON_OUT/aux_info/exp_gc.gz" ]; then
    echo "✓ Salmon created exp_gc.gz"
    gunzip -c "$SALMON_OUT/aux_info/exp_gc.gz" > "$SALMON_OUT/exp_gc.bin"
    echo "  (Note: Salmon uses binary GCFragModel format, we'll need to parse it)"
else
    echo "✗ Salmon did not create exp_gc.gz"
    echo "  This might mean GC bias wasn't computed or written"
fi

echo ""
echo "Step 2: Running our compute_expected_gc tool..."
echo "----------------------------------------"
./compute_expected_gc \
    --transcriptome "$TRANSCRIPTOME" \
    --output "$OUR_OUT/expected_gc.tsv" 2>&1

if [ -f "$OUR_OUT/expected_gc.tsv" ]; then
    echo "✓ Our tool created expected_gc.tsv"
    echo ""
    echo "First 10 lines of our output:"
    head -10 "$OUR_OUT/expected_gc.tsv"
    echo ""
    echo "Last 10 lines:"
    tail -10 "$OUR_OUT/expected_gc.tsv"
    echo ""
    echo "Total lines (should be 101):"
    wc -l "$OUR_OUT/expected_gc.tsv"
else
    echo "✗ Our tool did not create expected_gc.tsv"
    exit 1
fi

echo ""
echo "Step 3: Comparison"
echo "----------------------------------------"

# Check if probabilities sum to ~1.0
SUM=$(awk '{sum+=$2} END {print sum}' "$OUR_OUT/expected_gc.tsv")
echo "Sum of our probabilities: $SUM (should be ~1.0)"

if (( $(echo "$SUM > 0.99 && $SUM < 1.01" | bc -l) )); then
    echo "✓ Our probabilities sum correctly"
else
    echo "✗ Our probabilities don't sum to ~1.0"
fi

# Compare with Salmon if available
if [ -f "$SALMON_OUT/aux_info/exp_gc.gz" ]; then
    echo ""
    echo "Comparing with Salmon's expected GC..."
    python3 "$(dirname "$0")/compare_with_salmon.py" \
        "$SALMON_OUT/aux_info/exp_gc.gz" \
        "$OUR_OUT/expected_gc.tsv"
    
    if [ $? -eq 0 ]; then
        echo ""
        echo "✓ Parity check PASSED"
    else
        echo ""
        echo "✗ Parity check FAILED - distributions differ"
    fi
else
    echo ""
    echo "Note: Salmon's exp_gc.gz not found - skipping comparison"
    echo "      (Salmon may not have computed expected GC if --gcBias wasn't enabled)"
fi

echo ""
echo "=== Test Complete ==="
echo "Output files:"
echo "  Salmon: $SALMON_OUT/aux_info/exp_gc.gz"
echo "  Ours:   $OUR_OUT/expected_gc.tsv"
