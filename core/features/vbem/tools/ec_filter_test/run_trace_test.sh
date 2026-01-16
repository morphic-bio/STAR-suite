#!/bin/bash
set -e

BAM="/mnt/pikachu/test_jax/small/filtered/TWIST2_PTC_H09_GT24-03507_ATTGTGAA-GCAATGCA_S20_L001_R1_001_val/Aligned.toTranscriptome.out.bam"
TRANSCRIPTOME=$(find /mnt/pikachu -name "transcriptome.fa" -type f 2>/dev/null | head -1)
SALMON_BIN="/mnt/pikachu/salmon/build_trace/src/salmon"

if [ ! -f "$BAM" ]; then
    echo "BAM file not found: $BAM"
    exit 1
fi

if [ ! -f "$TRANSCRIPTOME" ]; then
    echo "Transcriptome file not found"
    exit 1
fi

if [ ! -f "$SALMON_BIN" ]; then
    echo "Salmon binary not found: $SALMON_BIN"
    exit 1
fi

echo "=== Running Salmon with trace ==="
echo "BAM: $BAM"
echo "Transcriptome: $TRANSCRIPTOME"
echo "Salmon: $SALMON_BIN"

rm -f /tmp/salmon_trace_run.txt
export SALMON_TRACE_FILE="/tmp/salmon_trace_run.txt"
export SALMON_TRACE_LIMIT="10000"

"$SALMON_BIN" quant \
    -t "$TRANSCRIPTOME" \
    -l A \
    -a "$BAM" \
    --dumpEqWeights \
    --noLengthCorrection \
    --noEffectiveLengthCorrection \
    --noFragLengthDist \
    --incompatPrior 0 \
    --scoreExp 1.0 \
    -o /tmp/salmon_trace_out_run 2>&1 | tail -10

echo ""
echo "=== Checking trace file ==="
if [ -f /tmp/salmon_trace_run.txt ]; then
    echo "Trace file created:"
    ls -lh /tmp/salmon_trace_run.txt
    echo ""
    echo "Line count: $(wc -l < /tmp/salmon_trace_run.txt)"
    echo ""
    echo "First 3 lines:"
    head -3 /tmp/salmon_trace_run.txt
    echo ""
    echo "Lines matching our format:"
    grep -E "^LH00341.*txpIDs=" /tmp/salmon_trace_run.txt | head -2 || echo "No matching lines found"
else
    echo "Trace file not created!"
fi
