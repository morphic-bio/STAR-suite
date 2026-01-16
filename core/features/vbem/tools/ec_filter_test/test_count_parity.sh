#!/bin/bash
# Test script to isolate count parity issues

set -e

BAM="/tmp/nfcore_ec_parity_test/star_align/Aligned.toTranscriptome.out.bam"
FASTA="/tmp/nfcore_ec_parity_test/transcriptome.fasta"
OUT_DIR="/tmp/count_parity_test"

mkdir -p "$OUT_DIR"

echo "=== Test 1: Baseline (current state) ==="
./ec_filter_cli --input "$BAM" --transcripts "$FASTA" \
    --score-exp 1.0 --range-factorization-bins 0 \
    --no-local-pruning --no-global-pruning \
    -o "$OUT_DIR/our_baseline.txt" 2>&1 | grep -E "Built|Zero-prob|reads"

echo ""
echo "=== Comparing with Salmon ==="
python3 compare_ecs.py \
    --salmon /tmp/salmon_no_err_1t_eq.txt \
    --ours "$OUT_DIR/our_baseline.txt" \
    --tolerance 1e-6 2>&1 | tail -5
