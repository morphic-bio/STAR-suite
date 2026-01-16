#!/bin/bash
# End-to-end GC bias parity test
# Compares our GC bias pipeline with Salmon's --gcBias output

set -e

TEST_DIR="${1:-/tmp/gc_bias_parity_test}"
BAM_FILE="${2:-/tmp/nfcore_ec_parity_test/star_align/Aligned.toTranscriptome.out.bam}"
TRANSCRIPTOME="${3:-/tmp/nfcore_ec_parity_test/transcriptome.fasta}"

if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file not found: $BAM_FILE"
    exit 1
fi

if [ ! -f "$TRANSCRIPTOME" ]; then
    echo "Error: Transcriptome not found: $TRANSCRIPTOME"
    exit 1
fi

echo "GC Bias Parity Test"
echo "=================="
echo "Test directory: $TEST_DIR"
echo "BAM file: $BAM_FILE"
echo "Transcriptome: $TRANSCRIPTOME"
echo ""

mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

# Step 1: Run Salmon with GC bias
echo "Step 1: Running Salmon with --gcBias..."
SALMON_OUT="$TEST_DIR/salmon_gc_out"
if [ ! -d "$SALMON_OUT" ] || [ ! -f "$SALMON_OUT/quant.sf" ]; then
    /mnt/pikachu/salmon/build/src/salmon quant \
        -t "$TRANSCRIPTOME" \
        -l A \
        -a "$BAM_FILE" \
        --gcBias \
        -p 1 \
        -o "$SALMON_OUT" 2>&1 | tail -10
else
    echo "  Salmon output already exists, skipping..."
fi

# Step 2: Sample FLD
echo ""
echo "Step 2: Sampling fragment length distribution..."
FLD_FILE="$TEST_DIR/fld.tsv"
if [ ! -f "$FLD_FILE" ]; then
    /mnt/pikachu/STAR-Flex/tools/sample_fld/sample_fld \
        --bam "$BAM_FILE" \
        --output "$FLD_FILE" 2>&1 | tail -3
else
    echo "  FLD file already exists, skipping..."
fi

# Step 3: Compute expected GC
echo ""
echo "Step 3: Computing expected GC distribution..."
EXPECTED_GC="$TEST_DIR/expected_gc.tsv"
if [ ! -f "$EXPECTED_GC" ]; then
    /mnt/pikachu/STAR-Flex/tools/compute_expected_gc/compute_expected_gc \
        --transcriptome "$TRANSCRIPTOME" \
        --fld "$FLD_FILE" \
        --output "$EXPECTED_GC" 2>&1 | tail -3
else
    echo "  Expected GC file already exists, skipping..."
fi

# Step 4: Collect observed GC
echo ""
echo "Step 4: Collecting observed GC from alignments..."
OBSERVED_GC="$TEST_DIR/observed_gc.tsv"
EC_FILE="$TEST_DIR/eq_classes.txt"
if [ ! -f "$OBSERVED_GC" ]; then
    /mnt/pikachu/STAR-Flex/tools/ec_filter_test/ec_filter_cli \
        --input "$BAM_FILE" \
        --transcripts "$TRANSCRIPTOME" \
        --gc-bias \
        --expected-gc "$EXPECTED_GC" \
        --observed-gc-out "$OBSERVED_GC" \
        --output "$EC_FILE" 2>&1 | tail -5
else
    echo "  Observed GC file already exists, skipping..."
fi

# Step 5: Compute effective lengths
echo ""
echo "Step 5: Computing GC-corrected effective lengths..."
EFF_LENGTHS="$TEST_DIR/effective_lengths.tsv"
if [ ! -f "$EFF_LENGTHS" ]; then
    /mnt/pikachu/STAR-Flex/tools/compute_gc_bias/compute_gc_bias \
        --expected-gc "$EXPECTED_GC" \
        --observed-gc "$OBSERVED_GC" \
        --fld "$FLD_FILE" \
        --transcriptome "$TRANSCRIPTOME" \
        --output "$EFF_LENGTHS" 2>&1 | tail -3
else
    echo "  Effective lengths file already exists, skipping..."
fi

# Step 6: Run EM/VB with GC-corrected effective lengths
echo ""
echo "Step 6: Running VB quantification with GC-corrected effective lengths..."
OUR_QUANT="$TEST_DIR/our_quant.sf"
if [ ! -f "$OUR_QUANT" ]; then
    /mnt/pikachu/STAR-Flex/tools/em_quant/em_quant \
        --ec "$EC_FILE" \
        --effective-lengths "$EFF_LENGTHS" \
        --vb \
        --output "$OUR_QUANT" 2>&1 | tail -5
else
    echo "  Quantification output already exists, skipping..."
fi

# Step 7: Compare results
echo ""
echo "Step 7: Comparing results..."
echo ""

# Create comparison script
COMPARE_SCRIPT="$TEST_DIR/compare_quant.py"
cat > "$COMPARE_SCRIPT" << 'PYEOF'
#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np

if len(sys.argv) < 3:
    print("Usage: compare_quant.py <salmon_quant.sf> <our_quant.sf>")
    sys.exit(1)

salmon_file = sys.argv[1]
our_file = sys.argv[2]

# Load Salmon output
salmon_df = pd.read_csv(salmon_file, sep='\t')
salmon_df = salmon_df.set_index('Name')

# Load our output
our_df = pd.read_csv(our_file, sep='\t')
our_df = our_df.set_index('Name')

# Merge on transcript name
merged = salmon_df.join(our_df, how='inner', rsuffix='_our')

print("Comparison Results:")
print("=" * 60)
print(f"Transcripts in Salmon: {len(salmon_df)}")
print(f"Transcripts in Ours: {len(our_df)}")
print(f"Transcripts in both: {len(merged)}")
print("")

# Compare effective lengths
if 'EffectiveLength' in merged.columns and 'EffectiveLength_our' in merged.columns:
    eff_len_corr = merged['EffectiveLength'].corr(merged['EffectiveLength_our'])
    eff_len_rmse = np.sqrt(((merged['EffectiveLength'] - merged['EffectiveLength_our'])**2).mean())
    eff_len_mean = merged['EffectiveLength'].mean()
    eff_len_rel_rmse = eff_len_rmse / eff_len_mean if eff_len_mean > 0 else 0
    
    print("Effective Length Comparison:")
    print(f"  Correlation: {eff_len_corr:.6f}")
    print(f"  RMSE: {eff_len_rmse:.2f}")
    print(f"  Relative RMSE: {eff_len_rel_rmse*100:.2f}%")
    print("")

# Compare TPM
if 'TPM' in merged.columns and 'TPM_our' in merged.columns:
    tpm_corr = merged['TPM'].corr(merged['TPM_our'])
    tpm_rmse = np.sqrt(((merged['TPM'] - merged['TPM_our'])**2).mean())
    
    print("TPM Comparison:")
    print(f"  Correlation: {tpm_corr:.6f}")
    print(f"  RMSE: {tpm_rmse:.6f}")
    print("")

# Compare NumReads
if 'NumReads' in merged.columns and 'NumReads_our' in merged.columns:
    numreads_corr = merged['NumReads'].corr(merged['NumReads_our'])
    numreads_rmse = np.sqrt(((merged['NumReads'] - merged['NumReads_our'])**2).mean())
    numreads_mean = merged['NumReads'].mean()
    numreads_rel_rmse = numreads_rmse / numreads_mean if numreads_mean > 0 else 0
    
    print("NumReads Comparison:")
    print(f"  Correlation: {numreads_corr:.6f}")
    print(f"  RMSE: {numreads_rmse:.6f}")
    print(f"  Relative RMSE: {numreads_rel_rmse*100:.2f}%")
    print("")

# Overall assessment
if eff_len_corr > 0.99 and tpm_corr > 0.95 and numreads_rel_rmse < 0.05:
    print("✅ PASS: Results match within tolerance")
    sys.exit(0)
else:
    print("⚠️  WARN: Some metrics outside tolerance")
    sys.exit(1)
PYEOF
chmod +x "$COMPARE_SCRIPT"

python3 "$COMPARE_SCRIPT" "$SALMON_OUT/quant.sf" "$OUR_QUANT"

echo ""
echo "Test complete. Results saved in: $TEST_DIR"
