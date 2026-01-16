# Error Model Harmonization Testing Guide

## Overview

This guide describes how to test and verify error model harmonization between our CLI and Salmon.

## Prerequisites

1. **Our CLI**: Built with `--use-error-model` support
2. **Salmon**: Modified with trace/matrix dump support (see `SALMON_INSTRUMENTATION_GUIDE.md`)
3. **Test data**: BAM file + transcriptome FASTA

## Test Workflow

### Stage 1: Baseline Initialization Test

**Goal**: Verify both models start with identical matrices.

```bash
# Our CLI - dump initial matrices (before processing any reads)
# Note: This requires adding a dump call before any reads are processed
./ec_filter_cli --input test.bam --transcripts transcriptome.fa \
    --use-error-model --dump-matrices /tmp/our_init \
    --output /tmp/our_ec.txt

# Compare initial matrices (should be identical)
python3 compare_matrices.py \
    --ours /tmp/our_init_final_left_bin0.tsv \
    --salmon /tmp/salmon_init_final_left_bin0.tsv \
    --tolerance 1e-10
```

**Expected**: 100% match (all cells within 1e-10)

### Stage 2: Single-Read Update Test

**Goal**: Process exactly 1 read and compare matrices.

```bash
# Create a BAM with only 1 read (use samtools view -s)
samtools view -s 0.00001 test.bam | samtools view -Sb - > single_read.bam

# Our CLI
./ec_filter_cli --input single_read.bam --transcripts transcriptome.fa \
    --use-error-model --dump-matrices /tmp/our_single \
    --output /tmp/our_ec.txt

# Salmon (single-threaded)
SALMON_DUMP_MATRICES=/tmp/salmon_single \
salmon quant -t transcriptome.fa -l A -a single_read.bam \
    --dumpEqWeights --noLengthCorrection --noEffectiveLengthCorrection \
    --noFragLengthDist -p 1 -o /tmp/salmon_out

# Compare matrices
python3 compare_matrices.py \
    --ours /tmp/our_single_final_left_bin0.tsv \
    --salmon /tmp/salmon_single_final_left_bin0.tsv
```

**Expected**: Identical matrices if CIGAR parsing and state indexing match

### Stage 3: Pre-Burnin Phase Test (0-5000 reads)

**Goal**: Verify matrices evolve identically during training phase.

```bash
# Limit to first 5000 reads
head -n 10000 test.bam | samtools view -Sb - > first_5k.bam

# Our CLI with checkpoint dumps
./ec_filter_cli --input first_5k.bam --transcripts transcriptome.fa \
    --use-error-model --dump-matrices /tmp/our_pre5k \
    --output /tmp/our_ec.txt

# Salmon with checkpoint dumps
SALMON_DUMP_MATRICES=/tmp/salmon_pre5k \
salmon quant -t transcriptome.fa -l A -a first_5k.bam \
    --dumpEqWeights --noLengthCorrection --noEffectiveLengthCorrection \
    --noFragLengthDist -p 1 -o /tmp/salmon_out

# Compare at each checkpoint
for checkpoint in 1000 2000 3000 4000 5000; do
    echo "Comparing checkpoint $checkpoint:"
    python3 compare_matrices.py \
        --ours /tmp/our_pre5k_pre${checkpoint}_left_bin0.tsv \
        --salmon /tmp/salmon_pre5k_pre${checkpoint}_left_bin0.tsv \
        --tolerance 1e-10
done
```

**Expected**: Matrices should match at each checkpoint (within 1e-10)

### Stage 4: Post Pre-Burnin Test (5000+ reads)

**Goal**: Verify errLike computation matches after model starts being used.

```bash
# Use first 10000 reads
head -n 20000 test.bam | samtools view -Sb - > first_10k.bam

# Our CLI with tracing
./ec_filter_cli --input first_10k.bam --transcripts transcriptome.fa \
    --use-error-model --error-model-trace /tmp/our_trace.txt \
    --trace-level 2 --output /tmp/our_ec.txt

# Salmon with tracing
SALMON_ERROR_MODEL_TRACE=/tmp/salmon_trace.txt \
SALMON_TRACE_LEVEL=2 \
salmon quant -t transcriptome.fa -l A -a first_10k.bam \
    --dumpEqWeights --noLengthCorrection --noEffectiveLengthCorrection \
    --noFragLengthDist -p 1 -o /tmp/salmon_out

# Compare traces
python3 compare_error_model.py \
    --ours /tmp/our_trace.txt \
    --salmon /tmp/salmon_trace.txt \
    --tolerance 1e-6 \
    --report /tmp/trace_comparison.txt
```

**Expected**: 
- Read 5001+ should have identical errLike values (within 1e-6)
- Reads < 5000 should have errLike=0.0 in both

### Stage 5: Full Run Comparison

**Goal**: Complete comparison with all instrumentation.

```bash
# Our CLI
./ec_filter_cli --input test.bam --transcripts transcriptome.fa \
    --use-error-model \
    --error-model-trace /tmp/our_full_trace.txt \
    --trace-level 2 \
    --dump-matrices /tmp/our_full \
    --output /tmp/our_ec.txt

# Salmon
SALMON_ERROR_MODEL_TRACE=/tmp/salmon_full_trace.txt \
SALMON_TRACE_LEVEL=2 \
SALMON_DUMP_MATRICES=/tmp/salmon_full \
salmon quant -t transcriptome.fa -l A -a test.bam \
    --dumpEqWeights --noLengthCorrection --noEffectiveLengthCorrection \
    --noFragLengthDist -p 1 -o /tmp/salmon_out

# Compare everything
python3 compare_error_model.py \
    --ours /tmp/our_full_trace.txt \
    --salmon /tmp/salmon_full_trace.txt \
    --tolerance 1e-6 \
    --report /tmp/full_comparison.txt

# Compare final matrices
python3 compare_matrices.py \
    --ours /tmp/our_full_final_left_bin0.tsv \
    --salmon /tmp/salmon_full_final_left_bin0.tsv \
    --tolerance 1e-10
```

## Known Differences to Verify

1. **Left/Right Read Determination**: We use `<=`, Salmon uses `<` when positions are equal
2. **Alpha Default**: Verify Salmon's default alpha value
3. **Read Bins**: Verify Salmon uses 6 bins by default (we use 6 in CLI)

## Troubleshooting

### Matrices diverge early (< 100 reads)
- Check CIGAR parsing logic
- Verify state indexing formula: `curStateIdx = curRefBase * 9 + curReadBase`
- Check two-bit encoding matches

### Matrices match but errLike differs
- Check likelihood computation (fg - bg)
- Verify transition probability lookup
- Check normalization (row sums)

### errLike matches but weights differ
- Check auxProb combination: `logFragProb + errLike + logAlignCompatProb`
- Verify normalization (auxDenom computation)
- Check weight assignment to ECs

## Success Criteria

- ✅ Transition matrices match within 1e-10 at all checkpoints
- ✅ errLike values match within 1e-6 for all alignments
- ✅ Final EC weights achieve 99%+ parity with error model enabled
