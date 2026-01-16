# Debug Tracing Guide

This guide explains how to use the debug tracing feature to compare em_quant with Salmon and identify where divergence occurs.

## Overview

Debug tracing logs detailed per-iteration and per-EC values for selected transcripts, allowing you to pinpoint exactly where em_quant and Salmon diverge.

## Enabling Tracing

### 1. Set Transcript IDs

Set the `EM_QUANT_DEBUG_TXPS` environment variable with comma-separated transcript IDs:

```bash
export EM_QUANT_DEBUG_TXPS=ENST00000465752,ENST00000464034,ENST00000488272
```

### 2. Run em_quant with Tracing

```bash
em_quant --vb --threads 1 \
    -e eq_classes.txt \
    -l quant.sf \
    -o output.tsv \
    --debug-trace /tmp/em_trace.txt \
    -v
```

## Trace File Format

### Header
```
iter	transcript	alpha	logNorm	expTheta	expected_count
# EC-level trace: EC	iter	ec_id	transcript	denom	expTheta	aux	contribution
```

### Per-Iteration Lines
```
-1	ENST00000465752	1.01	7.51817	0.000309963	0.0
0	ENST00000465752	2.88659	8.51603	0.000382601	2.87659
1	ENST00000465752	2.27027	8.51603	0.000358948	2.26027
```

Fields:
- `iter`: Iteration number (-1 for initial state)
- `transcript`: Transcript name
- `alpha`: Current alpha value (prior + expected_counts)
- `logNorm`: digamma(sum of all alphas)
- `expTheta`: exp(digamma(alpha) - logNorm)
- `expected_count`: Expected count from E-step

### EC-Level Lines
```
EC	0	45	ENST00000465752	0.0104374	0.000309963	0.333333	0.00989909
EC	0	71	ENST00000465752	0.00138222	0.000309963	0.5	1.45763
EC	0	208	ENST00000465752	SINGLE	1	1	1
```

Fields:
- `EC`: Marker for EC-level trace
- `iter`: Iteration number
- `ec_id`: EC index (0-based)
- `transcript`: Transcript name
- `denom`: Denominator (Σ expTheta[tid] * aux[i])
- `expTheta`: expTheta value for this transcript
- `aux`: Auxiliary weight (from EC weights or 1/effLen)
- `contribution`: Per-transcript contribution to expected_counts
- For single-transcript ECs: Shows "SINGLE" with full count

## Helper Scripts

### Validate Trace File

```bash
python3 validate_trace.py /tmp/em_trace.txt
```

Checks:
- Format correctness
- Header presence
- Consistency of logNorm across transcripts
- EC contributions sum to expected_count

### Extract Specific Data

```bash
# Extract iteration 0 for a specific transcript
python3 extract_trace.py /tmp/em_trace.txt --iter 0 --transcript ENST00000465752

# Extract all ECs for a specific transcript
python3 extract_trace.py /tmp/em_trace.txt --transcript ENST00000465752

# Extract a specific EC across all iterations
python3 extract_trace.py /tmp/em_trace.txt --ec 45
```

### Compare with Salmon Trace

```bash
python3 compare_traces.py /tmp/em_trace.txt /tmp/salmon_trace.txt

# Focus on specific transcript
python3 compare_traces.py /tmp/em_trace.txt /tmp/salmon_trace.txt \
    --transcript ENST00000465752

# Focus on specific iteration
python3 compare_traces.py /tmp/em_trace.txt /tmp/salmon_trace.txt \
    --iter 0

# Adjust tolerance
python3 compare_traces.py /tmp/em_trace.txt /tmp/salmon_trace.txt \
    --tolerance 1e-12
```

## Adding Salmon Instrumentation

To compare traces, you need to add similar logging to Salmon. Add instrumentation in:

**File**: `/mnt/pikachu/salmon/src/CollapsedEMOptimizer.cpp`  
**Location**: Line 133 (VBEMUpdate_ loop)

### What to Log

1. **Per-iteration** (after M-step, before convergence check):
   ```cpp
   // For each debug transcript:
   debug_stream << iter << "\t" << transcript_name << "\t"
                << alpha[i] << "\t" << logNorm << "\t"
                << expTheta[i] << "\t" << expected_counts[i] << "\n";
   ```

2. **Per-EC** (inside VBEMUpdate_ loop):
   ```cpp
   // For each EC containing debug transcripts:
   debug_stream << "EC\t" << iter << "\t" << ec_idx << "\t"
                << transcript_name << "\t" << denom << "\t"
                << expTheta[tid] << "\t" << aux << "\t"
                << contribution << "\n";
   ```

### Salmon Command (after adding flags)

```bash
salmon quant \
    -i index \
    -l A \
    -1 reads_1.fq -2 reads_2.fq \
    --dumpEqWeights \
    --threads 1 \
    -o salmon_out \
    --debugTrace /tmp/salmon_trace.txt \
    --debugTranscripts ENST00000465752,ENST00000464034
```

## Finding Divergence

### Step 1: Generate Both Traces

```bash
# em_quant
export EM_QUANT_DEBUG_TXPS=ENST00000465752
em_quant --vb --threads 1 \
    -e eq_classes.txt -l quant.sf \
    -o em_output.tsv \
    --debug-trace /tmp/em_trace.txt

# Salmon (after adding instrumentation)
salmon quant -i index -l A -1 r1.fq -2 r2.fq \
    --dumpEqWeights --threads 1 \
    -o salmon_out \
    --debugTrace /tmp/salmon_trace.txt \
    --debugTranscripts ENST00000465752
```

### Step 2: Compare Traces

```bash
python3 compare_traces.py /tmp/em_trace.txt /tmp/salmon_trace.txt
```

This will show:
- First divergence in per-iteration values
- First divergence in EC-level details
- Field-by-field differences

### Step 3: Analyze Divergence

The comparison will identify:
- **Weights parsing/order**: If `aux` values differ
- **Digamma**: If `logNorm` or `expTheta` differ
- **Truncation**: If values differ after many iterations
- **Convergence**: If final values differ

### Step 4: Inspect Specific Iteration/EC

```bash
# Extract iteration 0 from both traces
python3 extract_trace.py /tmp/em_trace.txt --iter 0 > em_iter0.txt
python3 extract_trace.py /tmp/salmon_trace.txt --iter 0 > salmon_iter0.txt

# Compare side-by-side
diff -y em_iter0.txt salmon_iter0.txt
```

## Example: Tracing an Outlier

For transcript `ENST00000465752` (max outlier, 10.11% difference):

```bash
export EM_QUANT_DEBUG_TXPS=ENST00000465752
em_quant --vb --threads 1 \
    -e ../../test/fixtures/salmon_eq/eq_classes.txt \
    -l ../../test/fixtures/salmon_eq/quant.sf \
    -o /tmp/test.tsv \
    --debug-trace /tmp/trace.txt

# Check initial state
python3 extract_trace.py /tmp/trace.txt --iter -1

# Check first iteration
python3 extract_trace.py /tmp/trace.txt --iter 0

# Check all ECs for this transcript
python3 extract_trace.py /tmp/trace.txt --transcript ENST00000465752 | grep "^EC"
```

The trace shows:
- Initial: alpha=1.01 (prior + unique_count=1)
- Iteration 0: 5 ECs contribute (4 multi-transcript, 1 single-transcript)
- Final: alpha=1.77278, expected_count=1.76278
- Salmon: NumReads=1.961

Comparing EC-level details will show where the divergence occurs.

## Troubleshooting

### No EC lines in trace

- Check that transcript IDs match exactly (case-sensitive)
- Verify transcript appears in multi-transcript ECs (single-transcript ECs are logged differently)

### Inconsistent logNorm values

- logNorm should be the same for all transcripts in the same iteration
- If different, check digamma computation

### EC contributions don't sum to expected_count

- May be due to rounding or ECs skipped due to minEQClassWeight
- Check if difference is within floating-point precision
