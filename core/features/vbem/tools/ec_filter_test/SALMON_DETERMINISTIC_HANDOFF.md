# Salmon Agent: Implement Deterministic Error Model Updates

## Background

We're debugging error model parity between our CLI and Salmon. We've identified that **stochastic sampling** during error model updates is the primary remaining source of divergence.

Both implementations randomly sample which alignments to use for updating the transition matrices. Even with identical read order, different random numbers select different alignments, causing matrix divergence.

## What We've Done (STAR-Flex CLI)

1. ✅ Added `--read-order-file` to process reads in Salmon's order
2. ✅ Added `--deterministic-updates` to update ALL alignments (no sampling)
3. ✅ Verified pre-burnin reads (first 5000) now match exactly

## What Salmon Needs

Add a flag or environment variable to disable stochastic sampling in the error model update.

### Location

File: `/mnt/pikachu/salmon/src/SalmonQuantifyAlignments.cpp`

Around line ~1200 (in the `batchReads` processing section), find:

```cpp
double r = uni(eng);
if (!burnedIn and r < std::exp(aln->logProb)) {
    /**
     * Update the bias sequence-specific bias model
     **/

    // Update the error model
    if (!useASWithoutCIGAR and salmonOpts.useErrorModel) {
        auto alignerScore = getAlignerAssignedScore(aln);
        alnMod.update(*aln, *alnPrimary, transcript, alignerScore, logForgettingMass);
    }
}
```

### Proposed Change

Option A: Environment variable (minimal change):

```cpp
static bool deterministicEM = (getenv("SALMON_DETERMINISTIC_EM") != nullptr);
double r = uni(eng);
if (!burnedIn and (deterministicEM || r < std::exp(aln->logProb))) {
    /**
     * Update the bias sequence-specific bias model
     **/

    // Update the error model
    if (!useASWithoutCIGAR and salmonOpts.useErrorModel) {
        auto alignerScore = getAlignerAssignedScore(aln);
        alnMod.update(*aln, *alnPrimary, transcript, alignerScore, logForgettingMass);
    }
}
```

Option B: Add to SalmonOpts (more robust):

1. In `SalmonOpts.hpp`, add:
```cpp
bool deterministicErrorModel{false};
```

2. In `ProgramOptionsGenerator.cpp`, add:
```cpp
("deterministicErrorModel",
 po::bool_switch(&(sopt.deterministicErrorModel))->default_value(false),
 "Update error model with all alignments (disable stochastic sampling)")
```

3. In `SalmonQuantifyAlignments.cpp`:
```cpp
double r = uni(eng);
if (!burnedIn and (salmonOpts.deterministicErrorModel || r < std::exp(aln->logProb))) {
    // ... error model update ...
}
```

## Testing

### Build and Run

```bash
cd /mnt/pikachu/salmon/build
make -j4

# Test with deterministic mode
SALMON_DETERMINISTIC_EM=1 \
SALMON_ERROR_MODEL_TRACE=/tmp/salmon_determ_trace.txt \
SALMON_TRACE_LEVEL=1 \
SALMON_DUMP_MATRICES=/tmp/salmon_determ_matrices \
/mnt/pikachu/salmon/build/src/salmon quant \
    -t /tmp/nfcore_ec_parity_test/transcriptome.fasta \
    -l A \
    -a /tmp/nfcore_ec_parity_test/star_align/Aligned.toTranscriptome.out.bam \
    --dumpEqWeights --noLengthCorrection --noEffectiveLengthCorrection \
    --noFragLengthDist -p 1 -o /tmp/salmon_determ_out
```

### Expected Outcome

1. Trace file should have same number of reads (29726)
2. All reads should have `modelUsed=1` after position 5000 (the 101 errLikeSum=0 bug should also be fixed if not already)
3. Matrix files should be generated at checkpoints

### Verification

After implementing, we can compare:

```bash
# Compare traces
diff /tmp/salmon_determ_trace.txt /tmp/deterministic_test/our_trace.txt

# If traces match, matrices should also match
```

## Current State

| Metric | Our CLI (deterministic) | Salmon (stochastic) | Expected After Fix |
|--------|------------------------|--------------------|--------------------|
| Pre-burnin errLike | 0 (correct) | 0 (correct) | Match |
| Post-burnin errLike | Computed | Computed | Should match |
| Matrix mean diff | 0.086 | - | < 0.001 |

## Files for Reference

- Our deterministic trace: `/tmp/deterministic_test/our_trace.txt`
- Our deterministic matrices: `/tmp/deterministic_test/our_matrices_*.tsv`
- Salmon's current trace: `/tmp/salmon_debug_trace.txt`
- Read order file: `/tmp/salmon_read_order.txt`

## Summary

1. Add environment variable `SALMON_DETERMINISTIC_EM` check
2. If set, always update error model (skip `r < exp(logProb)` check)
3. Rebuild and test
4. Compare traces - they should now match exactly

This is the final piece needed for error model parity verification.
