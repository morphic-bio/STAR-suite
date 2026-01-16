# Error Model Parity Status

## Summary
We've identified the root causes of errLike differences between our CLI and Salmon. Full parity requires changes on BOTH sides.

## What We've Implemented

### 1. Read Order Matching (✓ Implemented)
```bash
--read-order-file /tmp/salmon_read_order.txt
```
This processes reads in the same order as Salmon.

**Impact**: Pre-burnin reads (first 5000) now match exactly.

### 2. Deterministic Updates (✓ Implemented)
```bash
--deterministic-updates
```
This disables stochastic sampling and updates ALL alignments.

**Impact**: Makes our matrices trained MORE than Salmon's (worse parity).

## What Salmon Needs

For full parity, Salmon must ALSO disable stochastic sampling:

### Option 1: Add Salmon flag `--deterministic-error-model`
Modify `SalmonQuantifyAlignments.cpp`:
```cpp
// Current code:
double r = uni(eng);
if (!burnedIn and r < std::exp(aln->logProb)) {
    // Update error model
}

// Change to:
if (!burnedIn) {  // Always update when deterministic
    // Update error model
}
```

### Option 2: Use Environment Variable
```cpp
bool deterministicEM = (getenv("SALMON_DETERMINISTIC_EM") != nullptr);
if (!burnedIn and (deterministicEM || r < std::exp(aln->logProb))) {
    // Update error model
}
```

## Current Test Results

| Configuration | Mean errLike diff | Notes |
|--------------|------------------|-------|
| Our CLI (original order) | 0.82 | Different read order + stochastic |
| Our CLI (Salmon order) | 0.83 | Same order, still stochastic |
| Our CLI (Salmon order + deterministic) | 3.13 | We update more, Salmon still samples |

## Path to Parity

1. **Step 1** ✓: Match read order → Pre-burnin matches
2. **Step 2** ✓: Disable our stochastic sampling → Updates all alignments
3. **Step 3** (TODO): Salmon disables stochastic sampling → Should match

## Files Needed for Salmon Agent

To implement deterministic mode in Salmon:

1. **Modify**: `/mnt/pikachu/salmon/src/SalmonQuantifyAlignments.cpp`
   - Line ~1200: Add deterministic flag check around error model update

2. **Test Command**:
```bash
SALMON_DETERMINISTIC_EM=1 salmon quant \
    -t transcriptome.fasta \
    -l A \
    -a aligned.bam \
    --dumpEqWeights --noLengthCorrection --noEffectiveLengthCorrection \
    --noFragLengthDist -p 1 -o salmon_out
```

## Verification Steps

Once Salmon has deterministic mode:

1. Run both with deterministic updates:
```bash
# Our CLI
ec_filter_cli --read-order-file salmon_order.txt --deterministic-updates ...

# Salmon
SALMON_DETERMINISTIC_EM=1 salmon quant ...
```

2. Compare matrices at checkpoint 5000:
   - Should have < 1e-6 differences if implementations match

3. Compare errLike values:
   - Should match exactly for same reads
