# STAR VB/EM Harmonization with Salmon's VBEM Optimizer

**Date:** January 1, 2026  
**Goal:** Harmonize STAR's VB/EM quantification behavior with Salmon's VBEM optimizer

---

## Changes Implemented

### 1. VB Defaults in STAR.cpp

**File:** `source/STAR.cpp`

**Change:** Removed VB-specific parameter overrides. VB now uses the same defaults as EM (defined in `EMParams`):
- `max_iters = 10000` (was 200)
- `min_iters = 100` (now enforced)
- `tolerance = 0.01` (was 1e-8)
- `threads = 0` (use OMP default, multi-thread capable)

**Before:**
```cpp
if (params.use_vb) {
    params.max_iters = 200;
    params.tolerance = 1e-8;
}
params.threads = 1;  // Single-threaded
```

**After:**
```cpp
// Use defaults from EMParams (max_iters=10000, min_iters=100, tolerance=0.01)
// Do NOT override for VB - let VB use same defaults as EM for Salmon parity
params.threads = 0;  // 0 = use OMP default (multi-thread capable)
```

---

### 2. Salmon's VB Initialization Mix (fracObserved + uniformPrior)

**File:** `source/libem/vb_engine.cpp`

**Change:** Replaced uniform initialization with Salmon's VBEMOptimizer initialization:

**Salmon's approach:**
- Compute `projectedCounts[i]` from EC counts (includes multi-mappers, weighted by EC weights)
- Compute `fracObserved = min(0.999, totalWeight / numRequiredFragments)`
- Compute `uniformPrior = (1.0 - fracObserved) * (totalWeight / numActive)`
- Initialize: `alpha[i] = priorAlphas[i] + projectedCounts[i] * fracObserved + uniformPrior`
- If `!per_transcript_prior`, scale both by `eff_len` (per-nucleotide prior mode)

**Key improvements:**
- Uses `projectedCounts` (includes multi-mappers) instead of just `unique_counts`
- Mixes observed counts with uniform prior based on `fracObserved`
- Respects `num_required_fragments` parameter (default: 5e6, matching Salmon)

---

### 3. Enforce min_iters Before Convergence Check

**File:** `source/libem/vb_engine.cpp`

**Change:** Convergence is only checked after `min_iters` iterations (default: 100).

**Before:**
```cpp
// Check convergence every iteration
bool converged = true;
for (size_t i = 0; i < state.n; ++i) {
    if (alpha[i] > alphaCheckCutoff) {
        // ... check convergence
    }
}
```

**After:**
```cpp
// Only check convergence after min_iters (Salmon's approach)
bool converged = false;
if (iter + 1 >= params.min_iters) {
    // Check convergence using alpha relative difference
    converged = true;
    // ... convergence check
}
```

---

### 4. Configurable Alpha Cutoff and Zero Threshold

**File:** `source/libem/vb_engine.cpp`

**Change:** Replaced hardcoded constants with configurable values from `EMParams`:
- `alphaCheckCutoff` → `params.alpha_check_cutoff` (default: 1e-2)
- `minAlpha` → `params.zero_threshold` (default: 1e-8)

**Before:**
```cpp
constexpr double minAlpha = 1e-8;
constexpr double alphaCheckCutoff = 1e-2;
// ...
if (alpha[i] > alphaCheckCutoff) { ... }
if (result.counts[i] <= minAlpha) { ... }
```

**After:**
```cpp
double alphaCheckCutoff = params.alpha_check_cutoff;  // Configurable
double zero_threshold = params.zero_threshold;  // Configurable
// ...
if (alpha[i] > alphaCheckCutoff) { ... }
if (result.counts[i] <= zero_threshold) { ... }
```

---

### 5. ELBO Frequency Optimization

**File:** `source/libem/vb_engine.cpp`

**Change:** ELBO is only computed when needed (at convergence, max iterations, or debug mode), not every iteration.

**Before:**
```cpp
// Compute ELBO every iteration (expensive)
double curr_elbo = compute_elbo(...);
result.final_ll = curr_elbo;
```

**After:**
```cpp
// Only compute ELBO if converged, at max iterations, or debug mode
if (converged || iter + 1 >= params.max_iters || params.debug_trace) {
    double curr_elbo = compute_elbo(...);
    result.final_ll = curr_elbo;
}
// ... at end of function ...
// Compute ELBO at end if not already computed
if (result.final_ll == 0.0) {
    result.final_ll = compute_elbo(...);
}
```

---

### 6. Deterministic Path (Multi-Thread Capable)

**No changes needed** - The code already supports deterministic results when `params.threads == 1`:
- Thread-local storage with deterministic reduction
- Fixed iteration order
- No race conditions

For parity testing, run with `--runThreadN 1` externally.

---

## Parameter Defaults (EMParams)

All defaults match Salmon's VBEM optimizer:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_iters` | 10000 | Maximum iterations |
| `min_iters` | 100 | Minimum iterations before convergence check |
| `tolerance` | 0.01 | Relative change tolerance for convergence |
| `alpha_check_cutoff` | 1e-2 | Only check convergence for transcripts with alpha > cutoff |
| `zero_threshold` | 1e-8 | Truncation threshold (minAlpha) |
| `vb_prior` | 0.01 | Dirichlet concentration parameter |
| `num_required_fragments` | 5e6 | Used for fracObserved calculation |

---

## Expected Improvements

1. **Better Convergence:** VB will now run for at least `min_iters` iterations before checking convergence
2. **Salmon Parity:** Initialization matches Salmon's VBEMOptimizer exactly
3. **Performance:** ELBO computation optimized (only when needed)
4. **Flexibility:** All thresholds configurable via `EMParams`
5. **Multi-Thread:** Remains multi-thread capable (deterministic when `threads=1`)

---

## Validation

To validate these changes:

1. **Build STAR:**
   ```bash
   cd /mnt/pikachu/STAR-Flex/source
   make clean && make
   ```

2. **Run VB parity test:**
   ```bash
   STAR --runThreadN 1 \
        --quantMode TranscriptVB TranscriptomeSAM GeneCounts \
        --quantVBLibType ISF \
        --quantVBem 0 \
        ...
   ```

3. **Compare with Salmon VB:**
   ```bash
   salmon quant -l ISF --useVBOpt ...
   ```

4. **Expected results:**
   - STAR VB should converge (not stop at 200 iterations)
   - Initialization should match Salmon's projectedCounts + uniformPrior mix
   - Correlation should improve (especially for VB mode)

---

## Files Modified

1. `source/STAR.cpp` - Removed VB-specific parameter overrides
2. `source/libem/vb_engine.cpp` - Implemented Salmon's initialization, min_iters enforcement, configurable thresholds, ELBO optimization
3. `source/libem/em_types.h` - Already had correct defaults (no changes needed)

---

## Notes

- **Backward Compatibility:** The changes maintain backward compatibility for EM mode
- **Thread Safety:** Code remains multi-thread capable and deterministic when `threads=1`
- **Salmon Parity:** Initialization and convergence behavior now match Salmon's VBEM optimizer

