# Error Model Sample Size Recommendations

## Background

During parity testing, we observed that with only ~30,000 reads:
- **21% matrix variation** due to read ordering and stochastic sampling
- Training was highly sensitive to which reads were selected
- Required deterministic mode + fixed ordering to achieve parity

This highlights that the error model training is **unstable at small sample sizes**.

## Why Small Samples Are Problematic

The CIGAR-based error model has:
- **82 × 82 = 6,724 transition cells** per matrix
- **12 matrices total** (6 quality bins × 2 read orientations)
- **~80,000 cells** to populate with training data

With only 5,000 pre-burnin reads:
- Most cells receive **< 1 training observation** on average
- Stochastic sampling further reduces effective training data
- Small perturbations (ordering, sampling) cause large matrix changes

## Recommendations

### Minimum Sample Sizes

| Sample Size | Status | Expected Stability |
|-------------|--------|-------------------|
| < 50,000 reads | ⚠️ Warning | Highly unstable, results may vary significantly |
| 50,000 - 200,000 reads | 🔶 Caution | Moderately stable, some variation expected |
| 200,000 - 1,000,000 reads | ✅ Good | Stable for most applications |
| > 1,000,000 reads | ✅ Excellent | Highly stable, recommended for production |

### Implementation: Add Warning Guards

```cpp
// In error model initialization
const size_t MIN_RECOMMENDED_READS = 200000;
const size_t MIN_WARNING_READS = 50000;

void checkSampleSize(size_t totalReads) {
    if (totalReads < MIN_WARNING_READS) {
        std::cerr << "WARNING: Only " << totalReads << " reads available for error model training.\n"
                  << "Error model may be unstable. Recommend >= " << MIN_RECOMMENDED_READS << " reads.\n"
                  << "Consider using --no-error-model for small datasets.\n";
    } else if (totalReads < MIN_RECOMMENDED_READS) {
        std::cerr << "NOTE: " << totalReads << " reads is below recommended " << MIN_RECOMMENDED_READS 
                  << " for optimal error model stability.\n";
    }
}
```

### Documentation Note

Add to user documentation:

> **Error Model Sample Size**: The CIGAR-based error model requires sufficient training data 
> for stable probability estimates. We recommend:
> - **Minimum**: 50,000 aligned reads
> - **Recommended**: 200,000+ aligned reads
> - **Optimal**: 1,000,000+ aligned reads
>
> For datasets with fewer reads, consider using `--no-error-model` to skip error model 
> training and use uniform alignment probabilities instead.

## Salmon's Defaults (for reference)

```cpp
// From Salmon source
static constexpr size_t PRE_BURNIN_FRAGS = 5000;      // Training only
static constexpr size_t BURNIN_FRAGS = 5000000;       // Training + use
```

Salmon expects 5M reads for full model convergence, but uses the model after just 5K reads. This works because:
1. Large datasets (millions of reads) are typical for RNA-seq
2. The model continues updating during burn-in phase
3. Stochastic variation averages out with large samples

## Test Results Supporting These Recommendations

| Total Reads | Matrix Variation (stochastic) | Matrix Variation (deterministic) |
|-------------|------------------------------|----------------------------------|
| ~30,000 | 21% | 8% |

Even with deterministic mode, 8% matrix variation remained due to the 307-read pre-burnin threshold timing difference. This variation would be negligible with larger samples.

## Action Items

1. [ ] Add sample size warning to CLI when < 50K reads
2. [ ] Add `--no-error-model` fallback option for small datasets
3. [ ] Document minimum sample sizes in user guide
4. [ ] Consider adaptive pre-burnin threshold based on dataset size
