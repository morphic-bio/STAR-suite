# 100K Regression Test

This regression test compares MEX outputs between commit `2f89f17` (baseline) and current HEAD to ensure parity.

## Quick Start

```bash
cd /mnt/pikachu/STAR-Flex/tests
./run_100K_regression_test.sh
```

## What It Does

1. **Builds baseline**: Checks out commit `2f89f17`, builds STAR, runs the 100K test
2. **Builds current**: Checks out HEAD, builds STAR, runs the same test
3. **Compares outputs**: Compares all MEX files (matrix.mtx, barcodes.tsv, features.tsv) for:
   - Pooled raw output: `Solo.out/Gene/raw/`
   - Per-sample filtered output: `per_sample/{BC004,BC006,BC007,BC008}/Gene/filtered/`

## Output Locations

- **Baseline outputs**: `tests/regress_100K_baseline/`
- **Current outputs**: `tests/regress_100K_current/`

## Manual Comparison

To compare outputs manually:

```bash
./compare_mex_outputs.sh tests/regress_100K_baseline tests/regress_100K_current
```

## Test Parameters

Uses the same dataset and parameters as `run_flex_multisample_test.sh`:
- Dataset: `/storage/downsampled_100K/SC2300771/...` (8 lanes)
- Reference: `/storage/flex_filtered_reference/star_index`
- Flex mode enabled with multi-sample output

## Expected Runtime

- Baseline build + run: ~5-10 minutes
- Current build + run: ~5-10 minutes
- Comparison: < 1 minute
- **Total: ~10-20 minutes**

## Notes

- The script will checkout different commits, so ensure you have no uncommitted changes
- Original commit is restored at the end
- Output directories are preserved for manual inspection

