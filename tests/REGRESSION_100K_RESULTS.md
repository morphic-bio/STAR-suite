# 100K Regression Test Results

**Date:** December 23, 2025  
**Baseline Commit:** `2f89f17` (cleaned large files)  
**Current Commit:** `473d767` (HEAD)

## Test Configuration

- **Dataset:** 100K downsampled SC2300771 (8 lanes)
- **Reference:** `/storage/flex_filtered_reference/star_index`
- **Parameters:** Identical to `run_flex_multisample_test.sh`
- **Flex Mode:** Enabled with multi-sample output

## Test Results

### Status: ✅ PASSED

### Cell Counts: IDENTICAL

| Sample | Baseline | Current | Status |
|--------|----------|--------|--------|
| BC004  | 2,291    | 2,291  | ✓      |
| BC006  | 1,867    | 1,867  | ✓      |
| BC007  | 4,142    | 4,142  | ✓      |
| BC008  | 4,225    | 4,225  | ✓      |
| **TOTAL** | **12,525** | **12,525** | **✓** |
| **Total UMIs** | **544,439** | **544,439** | **✓** |

### MEX File Comparison

#### Pooled Raw Outputs (`Solo.out/Gene/raw/`)

| File | Status | Notes |
|------|--------|-------|
| `matrix.mtx` | ⚠ Order difference | Same entries, different order (619,225 entries) |
| `barcodes.tsv` | ✓ IDENTICAL | Exact match |
| `features.tsv` | ✓ IDENTICAL | Exact match |

#### Per-Sample Filtered Outputs (`per_sample/{SAMPLE}/Gene/filtered/`)

| Sample | matrix.mtx | barcodes.tsv | features.tsv |
|--------|------------|--------------|-------------|
| BC004  | ⚠ Order diff | ✓ IDENTICAL | ✓ IDENTICAL |
| BC006  | ⚠ Order diff | ✓ IDENTICAL | ✓ IDENTICAL |
| BC007  | ⚠ Order diff | ✓ IDENTICAL | ✓ IDENTICAL |
| BC008  | ⚠ Order diff | ✓ IDENTICAL | ✓ IDENTICAL |

### Analysis

**Matrix Entry Ordering:**
- All `matrix.mtx` files differ only in the order of entries
- When sorted, baseline and current matrices are **identical**
- Matrix Market format allows entries in any order
- This is **expected behavior** and does not affect downstream analysis

**Verification:**
```bash
# Sorted comparison shows zero differences
tail -n +3 matrix.mtx | sort | diff baseline current
# Result: No differences found
```

## Conclusion

**✅ MEX PARITY CONFIRMED**

All differences are due to entry ordering in sparse matrices, which is acceptable per Matrix Market format specification. All actual data (cell counts, barcodes, features, UMI counts) are **identical** between baseline and current.

The regression test confirms that changes between `2f89f17` and `473d767` do not affect the correctness of MEX outputs.

## Output Locations

- **Baseline outputs:** `tests/regress_100K_baseline/`
- **Current outputs:** `tests/regress_100K_current/`

## Test Scripts

- **Main test:** `tests/run_100K_regression_test.sh`
- **Comparison:** `tests/compare_mex_outputs.sh`
- **Documentation:** `tests/REGression_100K_README.md`

## Running the Test

```bash
cd /mnt/pikachu/STAR-Flex/tests
SKIP_CONFIRM=yes ./run_100K_regression_test.sh
```

Expected runtime: ~10-20 minutes (builds STAR twice, runs test twice)

