# Modularize assignBarcodes: Stage 1-2 Implementation Summary

**Date:** January 2026  
**Status:** Complete  
**Based on:** `plans/modularize_assignBarcodes_plan.md`

## Overview

This document summarizes the implementation of Stages 1-2 of the assignBarcodes modularization plan. The goal was to split the monolithic `assignBarcodes.c` into explicit stages while maintaining byte-identical output.

## What Was Implemented

### Stage 1: Extract In-Memory Counts Result

Created a new `pf_counts_result` structure to encapsulate deduped counts, replacing ad-hoc local variables in `finalize_processing()`.

**New Files:**
- `core/features/process_features/include/pf_counts.h` - Header with struct definition and API
- `core/features/process_features/src/pf_counts.c` - Implementation

**Key APIs:**
```c
// Initialize a counts result for n features
pf_counts_result* pf_counts_result_init(int n_features);

// Build deduped counts from hashes (calls find_deduped_counts internally)
pf_counts_result* pf_build_deduped_counts(
    data_structures *hashes,
    int n_features,
    uint16_t stringency,
    uint16_t min_counts
);

// Free all resources
void pf_counts_result_free(pf_counts_result *result);

// Reset histograms for reuse between passes
void pf_counts_result_reset_histograms(pf_counts_result *result);
```

**Structure:**
```c
typedef struct pf_counts_result {
    khash_t(u32ptr) *barcode_to_deduped_hash;  // barcode -> feature counts
    vec_u32_t **feature_hist;                   // per-feature UMI histograms
    int n_features;
    int *total_deduped_counts;                  // per-feature totals
    int *total_barcoded_counts;                 // per-feature raw counts
} pf_counts_result;
```

### Stage 2: Extract MEX Writer

Created a dedicated module for all output file generation, replacing the monolithic `printFeatureCounts()` function.

**New Files:**
- `core/features/process_features/include/mex_writer.h` - Header with config and API
- `core/features/process_features/src/mex_writer.c` - Implementation

**Key APIs:**
```c
// Write all MEX outputs (barcodes.txt, features.txt, matrix.mtx, stats.txt,
// feature_per_cell.csv, heatmaps, histograms)
int mex_write_all(const mex_writer_config *config, pf_counts_result *counts);

// Write only core files (barcodes.txt, features.txt, matrix.mtx)
int mex_write_core(const mex_writer_config *config, pf_counts_result *counts);

// Individual writers for specific outputs
int mex_write_stats(...);
int mex_write_feature_per_cell(...);
int mex_write_heatmaps(...);
```

**Configuration:**
```c
typedef struct mex_writer_config {
    const char *output_dir;
    feature_arrays *features;
    data_structures *hashes;
    statistics *stats;
    khash_t(strptr) *filtered_barcodes_hash;  // NULL for unfiltered
    int min_heatmap_counts;
} mex_writer_config;
```

### Wiring into finalize_processing()

`finalize_processing()` in `assignBarcodes.c` was refactored to use the new modular code:

**Before:**
```c
void finalize_processing(...) {
    int total_deduped_counts[...];
    int total_barcoded_counts[...];
    vec_u32_t **feature_hist = calloc(...);
    
    khash_t(u32ptr)* barcode_to_deduped_hash = kh_init(u32ptr);
    find_deduped_counts(hashes, barcode_to_deduped_hash, ...);
    
    printFeatureCounts(..., barcode_to_deduped_hash, NULL);
    if (filtered_barcodes_hash) {
        // reset feature_hist
        printFeatureCounts(..., barcode_to_deduped_hash, filtered_barcodes_hash);
    }
    
    // manual cleanup of barcode_to_deduped_hash
    // use local arrays for downstream operations
}
```

**After:**
```c
void finalize_processing(...) {
    pf_counts_result *counts_result = pf_build_deduped_counts(hashes, ...);
    
    mex_writer_config config = {
        .output_dir = directory,
        .features = features,
        .hashes = hashes,
        .stats = stats,
        .filtered_barcodes_hash = NULL,
        .min_heatmap_counts = min_heatmap
    };
    
    mex_write_all(&config, counts_result);  // unfiltered
    if (filtered_barcodes_hash) {
        config.filtered_barcodes_hash = filtered_barcodes_hash;
        mex_write_all(&config, counts_result);  // filtered
    }
    
    // use counts_result->feature_hist and counts_result->total_* for downstream
    pf_counts_result_free(counts_result);
}
```

## Output Files Produced

The mex_writer produces all the same outputs as the original `printFeatureCounts()`:

| File | Description |
|------|-------------|
| `barcodes.txt` | One barcode per line |
| `features.txt` | One feature name per line |
| `matrix.mtx` | Matrix Market coordinate format (feature × barcode) |
| `stats.txt` | Summary statistics |
| `feature_per_cell.csv` | Per-barcode: num_features, top_feature, total_umi |
| `feature_richness_histogram.html` | Features per barcode distribution |
| `feature_multiplicity_histogram.html` | UMI count distribution |
| `Feature_counts_heatmap.png` | Co-expression heatmap |
| `Feature_types_heatmap.png` | Feature richness heatmap |
| `deduped_counts_histograms.txt` | Per-feature UMI histograms (written by finalize_processing) |
| `umi_counts_histogram` | Cumulative UMI histogram (written by finalize_processing) |

For filtered output, files go to `{output_dir}/filtered/`.

## Key Implementation Details

### Multi-Pass Safety

The mex_writer properly handles being called twice (raw + filtered):
- Clears `counts->total_deduped_counts` and `counts->total_barcoded_counts` at start
- Destroys and recreates `counts->feature_hist[]` entries
- Clears cumulative histogram (`feature_hist[0]`) before rebuilding

### Preserved Behaviors

- Matrix Market header format: `%%MatrixMarket matrix coordinate real general`
- Metadata line: `%metadata_json: {"software_version": "assignBarcodes-0.1", "format_version": 1}`
- Barcode iteration order (khash order, not sorted)
- `translate_NXT` applied to barcode strings when enabled
- All stderr log lines preserved ("Skipped barcodes", "Processed barcodes", etc.)

## Testing

### Unit Tests

**`tests/test_pf_counts.c`** (4 tests):
- `pf_counts_result_init` - Verifies initialization
- `pf_counts_result_free_null` - Safe NULL handling
- `pf_counts_result_reset_histograms` - Histogram clearing
- `pf_counts_result_with_hash` - Manual hash population and retrieval

**`tests/test_mex_writer.c`** (5 tests):
- `mex_write_core_creates_files` - Basic file creation
- `matrix_header_format` - Matrix Market header validation
- `features_line_count` - Correct feature count
- `filtered_subdir` - Filtered output goes to `filtered/`
- `counts_cleared_between_passes` - Arrays properly cleared

### Automated Parity Test

**`tests/test_mex_writer_parity.c`** (2 tests):
Compares `mex_write_all()` output against `printFeatureCounts()` output byte-by-byte.

- `parity_core_files` - Compares unfiltered outputs:
  - `barcodes.txt` ✓ byte-identical
  - `features.txt` ✓ byte-identical
  - `matrix.mtx` ✓ byte-identical
  - `stats.txt` ✓ byte-identical
  - `feature_per_cell.csv` ✓ byte-identical

- `parity_filtered_output` - Compares filtered outputs:
  - `filtered/barcodes.txt` ✓ byte-identical
  - `filtered/features.txt` ✓ byte-identical
  - `filtered/matrix.mtx` ✓ byte-identical

Run parity tests:
```bash
cd core/features/process_features
make tests/test_mex_writer_parity
./tests/test_mex_writer_parity
```

### Regression Test

**`tests/test_assignbarcodes_regression.sh`** - Compares `assignBarcodes` output against a saved baseline fixture:

- Runs `assignBarcodes` on synthetic FASTQ fixture
- Compares output files byte-by-byte against saved baseline
- Files compared: `barcodes.txt`, `features.txt`, `matrix.mtx`, `stats.txt`, `feature_per_cell.csv`, `deduped_counts_histograms.txt`

Run regression test:
```bash
cd core/features/process_features
./tests/test_assignbarcodes_regression.sh

# To update baseline (after intentional changes):
./tests/test_assignbarcodes_regression.sh --update-baseline
```

### Smoke Test

Ran with real A375 CRISPR data:
- **Input:** 1000 reads, 11 CRISPR guide features, 10x whitelist
- **Results:** 116 barcodes, 942 raw counts, 131 deduped counts, 96% assignment rate
- **Output:** All files created with correct format

```bash
cd core/features/process_features
./assignBarcodes \
  -w /storage/A375/3M-5pgex-jan-2023.txt \
  -f features.csv \
  -d output \
  /path/to/fastqs \
  -b 16 -u 12
```

## Build Changes

**Makefile updates:**
```makefile
# Added to CORE_SRCS
CORE_SRCS = ... pf_counts.c mex_writer.c

# New test targets
TEST_PF_COUNTS_TARGET=tests/test_pf_counts
TEST_MEX_WRITER_TARGET=tests/test_mex_writer

test: $(TEST_API_TARGET) $(TEST_PF_COUNTS_TARGET) $(TEST_MEX_WRITER_TARGET)
```

## Files Changed

| File | Change |
|------|--------|
| `include/pf_counts.h` | **New** - Counts result struct and API |
| `src/pf_counts.c` | **New** - Counts result implementation |
| `include/mex_writer.h` | **New** - MEX writer API |
| `src/mex_writer.c` | **New** - MEX writer implementation |
| `src/assignBarcodes.c` | Modified - Wired new modules into finalize_processing() |
| `Makefile` | Modified - Added new sources and test targets |
| `tests/test_pf_counts.c` | **New** - Unit tests for pf_counts |
| `tests/test_mex_writer.c` | **New** - Unit tests for mex_writer |
| `tests/test_mex_writer_parity.c` | **New** - Byte-level parity test vs printFeatureCounts() |
| `tests/test_assignbarcodes_regression.sh` | **New** - Regression test against saved baseline |
| `tests/fixtures/assignbarcodes_baseline/` | **New** - Fixture data and expected outputs |

## What's NOT Changed

- `printFeatureCounts()` still exists (not deleted, just no longer called from pipeline)
- CLI interface unchanged
- Output file formats identical
- All existing functionality preserved

## Next Steps: Stage 3

Stage 3 will add the filter stage with EmptyDrops integration:

1. Create `barcode_filter.c` / `barcode_filter.h`
2. Create `emptydrops_bridge.c` to convert `pf_counts_result` to `scrna_matrix_input`
3. Wire `pf_run_emptydrops_premex()` into the pipeline
4. Add `--filtered_barcodes` list path as alternative to EmptyDrops

The modular structure now makes this straightforward:
```
Assignment → pf_build_deduped_counts()
         ↓
Filter    → pf_filter_barcodes() or EmptyDrops  [Stage 3]
         ↓
Output    → mex_write_all()
```

## Reference Documents

- `plans/modularize_assignBarcodes_plan.md` - Original plan
- `plans/emptydrops_wiring_handoff.md` - EmptyDrops integration context
- `plans/emptydrops_refactor_plan.md` - EmptyDrops defaults and behavior
