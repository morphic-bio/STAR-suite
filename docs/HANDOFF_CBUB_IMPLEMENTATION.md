# CB/UB Tag Injection Implementation Handoff

**Date**: 2026-01-26  
**Branch**: `bamcbub`  
**Status**: Core implementation complete, all tests passing (buffered path re-enabled, two-pass deprecated)

---

## Summary

This document provides context for continuing development on CB/UB (Cell Barcode / UMI Barcode) tag injection for STARsolo BAM output. The implementation enables CB and UB tags to be written to both sorted and unsorted BAM files for Flex and non-Flex modes.

---

## What Has Been Implemented

### 1. CB/UB Independence Policy
- **CB and UB tags are now independent** - each can be present/absent based on read validity
- Removed legacy `requireCbUbTogether` enforcement
- Removed ZI tag injection (was never needed)
- See: `docs/CB_UB_INDEPENDENCE_IMPLEMENTATION_SUMMARY.md`

### 2. Sorted BAM CB/UB Injection (Both Flex and Non-Flex)
- Works for gene-mapped reads (reads that go through Solo counting)
- Uses `packedReadInfo` array populated during Solo counting
- Flex mode: Uses parallel `trackReadIdsForTags` hash to maintain readId→CB/UB mapping while keeping `inlineHashMode` active
- Non-Flex mode: Uses standard `packedReadInfo` population via stream files
- **Tag injection happens during `bamSortByCoordinate()` merge phase**

### 3. Unsorted BAM CB/UB Injection
- Uses `SamtoolsSorter` in "noSort" (deterministic order) mode as `g_unsortedTagBuffer`
- Records buffered during mapping, tags injected post-Solo counting via `bamUnsortedWithTags()`
- Transcriptome BAM protected via `skipGlobalBuffer_` flag
- **Two-pass path deprecated**: `outBAMunsortedUseSoloTmp` forced to `false`, fields kept for ABI compatibility
- **Auto-skipProcessing removed**: No longer auto-enables `--soloSkipProcessing` for unsorted CB/UB
- See: `docs/UNSORTED_BAM_CBUB_IMPLEMENTATION_SUMMARY.md`

### 4. Bug Fix: Non-Flex STARsolo Crash
- **Root cause**: `streamReads` was uninitialized (garbage pointer) in `SoloReadFeature` when `iChunk < 0`
- **Fixes applied**:
  1. Initialize `streamReads` to `nullptr` in constructor (`SoloReadFeature.cpp`)
  2. Add null guard in `SoloReadInfoLoader::load()` with diagnostic message
  3. Add defensive checks before `loader.load()` in `SoloFeature_countCBgeneUMI.cpp`
- See: `docs/BUG_REPORT_NONFLEX_SOLO_CRASH.md`

---

## Key Files Modified

| File | Purpose |
|------|---------|
| `core/legacy/source/SoloReadFeature.cpp` | Initialize `streamReads` to `nullptr` |
| `flex/source/SoloReadInfoLoader.cpp` | Add null guard for `streamReads` |
| `core/legacy/source/SoloFeature_countCBgeneUMI.cpp` | Defensive checks, `trackReadIdsForTags` logic |
| `core/legacy/source/ParametersSolo.cpp` | `trackReadIdsForTags` activation |
| `core/legacy/source/BAMoutput.cpp` | Unsorted routing to `g_unsortedTagBuffer` |
| `core/legacy/source/BAMoutput.h` | `skipGlobalBuffer_` flag |
| `core/legacy/source/bamSortByCoordinate.cpp` | Tag injection during sort merge, `bamUnsortedWithTags()` |
| `core/legacy/source/GlobalVariables.cpp/h` | `g_unsortedTagBuffer` declaration |
| `core/legacy/source/STAR.cpp` | Initialize `g_unsortedTagBuffer`, call `bamUnsortedWithTags()` |
| `core/features/bamsort/source/SamtoolsSorter.cpp/h` | `noSort` mode for deterministic order buffering |
| `flex/source/SoloReadFeature.cpp` | `readIdTracker_` hash for Flex tag injection |

---

## Test Results

### Flex Mode (100K dataset, 2 lanes quick test)
After the buffered path re-enablement (2026-01-26), sorted and unsorted BAM now have identical CB/UB tag rates:
```
Test 1: Sorted BAM - PASS (CB: 79.38%, UB: 40.43%)
Test 2: Unsorted BAM - PASS (CB: 79.38%, UB: 40.43%)  
Test 3: Sorted vs Unsorted MEX - PASS (identical)
```

**Key improvement**: Unsorted BAM now matches sorted BAM tag rates (previously unsorted had higher rates due to using inline extraction before Solo correction).

### Non-Flex Mode (A375 100K GEX)
```
Test 1: Sorted BAM - PASS (CB: 120289, UB: 120241)
Test 2: Unsorted BAM - PASS (CB: 120289, UB: 120241)
Test 3: Tag comparison - PASS (0 mismatches)
Test 4: TranscriptomeSAM separation - PASS
```

### Regression Test Scripts
- `tests/run_cbub_regression_test.sh` - Flex sorted/unsorted test
- `tests/run_unsorted_cbub_a375_nonflex.sh` - Non-Flex A375 test

---

## Known Limitations

1. **Tag injection only for gene-mapped reads**: Reads that don't go through Solo counting (unmapped, filtered) won't have CB/UB tags. This is by design (same as non-Flex behavior).

2. **Memory overhead for Flex**: `trackReadIdsForTags` adds a parallel hash (~16 bytes per unique read) when CB/UB injection is enabled.

3. **Unsorted output order**: Deterministic order (spill files in creation order, then in-memory records). Not true FIFO, but consistent across runs.

---

## Architecture Overview

### Sorted BAM Path
```
Mapping Phase:
  ReadAlign → Solo streams/hash → packedReadInfo populated
  
Sort Phase (bamSortByCoordinate):
  For each sorted record:
    1. Extract readId from BAM
    2. Lookup CB/UB in packedReadInfo (or readIdTracker_ for Flex)
    3. Inject tags via CbUbTagInjector
    4. Write to sorted BAM
```

### Unsorted BAM Path
```
Mapping Phase:
  BAMoutput::unsortedOneAlign → g_unsortedTagBuffer (if not skipGlobalBuffer_)
  
Post-Solo Phase (bamUnsortedWithTags):
  1. Finalize g_unsortedTagBuffer
  2. Stream records in deterministic order (spills first, then in-memory)
  3. Lookup CB/UB, inject tags
  4. Write to Aligned.out.bam
```

---

## Future Work (from plans/unsorted_bam_cbub_insertion_plan.md)

The unsorted BAM implementation is complete. Potential enhancements:

1. **Memory optimization**: The `trackReadIdsForTags` hash could potentially be merged into the main inline hash to reduce memory overhead for very large datasets.

2. **Extended tag coverage**: Currently only gene-mapped reads get tags. Could extend to include reads that match whitelist but don't map to genes.

---

## How to Test

### Quick Flex test (2 lanes, ~90s)
```bash
cd /mnt/pikachu/STAR-suite
bash tests/run_cbub_regression_test.sh --quick
```

### Full Flex test (8 lanes, ~60s)
```bash
bash tests/run_cbub_regression_test.sh
```

### Non-Flex A375 test (~5min)
```bash
bash tests/run_unsorted_cbub_a375_nonflex.sh
```

---

## Build Instructions

```bash
cd /mnt/pikachu/STAR-suite/core/legacy/source
make STAR -j8
```

Binary location: `/mnt/pikachu/STAR-suite/core/legacy/source/STAR`

---

## Key Parameters

| Parameter | Effect |
|-----------|--------|
| `--outSAMattributes ... CB UB` | Enable CB/UB tag output |
| `--outSAMtype BAM SortedByCoordinate` | Sorted BAM with tags |
| `--outSAMtype BAM Unsorted` | Unsorted BAM with tags (automatic) |
| `--flex yes` | Enable Flex mode (uses inline hash) |

Note: Unsorted CB/UB tagging is now automatic when CB/UB are in `--outSAMattributes` (the legacy `--soloAddTagsToUnsorted` flag has been removed).

---

## Contact / References

- Implementation summaries: `docs/CB_UB_INDEPENDENCE_IMPLEMENTATION_SUMMARY.md`, `docs/UNSORTED_BAM_CBUB_IMPLEMENTATION_SUMMARY.md`
- Bug report: `docs/BUG_REPORT_NONFLEX_SOLO_CRASH.md`
- Original plan: `plans/unsorted_bam_cbub_insertion_plan.md`
