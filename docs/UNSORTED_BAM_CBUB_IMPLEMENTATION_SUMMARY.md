# Unsorted BAM CB/UB Tag Injection Implementation Summary

**Date**: 2026-01-26  
**Branch**: `bamcbub`  
**Status**: Complete  
**Regression Test**: `tests/run_cbub_regression_test.sh`

## Overview

This implementation adds CB/UB tag injection support for unsorted BAM output by reusing the existing sorted BAM tag injection path. The key insight is that both output modes can share the same `SoloFeature::addBAMtags()` hook, with the only difference being whether records are sorted by coordinate before output.

**Key Changes (2026-01-26)**:
- **Automatic CB/UB injection**: CB/UB tags are automatically injected when requested via `--outSAMattributes`
- **Deterministic order**: Records are emitted in a deterministic order (spill files in creation order, then in-memory records), not FIFO
- **Two-pass path removed**: The legacy `--soloAddTagsToUnsorted` flag and two-pass path have been completely removed
- **Auto-skipProcessing removed**: No longer auto-enables `--soloSkipProcessing` for unsorted BAM

## Problem Statement

Previously, CB/UB tag injection for unsorted BAM used a separate inline path (`CbUbTagInjector` in `BAMoutput::unsortedOneAlign()`), which:
1. Injected tags during mapping, before Solo counting completed
2. Could not use the corrected CB/UMI values from Solo processing
3. Was a separate code path from sorted BAM, risking drift

## Solution

Route unsorted BAM with CB/UB tags through the same pipeline as sorted BAM, using a "no-sort" mode that maintains **deterministic order** instead of coordinate sorting.

### Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        MAPPING PHASE                            │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ReadAlign::outputAlignments()                                  │
│        │                                                        │
│        ▼                                                        │
│  ┌─────────────────────────────────────────────────────────┐   │
│  │ BAMoutput::unsortedOneAlign()                           │   │
│  │                                                         │   │
│  │ if (g_unsortedTagBuffer != nullptr)                     │   │
│  │     g_unsortedTagBuffer->addRecord(bam, readId)         │   │
│  │     return; // Buffer, don't write directly             │   │
│  │ else                                                    │   │
│  │     // Legacy direct write (no CB/UB tags)              │   │
│  └─────────────────────────────────────────────────────────┘   │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                      SOLO COUNTING PHASE                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  - Inline hash deduplication                                    │
│  - CB rescue / correction                                       │
│  - readIdTracker_ → packedReadInfo population                   │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                       OUTPUT PHASE                              │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  bamUnsortedWithTags()                                          │
│        │                                                        │
│        ▼                                                        │
│  while (g_unsortedTagBuffer->nextRecord(&bam, &readId))         │
│      soloFeat->addBAMtags(bam, readId)  // Inject CB/UB         │
│      bgzf_write(bam)                     // Write to output     │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

## Implementation Details

### 1. SamtoolsSorter `noSort` Mode

Added a `noSort` parameter to `SamtoolsSorter` that:
- Buffers records to memory/disk in insertion order
- Skips coordinate sorting in `sortAndSpill()` and `finalize()`
- Uses deterministic iteration in `nextRecord()`: spill files first (in creation order), then in-memory records

```cpp
// core/features/bamsort/source/SamtoolsSorter.h
SamtoolsSorter(uint64_t maxRAM, int nThreads, const string& tmpDir, 
               Parameters& P, bool noSort = false);
```

### 2. Global Unsorted Buffer

Added `g_unsortedTagBuffer` to store records during mapping:

```cpp
// core/legacy/source/GlobalVariables.h
extern SamtoolsSorter* g_unsortedTagBuffer;
```

Created when unsorted BAM + CB/UB tags requested:
```cpp
// core/legacy/source/STAR.cpp
if (P.outBAMunsorted && P.pSolo.samAttrYes) {
    g_unsortedTagBuffer = new SamtoolsSorter(P.limitBAMsortRAM,
                                              P.outBAMsortingThreadNactual,
                                              P.outBAMsortTmpDir,
                                              P,
                                              true);  // noSort = true
}
```

### 3. Record Routing

Modified `BAMoutput::unsortedOneAlign()` to route through buffer:

```cpp
// core/legacy/source/BAMoutput.cpp
if (g_unsortedTagBuffer != nullptr) {
    uint32_t readId = static_cast<uint32_t>(iReadAll);
    g_unsortedTagBuffer->addRecord(bamIn, bamSize, readId, hasY);
    return;  // Don't write directly
}
```

### 4. Output with Tag Injection

New `bamUnsortedWithTags()` function streams records with tags:

```cpp
// core/legacy/source/bamSortByCoordinate.cpp
void bamUnsortedWithTags(Parameters &P, Genome &genome, Solo &solo) {
    g_unsortedTagBuffer->finalize();
    
    while (g_unsortedTagBuffer->nextRecord(&bamData, &bamSize, &readId, &hasY)) {
        // Inject CB/UB tags using same path as sorted BAM
        solo.soloFeat[...]->addBAMtags(bam0, size0, bam1, readId);
        bgzf_write(bgzfOut, bam0, size0);
    }
}
```

### 5. ReadId Tracking for Flex

Extended `trackReadIdsForTags` to also enable for unsorted BAM:

```cpp
// core/legacy/source/ParametersSolo.cpp
if ((pP->outBAMcoord || pP->outBAMunsorted) && inlineHashMode) {
    trackReadIdsForTags = true;
}
```

## Files Modified

| File | Changes |
|------|---------|
| `core/features/bamsort/source/SamtoolsSorter.h` | Added `noSort_` flag, deterministic iteration state |
| `core/features/bamsort/source/SamtoolsSorter.cpp` | Implemented noSort mode, deterministic iteration |
| `core/legacy/source/GlobalVariables.h` | Added `g_unsortedTagBuffer` declaration |
| `core/legacy/source/GlobalVariables.cpp` | Added `g_unsortedTagBuffer` definition |
| `core/legacy/source/STAR.cpp` | Initialize buffer, call `bamUnsortedWithTags()` |
| `core/legacy/source/BAMoutput.h` | Added `skipGlobalBuffer_` flag and setter |
| `core/legacy/source/BAMoutput.cpp` | Route through buffer when active (respecting skipGlobalBuffer_) |
| `core/legacy/source/ReadAlignChunk.cpp` | Set skipGlobalBuffer=true for transcriptome BAM |
| `core/legacy/source/bamSortByCoordinate.h` | Declare `bamUnsortedWithTags()` |
| `core/legacy/source/bamSortByCoordinate.cpp` | Implement `bamUnsortedWithTags()` |
| `core/legacy/source/ParametersSolo.cpp` | Enable trackReadIdsForTags for unsorted |

## Test Results

Both sorted and unsorted BAM output now produce identical tag injection results:

| Metric | Sorted BAM | Unsorted BAM |
|--------|------------|--------------|
| Total reads | 378,186 | 378,186 |
| Reads with CB tag | 300,217 (79.38%) | 300,217 (79.38%) |
| Reads with UB tag | 152,905 (40.43%) | 152,905 (40.43%) |
| Reads with both | 152,905 | 152,905 |

Test scripts:
- `tests/run_sorted_bam_cbub_test.sh`
- `tests/run_unsorted_bam_cbub_test.sh`

## Usage

CB/UB tag injection for unsorted BAM is **automatic** when:
1. `--outSAMtype BAM Unsorted` is specified
2. `--outSAMattributes` includes `CB` and/or `UB`
3. Solo counting is enabled (e.g., `--soloType CB_UMI_Simple`)

Unsorted CB/UB tags are automatically enabled via the buffered mode when the above conditions are met.

Example:
```bash
STAR \
  --outSAMtype BAM Unsorted \
  --outSAMattributes NH HI AS nM NM GX GN CB UB \
  --soloType CB_UMI_Simple \
  ...
```

The log will show:
```
NOTE: Using buffered mode for unsorted BAM CB/UB tag injection.
...
Jan 25 19:52:08 ..... writing unsorted BAM with CB/UB tags
```

## Benefits

1. **Code Unification**: Sorted and unsorted paths share `addBAMtags()` hook
2. **Consistency**: Both paths produce identical tag values for same reads
3. **Maintainability**: Single tag injection logic to maintain
4. **Correctness**: Tags use corrected CB/UMI values from Solo counting

## Important Fix: Transcriptome BAM Separation

When `--quantMode TranscriptomeSAM` is used alongside unsorted BAM, the transcriptome 
output (`Aligned.toTranscriptome.out.bam`) must go to a separate file, not through 
`g_unsortedTagBuffer`.

**Fix implemented**: Added `skipGlobalBuffer_` flag to `BAMoutput` class:
- Set to `true` for `chunkOutBAMquant` (transcriptome output)
- `unsortedOneAlign()` checks this flag before routing to `g_unsortedTagBuffer`
- Ensures transcriptome BAM writes directly to its own file

Without this fix, transcriptome records would be incorrectly merged into the main 
unsorted BAM output.

## Regression Test Results (8-lane 100K dataset)

Run with: `tests/run_cbub_regression_test.sh`

| Test | Status | Details |
|------|--------|---------|
| Sorted BAM tags | PASS | 78.73% CB, 40.28% UB |
| Sorted MEX vs gold | PASS | Cells: 12525 (exact), UMIs: +4 |
| Unsorted BAM tags | PASS | 78.73% CB, 40.28% UB |
| Unsorted MEX vs gold | PASS | Cells: 12525 (exact), UMIs: +4 |
| Sorted vs Unsorted | PASS | Identical outputs |

The +4 UMI difference is a pre-existing known difference, not caused by this implementation.

## Non-Flex Mode Support

Non-Flex STARsolo (standard GEX without `--flex yes`) is fully supported for CB/UB tag injection. A bug discovered during testing (null `streamReads` pointer crash) was fixed by:

1. Initializing `streamReads` to `nullptr` in constructor (`SoloReadFeature.cpp`)
2. Adding null guard in `SoloReadInfoLoader::load()` 
3. Adding defensive checks before `loader.load()` calls

Test results (A375 100K GEX):
- Total reads: 143,229
- CB tags: 60,447 (42%)
- UB tags: 60,430 (42%)

See `docs/BUG_REPORT_NONFLEX_SOLO_CRASH.md` for details.

## Regression Test (2026-01-26)

After the buffered path re-enablement, sorted and unsorted BAM CB/UB tag rates now match:

```
==============================================
CB/UB Tag Injection Regression Test Suite
==============================================
Mode: Quick (2 lanes)

Test 1: Sorted BAM CB/UB Tag Injection - PASS
        Reads: 378186, CB: 300217 (79.38%), UB: 152905 (40.43%)
        Sorted MEX output generated - PASS
        Cells: 13058, UMIs: 140877

Test 2: Unsorted BAM CB/UB Tag Injection (automatic) - PASS
        Reads: 378186, CB: 300217 (79.38%), UB: 152905 (40.43%)
        Unsorted MEX output generated - PASS
        Cells: 13058, UMIs: 140877

Test 3: Sorted vs Unsorted Consistency - PASS
        Sorted and Unsorted MEX outputs are identical
        Cells: 13058, UMIs: 140877

Summary: 5/5 tests passed
```

**Key Result**: Unsorted BAM now has the same CB/UB tag rates as sorted BAM (79.38% CB, 40.43% UB),
confirming the buffered path is working correctly.

## Memory Considerations

The buffered approach adds memory overhead:
- Records are buffered in memory up to `--limitBAMsortRAM`
- Excess records spill to disk in `--outTmpDir`
- Per-read tracking via `readIdTracker_` adds ~12 bytes/read

For large datasets, ensure adequate `--limitBAMsortRAM` and temp disk space.
