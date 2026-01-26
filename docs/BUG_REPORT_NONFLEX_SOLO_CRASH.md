# Bug Report: Non-Flex STARsolo Crash with CB/UB Tags

**Date**: 2026-01-25  
**Branch**: `bamcbub`  
**Commit**: `71a4cdb` (Add CB/UB tag injection for sorted BAM with Flex support)  
**Severity**: High (crash)  
**Affected**: Non-Flex STARsolo with `--outSAMattributes ... CB UB`  
**Status**: **FIXED** (2026-01-25)

## Summary

Non-Flex STARsolo crashes with a segmentation fault during Solo counting when CB/UB tags are requested in `--outSAMattributes`. The crash occurs regardless of BAM output type (sorted or unsorted). Flex mode (`--flex yes`) is NOT affected.

## Reproduction

```bash
# Minimal reproduction (crashes)
STAR \
  --runThreadN 4 \
  --genomeDir /storage/autoindex_110_44/bulk_index \
  --readFilesIn \
    /storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/gex/downsampled_100000_v2/1k_CRISPR_5p_gemx_gex_S2_L001_R2_001.fastq.gz \
    /storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/gex/downsampled_100000_v2/1k_CRISPR_5p_gemx_gex_S2_L001_R1_001.fastq.gz \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloCBwhitelist /storage/A375/3M-5pgex-jan-2023.txt \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloMultiMappers Rescue \
  --soloCellFilter None \
  --clipAdapterType CellRanger4 \
  --soloFeatures Gene \
  --soloStrand Unstranded \
  --alignEndsType Local \
  --chimSegmentMin 1000000 \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix /tmp/test/
```

**Output:**
```
Jan 25 21:22:10 ..... started Solo counting
Jan 25 21:22:10 ... Starting Solo post-map for Gene
Jan 25 21:22:10 ... Allocated and initialized readInfo array, nReadsInput = 100001
Segmentation fault (core dumped)
```

## Stack Trace (GDB)

```
#0  soloInputFeatureUMI(std::fstream*, int, bool, ...) ()
#1  SoloReadInfoLoader::load(SoloReadFeature&, SoloReadInfoMode, ...) ()
#2  SoloFeature::countCBgeneUMI() ()
#3  SoloFeature::processRecords() ()
#4  Solo::processAndOutput() ()
#5  main ()
```

## Root Cause Analysis

### Crash Location

The crash occurs in `flex/source/SoloReadInfoLoader.cpp` at lines 41-42:

```cpp
void SoloReadInfoLoader::load(SoloReadFeature &rf, ...) {
    rf.streamReads->flush();      // <-- CRASH: streamReads is invalid
    rf.streamReads->seekg(0,std::ios::beg);
    ...
}
```

### Why `streamReads` is Invalid

In `SoloReadFeature` constructor (`core/legacy/source/SoloReadFeature.cpp`):

```cpp
SoloReadFeature::SoloReadFeature(int32 feTy, Parameters &Pin, int iChunk)
             : featureType(feTy), P(Pin), pSolo(P.pSolo), inlineHash_(nullptr), readIdTracker_(nullptr)
{
    ...
    if (pSolo.inlineHashMode) {
        inlineHash_ = kh_init(cg_agg);
        streamReads = nullptr;  // Explicitly set to nullptr
        ...
    } else if (iChunk>=0) {
        streamReads = &fstrOpen(...);  // Opened for per-thread objects
    }
    // NOTE: streamReads is UNINITIALIZED if inlineHashMode=false AND iChunk<0
}
```

For non-Flex mode:
- `inlineHashMode` = false
- Per-thread `SoloReadFeature` objects have `iChunk >= 0` → `streamReads` is opened
- But the code path reaches `SoloReadInfoLoader::load` with potentially invalid `streamReads`

### Suspected Issue

The changes in commit `71a4cdb` modified `SoloFeature_countCBgeneUMI.cpp`:

**Before (perturb branch):**
```cpp
if (pSolo.readInfoYes[featureType] && !(pSolo.soloFlexMinimalMemory && pSolo.inlineHashMode)) {
    resetPackedStorage(nReadsInput);
    ...
}
```

**After (bamcbub branch):**
```cpp
bool needPackedReadInfo = pSolo.readInfoYes[featureType] || pSolo.trackReadIdsForTags;
bool skipForMinimalMemory = pSolo.soloFlexMinimalMemory && pSolo.inlineHashMode && !pSolo.trackReadIdsForTags;
if (needPackedReadInfo && !skipForMinimalMemory) {
    resetPackedStorage(nReadsInput);
    ...
}
```

The logic change may have altered which code paths are executed in non-Flex mode, exposing a latent bug where `SoloReadInfoLoader::load` is called without proper null-checking of `streamReads`.

## Verification

| Branch | Non-Flex A375 | Flex 100K |
|--------|---------------|-----------|
| `perturb` | PASS | PASS |
| `bamcbub` (71a4cdb) | CRASH | PASS |
| `bamcbub` + unsorted | CRASH | PASS |

The crash exists in the committed `bamcbub` branch **before** the unsorted BAM changes.

## Files Involved

1. **`flex/source/SoloReadInfoLoader.cpp`** - Missing null check for `streamReads`
2. **`core/legacy/source/SoloReadFeature.cpp`** - `streamReads` uninitialized in edge case
3. **`core/legacy/source/SoloFeature_countCBgeneUMI.cpp`** - Modified logic in bamcbub

## Fix Applied

Three fixes were implemented:

### 1. Initialize `streamReads` to `nullptr` in constructor

**File**: `core/legacy/source/SoloReadFeature.cpp`

```cpp
SoloReadFeature::SoloReadFeature(int32 feTy, Parameters &Pin, int iChunk)
             : featureType(feTy), P(Pin), pSolo(P.pSolo), 
               streamReads(nullptr),  // <-- Added explicit initialization
               inlineHash_(nullptr), readIdTracker_(nullptr)
{
```

### 2. Add guard in `SoloReadInfoLoader::load`

**File**: `flex/source/SoloReadInfoLoader.cpp`

```cpp
void SoloReadInfoLoader::load(SoloReadFeature &rf, ...) {
    // Guard against null streamReads
    if (rf.streamReads == nullptr) {
        std::cerr << "EXITING because of fatal ERROR: SoloReadInfoLoader::load called with null streamReads\n"
                  << "This typically happens when:\n"
                  << "  1. Inline hash mode is active but loader was called\n"
                  << "  2. SoloReadFeature was constructed with iChunk < 0\n"
                  << "  3. readFeatAll[ii] assignment mismatch\n"
                  << "featureType=" << rf.featureType << " inlineHashMode=" << rf.pSolo.inlineHashMode << "\n";
        exit(1);
    }
    rf.streamReads->flush();
    ...
}
```

### 3. Add defensive check before `loader.load()`

**File**: `core/legacy/source/SoloFeature_countCBgeneUMI.cpp`

```cpp
for (int ii=0; ii<P.runThreadN; ii++) {
    // Defensive check: verify readFeatAll[ii] and its streamReads are valid
    if (readFeatAll[ii] == nullptr) {
        // Exit with error about null readFeatAll
    }
    if (!pSolo.inlineHashMode && readFeatAll[ii]->streamReads == nullptr) {
        // Exit with error about null streamReads in non-inline-hash mode
    }
    loader.load(*readFeatAll[ii], ...);
}
```

## Verification

After fix:

| Test | Result |
|------|--------|
| A375 non-Flex (no CB/UB) | PASS |
| A375 non-Flex (CB/UB tags) | PASS |
| Flex sorted BAM (CB/UB) | PASS |
| Flex unsorted BAM (CB/UB) | PASS |

A375 output statistics:
- Total reads: 143,229
- CB tags: 60,447 (42%)
- UB tags: 60,430 (42%)

## Impact

- **Flex mode**: Not affected (works correctly)
- **Non-Flex mode**: Crashes when CB/UB tags requested
- **Workaround**: Do not request CB/UB tags in `--outSAMattributes` for non-Flex runs

## Test Scripts

- `tests/run_unsorted_cbub_a375_nonflex.sh` - Documents the issue
- `tests/run_cbub_regression_test.sh` - Flex tests (passing)

## Related

- Commit: `71a4cdb` - Add CB/UB tag injection for sorted BAM with Flex support
- Branch: `bamcbub`
- Prior working branch: `perturb`
