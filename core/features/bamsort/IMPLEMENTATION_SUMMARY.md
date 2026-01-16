# Internal Spill-to-Disk Sorter Implementation Summary

## Overview

This document summarizes the implementation of an internal spill-to-disk coordinate-sorting backend for STAR-Flex, selectable via `--outBAMsortMethod samtools`. The implementation uses raw BAM buffer sorting with bounded-memory spill-to-disk and k-way merge, implemented as an internal STAR component (no samtools code vendoring required).

**Key Design Principles:**
- Raw buffer sorting: Parse BAM core fields directly from raw buffers without converting to `bam1_t` structures
- Memory-bounded: Check memory usage and spill to disk when `--limitBAMsortRAM` is exceeded
- QNAME tie-breaking: Deterministic ordering via byte-level QNAME comparison when numeric keys are equal (REQUIRED)
- No BAM mutation: Sorting metadata stored separately, never modifies BAM record content
- Existing BGZF I/O: Reuses STAR's existing `bgzf_write()` and `outBAMwriteHeader()` paths

**Non-Goals:**
- No htslib upgrade: Uses existing vendored htslib BGZF APIs only
- No samtools code vendoring: Sorting logic implemented internally
- No samtools ordering parity: Deterministic ordering defined by STAR spec

---

## Files Created

### Core Implementation

1. **`source/SamtoolsSorter.h`**
   - `SortKey` struct: Numeric fields (tid, pos, flag, mtid, mpos, isize) for sorting
   - `BAMRecord` struct: Holds raw BAM data with SortKey, uses move semantics
   - `BAMRecordComparator`: Compares SortKey + QNAME bytes for deterministic ordering
   - `SamtoolsSorter` class: Main sorter with spill-to-disk and k-way merge
   - `HeapEntry` struct: For k-way merge with QNAME tie-breaking
   - `HeapLess` comparator: Min-heap comparator with REQUIRED QNAME tie-breaking
   - Spec comment header: Documents QNAME tie-breaking as hard requirement
   - Assert guards: `assert(bamPtr != nullptr)` in heap comparator

2. **`source/SamtoolsSorter.cpp`**
   - `computeSortKey()`: Parses BAM core fields from raw buffer (bam32[1]=tid, bam32[2]=pos, bam32[4]=flag, bam32[6]=mtid, bam32[7]=mpos, bam32[8]=isize)
   - `compareQName()`: Unified QNAME byte comparison helper (FNV-1a not used, direct byte comparison)
   - `addRecord()`: Thread-safe buffer with memory check, triggers spill when threshold exceeded
   - `sortAndSpill()`: Sorts in-memory buffer and writes to temp file with explicit key serialization
   - `finalize()`: Initializes k-way merge
   - `nextRecord()`: Streams records from k-way merge with QNAME tie-breaking
   - `SpillFileReader`: Reads sorted chunks from spill files for k-way merge

### Integration Points

3. **`source/GlobalVariables.h`**
   - Added: `extern SamtoolsSorter* g_samtoolsSorter;`

4. **`source/GlobalVariables.cpp`**
   - Added: `SamtoolsSorter* g_samtoolsSorter = nullptr;`

5. **`source/STAR.cpp`**
   - Initializes `g_samtoolsSorter` before alignment threads start
   - Cleans up `g_samtoolsSorter` after sorting completes

6. **`source/BAMoutput.cpp`**
   - Modified `coordOneAlign()`: Branches to `g_samtoolsSorter->addRecord()` when `P.outBAMsortMethod == "samtools"`

7. **`source/bamSortByCoordinate.cpp`**
   - Added `bamSortSamtoolsFinalize()`: Handles final merge and output streaming
   - Implements Y/noY splitting with explicit `keepBAM` semantics
   - Uses `genome.yTids` for Y-chromosome detection
   - Modified `bamSortByCoordinate()`: Branches to samtools finalization path

8. **`source/Parameters.h`**
   - Added: `string outBAMsortMethod;` // "star" (default) or "samtools"

9. **`source/Parameters.cpp`**
   - Registered `outBAMsortMethod` parameter
   - Added validation: accepts "star" or "samtools" (case-insensitive)

10. **`source/parametersDefault`**
    - Documented `outBAMsortMethod` parameter with default value "star"

11. **`source/Genome.h`** (pre-existing)
    - Uses: `std::unordered_set<int> yTids;` // tid indices for Y-chromosome contigs (pre-existing)

12. **`source/Genome.cpp`** (pre-existing)
    - Uses: `yTids` set populated during genome loading (pre-existing, used by `bamSortByCoordinate.cpp`)

13. **`source/Makefile`**
    - Added `SamtoolsSorter.o` to OBJECTS list
    - Added `test_SamtoolsSorter` target for unit testing

### Testing

14. **`test/SamtoolsSorter_test.cpp`**
    - Unit test harness with raw BAM record construction
    - Tests QNAME tie-break (readA before readB with same numeric keys)
    - Tests unmapped-last behavior (tid=-1 records sort last)
    - Tests spill handling with tiny maxRAM

15. **`tests/run_ychrom_spill_sorter.sh`**
    - Fixture-based smoke test using `tests/ychrom_test/`
    - Tests spill vs no-spill determinism
    - Validates Y/noY split outputs
    - Checks `SO:coordinate` header
    - Verifies record counts

16. **`tests/run_solo_smoke.sh`**
    - Self-contained scRNA-seq smoke test
    - Builds tiny genome, whitelist, and reads
    - Tests STARsolo + coordinate sorting
    - **Dependencies**: Requires `samtools` command and built `source/STAR` binary

---

## Key Implementation Details

### Sorting Specification (REQUIRED)

**QNAME tie-breaking is part of STAR's sort spec (non-optional).**

When numeric keys are equal (tid, pos, flag, mtid, mpos, isize), ordering MUST be deterministic by QNAME byte comparison. This is enforced via:
- Spec comment header in `SamtoolsSorter.h`
- `assert(bamPtr != nullptr)` guards in heap comparator
- Unified `compareQName()` helper used by both `BAMRecordComparator` and `HeapLess`

### Core Field Extraction

BAM core fields parsed directly from raw buffer:
- `bam32[1]=tid`, `bam32[2]=pos`
- `bam32[4]=flag_nc` (flag extracted from upper 16 bits)
- `bam32[6]=mtid`, `bam32[7]=mpos`, `bam32[8]=isize`

Unmapped reads (`tid=-1` or `pos=-1`) normalized to `INT32_MAX` to sort last.

### Spill-to-Disk Format

Explicit serialization to avoid padding/ABI issues:
```
[bamSize:uint32][hasY:uint8][key.tid:int32][key.pos:int32][key.flag:uint16][key.mtid:int32][key.mpos:int32][key.isize:int32][bamData:bytes]
```

Note: This is a temp-only format, not a stable ABI.

### K-Way Merge

- Uses min-heap (`std::priority_queue` with `HeapLess` comparator)
- Each `HeapEntry` contains `SortKey`, `sourceId`, `recordIdx`, and `bamPtr`
- QNAME comparison happens in comparator when numeric keys are equal
- Streams records without loading all spill files into memory simultaneously

### Y/noY Split Behavior

| emitNoYBAM | keepBAM | Output Files |
|------------|---------|--------------|
| no | - | `Aligned.sortedByCoord.out.bam` (primary) |
| yes | no | `_Y.bam` + `_noY.bam` only (no primary) |
| yes | yes | `_Y.bam` + `_noY.bam` + `Aligned.sortedByCoord.out.bam` |

Y-chromosome detection uses `genome.yTids` for consistency with existing STAR logic.

---

## Compilation Status

✅ All source files compile successfully with no errors or warnings
✅ `make test_SamtoolsSorter` target added (requires `libflex.a`, `libtrim.a`, `libem.a` to be built)
✅ Integration with existing STAR build system complete

**Note**: To run `make test_SamtoolsSorter`, ensure libraries are built:
```bash
cd source
make libflex libtrim libem
make test_SamtoolsSorter
```

**Test Script Dependencies**:
- `tests/run_ychrom_spill_sorter.sh`: Requires `samtools` command and built `source/STAR` binary
- `tests/run_solo_smoke.sh`: Requires `samtools` command and built `source/STAR` binary

---

## Testing Strategy

### Tier 1: Unit Test
- `make test_SamtoolsSorter` runs `test/SamtoolsSorter_test.cpp`
- Validates QNAME tie-break, unmapped-last, and spill handling

### Tier 2: Fixture-Based Smoke Test
- `tests/run_ychrom_spill_sorter.sh` uses `tests/ychrom_test/`
- Tests spill vs no-spill determinism and Y/noY routing

### Tier 3: scRNA-seq Smoke Test (optional)
- `tests/run_solo_smoke.sh` validates STARsolo + coordinate sorting

---

## Backward Compatibility

- Legacy STAR bin sorter remains default (`--outBAMsortMethod star`)
- No behavior change for users who don't specify the new parameter
- Legacy code path completely unchanged
- Rollback: Set `--outBAMsortMethod star` to use legacy sorter

---

## Notes

- **Thread-safety**: `addRecord()` uses mutex for thread-safe buffer access
- **Memory accounting**: Approximate (`bamSize + sizeof(BAMRecord) + 64` overhead per record)
- **Temp file cleanup**: Handles failures gracefully (cleanup in destructor)
- **BAM output**: Uses `bgzf_write()` + `outBAMwriteHeader()` (consistent with existing STAR code)
- **No htslib upgrade required**: Uses existing BGZF I/O APIs only

---

## Recent Fixes (December 20, 2025)

1. **limitBAMsortRAM initialization**: Moved fallback logic earlier in `STAR.cpp` so that `limitBAMsortRAM == 0` is resolved before `SamtoolsSorter` initialization, ensuring spill-to-disk triggers correctly.

2. **Test script hardening**: Added guards to `run_ychrom_spill_sorter.sh` and `run_solo_smoke.sh`:
   - Check for `samtools` command availability
   - Verify STAR binary exists
   - Verify STAR index exists (for ychrom test)

3. **Unit test enhancement**: Added spill file check in `SamtoolsSorter_test.cpp` to verify spill-to-disk is triggered with tiny maxRAM. Made spill check assertive (fails if no spills) and added temp dir cleanup.

4. **Thread-safety fix**: Added `spillFileCounterMutex_` to protect `spillFileCounter_` from concurrent access in `sortAndSpill()`, preventing filename collisions.

5. **Error handling**: Changed silent failure to `exitWithError()` in `bamSortSamtoolsFinalize()` when sorter is null.

6. **Makefile test target**: Updated `test_SamtoolsSorter` to link against `libflex.a`, `libtrim.a`, `libem.a` (required for linking).

---

## Implementation Date

December 20, 2025

---

## Related Documentation

- Plan: `.cursor/plans/samtools_bam_sort_backend_f3c39581.plan.md`
- Original request: `plans/bam-sort-samtools-backend.md`

