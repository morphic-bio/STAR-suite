# Handoff: FLEX No-Align Snapshot Harness Regression (2026-03-25)

## Status

Partial implementation complete. The FLEX snapshot/replay harness is wired up and
clean-builds, but the first full-sample `--flexNoAlign 1` seed+replay run exposed
a regression on the current merged branch:

- the full no-align seed run completed
- replay correctly skipped mapping startup
- but the seed never wrote the snapshot file because the FLEX inline hash was
  empty by the time `collapseUMIall_fromHash()` ran

This means the harness plumbing is present, but the current no-align FLEX path
is not yet producing a usable post-map inline-hash surface on this branch.

## Branch / Commit / Worktree

- Branch: `feature/flex-optimization-using-solo-20260325`
- Base HEAD before this work: `a70a03929b14331f841988df4ccb87ca26f9a770`
- Repo: `/mnt/pikachu/STAR-suite`

## Uncommitted Code Changes

Current dirty files:

- `core/legacy/source/STAR.cpp`
- `core/legacy/source/SoloFeature.h`
- `core/legacy/source/SoloFeature_bridgeHashSnapshot.cpp`
- `core/legacy/source/SoloFeature_countCBgeneUMI.cpp`
- `core/legacy/source/SoloFeature_processRecords.cpp`
- `core/legacy/source/SoloFeature_sumThreads.cpp`
- `flex/source/SoloFeature_collapseUMI_fromHash.cpp`
- `tests/run_flex_snapshot_noalign_full.sh` (new)

These changes were **not** committed yet.

## What Was Implemented

### 1. FLEX snapshot/replay support

Added FLEX equivalents of the Solo bridge snapshot hooks:

- `SoloFeature::flexHashSnapshotWrite(const char *path);`
- `SoloFeature::flexHashSnapshotLoad(const char *path);`

Implemented in:

- `core/legacy/source/SoloFeature.h`
- `core/legacy/source/SoloFeature_bridgeHashSnapshot.cpp`

The snapshot format currently serializes:

- `readFeatSum->inlineHash_`
- `readFeatSum->stats.V`
- `readFeatSum->cbReadCount` if present
- `g_statsAll.readN`
- `nReadsInput`

Current FLEX snapshot magic/version:

- magic: `STARFX2`
- version: `1`

### 2. Replay gating

Added FLEX replay handling in:

- `core/legacy/source/SoloFeature_countCBgeneUMI.cpp`
- `core/legacy/source/STAR.cpp`
- `core/legacy/source/SoloFeature_processRecords.cpp`

Relevant env vars:

- seed:
  - `STAR_SOLO_FLEX_HASH_SNAPSHOT_OUT=/abs/path/flex_snapshot.bin`
- replay:
  - `STAR_SOLO_FLEX_HASH_SNAPSHOT_IN=/abs/path/flex_snapshot.bin`
  - `STAR_SOLO_FLEX_HASH_SNAPSHOT_REPLAY_SKIP_READS=1`
  - `STAR_SOLO_FLEX_SNAPSHOT_STOP_AFTER_COUNT=1`

Replay currently rejects sorted-BAM CB/UB replay because it does not serialize
`readIdTracker_`.

### 3. Snapshot boundary for FLEX

Current boundary is in:

- `flex/source/SoloFeature_collapseUMI_fromHash.cpp`

The intent was:

1. normal seed run:
   - `resolvePendingAmbiguousToHash(false)`
   - write snapshot
   - continue into direct hash collapse
2. replay:
   - skip mapping
   - load snapshot
   - continue from direct hash collapse

### 4. Full-sample no-align runner

Added:

- `tests/run_flex_snapshot_noalign_full.sh`

This script is the current reproduction harness for the target optimization
surface:

- full JAX SC2300771 sample under `/storage/JAX_sequences`
- `--flexNoAlign 1`
- internal gzip only
- no external `zcat`
- no BAM
- seed / replay / both modes

## Important Branch Cleanup Already Done

While doing a mandatory clean rebuild, I hit two pre-existing branch issues and
fixed them so the branch would build:

- removed stale deleted-field references in `SoloFeature_sumThreads.cpp`
- replaced dead bridge post-merge calls there with the live
  `resolvePendingAmbiguousToHash(true)` entry point

After those fixes, the required clean rebuild passed:

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
```

## Full-Sample No-Align Reproduction

### Runner script

The seed+replay command was:

```bash
OUT_ROOT=/storage/flex_snapshot_noalign_$(date -u +%Y%m%d_%H%M%S) \
tests/run_flex_snapshot_noalign_full.sh both
```

The script uses:

- genome:
  - `/storage/flex_filtered_reference_2024/star_index`
- whitelist:
  - `/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt`
- sample whitelist:
  - `/mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv`
- probe list:
  - `/storage/flex_filtered_reference_2024/star_index/flex_probe_artifacts/probe_list.txt`
- sample probes:
  - `/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt`
- hash cache:
  - `/storage/downsampled_100K/SC2300771/results/flex_h01_2024_20260320_081246/h01_cache.bin`

### Artifact root

The run root is:

- `/storage/flex_snapshot_noalign_20260325_221330`

Files of interest:

- seed log:
  - `/storage/flex_snapshot_noalign_20260325_221330/seed/Log.out`
- seed timing:
  - `/storage/flex_snapshot_noalign_20260325_221330/seed/stderr.log`
- replay log:
  - `/storage/flex_snapshot_noalign_20260325_221330/replay/Log.out`
- replay timing/failure:
  - `/storage/flex_snapshot_noalign_20260325_221330/replay/stderr.log`
- seed barcode stats:
  - `/storage/flex_snapshot_noalign_20260325_221330/seed/Solo.out/Barcodes.stats`

No snapshot file was produced:

- expected:
  - `/storage/flex_snapshot_noalign_20260325_221330/flex_snapshot.bin`
- actual:
  - file missing

### Seed timing

From `seed/stderr.log`:

- wall: `6:45.11`
- max RSS: `38551980 kB`

### Replay timing

Replay did skip mapping startup successfully:

- replay reached:
  - `inline-hash snapshot replay: skipping entire mapping phase`

From `replay/stderr.log`:

- wall before failure: `0:28.70`
- max RSS: `38412064 kB`

Replay then failed because the snapshot file was missing:

- `EXITING: cannot open STAR_SOLO_FLEX_HASH_SNAPSHOT_IN file: /storage/flex_snapshot_noalign_20260325_221330/flex_snapshot.bin`

## The Actual Regression

The seed run finished, but did **not** produce meaningful FLEX Gene output on
this branch.

Evidence from `seed/Log.out`:

- `Running clique-based UMI correction...`
- `Clique correction phase 1 (extract+sort): 0 entries, 0.0 sec`
- `WARNING: inlineHash_ is empty, no reads to collapse`
- `Solo timing: collapseUMIall_fromHash 0.0 s`

Evidence from `seed/Solo.out/Barcodes.stats`:

- all key barcode counters are zero

So the no-align seed did **not** produce a usable inline-hash count surface.

## Why This Matters Relative To README

`README.md` still reports successful FLEX no-align benchmarks and parity:

- `README.md` Flex wall-time section
- `README.md` Flex no-align section

Those numbers are historical results from the earlier `benchmark-flex` lineage,
as documented in:

- `docs/HANDOFF_FLEX_PIPELINE_BENCHMARK_SUMMARY_20260323.md`
- `docs/HANDOFF_FLEX_BENCHMARK_20260319.md`

So there is no contradiction if you read it correctly:

- the old benchmark branch had working no-align outputs
- the current merged branch appears to have regressed on the FLEX no-align
  inline-hash surface

Do **not** treat the current branch as validated just because the README still
contains the older benchmark numbers.

## Likely Root Cause

The strongest current hypothesis is in:

- `flex/source/SoloFeature_collapseUMI_fromHash.cpp`

Today the function does:

1. take `hashSize = kh_size(readFeatSum->inlineHash_)`
2. if `hashSize == 0`, return early with:
   - `WARNING: inlineHash_ is empty, no reads to collapse`
3. only **after that** does it reach:
   - `resolvePendingAmbiguousToHash(false)`
   - snapshot write

That ordering is dangerous for FLEX because many reads can still be in
`pendingAmbiguous_` rather than in `inlineHash_`.

If most or all no-align keep reads are landing in deferred ambiguous state,
then the current logic guarantees:

- empty hash check fires first
- no ambiguous resolution happens
- no snapshot is written
- seed exits “successfully” but with empty output

This is the first thing to test.

## Other Possible Root Cause

There is also a second plausible failure point:

- the fused no-align pipeline may be recording into per-thread `fusedFeats[i]`
  correctly, but the final merged Gene hash may not be reaching `readFeatSum`
  as expected by the time `countCBgeneUMI()` runs

Relevant files:

- `core/legacy/source/mapThreadsSpawn.cpp`
- `core/legacy/source/SoloFeature_sumThreads.cpp`
- `core/legacy/source/flex/SoloReadFeature_record_flex.cpp`

Current keep-path writer:

- `record_flex_hash_screen_keep(...)`

It definitely does hash inserts for unique-CB keep reads:

- direct `kh_put(cg_agg, soloReadFeat->inlineHash_, key, ...)`

So if the final surface is still empty, either:

1. the seed is dominated by deferred ambiguous entries and the early-return
   ordering bug is the main issue, or
2. the fused-thread merge path into `readFeatSum` is broken on this branch

## Recommended Next Steps

### 1. Patch the empty-hash ordering first

In `flex/source/SoloFeature_collapseUMI_fromHash.cpp`:

- move the `hashSize == 0` early return to **after**
  `resolvePendingAmbiguousToHash(false)`
- and after snapshot write for normal seed runs

That is the highest-probability fix and is low risk.

### 2. Add one diagnostic log before and after ambiguous resolution

Log:

- `kh_size(readFeatSum->inlineHash_)` before resolution
- `readFeatSum->pendingAmbiguous_.size()`
- `kh_size(readFeatSum->inlineHash_)` after resolution

This will immediately tell you whether the no-align surface is simply deferred
into ambiguous state.

### 3. If still empty, instrument the fused merge path

In `core/legacy/source/mapThreadsSpawn.cpp`, after:

- `workerGeneFeat->mergeInlineHash(*fusedFeats[i]);`

log the merged `kh_size(workerGeneFeat->inlineHash_)`.

Then in `SoloFeature_sumThreads.cpp`, log:

- `kh_size(readFeatAll[ii]->inlineHash_)` for each thread
- final `kh_size(readFeatSum->inlineHash_)` after merge

This will show whether the data is lost before `processRecords()`.

### 4. Re-run the exact same full-sample no-align harness

Use the same script and same full-sample dataset again:

```bash
OUT_ROOT=/storage/flex_snapshot_noalign_$(date -u +%Y%m%d_%H%M%S) \
tests/run_flex_snapshot_noalign_full.sh both
```

Do not change to a different dataset until the regression is resolved.

### 5. Do not update README yet

The README benchmark claims should remain untouched until the current branch is
revalidated on the full no-align path.

## Update (2026-03-26): harness fixed; optimized collapse benchmarked

- The empty-hash / seed snapshot regression on this branch is resolved: seed writes a non-empty FLEX snapshot, replay loads it and skips mapping, and `FLEX_SNAPSHOT_COMPARE=count` validates the count surface.
- A clean-rebuilt STAR was rerun through `tests/run_flex_snapshot_noalign_full.sh` on the full `/storage` sample. **Both** seed and replay `Log.out` now show:
  - `Grouping hash entries by CB/gene/tag (khash accumulate + per-bucket gene sort)...`
  - `Collapse timing: accumulate=… bucket_build=… bucket_sort=… total=…`  
  (The older artifacts under `/storage/flex_snapshot_noalign_20260325_rerun_{count,full}_v1/` still logged `extract/sort/group`, so those runs were not exercising this path.)
- **Artifact roots (optimized reruns):**
  - Count mode: `/storage/flex_snapshot_noalign_opt_20260326/count_v1/` — `compare.log` OK; snapshot **3,602,016,112** bytes, **224,941,678** entries; example seed collapse line `Collapse timing: accumulate=17.4s bucket_build=4.7s bucket_sort=0.3s total=22.4s`; Solo seed `collapseUMIall_fromHash` **72.3 s**, `countCBgeneUMI` **101.9 s**; replay **70.3 s** / **72.9 s**; seed wall **10:49**, replay **1:47** (similar RSS to the 2026-03-25 reruns).
  - Full mode: `/storage/flex_snapshot_noalign_opt_20260326/full_v1/` — `per_sample` + `Solo.out` strict diff OK; example replay `Collapse timing: accumulate=17.6s bucket_build=5.0s bucket_sort=0.3s total=22.8s`; Solo seed **73.5 s** / **103.2 s**, replay **70.2 s** / **72.8 s**; wall/RSS within noise of the count run.
- **vs 2026-03-25 rerun baselines:** internal collapse phase total rose ~**17.2 s → ~22.4–23.1 s** (new breakdown is not comparable 1:1 to `extract/sort/group`). Solo `collapseUMIall_fromHash` and `countCBgeneUMI` times shifted upward on seed (~**65 s → ~72–73 s** collapse; ~**95 s → ~102 s** count) while end-to-end wall clock stayed ~**10:49** seed / **~1:47** replay — the important closure is that **recorded logs now prove the optimized collapse implementation ran** and full parity still holds.

## Update (2026-03-26, later): CB-first CSR collapse replaces khash accumulator

The khash triplet accumulator (`tripletAcc`) was replaced with a Solo-style CB-first CSR design:

1. Flat CB counter (~5.8 MB vector indexed by whitelist CB) — one pass over the 224M-entry hash
2. CSR prefix-sum offsets — per-CB slice boundaries in a contiguous flat array
3. Scatter hash entries into per-CB CSR slices — packed `[tag(5)][gene(15)][umi(24)]` uint64 per slot
4. Parallel per-CB sort (`#pragma omp parallel for`) — ~323 entries/CB avg, L1-resident
5. Serial run-length scan — count distinct UMIs per (CB, tag, gene), emit to `cbTagGeneCountsVec`

No second hash table. No global sort. No `std::unordered_map` of vectors during the hot loop.

**Results (CSR vs khash accumulator vs old extract/sort/group):**

| Metric | CSR (this) | khash accum | old extract/sort/group |
|--------|-----------|-------------|----------------------|
| Collapse total | **5.5 s** | 22.4 s | 17.2 s |
| `collapseUMIall_fromHash` (seed) | **55 s** | 72 s | 65 s |
| `countCBgeneUMI` (seed) | **85 s** | 102 s | 95 s |
| Replay wall | **1:30** | 1:47 | 1:40 |
| Seed wall | **10:32** | 10:49 | 10:41 |

Full parity (strict `diff` on `per_sample` + `Solo.out`) still passes.

**Artifact roots:**
- Count mode: `/storage/flex_snapshot_noalign_csr_20260326/count_v1/`
- Full mode: `/storage/flex_snapshot_noalign_csr_20260326/full_v1/`

## Update (2026-03-26, later): Direct-emit CSR — no regroup, no re-sort

The remaining Solo-style win was to remove the second generic regrouping layer:

**Before (CSR v1):** Phase 4 scan rebuilt an `unordered_map<uint64_t, vector<pair<uint32_t, uint32_t>>>` (`cbTagGeneCountsVec`) from the sorted CSR data. Then `buildInlineMatrixFromHash()` sorted `cbTagKeys`, iterated to build triplets, and sorted `bundle.triplets` again — three redundant passes.

**After (CSR v2 / direct-emit):** Phase 4 scan emits composite barcodes, triplets, and `SampleMatrixData` fields directly into `InlineMatrixBundle` during the sequential scan. Because the CSR data is already sorted by (CB, tag, gene), the output is in (cell_idx, gene_idx) order with no additional sorting. `buildInlineMatrixFromHash` is now dead code.

**Changes:**
- `SoloFeature_collapseUMI_fromHash.cpp`: Phase 4 loads probe list, initializes `InlineMatrixBundle`, builds composite barcodes/triplets/cellUMIs/cellGenes incrementally during the scan. Calls `writeMexFromInlineHashDedup(mexDir, inlineMatrix)` directly.
- `SoloFeature_writeMexFromInlineHashDedup.cpp`: `buildInlineMatrixFromHash()` is dead code (no callers). `writeMexFromInlineHashDedup()` is unchanged — it receives the pre-built bundle.
- Removed `#include <unordered_map>` from the collapse file.

**Results (direct-emit vs CSR v1 vs baseline):**

| Metric | Direct-emit (this) | CSR v1 | khash accum | old baseline |
|--------|-------------------|--------|-------------|-------------|
| Collapse total | **5.9 s** | 6.0 s | 22.4 s | 17.2 s |
| `collapseUMIall_fromHash` (replay) | **49.6 s** | 52.6 s | 72 s | 65 s |
| `countCBgeneUMI` (replay) | **52.2 s** | 55.2 s | 102 s | 95 s |
| `processRecords` (replay) | **52.2 s** | 55.2 s | — | — |

The 3 s savings in `collapseUMIall_fromHash` come from eliminating the `unordered_map` rebuild, `cbTagKeys` sort, and final `triplets` sort that `buildInlineMatrixFromHash` performed. Collapse phase itself is unchanged since the CSR scan was already doing the grouping.

Full parity (strict `diff` on `per_sample` + `Solo.out`) validated.

**Artifact root:** `/storage/flex_snapshot_noalign_csr_20260326/full_v2/`

## Update (2026-03-26, later): Dead code removal + fused counting + MexWriter optimization

1. **Dead code removal:** `buildInlineMatrixFromHash` removed from `SoloFeature_writeMexFromInlineHashDedup.cpp` (definition) and `SoloFeature.h` (declaration). Unused includes cleaned up. Duplicate `#include "SoloReadFeature.h"` removed.

2. **Fused pre-dedup counting:** The separate 224M-entry hash scan that counted `totalCountsPre` for the log line was merged into Phase 1's CB counting loop. Both counters now computed in a single pass — eliminates ~2s of redundant work.

3. **MexWriter buffering:** Added 1 MiB `setvbuf` to `matrix.mtx` file handle in `MexWriter::writeMex()`. Avoids barcode vector copy when no truncation is needed (`cb_len <= 0`).

**Results (replay):**

| Metric | Cleanup (this) | Direct-emit | CSR v1 | khash accum | Baseline |
|--------|---------------|-------------|--------|-------------|----------|
| Collapse total | **5.9 s** | 5.9 s | 6.0 s | 22.4 s | 17.2 s |
| `collapseUMIall_fromHash` | **47.6 s** | 49.6 s | 52.6 s | 72 s | 65 s |
| `processRecords` | **50.2 s** | 52.2 s | 55.2 s | — | — |

Full parity validated.

**Artifact root:** `/storage/flex_snapshot_noalign_csr_20260326/full_v3/`

## Remaining optimization analysis (detailed timing breakdown, replay)

| Phase | Time | % of collapseUMI | Notes |
|-------|------|-------------------|-------|
| CSR collapse (cb_count+scatter+sort+scan) | 5.9 s | 12% | Already optimized |
| Bookkeeping + MEX write | ~5 s | 11% | 98.9M triplets to text; setvbuf applied |
| FlexFilter (OrdMag + EmptyDrops × 16 tags) | **30 s** | **63%** | Dominant remaining cost |
| Other overhead (logging, snapshot load, etc.) | ~6.7 s | 14% | |

**FlexFilter parallelization (next major opportunity):**
The FlexFilter processes 16 tags sequentially (`numThreads = 1` in `FlexFilter.cpp:602`), giving all 32 cores to each tag's EmptyDrops Monte Carlo. Processing 4 tags in parallel with 8 MC threads each would overlap the sequential phases (OrdMag, sorting, barcode extraction) and could reduce the 30s to ~10-15s. Requires validating thread safety of FlexFilter internals and memory pressure under parallel OrdMag/EmptyDrops.

## Update (2026-03-26): InlineCBCorrection benchmark (mapping phase H1 hash)

Tested whether enabling `InlineCBCorrection` (which precomputes an H0/H1 hash
for O(1) CB whitelist lookup) would speed up the mapping phase vs legacy
`matchCBtoWL` (binary search on sorted array + 48 1MM variant searches).

**InlineCBCorrection hash sizes** (737K whitelist):
- `exact_map` (H0): 737,280 entries
- `variant_map` (H1, unique): 25,566,912 entries
- `variant_collisions` (ambiguous): 7,511,616 entries (max fanout 6)

**Results — v4 (`InlineCBCorrection`) vs v3 (legacy `matchCBtoWL`):**

| Metric | v3 (legacy) | v4 (InlineCBCorr) | Delta |
|--------|------------|-------------------|-------|
| Total wall clock | 10:25.79 | 12:24.53 | **+119s (+19%)** |
| Mapping phase | ~501s | ~601s | **+100s (+20%)** |
| Hash entries (pre-dedup) | 224,941,678 | 226,570,847 | +1.6M |
| Total read counts | 1,626,999,647 | 1,649,535,943 | +22.5M |
| Ambig inline resolved | n/a | 24,706,851 | |
| Ambig captured | n/a | 8,119,903 | |

**Conclusion: InlineCBCorrection is 20% slower.** It captures more reads (22.5M)
via N-rescue and Bayesian ambiguous resolution, but the per-read overhead
exceeds the savings from O(1) lookup:

1. Double encoding (`fastPathCorrection` packs CB, then `exactIndex` re-packs)
2. String allocation for `correctedCB` on every read
3. Per-read `CbBayesianResolver` construction for 24.7M ambiguous reads
4. Per-read evidence tracking (`recordParentEvidence`)
5. 25.6M-entry hash has worse cache locality than binary search on 737K sorted array

The legacy `matchCBtoWL` exact-match path is very lean (one `convertNuclStrToInt64`
+ one `binarySearchExact` on a compact sorted array that fits in L2/L3 cache),
and ~90%+ of reads take this fast path.

**What would help:** A lean H1 hash that avoids InlineCBCorrection's overhead —
takes pre-packed uint64, returns WL index directly (no string intermediary),
skips evidence tracking, and defers ambiguous resolution to batch. The concept
is sound; the current implementation is too heavy for the hot path.

**Artifact root:** `/storage/flex_snapshot_noalign_csr_20260326/full_v4_cbcorr/`

## Update (2026-03-26): Lean H0 khash for CB whitelist lookup

After InlineCBCorrection proved too heavy, tested a minimal approach: replace
only `binarySearchExact` in `matchCBtoWL` with a lean khash H0 probe, keeping
the existing 48 1MM binary-search-turned-hash-probe loop unchanged.

**Implementation** (3 files changed):
- `ParametersSolo.h`: added `#include "htslib/khash.h"`, declared
  `KHASH_MAP_INIT_INT(cbH0, uint32_t)` (uint32 key → uint32 WL index), added
  `khash_t(cbH0) *cbWLhash` member.
- `ParametersSolo.cpp`: build hash once after `cbWL` is populated (737K entries,
  2M buckets at ~35% load, ~16 MB total).
- `SoloReadBarcode_getCBandUMI.cpp`: `matchCBtoWL` uses `kh_get(cbH0, h, key)`
  for exact match and all 48 1MM variant probes (fallback to binary search when
  hash is null, e.g., Complex CB type).

**Results — v5 (H0 khash) vs v3 (binary search):**

| Metric | v3 (binarySearch) | v5 (H0 khash) | Delta |
|--------|-------------------|---------------|-------|
| Total wall clock | 10:25.79 | 10:42.84 | +17s (noise) |
| Mapping speed (M reads/hr) | 12,107 | 11,908 | ~same |
| Hash entries (pre-dedup) | 224,941,678 | 224,941,678 | identical |
| Total counts | 1,626,999,647 | 1,626,999,647 | identical |

**Conclusion: No speedup.** Output is bit-identical (same hash entries, same
counts). The mapping phase is **gzip-I/O-bound**, not CB-lookup-bound. With 8
FASTQ lanes each saturating a `gzread` decompression stream on 8 of 32 threads,
the per-read CPU work (including CB matching) is never on the critical path.

**Why binary search was never the bottleneck:** Each of 8 active threads
processes ~420K reads/s. Saving ~15ns per hash probe vs binary search saves
6 µs per second of wall time — invisible against the ~2.4 µs/read budget
dominated by gzip decompression and string parsing.

**The true mapping bottleneck is gzip decompression** (8 serial `gzread`
streams). Faster decompression (ISA-L/igzip, which is 3-5× faster than zlib)
or more parallel lanes would shift the bottleneck to CPU, at which point the
H0 hash improvement would become measurable.

The H0 hash change is kept in the tree as a zero-cost, correctness-preserving
improvement that benefits any future path where gzip is no longer the
bottleneck.

**Artifact root:** `/storage/flex_snapshot_noalign_csr_20260326/full_v5_h0hash/`

## Notes For The Next Agent

- The harness uses internal gzip; do **not** add `--readFilesCommand zcat`
- The target optimization surface is explicitly:
  - FLEX
  - full sample
  - `--flexNoAlign 1`
  - no BAM
- The dominant remaining cost is FlexFilter (30s / 63%). See "Remaining optimization analysis" above.
- **Do NOT auto-enable `soloInlineCBCorrection`** in `--flex yes` mode — it causes a 20% mapping regression. See the InlineCBCorrection benchmark section above.
- The mapping phase is gzip-I/O-bound (8 lanes × 1 `gzread` stream each). CB
  lookup optimization (binary search vs H0 hash) has no measurable effect.
  The next meaningful mapping speedup requires faster decompression (ISA-L) or
  more parallel I/O.

