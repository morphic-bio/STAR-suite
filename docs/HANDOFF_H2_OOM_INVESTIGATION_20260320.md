# Handoff: H2 Hash Cache Generation OOM — Memory Budget Mismatch

**Date**: 2026-03-20
**Branch**: `benchmark-flex`
**Status**: Streaming workaround in place; dedup container rewritten to packed khash, but full H2 OOM root cause is still open

## Problem

During H2 (double-mismatch) hash cache generation, the STAR process was
OOM-killed at **~67 GB RSS** after ~15 wall-minutes of H2 work (30 threads).

The expected memory budget does **not** justify 67 GB:

| Component | Size |
|---|---|
| Genome (loaded, shared) | ~35 GB |
| H0+H1 records in `threadRecords` (8.8M × 24 B) | ~0.2 GB |
| H2 KEEP records (worst case 589M × 24 B) | ~14 GB |
| `h01SeqSet` (8.8M entries × ~64 B/bucket) | ~0.6 GB |
| Vector doubling overhead (30 independent vectors) | up to 2× records |
| **Max expected total** | **~64 GB** |

Even at 100% KEEP rate with worst-case vector doubling, this is borderline.
But realistically the H2 KEEP rate should be well below 100% — H1 was 89% and
H2 (2 mismatches) should be lower. At a more realistic 50% KEEP rate:

| Component | Size |
|---|---|
| Genome | ~35 GB |
| H2 KEEP records (295M × 24 B) | ~7 GB |
| Vector overhead | ~7 GB |
| Other (H0/H1 records, h01SeqSet) | ~1 GB |
| **Expected total** | **~50 GB** |

**The 67 GB RSS implies ~32 GB above the genome baseline, which cannot be
explained by 24-byte records alone at any plausible KEEP rate.** Something
else is consuming memory.

## What Was Ruled Out

1. **`soloRead->record()`** — Correctly gated by `hashCacheSynthProbe_` flag
   (`ReadAlign_outputAlignments.cpp:998`). Solo pipeline does not accumulate
   data during synthetic validation.

2. **`storeQnameMapping` / `g_iReadToQname`** — Gated by env var
   `STAR_INLINE_TRACE_QNAME` which defaults off. Returns immediately on line
   151 of `flex/SoloReadFeature_record_flex.cpp`. No map growth.

3. **`pendingAmbiguous_` CB map** — Synthetic reads have `hammingDist=255`
   (no match), so `match.ambiguous` is false. The ambiguous accumulation
   path at `ReadAlign_outputAlignments.cpp:709` is never entered.

4. **BAM output buffers** — `--outSAMtype None` means `P.outSAMbool`,
   `P.outBAMunsorted`, `P.outBAMcoord` are all false. `writeSAM()` only
   sets `mateMapped` and returns.

5. **Record struct padding** — `FlexHashScreenCache::Record` is exactly 24
   bytes: `uint64_t(8) + uint64_t(8) + uint32_t(4) + uint8_t(1) +
   uint8_t(1) + uint16_t(2) = 24`, naturally aligned, no padding.

## Suspects Still Open

### A. Per-call memory growth in the alignment pipeline

Each of 589M H2 variants calls `flexHashCacheValidateSyntheticPair`, which
runs the full alignment pipeline:

```
mapOneRead() → peOverlapMergeMap() → multMapSelect() → mappedFilter()
→ transformGenome() → chimericDetection() → waspMap() → outputAlignments()
```

Key file: `core/legacy/source/ReadAlign_hashCacheSynth.cpp`

**ReadAlign** objects are allocated once per thread (30 total) in
`ReadAlignChunk` constructors before the loop. Their internal buffers
(suffix arrays, seed stores, transcript arrays) should be reused across
calls. However, over ~20M calls per thread, even a tiny per-call leak
would compound:

- 1 byte/call × 20M calls × 30 threads = **600 MB**
- 16 bytes/call × 20M calls × 30 threads = **9.6 GB** ← this alone would explain the gap

Specific areas to check for growth:

- **`trMult` / `trAll`** — Transcript arrays. Are they bounded, or do they
  grow with `nTr`?
- **`readAnnot.annotFeatures`** — `readAnnot.reset()` is called at the top
  of each synthetic validation. Verify it truly deallocates, not just clears
  size while retaining capacity.
- **`mapOneRead()` seed windows** — Internal `winBin`, `nWA`, `nWAmax`
  arrays. Check if `nWAmax` ratchets up and never shrinks.
- **`chimericDetection()`** — Chimeric transcript arrays may grow.
- **`cigarString` / `chrName`** — `generateCigarP()` creates `std::string`
  stored in `Transcript::cigarString`. Should be overwritten each call, but
  verify the old string's heap allocation is freed (not leaked by a
  pointer-to-Transcript indirection).

### B. `outputAlignments()` side effects

Even though Solo recording is skipped, `outputAlignments()` still runs:

- CIGAR string generation (line 984)
- chrName assignment (line 986)
- `storeQnameMapping` call (no-op, but called)
- `writeSAM()` call (mostly no-op with `--outSAMtype None`)
- `outReadsUnmapped()` for `unmapType >= 0`

`outputAlignments()` is called AFTER we check `unmapType` and `nTr` in
`flexHashCacheValidateSyntheticPair`, but for KEEP/DENY reads (which did
map), it runs the full output path. With `--outSAMtype None`, most
branches are dead, but there may be allocations in the setup code that
persist.

### C. Genome-internal growth

`transformGenome()` at line 99 of `ReadAlign_hashCacheSynth.cpp` transforms
coordinates. Check whether this path appends to any genome-level shared
structure (e.g., splice junction arrays via `recordSJ`).

## Current Workaround

The streaming fix (already committed) avoids in-memory H2 record
accumulation by writing per-thread binary temp files
(`h2_tmp_t<tid>.bin`). This eliminates the vector-doubling overhead
and caps record memory at near-zero. However, **if the leak is in the
alignment pipeline itself, the streaming fix does not prevent it** — the
process will still grow over 589M alignment calls.

## Follow-up: Dedup Container Rewrite

The next memory-shape fix is now implemented in
[`FlexHashCacheGenerate.cpp`](/mnt/pikachu/STAR-suite/core/legacy/source/FlexHashCacheGenerate.cpp):

- H0/H1 records are merged immediately into a packed `khash` dedup table
  keyed by `(seqLo, seqHi, sampleIdx)`.
- H2 temp records are streamed back into the same dedup table and discarded
  immediately after merge.
- The final `vector<Record>` is only materialized once, just before
  `writeHashCacheFile()`.

This removes the previous

- retained `threadRecords` payload after H0/H1 generation,
- heavyweight `std::unordered_map<DedupKey, DedupBucket>` node overhead,
- and the triple-lifetime pattern of `records + dedup map + final vector`.

### 100-probe sanity check after khash rewrite

Replayed the reduced `--hashCacheParentLimit 100` sample after a clean rebuild.

`H0,H1,H2`:
- old: `937,819` final records, `34,735,780 kB` peak RSS
- patched: `937,819` final records, `34,707,140 kB` peak RSS

`H0,H1` control:
- old reference sample: `16,591` final records, `34,618,664 kB` peak RSS
- patched rerun A: `16,590` final records, `34,618,580 kB` peak RSS
- patched rerun B: `16,591` final records, `34,619,104 kB` peak RSS

Interpretation:
- the khash rewrite preserves the full H0/H1/H2 100-probe count exactly
- the small H0/H1 control is not bit-stable run-to-run even on the patched
  binary, so the 1-record delta is current synthetic-validation
  nondeterminism, not a clean deterministic regression from the container rewrite
- at 100 probes the RSS effect is small because the genome dominates; the
  rewrite is aimed at removing the large dedup-container overhead at full H2 scale

## Follow-up: Synthetic R1 Fix + ASAN Check

The synthetic R1 generator had a real off-by-one bug:

- [`buildR1FromParams()`](/mnt/pikachu/STAR-suite/core/legacy/source/FlexHashCacheGenerate.cpp)
  was writing CB/UMI at `cbS` / `umiS`
- while extraction reads at `cbS-1` / `umiS-1` in
  [`SoloReadBarcode_getCBandUMI.cpp`](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadBarcode_getCBandUMI.cpp)

That is now fixed. After the patch, the synthetic validation logs moved from:

- `whitelistIdx=0, hammingDist=255, extractedCbIdxPlus1_=0`

to:

- `whitelistIdx=1, hammingDist=0, extractedCbIdxPlus1_=1`

so the synthetic generator is now actually exercising the intended CB surface.

### Reduced ASAN run

Ran a clean sanitizer build (`ASAN=1`) on a reduced sample:

- `H1,H2`
- `--hashCacheParentLimit 20`
- `--runThreadN 1`

Artifacts:
- `/tmp/h2_asan_20260320/`

Observed result:
- no AddressSanitizer invalid-access report in the synthetic validation path
- LeakSanitizer reported only **44,538 bytes in 728 allocations**
- dominant direct leaks are tiny allocations rooted in `Parameters::Parameters()`
  during `Genome::genomeLoad()`

Interpretation:
- ASAN does **not** show a large leak in the H1/H2 synthetic validation path on
  this reduced fixture
- the reported leaks are far too small to explain the full-run OOM
- the remaining likely causes are still:
  - very large steady-state container/accounting overhead at full H2 scale
  - or per-call retained capacity / ratcheting buffers that do not appear as
    classic unreachable leaks at process exit

## Files to Investigate

| File | Relevance |
|---|---|
| `core/legacy/source/ReadAlign_hashCacheSynth.cpp` | Entry point for synthetic validation; full pipeline invocation |
| `core/legacy/source/FlexHashCacheGenerate.cpp` | H2 generation loop, streaming I/O |
| `core/legacy/source/ReadAlign_mapOneRead.cpp` | Seed/window arrays, `nWAmax` growth |
| `core/legacy/source/ReadAlign.h` | Buffer declarations, array sizes |
| `core/legacy/source/ReadAlignChunk.cpp` | Per-thread ReadAlign allocation |
| `core/legacy/source/ReadAlign_chimericDetection.cpp` | Chimeric arrays |
| `core/legacy/source/ReadAlign_outputAlignments.cpp` | Output path side effects |
| `core/legacy/source/Transcript.h` | `cigarString`, `chrName` members |

## Reproduction

```bash
# Build
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR

# Run (will OOM at ~67 GB with 30 threads on 128 GB machine)
core/legacy/source/STAR \
  --runMode hashCacheGenerate \
  --runThreadN 30 \
  --genomeDir /storage/flex_filtered_reference_2024/star_index \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 \
  --soloBarcodeReadLength 0 \
  --soloCBwhitelist /storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt \
  --flex yes \
  --soloProbeList /storage/flex_filtered_reference_2024/star_index/probe_gene_list.txt \
  --soloFlexExpectedCellsPerTag 3000 \
  --soloSampleWhitelist /mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv \
  --soloSampleProbes /mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt \
  --soloSampleProbeOffset 68 \
  --soloFlexAllowedTags /mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloCellFilter None \
  --soloFeatures Gene \
  --alignEndsType Local \
  --soloStrand Unstranded \
  --soloKeysCompat cr \
  --outSAMtype None \
  --hashCacheTiers "H0,H1,H2" \
  --hashCacheOutput /tmp/test_h012.bin \
  --outFileNamePrefix /tmp/test_h012/
```

Monitor with: `while sleep 10; do ps -o rss= -p $PID | awk '{printf "%.1f GB\n", $1/1048576}'; done`

## Debugging Approach

1. **Instrument per-thread high-water-mark**: Add a periodic (every 100K
   calls) RSS check inside the H2 loop to measure growth rate per alignment
   call. Use `/proc/self/statm` or `getrusage(RUSAGE_SELF)`.

2. **Isolate the pipeline**: Run H2 generation but replace
   `flexHashCacheValidateSyntheticPair` with a no-op (return -1
   immediately). If RSS stays flat, the leak is in the alignment path. If
   it still grows, it's in the loop scaffolding.

3. **Binary search the pipeline**: Progressively uncomment pipeline stages
   in `flexHashCacheValidateSyntheticPair` (mapOneRead only, then add
   multMapSelect, then mappedFilter, etc.) to find which stage introduces
   growth.

4. **Heap profiler**: Run under `heaptrack` or `MALLOC_TRACE` with a
   reduced probe set (`--hashCacheParentLimit 100`) to get an allocation
   backtrace.

## Related Documents

- `docs/HASH_CACHE_TIERED_DESIGN.md` — Tiered cache architecture
- `docs/RUNBOOK_HASH_CACHE_PRUNING_AND_H2_20260319.md` — H1/H2 generation runbook
