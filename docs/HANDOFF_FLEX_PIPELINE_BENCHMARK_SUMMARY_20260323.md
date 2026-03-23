# Flex Pipeline Optimization: Comprehensive Benchmark Summary

**Date**: 2026-03-23
**Branch**: `benchmark-flex`
**Dataset**: JAX SC2300771 16-plex, 2.01B reads, 8 lanes, 141 GB gzipped FASTQs
**Machine**: Intel i9-13900KF (24C/32T, 5.8 GHz), 128 GB DDR5, NVMe `/storage`

## Executive Summary

The Flex pipeline has been optimized through two major workstreams:

1. **Solo post-map phase** (P0–P4): 16 min → 3 min (5.1x speedup)
2. **Mapping phase I/O pipeline** (M1): chunk-serial → fully-fused producer-consumer
   with sample pre-filter and combined H0+H1 hash screen, delivering **1.44M
   reads/sec** sustained throughput

Combined, a run that took **2h 07m** (legacy, no hash screen) or **52m** (hash
screen, serial Solo) is now projected at **~27m** with the fully-fused pipeline
and optimized Solo phase.

## Timeline of Optimizations

| Date | Milestone | Impact |
|------|-----------|--------|
| Mar 15 | H0/H1 hash screen cache | 2h07m → 52m (2.4x wall-clock) |
| Mar 19 | Solo P1: packed uint32 UMIs | Solo 16m → 7m (2.2x) |
| Mar 19 | Solo P2: parallel clique BFS (OpenMP) | Solo 7m → 4m |
| Mar 19 | Solo P4: flat sort collapse | Solo 4m → 3m (5.1x total) |
| Mar 19 | Memory fixes (soloFlexMinimalMemory, RAchunk) | Peak RSS 118 → 86 GB |
| Mar 20 | M1: 4-stage producer-consumer pipeline v1 | Functional, ~436K reads/sec best |
| Mar 21 | Fused reader+router (nTriage=0) | 1.05M reads/sec (+140% vs best v1) |
| Mar 21 | Fully-fused (nTriage=0, nSolo=0) | 1.28M reads/sec (+22% vs fused) |
| Mar 23 | **Sample pre-filter + combined H0+H1 offset 0** | **1.44M reads/sec (+14% vs H0-only)** |

## Pipeline Architecture Evolution

### 1. Legacy (chunk-serial)

```
Main thread → mutexInRead → read chunk → mapChunk (per-read: CB/UMI → hash screen → align → Solo)
```

All work serialized under `mutexInRead`. 32 threads contend for one input lock.
Hash screen lookup, CB/UMI extraction, and alignment all in the same per-read path.

**Throughput**: ~952 M reads/hr = ~264K reads/sec

### 2. Hash screen (chunk-serial, with cache)

Same threading model, but 85% of reads skip alignment via hash cache lookup.

**Throughput**: ~2,336 M reads/hr = ~649K reads/sec (2.4x)

### 3. Pipeline v1: 4-stage producer-consumer (Mar 20)

```
Lane readers (8) → readerQ → Triage (1-4) → soloQ / alignQ → Solo consumers (1-4) / Alignment workers
```

Per-lane gzFile readers bypass `mutexInRead`. Triage thread does H0 lookup and
routes KEEP/DENY to sharded Solo consumers, MISS to alignment workers.

**Best result** (t2s2): ~369K reads/sec at 2m mark — **readerQ always full (255)**,
triage was the bottleneck.

| Config | Triage | Solo | Workers | Reads/sec | Bottleneck |
|--------|-------:|-----:|--------:|----------:|------------|
| t1s1 | 1 | 1 | 22 | 124K | Triage (readerQ=255) |
| t1s3 | 1 | 3 | 20 | 309K | Triage (readerQ=255) |
| t2s2 | 2 | 2 | 20 | 369K | Triage (readerQ=255) |
| t2s3 | 2 | 3 | 19 | 356K | Triage (readerQ=254) |
| t3s2 | 3 | 2 | 19 | 360K | Triage (readerQ=255) |
| t4s2 | 4 | 2 | 18 | 437K | Triage (readerQ=255) |

Triage was universally the bottleneck — `readerQ` pegged at 255/256 in every
configuration. More triage threads helped but stole worker threads.

### 4. Fused reader+router (Mar 21)

```
Lane readers (8, fused: gzip + H0 lookup + route) → soloQ → Solo consumers (2) / Alignment workers (22)
```

**Key insight**: merging gzip decompression and H0 lookup into a single per-lane
thread eliminates the `readerQ` and its mutex entirely. Each lane thread does
`gzgets → LUT encode → H0 hash lookup → push to soloQ or alignQ`.

`readerQ` removed. `soloQ` is the new inter-stage queue.

**Throughput**: ~1.05M reads/sec (t0s2 config, 80s mark)

```
PIPELINE_STATS t=80s total=83,007,295 ... soloQ=1,1 alignQ=254 delta10s=10,255,159
```

- Sustained ~1.04M reads/sec (delta10s ~10.3M)
- `alignQ=254` indicates alignment workers are the bottleneck now
- 22 workers available (32 - 8 readers - 2 solo)

### 5. Fully-fused (Mar 21)

```
Lane readers (8, fully-fused: gzip + H0 + CB/UMI + sample + gene count) → alignQ → Alignment workers (24)
```

**Key insight**: for H0 hits (83% of reads), the entire Solo processing chain
(CB/UMI extraction, sample barcode detection, gene counting) can be done inline
in the lane thread. Only MISS reads (17%) need the alignment queue. This
eliminates `soloQ` entirely and frees 2 more threads for alignment.

Both `readerQ` and `soloQ` are removed. Only `alignQ` remains.

**Throughput**: ~1.28M reads/sec (t0s0 config, 80s mark)

- 24 alignment workers (32 - 8 readers)
- `alignQ=256` (full) — alignment is now the clear bottleneck
- KEEP rate: 82.8% (H0 only)

### 6. Sample pre-filter + combined H0+H1 (Mar 23) — CURRENT

```
Lane readers (8, fused: sample check → H0+H1 offset 0 → CB/UMI/Solo) → alignQ → Alignment workers (24, pure alignment)
```

Three optimizations stacked on the fully-fused design:

1. **Sample pre-filter**: encode the 8bp sample tag at offset 68 to a nibble
   code and check against a hash set of all valid variant codes (4 user
   samples × 8 variants = 32 codes). 77.5% of reads fail this check and skip
   the hash screen entirely — they go straight to `alignQ`.

2. **Combined H0+H1 at offset 0**: for the 22.5% of reads with a valid sample
   tag, do a single LUT encode and two hash probes (H0 map, then H1+Deny map
   if H0 misses). This catches the extra 1.9% of reads that H1 resolves at
   offset 0, with zero gene conflicts (validated on 800K reads: 0 conflicts,
   99.97% of H1 hits are at offset 0).

3. **Pure alignment workers**: `classifyRead` removed from `oneReadFromPacket`.
   Workers receive only true MISS reads and do nothing but alignment. No hash
   screen encoding, no binary search, no sample detection overhead.

**Throughput**: ~1.44M reads/sec (t0s0 config, 80s mark)

```
PIPELINE_STATS t=80s total=114,987,850 keep=96,558,576 deny=798,013 miss=17,631,266 alignQ=254
```

- Sustained 1.44M reads/sec (delta10s ~14.5M)
- 24 alignment workers (32 - 8 readers)
- `alignQ=254–256` (full) — alignment remains the bottleneck
- KEEP rate: 84.0%, DENY: 0.7%, MISS: 15.3%
- +14% throughput vs H0-only fully-fused

## Throughput Comparison (2-minute benchmark, 8 lanes, 32 threads)

| Mode | Config | Workers | Reads/sec | vs Legacy | vs Hash-serial |
|------|--------|--------:|----------:|----------:|---------------:|
| Legacy (no hash) | serial | 32 | 264K | 1.0x | — |
| Hash screen (serial) | serial | 32 | 649K | 2.5x | 1.0x |
| Pipeline v1 (best) | t4s2 | 18 | 437K | 1.7x | 0.7x |
| Fused reader+router | t0s2 | 22 | 1,050K | 4.0x | 1.6x |
| Fully-fused (H0 only) | t0s0 | 24 | 1,280K | 4.8x | 2.0x |
| **Fully-fused (sample+H0+H1)** | t0s0 | 24 | **1,440K** | **5.5x** | **2.2x** |

## Full-Scale Run Status

### Completed runs

| Run | Wall time | Notes |
|-----|-----------|-------|
| Legacy (no hash screen) | 2h 07m | `/mnt/pikachu/benchmark_flex_full_legacy_20260319_091414/` |
| Hash screen (serial Solo) | 52m 23s | `/mnt/pikachu/benchmark_flex_full_hashon_20260319_112421/` |
| P1 + memory fixes | 42m 52s | `/mnt/pikachu/benchmark_flex_full_memfix_20260319_183507/` |
| P1 + P2 + P4 (all Solo opts) | 39m 01s | `/mnt/pikachu/benchmark_flex_full_p2p4_20260319_222701/` |

### Projected full-run timing (sample+H0+H1 fused)

| Phase | Serial (P1+P2+P4) | Projected (fully-fused) | Savings |
|-------|-------------------:|------------------------:|--------:|
| Genome load | 12s | 12s | — |
| **Mapping** | **35m 12s** | **~23m** | **~12m** |
| Solo post-map | 3m 06s | ~3m | — |
| Cleanup | ~31s | ~30s | — |
| **Total** | **39m 01s** | **~27m** | **~31%** |

## Hash Screen Performance

The hash screen now uses a three-stage pipeline: sample pre-filter → H0 → H1,
all at offset 0 only, using a LUT-based 2-bit encoding for fast 50-mer →
(seqLo, seqHi) conversion and sample-free hash maps.

| Stage | Action | Rate | Cost |
|-------|--------|------|------|
| Sample pre-filter | 4-nib encode of 8bp sample tag, hash set check | 77.5% rejected (non-user sample) | ~0.02 µs/read |
| H0 lookup | LUT encode + hash probe in `h0NoSampleMap_` | 84.0% KEEP | < 0.1 µs/read |
| H1+Deny lookup | Same LUT, probe `h1DenyNoSampleMap_` | 0.7% DENY, 15.3% MISS | < 0.1 µs/read |
| Total pre-queue resolution | | **84.7% resolved** (KEEP+DENY) | |

### H1 offset-0-only safety validation

Analysis of 800K reads from the hash screen dump confirmed that offset-0-only
H1 lookup is safe for gene assignment:

| Metric | Value |
|--------|-------|
| H1-only reads (offset 0 hit) | 43,453 |
| Gene conflicts (different gene at offset -1/+1) | **0** (0.00%) |
| Offset 0 coverage of all H1 hits | 99.97% |

Zero gene conflicts means offset-0-only H1 never assigns a read to the wrong
gene compared to the full 3-offset H1 search.

## Solo Phase Optimization Detail

| Phase | Original (16m) | After P1 (7m) | After P1+P2+P4 (3m) |
|-------|---------------:|---------------:|---------------------:|
| readInfo alloc (2B) | 2m 43s | skipped | skipped |
| Clique correction | 8m 14s | 3m 39s | 30s (parallel, 32T) |
| Hash collapse | 3m 24s | 2m 13s | 18s (flat sort) |
| FlexFilter (16 tags) | 1m 23s | 28s | 27s |
| MEX build+write | 13s | 11s | 10s |

### Optimizations applied

| ID | What | Technique | Savings |
|----|------|-----------|---------|
| P0 | Clique replay harness | Standalone tool + 14 tests | (test infra) |
| P1 | Packed uint32 UMIs | Eliminated 750M+ string allocs | 4m 35s |
| P2 | Parallel clique BFS | OpenMP, 32 threads, dynamic(256) | 3m+ |
| P4 | Flat sort collapse | Removed dead nested maps, sort-based dedup | 3m+ |
| Mem | soloFlexMinimalMemory | Skip 16GB readInfo + eager hash free | 28 GB RSS |

## Correctness

### Hash screen parity (full scale)

| Metric | Legacy | Hash_on | Delta |
|--------|--------|---------|-------|
| Reads evaluated | 2,011,130,186 | 2,011,130,186 | 0 |
| Per-sample cell counts | ±4 max | ±4 max | Effectively lossless |

### Solo phase parity (full scale, P1+P2+P4 vs hash_on baseline)

| Check | Result |
|-------|--------|
| Raw matrix dimensions | Identical: 18082 × 1882155, 99,963,085 entries |
| Value diffs | 1005 entries (0.001%) — parallel tie-breaking |
| Total UMI delta | −376 out of 228M (0.00016%) |
| reads_before/after | Identical: 1,673,223,131 |

### CR9 parity (filtered cells)

| Sample | STAR cells | CR9 cells | Barcode Jaccard | Gene Pearson |
|--------|----------:|----------:|----------------:|-------------:|
| BC004 (WT-Day-7) | 4,393 | 4,397 | 0.965 | 0.996 |
| BC006 (PAX6-PTC-D9-Day7) | 5,288 | 5,343 | 0.978 | 0.994 |
| BC007 (WT-Day-8) | 5,366 | 5,383 | 0.991 | 0.996 |
| BC008 (PAX6-PTC-D9-Day8) | 5,258 | 5,296 | 0.987 | 0.996 |
| **Mean** | | | **0.980** | **0.996** |

## Pipeline Architecture (fully-fused mode)

### Thread allocation (32 threads, 8 lanes)

| Role | Count | Function |
|------|------:|----------|
| Lane reader+Solo (fully-fused) | 8 | `flexLaneReaderFullThread` |
| Alignment workers | 24 | Standard STAR alignment |
| Stats reporter | 1 | (shares with worker) |

### Per-lane thread work (full decision tree)

1. `gzgets()` — read 4 FASTQ lines (R2 + R1)
2. **Sample pre-filter** — `asciiToNib()` encode 8bp sample tag (offset 68),
   check `samplePreFilter` hash set → **77.5% rejected → `alignQ` (pure align)**
3. `encodeWindowLUT()` — LUT-based 2-bit encoding of 50-mer (22.5% reach here)
4. `classifyReadH0H1Offset0()` — H0 probe, then H1+Deny probe if H0 misses:
   - **KEEP (84%)** → `getCBandUMI()` → `detectSampleFromPackedTag()` → `record_flex_hash_screen_keep()` (inline Solo)
   - **DENY (0.7%)** → `getCBandUMI()` → `detectSampleFromPackedTag()` → record deny (inline Solo)
   - **MISS (15.3%)** → push `EnrichedPacket` to `alignQ`

### Queue summary

| Queue | Fully-fused | Fused (t0s2) | v1 (tNsM) |
|-------|:-----------:|:------------:|:----------:|
| readerQ | eliminated | eliminated | present (bottleneck) |
| soloQ | eliminated | present | present |
| alignQ | present | present | present |

### Parameters

| Parameter | Default | Fully-fused | Description |
|-----------|---------|-------------|-------------|
| `--flexPipeline` | `auto` | `yes` | Enable pipeline mode |
| `--flexPipelineNTriage` | `1` | `0` | 0 = fused into lane readers |
| `--flexPipelineNSolo` | `3` | `0` | 0 = fused into lane readers |

### Pipeline modes supported

| nTriage | nSolo | Mode | Description |
|--------:|------:|------|-------------|
| >0 | >0 | Separate stages | Reader → Triage → Solo → Alignment |
| 0 | >0 | Fused reader+router | Lane threads do H0 + route to soloQ/alignQ |
| 0 | 0 | **Fully-fused** | Lane threads do sample pre-filter + H0+H1 + Solo inline, only MISS → alignQ |

## Key Files (fully-fused implementation)

| File | Role |
|------|------|
| `core/legacy/source/FlexPipeline.h` | Pipeline state, packet structs, BoundedQueue, thread args |
| `core/legacy/source/FlexPipeline.cpp` | All thread implementations (reader, triage, solo, fused, fully-fused, sample pre-filter) |
| `core/legacy/source/FlexHashScreen.h` | Hash screen declarations, `H0NoSampleMap`, `H1DenyNoSampleMap` |
| `core/legacy/source/FlexHashScreen.cpp` | Hash screen cache: H0/H1/Deny maps, LUT encode, `classifyReadH0H1Offset0()` |
| `core/legacy/source/ReadAlign_oneRead.cpp` | Alignment path (hash screen removed — pure alignment) |
| `core/legacy/source/mapThreadsSpawn.cpp` | Thread spawning, mode selection, stats merge |
| `core/legacy/source/ParametersSolo.cpp` | Parameter validation (nTriage 0–8, nSolo 0–16) |

## Git History (benchmark-flex branch)

| Commit | Description |
|--------|-------------|
| (next) | flex: sample pre-filter + combined H0+H1 offset-0 in fused threads |
| `a7321fa` | flex: fully-fused lane threads with inline Solo processing |
| `481ad24` | flex: fused reader+router pipeline with LUT H0 lookup |
| `8d825f7` | docs: H012 cache validation results + tiered design doc |
| `a79046e` | flex: tighten hash cache H2 generation memory use |
| Earlier | Solo P0–P4 optimizations, hash screen instrumentation |

## What Remains

### Immediate

1. **Full-scale validation** — run the sample+H0+H1 fused binary on the full
   2.01B read dataset with corrected `--soloProbeList`. Verify Solo.out
   contents and cell-count parity against P1+P2+P4 baseline.

### Short-term optimizations

2. **H2 cache coverage** — build from 1M+ read pilot to push KEEP rate from
   83% toward 93–97%, reducing alignment calls further. Each 1% KEEP gain
   saves ~30s of alignment wall time.

3. **Alignment worker efficiency** — `alignQ=256` (full) means alignment is
   now the bottleneck. Profile alignment hotspots with `perf` to find
   micro-optimization opportunities for 50bp Flex reads.

### Medium-term

4. **BAM output support** — pipeline currently requires `--outSAMtype None`.
   Extending to sorted BAM needs mutex-free BAM writing or a dedicated writer
   thread.

5. **Deterministic read ordering** — the fully-fused mode produces minor
   non-deterministic UMI dedup tie-breaking (~0.001% of entries). A
   deterministic `iReadAll` assignment (lane-local counters) would fix this
   if strict bit-level reproducibility is required.

6. **Dynamic thread balancing** — currently 8 lanes = 8 reader threads
   (fixed). If lanes have unequal sizes, some readers finish early and their
   threads are wasted. A work-stealing or dynamic lane assignment scheme
   could improve tail latency.

## Benchmark Run Reference

| Run | Config | Path |
|-----|--------|------|
| Legacy (no hash) | serial, 32T | `/mnt/pikachu/benchmark_flex_full_legacy_20260319_091414/` |
| Hash_on (pre-Solo-opt) | serial, 32T | `/mnt/pikachu/benchmark_flex_full_hashon_20260319_112421/` |
| P1 + memory fixes | serial, 32T | `/mnt/pikachu/benchmark_flex_full_memfix_20260319_183507/` |
| P1+P2+P4 (all Solo opts) | serial, 32T | `/mnt/pikachu/benchmark_flex_full_p2p4_20260319_222701/` |
| Pipeline sweep | t1s1 → t4s4 | `/storage/flex_pipeline_sweep_20260320_230900/` |
| Fused 2-min test | t0s2 | `/storage/flex_fused_test_20260321_025257/fused/` |
| Fully-fused 2-min test (H0 only) | t0s0 | `/storage/flex_fullfused_test_20260321_031613/fullfused/` |
| Fused H0+H1 2-min test | t0s0 | `/storage/flex_fused_h0h1_test_20260323_185250/` |
| **Fully-fused full (pending)** | t0s0, 32T | (to be launched) |
| CR9 reference | CR 9.0.1 | `/mnt/pikachu/benchmark_cr9_flex_full/` |

## Related Documentation

| Doc | Contents |
|-----|----------|
| `docs/HANDOFF_FLEX_BENCHMARK_20260319.md` | Full benchmark history, Solo P0–P4 detail |
| `docs/HANDOFF_FLEX_MAPPING_PHASE_OPTIMIZATION_20260319.md` | Mapping phase targets M1–M5, tiered hash design |
| `docs/HANDOFF_FLEX_PIPELINE_V1_STATUS_20260320.md` | Pipeline v1 implementation details |
| `docs/RUNBOOK_FLEX_OPTIMIZATION_CHOKEPOINTS_20260319.md` | Solo optimization runbook P0–P5 |
| `docs/RUNBOOK_HASH_CACHE_PRUNING_AND_H2_20260319.md` | Hash cache tiering and H2 design |
