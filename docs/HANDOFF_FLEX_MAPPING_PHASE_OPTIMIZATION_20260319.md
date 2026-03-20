# Handoff: Flex Mapping Phase Optimization (2026-03-19)

## Status: IN-PROGRESS — Phase 1-3 implemented (dump + replay harness + unit tests), Phase 4+ pending

## Active Work Item

**Tiered hash screen lookup + cache pruning + H2 coverage expansion.**

Detailed runbook: `docs/RUNBOOK_HASH_CACHE_PRUNING_AND_H2_20260319.md`

### Decision summary from planning discussion

1. **Tiered lookup is the primary structural optimization.** The current
   `FlexHashScreen.cpp` stores all 8M cache records in a single sorted
   `std::vector<Record>` (185 MB, 5× larger than L3). Every read does 3-6
   binary searches (23 steps each) through this one array. Splitting into
   H0 (exact, cacheClass=0), H1 (1MM, cacheClass=1), and deny tiers means
   the dominant path — H0 exact match at offset 0, resolving ~70-80% of
   reads — searches a small array that fits in L2/L3 instead of the full 8M.

2. **Dead-weight pruning is free with tiering.** Records where
   cacheClass==2 + negativeCode!=ProbeAmbig + resolvedGeneIdx15==0 are
   invisible at runtime (`classifyHits` finds them, skips them via
   `continue`, and the read gets PASS — identical to a cache miss).
   The tiered loader simply doesn't load them into any tier.

3. **H2 cache only grows H1, not H0.** Probes are fixed sequences so
   H0 (exact matches) stays the same size regardless of pilot size. A
   larger pilot discovers more 1MM error variants, which go into H1.
   With tiered lookup the H0 fast path stays small and hot.

4. **Harness first, not Phase 0 audit.** Skip measuring the cache
   composition in isolation. Instead, build a standalone replay harness
   (following the clique replay pattern) that loads the cache, replays
   dumped read inputs, and compares flat vs tiered decisions. This gives
   both the composition stats and a rapid correctness/iteration tool.

### Next steps for the implementing agent

Follow the runbook phases in order:

1. **Phase 1**: Add env-gated dump (`STAR_DUMP_HASH_SCREEN=<path>`) in
   `ReadAlign_oneRead.cpp` at line 314 (after `classifyRead` call). Dump
   inputs (readSeq, readLen, sampleIdx) + decision (action, geneIdx15,
   cacheClass, offset) per read. One STAR build, one 100K harness run.

2. **Phase 2**: Build `flex/tools/hash_screen_replay/` — standalone tool
   with both flat (reference) and tiered lookup implementations. Shares
   `Record` type and `encodeWindow` with `FlexHashScreen.cpp` but no
   `ParametersSolo` dependency. Template: `flex/tools/clique_replay/`.

3. **Phase 3**: Unit tests (12 cases covering every decision path). Each
   test constructs a small in-memory cache, runs both flat and tiered,
   asserts identical decisions.

4. **Phase 4**: Validate flat replay matches STAR dump (zero mismatches),
   then validate tiered matches flat (zero mismatches). Sub-second runs.

5. **Phase 5**: Iterate on the tiered implementation (H0-offset-0 short-
   circuit, parameter sweeps, H2 cache files) — all in the harness, no
   STAR rebuilds. The `--stats` flag on the replay tool reports tier sizes
   and dead-weight count from the loaded cache.

6. **Phase 6-7**: Port validated tiered code to `FlexHashScreen.cpp` (two
   files: `.h` and `.cpp`), one STAR rebuild, 100K harness validation.

7. **Phase 8**: H2 cache (larger pilot), tested in harness first.

## Context

The Flex Solo post-map phase has been optimized from 16 minutes to 3 minutes
(5.1x speedup). The mapping phase at **~36 minutes is now 92% of the 39-minute
total wall time** and is the clear next optimization target.

Prior work:
- `docs/HANDOFF_FLEX_BENCHMARK_20260319.md` — full benchmark history and Solo optimization results
- `docs/RUNBOOK_FLEX_OPTIMIZATION_CHOKEPOINTS_20260319.md` — detailed Solo phase optimization runbook

## Machine Specs

- **CPU**: Intel i9-13900KF, 24 cores / 32 threads (16P + 8E), max 5.8 GHz
- **RAM**: 128 GB DDR5
- **Storage**: `/storage` on NVMe (1.8 TB), `/mnt/pikachu` on RAID (60 TB HDD)
- **OS**: Linux 6.8.0

## Current Pipeline Timing (best run: P1+P2+P4)

**Run**: `/mnt/pikachu/benchmark_flex_full_p2p4_20260319_222701/`

| Phase | Wall time | % of total |
|-------|-----------|------------|
| Genome load | 12s | 0.5% |
| **Mapping (32 threads)** | **35m 12s** | **90%** |
| Solo post-map | 3m 06s | 8% |
| Cleanup | ~31s | 1.3% |
| **Total** | **39m 01s** | |

### Mapping phase internal breakdown

The mapping phase processes 2.01B reads with 32 threads. The hash screen
diverts 85.2% of reads before alignment:

```
Total input reads:       2,011,130,186
  → KEEP (hash screen):  1,698,725,577  (84.47%)  — bypasses aligner
  → DENY (ambiguous):       14,727,817  ( 0.73%)  — bypasses aligner
  → PASS (to aligner):    297,676,792  (14.80%)  — aligned by STAR core
```

#### Per-stage cumulative thread-time

From `Log.final.out` Flex stage timing counters:

| Stage | Calls | Total thread-ms | Mean µs/call | Notes |
|-------|------:|---------:|-------------:|-------|
| Sample calling (pre-align) | 2,011,130,186 | 380,171 | 0.19 | Hash screen + sample barcode decode for ALL reads |
| Sample calling (output path) | 297,676,792 | 70,811 | 0.24 | Post-alignment output routing for PASS reads |
| Core alignment | 297,676,792 | 28,119,597 | 94.46 | STAR aligner for PASS reads |

#### Wall-time estimate for mapping sub-phases

With 32 threads and the cumulative thread-times above:

| Sub-phase | Thread-time (s) | Est. wall-time (32 threads) | Notes |
|-----------|------:|---------:|-------|
| Genome load (disk → RAM) | — | ~12s | 44 GB index (17G Genome + 26G SA + 1.6G SAindex) |
| Hash screen (all 2B reads) | 380s | ~12s | Well parallelized; 0.19 µs/read |
| Core alignment (298M reads) | 28,120s | ~14.5m | ~94 µs/read; 96% of mapping thread-time |
| Output routing (298M reads) | 71s | ~2.2s | Post-alignment per-read processing |
| I/O + chunk serialization | — | ~8m | Mutex waits, zcat decompression, chunk reads |

The **~8 minutes** of unaccounted wall time between measured thread-time and
actual mapping wall (~35m) represents I/O overhead: reading 141 GB of gzipped
FASTQs, chunk serialization under `mutexInRead`, and potential thread idle time.

### Hash screen architecture (current)

The cache is a single flat sorted array:

| Metric | Value |
|--------|-------|
| Records | 8,074,151 |
| Record size | 24 bytes |
| Total size | 185 MB |
| i9-13900KF L3 | 36 MB (cache is 5× larger) |
| Binary search depth | log2(8M) ≈ 23 steps |
| Lookups per read | 3–6 (3 offsets × 1–2 sample fallbacks) |

`FlexHashScreen.cpp::classifyRead` encodes the R2 50-mer at offsets
{0, +1, −1} and does a `std::lower_bound` through the entire 8M-record
vector for each. With the array 5× larger than L3, the later binary-search
steps hit main memory (~100 ns each).

**Key insight**: `classifyHits` only produces terminal decisions from three
record types:
- cacheClass 0/1, resolvedGeneIdx15 > 0 → KEEP
- cacheClass 2, negativeCode == ProbeAmbig → DENY
- cacheClass 0/1, sampleIdx > 0 (mismatch) → **unreachable**: `findRecord`
  only returns exact-sample or sample=0 fallback, so a wrong-sample record
  is never found. The read gets PASS (cache miss) instead.

Records with cacheClass==2, negativeCode!=ProbeAmbig, resolvedGeneIdx15==0
are found by binary search, examined, skipped (`continue`), and the read
gets PASS — identical to a cache miss. These are dead weight.

### Proposed tiered architecture

Split at load time into per-class sorted arrays:

| Tier | Contents | Est. records | Est. size | Fits in |
|------|----------|-------------:|----------:|---------|
| H0 | cacheClass 0, gene > 0 | ~50-300K | 1-7 MB | L3 / L2 |
| H1 | cacheClass 1, gene > 0 | ~1-5M | 24-120 MB | Partial L3 |
| Deny | cacheClass 2, ProbeAmbig | ~50-200K | 1-5 MB | L3 |
| Dropped | dead weight | ? | ? | not loaded |

Fast path: H0 at offset 0 → single binary search (log2(300K) ≈ 18 steps,
all L2/L3) → immediate KEEP for ~70-80% of reads.

### Comparison across runs

| Run | Reads aligned | Core align (ms) | µs/read | Mapping wall |
|-----|-------------:|----------------:|--------:|-------------:|
| Legacy (no hash screen) | 2,011,130,186 | 173,704,688 | 86.37 | ~1h 51m |
| Hash_on (current) | 297,676,792 | 28,119,597 | 94.46 | ~35m |

The hash screen reduces reads-to-align by 85%, but the PASS reads are slightly
harder to align (94 vs 86 µs/read) because the easy cache hits were removed.
Total alignment thread-time dropped 6.2x, but wall time only dropped 3.2x —
the I/O and serialization overhead (chunk reads, zcat) is now a larger fraction.

## All Optimization Targets (M1–M5)

### M1: I/O pipeline optimization (HIGH — ~8m of wall time)

**Problem**: ~8 minutes of the 35-minute mapping wall time is I/O overhead:
reading 141 GB of gzipped FASTQs through `zcat`, with `mutexInRead`
serializing chunk reads.

**Investigation areas**:

1. **Parallel decompression**: `zcat` spawns one decompression process per
   input file, but STAR reads chunks serially under a mutex. Consider
   `pigz` or `libdeflate` for faster decompression, or pre-decompress to
   tmpfs/NVMe.

2. **Chunk size tuning**: `limitIObufferSize` is set to 50M. Larger chunks
   reduce mutex contention but increase memory. Smaller chunks improve
   load balancing. Profile the actual mutex wait time.

3. **Double buffering**: Pre-read the next chunk while the current one is
   being mapped. This requires a producer-consumer pattern rather than the
   current "acquire mutex → read → release → map" serial pattern.

4. **Pre-decompressed input**: For repeated benchmarks on the same data,
   decompressing to NVMe once (`pigz -dk`) and using uncompressed input
   eliminates decompression overhead entirely.

**Status**: Not started. Lower priority than M2 (hash screen tiering) for
now; revisit if hash screen gains are insufficient.

### M2: Hash screen tiering + coverage (ACTIVE — see runbook)

**Tiered lookup** (structural): split cache into H0/H1/deny tiers, search
H0 first. Reduces per-read binary-search cost for the 85% of reads that
are exact matches.

**Cache pruning** (free): dead-weight records not loaded into any tier.

**H2 coverage** (follow-on): build cache from 1M+ read pilot to push
KEEP rate from 85% toward 93-97%, reducing alignment calls.

**Runbook**: `docs/RUNBOOK_HASH_CACHE_PRUNING_AND_H2_20260319.md`

### M3: Alignment speedup for PASS reads (MEDIUM-LOW)

**Problem**: PASS reads are harder to align (94 µs/read vs 86 µs for legacy
all-reads average). This is inherent — the easy cache-resolvable reads
were removed. But 94 µs × 298M reads × (1/32 threads) = ~14.5 minutes of
wall time.

**Investigation areas**:

1. **Profile alignment hotspots**: Use `perf record` during a short run to
   identify hot functions within `mapOneRead()`. Seed generation vs
   extension vs scoring vs splice detection.

2. **Alignment parameter tuning**: Current settings allow generous
   multi-mapping (`outFilterMultimapNmax 10000`). For Flex probes (50-mer
   reads), this may generate excessive candidate loci.

3. **Shorter read optimization**: Flex R2 reads are ~50 bp (average mapped
   length 47.4). STAR's seed/extend is designed for 100–150 bp RNA-seq.

4. **SIMD / vectorization**: Check compiler flags and whether hot inner
   loops in the aligner are auto-vectorized. AVX-512 is available on
   i9-13900KF.

**Status**: Not started.

### M4: Genome load optimization (LOW — 12s)

`--genomeLoad LoadAndKeep` for shared-memory reuse across runs.
Not useful for single-run benchmarks.

**Status**: Not started.

### M5: Thread utilization and load balancing (LOW-MEDIUM)

Profile per-thread alignment time vs mutex wait time; check
`--dynamicThreadInterface 1` for Flex.

**Status**: Not started.

## Estimated Impact Summary

| Target | Est. savings | Difficulty | Risk |
|--------|-------------|-----------|------|
| M1: I/O pipeline | 3–8m | Medium | Low (no alignment changes) |
| **M2: Hash screen tiering + H2** | **3–10m** | **Low-Medium** | **Low** |
| M3: Alignment speedup | 1–3m | High | Medium (parameter tuning risk) |
| M4: Genome load | 12s | Trivial | None |
| M5: Thread utilization | 1–3m | Medium | Low |

**Conservative target**: M2 + M1 could bring total wall time from 39m to
~25–30m without touching alignment code.

## Architecture of the Mapping Phase

### Thread model

STAR uses a chunk-based threading model:

1. **Main thread** spawns N-1 worker threads via `pthread_create`
2. Each thread runs `ReadAlignChunk::processChunks()` in a loop:
   - Acquire `mutexInRead` → read a chunk of FASTQ data → release
   - `mapPermitAcquire()` (throttle if dynamic threading enabled)
   - `mapChunk()` → loop over reads → `oneRead()` per read
   - `mapPermitRelease()`
3. Reads within each chunk are processed sequentially per thread

### Per-read flow (`ReadAlign_oneRead.cpp`)

For each read:

1. Parse CB + UMI from R1 (sample calling pre-align)
2. Hash screen check → KEEP/DENY/PASS
3. If PASS: `mapOneRead()` → seed + extend alignment
4. `multMapSelect()`, `mappedFilter()`, `outputAlignments()`
5. Record to Solo inline hash

### Key serialization points (mutexes)

| Mutex | Location | What it protects |
|-------|----------|-----------------|
| `mutexInRead` | `ReadAlignChunk_processChunks.cpp` | Chunk reading from input streams |
| `mutexOutSAM` | `ReadAlignChunk_mapChunk.cpp` / `BAMoutput.cpp` | SAM/BAM writes |
| `mutexOutBAM1` | `BAMoutput.cpp` | BAM sort bin updates |
| `mutexStats` | `ReadAlignChunk_mapChunk.cpp` | Stats aggregation |

With `--outSAMtype None`, `mutexOutSAM` and `mutexOutBAM1` are inactive. The
dominant serialization point is `mutexInRead` for chunk I/O.

### Key source files

| File | Role |
|------|------|
| `core/legacy/source/ReadAlignChunk_processChunks.cpp` | Main per-thread loop, chunk I/O |
| `core/legacy/source/ReadAlignChunk_mapChunk.cpp` | Per-chunk read loop, calls `oneRead()` |
| `core/legacy/source/ReadAlign_oneRead.cpp` | Per-read: hash screen, alignment, output |
| `core/legacy/source/ReadAlign_mapOneRead.cpp` | Core seed+extend alignment |
| `core/legacy/source/mapThreadsSpawn.cpp` | Thread creation |
| `core/legacy/source/ThreadControl.h` | Thread entry point |
| `core/legacy/source/Genome_genomeLoad.cpp` | Genome index loading |
| `core/legacy/source/FlexHashScreen.h` | Cache record type, classifyRead API |
| `core/legacy/source/FlexHashScreen.cpp` | Cache loader, flat binary search, classifyHits |
| `flex/source/SoloReadFeature_record_flex.cpp` | Per-read Solo recording (hash insert) |

## Hard Constraints

1. **Always clean-rebuild before timing** — `make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR`
2. **Serialize benchmark runs** — no parallel STAR jobs on this machine
3. **Raw-MEX parity is the primary correctness surface** — any change must produce
   identical (or near-identical with documented explanation) MEX output
4. **Do not change alignment semantics** without explicit approval
5. **Hash screen cache** is a fixed input file — do not modify it during mapping
   optimization. Cache improvements are a separate track (H2).

## Key Paths

| Item | Path |
|------|------|
| STAR binary | `/mnt/pikachu/STAR-suite/core/legacy/source/STAR` |
| Genome index | `/storage/flex_filtered_reference/star_index/` (44 GB total) |
| Hash cache | `/storage/downsampled_100K/SC2300771/results/flex_h01_full_cache_20260315_153914/reclassified/sequence_cache.bin` (185 MB, 8M records) |
| Full FASTQs | `/storage/JAX_sequences/SC2300771_GT23-14630_*` (141 GB, 8 lanes × 2 reads) |
| 100K harness FASTQs | `/storage/downsampled_100K/SC2300771/` |
| Current best run | `/mnt/pikachu/benchmark_flex_full_p2p4_20260319_222701/` |
| Baseline run (hash_on) | `/mnt/pikachu/benchmark_flex_full_hashon_20260319_112421/` |
| Legacy run (no hash) | `/mnt/pikachu/benchmark_flex_full_legacy_20260319_091414/` |
| CR9 reference run | `/mnt/pikachu/benchmark_cr9_flex_full/` |
| CB whitelist | `/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt` |
| Sample whitelist | `/mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv` |
| Probe list | `/storage/flex_filtered_reference/filtered_reference/probe_list.txt` |
| Sample probes | `/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt` |
| **Runbook** | `docs/RUNBOOK_HASH_CACHE_PRUNING_AND_H2_20260319.md` |
| **Clique replay (template)** | `flex/tools/clique_replay/` |
| **New replay tool (planned)** | `flex/tools/hash_screen_replay/` |

## Full Benchmark Command (current best)

```bash
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="/mnt/pikachu/benchmark_flex_mapping_opt_${TIMESTAMP}"
TMPDIR="/storage/tmp_flex_map/STARtmp"
rm -rf "${TMPDIR}"
mkdir -p "${OUTDIR}"

export STAR_FLEX_HASH_SCREEN_CACHE="/storage/downsampled_100K/SC2300771/results/flex_h01_full_cache_20260315_153914/reclassified/sequence_cache.bin"

/mnt/pikachu/STAR-suite/core/legacy/source/STAR \
  --runThreadN 32 \
  --outTmpDir "${TMPDIR}" \
  --genomeDir /storage/flex_filtered_reference/star_index \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
  --soloCBwhitelist /storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt \
  --flex yes \
  --soloFlexExpectedCellsPerTag 3000 \
  --soloSampleWhitelist /mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv \
  --soloProbeList /storage/flex_filtered_reference/filtered_reference/probe_list.txt \
  --soloSampleProbes /mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt \
  --soloSampleProbeOffset 68 \
  --soloFlexAllowedTags /mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv \
  --soloFlexOutputPrefix "${OUTDIR}/per_sample" \
  --limitIObufferSize 50000000 50000000 \
  --outSJtype None \
  --soloMultiMappers Rescue \
  --alignIntronMax 500000 \
  --outFilterMismatchNmax 6 \
  --outFilterMismatchNoverReadLmax 1.0 \
  --outFilterMatchNmin 25 \
  --outSAMunmapped None \
  --outFilterMatchNminOverLread 0 \
  --outFilterMultimapNmax 10000 \
  --outFilterMultimapScoreRange 4 \
  --outSAMmultNmax 10000 \
  --winAnchorMultimapNmax 200 \
  --outSAMprimaryFlag AllBestScore \
  --outFilterScoreMin 0 \
  --outFilterScoreMinOverLread 0 \
  --outSAMattributes None \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloCellFilter None \
  --clipAdapterType CellRanger4 \
  --soloFeatures Gene \
  --alignEndsType Local \
  --soloStrand Unstranded \
  --chimSegmentMin 1000000 \
  --soloKeysCompat cr \
  --outSAMtype None \
  --soloSampleSearchNearby no \
  --readFilesCommand zcat \
  --readFilesIn \
    /storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz,\
/storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R2_001.fastq.gz,\
/storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R2_001.fastq.gz,\
/storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R2_001.fastq.gz,\
/storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R2_001.fastq.gz,\
/storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R2_001.fastq.gz,\
/storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R2_001.fastq.gz,\
/storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R2_001.fastq.gz \
    /storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz,\
/storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R1_001.fastq.gz,\
/storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R1_001.fastq.gz,\
/storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R1_001.fastq.gz,\
/storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R1_001.fastq.gz,\
/storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R1_001.fastq.gz,\
/storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R1_001.fastq.gz,\
/storage/JAX_sequences/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R1_001.fastq.gz \
  --outFileNamePrefix "${OUTDIR}/" \
  > "${OUTDIR}/stdout.log" 2> "${OUTDIR}/stderr.log"
```

## Mistakes to Avoid

1. **Do NOT `rm -rf` benchmark output directories** without first moving them.
2. **Do NOT re-run benchmarks that already completed.** Check existing results first.
3. **Serialize benchmark runs** (per AGENTS.md) — don't run STAR and CR9 simultaneously.
4. **Always clean-rebuild before benchmarking.**
5. **Do NOT change alignment parameters and assume output is identical.** Always validate
   raw-MEX parity after any parameter change.
6. **The hash screen cache is read-only during mapping.** Cache building is a separate
   process that runs on downsampled data before the main benchmark.
7. **Do NOT start implementation without building the replay harness first.** The harness
   is the iteration surface. Follow the runbook phases in order.
8. **FlexHashScreen.cpp is in `core/legacy/source/`**, not `flex/source/`. The `flex/source/`
   reference in earlier docs was incorrect.
