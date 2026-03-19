# Handoff: Flex Mapping Phase Optimization (2026-03-19)

## Status: NOT STARTED — analysis and targets identified

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

### Comparison across runs

| Run | Reads aligned | Core align (ms) | µs/read | Mapping wall |
|-----|-------------:|----------------:|--------:|-------------:|
| Legacy (no hash screen) | 2,011,130,186 | 173,704,688 | 86.37 | ~1h 51m |
| Hash_on (current) | 297,676,792 | 28,119,597 | 94.46 | ~35m |

The hash screen reduces reads-to-align by 85%, but the PASS reads are slightly
harder to align (94 vs 86 µs/read) because the easy cache hits were removed.
Total alignment thread-time dropped 6.2x, but wall time only dropped 3.2x —
the I/O and serialization overhead (chunk reads, zcat) is now a larger fraction.

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
| `flex/source/FlexHashScreen.cpp` | H0/H1 hash screen implementation |
| `flex/source/SoloReadFeature_record_flex.cpp` | Per-read Solo recording (hash insert) |

## Optimization Targets

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

**Validation**: Exact same reads must be processed. Compare
`Log.final.out` stats, raw-MEX, and filtered-MEX against the baseline.

### M2: Hash screen coverage improvement (MEDIUM — potential ~4m savings)

**Problem**: 14.8% of reads (298M) fall through to the aligner because the
hash cache was built from a 100K-read downsampled pilot. Increasing cache
coverage reduces the most expensive per-read cost (94 µs alignment).

**Investigation areas**:

1. **Larger pilot**: Build cache from 1M or 10M reads instead of 100K. The
   current cache has 8M records (185 MB). A 10x larger pilot might cover
   95%+ of reads.

2. **Multi-round cache building**: Run a first pass with the current cache,
   collect PASS-read sequences, and build an augmented cache.

3. **Cache hit rate analysis**: Profile which reads miss the cache and why.
   Are they rare sequence variants? Error-containing reads? Novel splice
   junctions?

**Impact estimate**: If coverage goes from 85% to 95%, alignment calls drop
from 298M to ~100M, saving ~$\frac{198M \times 94\mu s}{32} \approx 10$
minutes of wall time. Even 90% coverage would save ~5 minutes.

**Validation**: Raw-MEX parity is the primary surface. Hash screen
KEEP/DENY/PASS counts will change. Alignment stats will change (fewer
reads aligned). Cell counts should remain near-identical.

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
   reads), this may generate excessive candidate loci. Profiling whether
   reducing `winAnchorMultimapNmax` (currently 200) or
   `outFilterMultimapScoreRange` (currently 4) improves speed without
   changing results.

3. **Shorter read optimization**: Flex R2 reads are ~50 bp (average mapped
   length 47.4). STAR's seed/extend is designed for 100–150 bp RNA-seq.
   There may be unnecessary work in the extension phase for these short
   reads.

4. **SIMD / vectorization**: Check compiler flags and whether hot inner
   loops in the aligner are auto-vectorized. AVX-512 is available on
   i9-13900KF.

**Validation**: Exact raw-MEX parity is required. Any parameter change must
be proven to produce identical output before benchmarking speed.

### M4: Genome load optimization (LOW — 12s)

**Problem**: 44 GB genome index is loaded from disk in 12 seconds. This is
fast on NVMe but could be near-zero with shared memory.

**Options**:

1. **`--genomeLoad LoadAndKeep`**: Loads genome into shared memory; subsequent
   runs use the shared copy. Saves 12s per run after the first.
   Not useful for single-run benchmarks, but valuable for production
   pipelines running multiple samples against the same reference.

2. **`--genomeLoad LoadAndRemove`**: Same as above but cleans up after the
   last run.

3. **Genome trimming**: The current genome index is 44 GB for a filtered
   Flex reference. Verify the reference only contains chromosomes with
   probe targets; removing unused contigs shrinks the index.

**Validation**: Transparent — identical output for any `genomeLoad` mode.

### M5: Thread utilization and load balancing (LOW-MEDIUM)

**Problem**: 32 threads with dynamic chunk-based scheduling. Thread idle
time during the last chunk and mutex contention may reduce parallel
efficiency.

**Investigation areas**:

1. **Thread utilization profile**: Run with `perf stat` or add per-thread
   timing to measure actual alignment time vs mutex wait time per thread.

2. **Dynamic thread interface**: `--dynamicThreadInterface 1` exists for
   CR-compat mode. Check if it helps for Flex.

3. **Chunk size vs thread count tradeoff**: With 32 threads and 50MB chunks,
   there are ~2800 chunks (141 GB ÷ 50 MB). That's ~88 chunks per thread,
   which should be sufficient for load balancing.

## Estimated Impact Summary

| Target | Est. savings | Difficulty | Risk |
|--------|-------------|-----------|------|
| M1: I/O pipeline | 3–8m | Medium | Low (no alignment changes) |
| M2: Hash screen coverage | 3–10m | Low-Medium | Low (more KEEP = fewer alignments) |
| M3: Alignment speedup | 1–3m | High | Medium (parameter tuning risk) |
| M4: Genome load | 12s | Trivial | None |
| M5: Thread utilization | 1–3m | Medium | Low |

**Conservative target**: M1 + M2 could bring total wall time from 39m to
~25–30m without touching alignment code.

**Aggressive target**: All five could approach 20m (2.6x total speedup from
legacy 52m, 6.3x from original 2h07m).

## Suggested Work Plan

### Phase 1: Instrument and profile

1. Add per-thread mapping timing (wall time in `mapChunk`, mutex wait time)
2. Run `perf record` during a 2M-read harness to identify alignment hotspots
3. Measure actual zcat decompression throughput vs alignment throughput
4. Profile chunk serialization overhead (time spent waiting on `mutexInRead`)

### Phase 2: Quick wins

1. Try pre-decompressed input (uncompressed FASTQs on NVMe) to measure I/O
   ceiling — this tells you the maximum possible I/O improvement
2. Try `--genomeLoad LoadAndKeep` to eliminate genome load
3. Build a larger hash cache (1M or 10M read pilot) and measure coverage

### Phase 3: Deeper optimization

1. Based on profiling, implement the most impactful change from M1–M5
2. Validate on 100K harness, then 2M harness, then full benchmark
3. Document timing deltas

## Hard Constraints

1. **Always clean-rebuild before timing** — `make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR`
2. **Serialize benchmark runs** — no parallel STAR jobs on this machine
3. **Raw-MEX parity is the primary correctness surface** — any change must produce
   identical (or near-identical with documented explanation) MEX output
4. **Do not change alignment semantics** without explicit approval
5. **Hash screen cache** is a fixed input file — do not modify it during mapping
   optimization. Cache improvements are a separate track.

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
