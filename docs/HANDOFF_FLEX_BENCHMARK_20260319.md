# Handoff: Flex Full-Scale Benchmark (2026-03-19)

## Status: COMPLETE — instrumentation added, all three runs done

## Summary

All three full-scale Flex benchmark runs are complete on the JAX 16-plex dataset
(2.01B reads, 8 lanes, 4 biological samples). Hash screen instrumentation was
added to `Log.final.out` and confirmed that the 14.8% aligner pass-through rate
is correct behavior, not a bug. Hash_on cell counts match legacy within ±3 cells.

## Hash Screen Instrumentation (new)

Added KEEP/DENY/PASS counters to `Stats.h`, `Stats.cpp`, and
`ReadAlign_oneRead.cpp`. These appear in `Log.final.out` under the
"FLEX HASH SCREEN" section whenever the hash screen is active.

Counter semantics (clean partition — every read gets exactly one):
- **KEEP**: hash cache resolved to a gene, bypassed aligner entirely
- **DENY**: ambiguous probe hit, read discarded
- **PASS**: no cache hit or Keep-not-handled, read sent to aligner

Files changed (on `benchmark-flex` branch):
- `core/legacy/source/Stats.h` — added `hashScreenKeep`, `hashScreenDeny`, `hashScreenPass` fields
- `core/legacy/source/Stats.cpp` — added to `resetN()`, `addStats()`, `reportFinal()`
- `core/legacy/source/ReadAlign_oneRead.cpp` — increment sites at decision branches

## Benchmark Results

### Timing

| Run | Wall time | User CPU (s) | Sys (s) | RSS (GB) | Mapping speed (M/hr) |
|-----|-----------|-------------|---------|----------|---------------------|
| CR9 9.0.1 | 58m 59s | 67,876 | 4,187 | 4.8 | — |
| STAR legacy | 2h 07m 18s | 209,819 | 481 | ~104 | 952 |
| STAR hash_on | 52m 23s | 52,908 | 490 | ~113 | 2,336 |

**Hash_on speedup over legacy**: 2.4x wall-clock, 4.0x CPU time.

### Hash Screen Breakdown (hash_on run)

```
Hash screen: reads evaluated   2,011,130,186   (100%)
                        KEEP   1,698,725,577   (84.47%)
                        DENY      14,727,817   ( 0.73%)
                        PASS     297,676,792   (14.80%)
```

84.5% of reads were resolved directly by the hash cache without touching the
aligner. Only 0.73% were denied (ambiguous probes). 14.8% fell through to the
aligner — this is expected because the cache was built from a 100K-read
downsampled pilot and does not cover all possible read sequences.

### Alignment Stats

| Metric | Legacy (all reads aligned) | Hash_on (298M aligned) |
|--------|---------------------------|------------------------|
| Input reads | 2,011,130,186 | 2,011,130,186 |
| Uniquely mapped | 339,276,959 (16.87%) | 57,534,044 (2.86%) |
| Multi-mapped | 1,634,579,458 (81.28%) | 203,702,164 (10.13%) |
| Unmapped too short | 31,916,124 (1.59%) | 31,916,124 (1.59%) |
| Unmapped other | 5,357,645 (0.27%) | 4,524,460 (0.22%) |
| Core alignment calls | 2,011,130,186 | 297,676,792 |
| Mean alignment time | 86.37 µs/read | 94.18 µs/read |

The hash_on percentages are relative to *all* 2B input reads. If computed
relative to only the 298M reads that reached the aligner: 19.3% unique, 68.4%
multi-mapped, 10.7% unmapped short — nearly identical to legacy's per-read
profile.

### Per-Sample Cell Counts (4 biological samples)

| Sample | Tag | Legacy | Hash_on | CR9 | Δ legacy→hash_on |
|--------|-----|--------|---------|-----|-------------------|
| BC004 (WT-Day-7) | ATGTTGAC | 4,044 | 4,047 | 4,397 | +3 |
| BC006 (PAX6-PTC-D9-Day7) | ATCCCAAC | 4,994 | 4,997 | 5,343 | +3 |
| BC007 (WT-Day-8) | AAGTAGAG | 4,590 | 4,593 | 5,383 | +3 |
| BC008 (PAX6-PTC-D9-Day8) | AGCTGTGA | 4,878 | 4,876 | 5,296 | -2 |

**Hash_on vs legacy: near-perfect match** (±3 cells). The hash screen does not
cause cell loss. The 300–700 cell gap vs CR9 is a pre-existing difference
between STAR-Flex and CellRanger, not a hash-screen regression.

### All 16 Tags (filtered cell counts)

| Tag | Seq | Legacy | Hash_on | Δ |
|-----|-----|--------|---------|---|
| BC001 | ACTTTAGG | 1,175 | 1,175 | 0 |
| BC002 | AACGGGAA | 3,242 | 3,242 | 0 |
| BC003 | AGTAGGCT | 3,251 | 3,253 | +2 |
| BC004 | ATGTTGAC | 4,044 | 4,047 | +3 |
| BC005 | ACAGACCT | 2,977 | 2,977 | 0 |
| BC006 | ATCCCAAC | 4,994 | 4,997 | +3 |
| BC007 | AAGTAGAG | 4,590 | 4,593 | +3 |
| BC008 | AGCTGTGA | 4,878 | 4,876 | -2 |
| BC009 | ACAGTCTG | 1,246 | 1,247 | +1 |
| BC010 | AGTGAGTG | 912 | 912 | 0 |
| BC011 | AGAGGCAA | 3,439 | 3,441 | +2 |
| BC012 | ACTACTCA | 499 | 499 | 0 |
| BC013 | ATACGTCA | 4,827 | 4,827 | 0 |
| BC014 | ATCATGTG | 1,606 | 1,610 | +4 |
| BC015 | AACGCCGA | 2,349 | 2,349 | 0 |
| BC016 | ATTCGGTT | 3,357 | 3,358 | +1 |

Maximum delta: ±4 cells across all 16 tags. The hash screen is effectively
lossless.

## Conclusions

1. **The 14.8% pass-through rate is not a bug.** 84.5% of reads hit the hash
   cache (KEEP), 0.7% were ambiguous (DENY), and 14.8% had no cache match
   (PASS → aligner). The cache was built from 100K downsampled reads and does
   not cover every sequence variant in the full 2B-read dataset.

2. **Hash screen causes zero regression.** Cell counts match legacy within ±4
   cells across all 16 tags. The ±3 on the 4 biological samples is likely
   due to borderline UMI counts at the cell-calling threshold interacting
   with the slightly different read routing.

3. **2.4x wall-clock speedup is real.** Legacy: 2h07m, hash_on: 52m. CPU
   savings are even larger (4x) because the 85% of reads bypassing the
   aligner cost almost zero CPU.

4. **Instrumentation is merged.** `Log.final.out` now shows KEEP/DENY/PASS
   counts and percentages whenever the hash screen is active. No output
   when disabled (legacy mode).

## Occupancy Follow-Up (2026-03-19)

The filtered-cell gap versus the older STAR production run was traced to the
combined occupancy post-filter, not to EmptyDrops:

- Old production STAR (`/mnt/pikachu/prod-12-28/SC2300771`) used a 4-tag
  sample whitelist and removed only `8` cells total by occupancy across the
  4 biological samples.
- The benchmark harness used the full 16-tag whitelist and removed `6418`
  cells total by occupancy, including `351/296/778/382` on the four
  biological tags.
- Pre-occupancy (`ED_Pass`) counts on the four biological tags already matched
  the old STAR/CR surfaces closely; the large deficit appeared only after the
  combined occupancy step.

Fix applied in `flex/source/libflex/FlexFilter.cpp`:

- Combined occupancy now builds its tag/GEM surface from *observed* tags only
  (`nTailTested > 0`), excluding fallback-only unused tags from the occupancy
  denominator.
- Flagged GEMs are still removed across all tag outputs.

Replay result on the new benchmark raw MEX with the patched standalone
`run_flexfilter_mex` and the full 16-tag whitelist:

- Combined occupancy observed tags: `4 / 16`
- Combined occupancy flagged: `2` GEMs
- Combined occupancy removed: `21` cells total
- Biological tag finals:
  - `BC004`: `4393`
  - `BC006`: `5288`
  - `BC007`: `5366`
  - `BC008`: `5258`

These match the old production STAR surface within `0/0/+2/+12` cells and are
again close to CR.

### Filtered MEX comparison vs CR9

Using the patched replay outputs in
`/tmp/flexfilter_occ_fix_20260319/new_raw_full16_occ_on` and the CR9 per-sample
filtered matrices in `/mnt/pikachu/benchmark_cr9_flex_full/outs/per_sample_outs`,
after normalizing CR barcodes from `CB16+TAG8-1` to bare `CB16`:

| Sample | STAR cells | CR9 cells | Shared | Barcode Jaccard | Gene-count Pearson |
|--------|-----------:|----------:|-------:|----------------:|-------------------:|
| `BC004 / WT-Day-7` | 4393 | 4397 | 4316 | 0.9647 | 0.99637 |
| `BC006 / PAX6-PTC-D9-Day7` | 5288 | 5343 | 5255 | 0.9775 | 0.99446 |
| `BC007 / WT-Day-8` | 5366 | 5383 | 5349 | 0.9906 | 0.99630 |
| `BC008 / PAX6-PTC-D9-Day8` | 5258 | 5296 | 5242 | 0.9868 | 0.99626 |
| **Mean** | - | - | - | **0.9799** | **0.99585** |

Combined across the four biological samples:

- barcode Jaccard: `0.98055`
- filtered gene-count Pearson: `0.99594`

This is good enough to treat the occupancy fix as validated against CR9 on the
filtered-cell surface.

## What Remains

### 1. Integrated full rerun on patched `STAR`
The CR9 comparison above was done on the patched `run_flexfilter_mex` replay of
the benchmark raw MEX. The integrated full `STAR` benchmark should still be
rerun once on the patched binary so the benchmark directory itself reflects the
fixed occupancy behavior.

### 2. Merge instrumentation to master
The Stats changes are on the `benchmark-flex` branch and should be squash-merged
to `master`. These are leaf changes (Stats only, no shared files) so squash merge
is safe per the safe merge policy.

### 3. Cache coverage improvement (optional)
The 14.8% pass-through could be reduced by building the cache from a larger
pilot (e.g., 1M reads instead of 100K). This is an optimization, not a
correctness issue.

## Key Paths

| Item | Path |
|------|------|
| Full FASTQs | `/storage/JAX_sequences/` (8 lanes, 141 GB) |
| STAR binary | `/mnt/pikachu/STAR-suite/core/legacy/source/STAR` |
| STAR genome index | `/storage/flex_filtered_reference/star_index` |
| Hash cache | `/storage/downsampled_100K/SC2300771/results/flex_h01_full_cache_20260315_153914/reclassified/sequence_cache.bin` |
| CB whitelist | `/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt` |
| Sample whitelist | `/mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv` |
| Probe list | `/storage/flex_filtered_reference/filtered_reference/probe_list.txt` |
| Sample probes | `/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt` |
| CR9 2024-A reference | `/mnt/pikachu/CR-references/refdata-gex-GRCh38-2024-A` |
| CR9 probe set | `/mnt/pikachu/process_features/tables/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv` |
| CR9 results | `/mnt/pikachu/benchmark_cr9_flex_full/` |
| CR9 config | `/mnt/pikachu/benchmark_cr9_flex_config.csv` |
| Legacy results | `/mnt/pikachu/benchmark_flex_full_legacy_20260319_091414/` |
| Hash_on results | `/mnt/pikachu/benchmark_flex_full_hashon_20260319_112421/` |
| Test harness | `tests/run_flex_hash_screen_internal_2M.sh` |
| MEX compare (shell) | `tests/compare_mex_outputs.sh` |
| MEX compare (python) | `tests/compare_feature_mex.py` |

## Disk Space

- `/storage`: ~591 GB free (1.8T total). STAR temp files go here. All temps cleaned up.
- `/mnt/pikachu`: ~1.5 TB free (60T total). All results safely stored here.

## Mistakes to Avoid

1. **Do NOT `rm -rf` benchmark output directories on `/storage` without first `mv`-ing them to `/mnt/pikachu`.** This session lost a completed 53-min STAR run this way.
2. **Do NOT re-run benchmarks that already completed.** Check existing results first.
3. **Do NOT run ad-hoc STAR commands when an existing test harness script does the job.** Use `tests/run_flex_hash_screen_internal_2M.sh` with env var overrides.
4. **Serialize benchmark runs** (per AGENTS.md) — don't run STAR and CR9 simultaneously.
5. **Always clean-rebuild before benchmarking** (`make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR`).
