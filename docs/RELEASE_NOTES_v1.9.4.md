# STAR Suite v1.9.4 Release Notes

Date: 2026-09-13

STAR Suite 1.9.4 makes the half-probe method the Flex default and retires the
other Flex probe-assignment routes as legacy, so that every Flex run uses one
method unless the user asks otherwise. It also cuts the memory needed for very
large Flex runs by more than half. It includes the merged changes through
`c56cd5f` and the half-probe default change.

`STAR --version` reports `1.9.4`. Debian source packaging uses `1.9.4-1`;
Ubuntu packages use `1.9.4-1~ubuntu22.04.1` and
`1.9.4-1~ubuntu24.04.1`. Upstream STAR remains `2.7.11b`, genome-index
compatibility remains `2.7.4a`, and legacy compatibility remains `2.7.1a`.
Existing indexes do not need rebuilding; Flex runs now need a half-probe cache.

## Flex: the half-probe method is the default

- **`--flex yes` assigns reads to probes from the half-probe (H1X2) cache and
  aligns nothing.** A read is assigned when it matches a probe exactly, or when
  one of its two 25-base halves identifies a single probe and the full 50 bases
  differ from that probe by at most ten positions. Halves that point to
  different probes are rejected as ambiguous. No genome is loaded.
- **The route defaults are set for you.** With `--flex yes`, `flexNoAlign`
  becomes 1, `flexPipelineNTriage` and `flexPipelineNSolo` become 0 (the fully
  fused pipeline), and `--outSAMtype` becomes `None`. Any of these given
  explicitly is kept, and the log records which were defaulted.
- **The cache is required.** A Flex alignment run without an H1X2 cache stops
  and prints how to build one: `--runMode hashCacheGenerate --hashCacheTiers
  H0,H1X2`. `--hashCacheTiers` now defaults to `H0,H1X2`. When
  `--soloHashScreenFile` is omitted, STAR looks next to the probe list for
  `flex_h01x2_cache.half.khash`, then `flex_h01x2_sequence_cache.bin`.
- **Results do not change.** On the JAX SC2300771 fixture (8 lanes of 100,000
  read pairs each), plain `--flex yes` with the cache produces output identical,
  file for file, to the explicit route flags and to STAR Suite 1.9.3 (139 files).

## Flex: other routes are legacy

- **New `--flexLegacy yes|no` (default `no`).** Each of the following now stops
  with an explanation unless `--flexLegacy yes` is given: an H0/H1 cache without
  the H1X2 tier, aligning cache misses (`--flexNoAlign 0`), `--no-hash-screen
  yes`, `--soloInlineHashMode no`, and SAM or BAM output. Because BAM output
  needs alignment, Flex BAM with CB/UB tags and Y-chromosome splitting are
  legacy too.
- **The H1 and H2 cache tiers are legacy.** Generating them logs a notice.
- **Spatial Flex is legacy until it is migrated.** `--soloSpatialFlexIntegrated
  yes` runs only with `--flexLegacy yes`.
- **Cache generation is unaffected.** `--runMode hashCacheGenerate` is exempt,
  so a half-probe cache can always be built.
- **To reproduce a result from an earlier release,** add `--flexLegacy yes` to
  the original command.

## Flex memory

- **A full 320k Flex run needs 57% less memory.** On all 7,303,142,230 read
  pairs (48 threads, no spill), peak resident memory is 75.3 GiB from CBQ input
  and 74.4 GiB from BGZF FASTQ, against 176.1 GiB for 1.9.3 on the same
  instance. CBQ wall time falls 9.7%. All 95 compared scientific output files
  and all 333,439 called cells are identical. These figures were measured on a
  48-thread cloud instance; smaller hosts have not yet been measured.
- **How.** Finished barcode buckets are read once and freed; the sorted merge
  flows straight through within-gene UMI correction; the caller writes raw
  matrices directly and shares one sparse count matrix across samples; gDNA
  diagnostics keep compact counts instead of a molecule ledger; barcode strings
  are borrowed rather than copied; packed records are merged without expanding
  them; and each in-memory record now takes eight bytes. The barcode is stored
  relative to its bucket (256 buckets leave 2,880 barcodes per bucket, which fit
  in 12 bits), and counts above 62 are kept exactly in a separate list.

## Fixes

- **Flex cell calling read one setting from uninitialized memory.** An old copy
  of `libflex` under `core/legacy/source/` shadowed the live `FlexFilter.h` when
  STAR's Flex caller was compiled, while `libflex.a` was compiled against the
  live header. The live header has a `useThreadPermits` field that the old copy
  lacks, so the library read that setting from an uninitialized byte. The
  observed effect was the tag-aware caller's thread scheduling: one test run
  used fixed scheduling (2 workers per sample group instead of sharing 32), and
  its matrices and cell calls were identical to the permit-scheduled runs. The
  defect is also present in 1.9.3; every archived 1.9.3 benchmark Flex run
  logged the intended permit scheduling. The stale copy is removed.

## Launchers, recipes and documentation

- **`scripts/run_flex_cr_config.sh`** now defaults to `--out-samtype none`,
  accepts `--hash-cache FILE`, adds `--flexLegacy yes` for the `bam-unsorted`
  and `bam-sorted` modes, and rejects unrecognised `--out-samtype` values
  instead of silently producing BAM.
- **The `star_flex_fixed_rna` and `star_flex_fixed_rna_cbq` recipes** require
  `solo_hash_screen_file`, accept an optional `flex_legacy`, and the CBQ recipe
  defaults `out_sam_kind` to `None`.
- **`README.md` and `flex/README_flex.md`** describe the half-probe route as the
  Flex method and mark the alignment routes as legacy.

## Tests

- New `tests/test_flex_v194_default_route.sh` checks the default route against
  the explicit flags (and, when `STAR_REF_BIN` is set, a reference binary), each
  legacy route's stop message, and that `--flexLegacy yes` and cache generation
  are not blocked. It also requires permit scheduling with the full thread
  budget in every sample group on both the default and legacy routes, and, in a
  build tree, that every `libflex` header STAR depends on resolves under
  `flex/source/libflex`.
- Existing Flex tests that exercise alignment, BAM output or H0/H1 caches now
  pass `--flexLegacy yes`, so they keep testing what they tested before.

## Known limitations

- Spatial Flex has not been migrated to the half-probe route.
- The fused Flex CBQ planner still ignores `--readMapNumber`
  (`docs/TODO_FLEX_CBQ_READ_LIMIT_20260912.md`); this is not new in 1.9.4.
- No command writes the compact `.half.khash` form. `--runMode
  hashCacheGenerate --hashCacheTiers H0,H1X2` writes an FH01SEQ1 cache with the
  H1X2 tier, which the default route accepts.
