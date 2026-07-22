# Runbook: Solo Overnight Benchmarks

Date: 2026-03-25  
Audience: future Solo/perturb benchmarking passes after the 2026-03-25 bridge-collapse optimization work

## Goal

Freeze the current Solo optimization state into a benchmark plan that can be
rerun later without reconstructing context.

This runbook defines:

- the current optimized STAR build to benchmark
- the no-BAM GEX-only comparison matrix for UCSF
- the perturb-seq benchmark matrix for UCSF, MSK, and A375
- the requirement to use internal gzip on STAR runs
- the 100K smoke-fixture family that should be carried forward

## Current Optimized Solo Build

The current optimized non-Flex Solo bridge candidate lives in the clean detached
worktree:

- `/tmp/star-suite-v10-redesign-20260325`

Current tip of that worktree:

- `004c36e9` `Use khash for Solo bridge CB counts`

Important preceding performance commits in that line:

- `cf517678` `Parallelize Solo bridge CB collapse`
- `a1bd4231` `Bucket Solo bridge collapse by gene`
- `183c2160` `Add Solo bridge replay harness cleanup`
- `f7bb16e4` `Optimize Solo bridge hash collapse`
- `6722f71` `Redesign Solo bridge tuple and ambiguous aggregates`

For overnight runs, treat `004c36e9` as the current optimized STARsolo
candidate unless a newer explicit replacement is nominated.

## Benchmark Hygiene

- Do not run benchmark jobs in parallel.
- Do not use the dirty `/mnt/pikachu/STAR-suite` checkout for production
  timing runs.
- Use a clean clone or clean worktree for overnight benchmarks.
- Always clean rebuild before benchmarking:

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
```

- Record all new run roots in `tests/ARTIFACTS.md` after the runs complete.

## README Metric Bundle

The overnight benchmark tables should collect the same top-level metrics that
are already summarized in `README.md`.

### For GEX-only runs

Collect:

- filtered cell count
- wall time
- max RSS
- speedup vs comparator
- top-line `Summary.csv` values
- if a CR reference is available:
  - barcode Jaccard
  - gene Pearson
  - cell Pearson

### For perturb-seq runs

Collect:

- filtered cell count
- wall time
- max RSS
- speedup vs CR9
- barcode Jaccard
- gene Pearson
- cell Pearson
- CRISPR match rate
- UMI Pearson if already available from the parity surface

Also collect a phase breakdown so the README table can be refreshed later:

- genome load
- feature assignment
- mapping
- Solo counting
- PfMulti merge + calling
- MEX writing / `writeCombinedMex` where relevant

If an overnight benchmark cannot produce one of these fields cleanly, note the
missing field explicitly in the benchmark summary rather than silently dropping
it.

## Mandatory STAR Benchmark Settings

### Internal gzip

For all STAR benchmark runs in this plan:

- use gzipped FASTQs directly
- do **not** pass `--readFilesCommand zcat`

Rationale:

- removes external decompressor process overhead
- gives cleaner timing
- matches the intended post-2026-03-25 benchmark surface

If an existing script hardcodes `--readFilesCommand zcat`, fork or override the
script for the overnight run rather than benchmarking the external decompressor
path.

### No-BAM default

For the primary overnight comparison matrix:

- use `--outSAMtype None`

Only the smoke/fixture track below should keep sorted-BAM or Y-removal variants.

### Canonical Solo settings for optimized STAR GEX-only runs

For the optimized Solo no-BAM GEX-only runs:

```bash
export STAR_SOLO_NONFLEX_HASH_BRIDGE=1
```

and include:

```bash
--soloInlineHashMode yes
--soloMultiMappers Unique
--soloUMIdedup 1MM_CR
--soloUMIfiltering MultiGeneUMI_CR
--soloCellFilter EmptyDrops_CR
--soloCrMultimapRescue yes
--soloFeatures GeneFull
--outSAMtype None
```

### Canonical perturb-seq settings

Use the shared perturb-seq methodology in:

- [PAPER_BENCHMARK_METHODOLOGY.md](/mnt/pikachu/STAR-suite/docs/PAPER_BENCHMARK_METHODOLOGY.md)

In particular, do not omit:

- `--soloCrMultimapRescue yes`
- `--dynamicThreadInterface 1`
- `--crAssignConsumerThreads -1`
- `--crAssignSearchThreads 1`
- `--soloFeatures GeneFull`

## Benchmark Track A: UCSF GEX-Only, No BAM

This track is the simplest direct Solo timing comparison.

### Inputs

- corrected UCSF full GEX root:
  - `/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2/GEX`
- genome:
  - `/storage/autoindex_110_44/bulk_index`
- Solo whitelist:
  - `/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-february-2018_TRU.txt`

### Runs

Run these arms serially. **Do the historical vanilla STARsolo arm last** in the
**full overnight plan** (after the perturb no-BAM matrix; see [Recommended
Overnight Order](#recommended-overnight-order) step 5). Do not wrap it in a
short tool timeout; allow multi-hour wall time. Logical order for the three
arms:

1. optimized STARsolo no-BAM (internal gzip;
   `scripts/paper/run_ucsf_gexonly_no_bam_benchmark.sh` `--modern-optimized`)
2. CR9 GEX-only no-BAM reference (only if you are timing a fresh CR9 run;
   otherwise reuse the stored reference run)
3. vanilla STARsolo no-BAM (`7a7fb08` / `--historical-vanilla`) — **last**

### Binary provenance

#### CR9

If timing a fresh Cell Ranger run:

- prefer `--create-bam=false` when the command surface supports it
- otherwise note the exception explicitly in the benchmark summary

If only parity/reference comparison is needed, reuse the existing UCSF CR run:

- `/storage/ucsf-full/bench_20260218_dynamic_first/cellranger_runs/cr_full_iPSC2_1_AALG2_crstar32_20260218_205804`

#### Vanilla STARsolo

Use the pre-optimization baseline already referenced in the Solo handoff:

- base commit: `7a7fb08`

If a fresh clean worktree is needed, build that exact commit and record the
worktree path in the benchmark notes.

#### Optimized STARsolo

Use:

- commit `004c36e9`
- clean worktree or clean clone built from that revision

### Suggested output roots

Use a fresh root such as:

- `/storage/solo_overnight_20260325/ucsf_gexonly_no_bam/`

with subdirectories:

- `cr9/`
- `star_vanilla_7a7fb08/`
- `star_optimized_004c36e9/`

### Metrics to collect

- wall time
- max RSS
- speedup vs comparator
- filtered cell count
- mapping elapsed
- `collapseUMIall` or `collapseUMIall_fromBridgeHash`
- `countCBgeneUMI`
- `processRecords`
- top-line `Summary.csv`
- barcode Jaccard vs CR9 where available
- gene Pearson vs CR9 where available
- cell Pearson vs CR9 where available
- parity vs CR9 filtered/raw GeneFull MEX where applicable

## Benchmark Track B: Perturb-Seq, No BAM

This track measures the integrated perturb workflow with the current STAR
optimizations and internal gzip.

### Shared requirements

- internal gzip only
- `--outSAMtype None`
- canonical perturb parameters from the paper methodology
- optimized STAR build at `004c36e9`

### B1. UCSF EBs2_2

Use the existing paper benchmark inputs:

- script reference:
  - [run_ucsf_ebs2_2_benchmark.sh](/mnt/pikachu/STAR-suite/scripts/paper/run_ucsf_ebs2_2_benchmark.sh)
- dataset root:
  - `/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2`
- feature ref:
  - `/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv`
- genome:
  - `/storage/autoindex_110_44/bulk_index`

Override the script behavior for the overnight run so that STAR uses internal
gzip rather than `--readFilesCommand zcat`.

Suggested root:

- `/storage/solo_overnight_20260325/ucsf_perturb_no_bam/`

### B2. MSK 30polyKO

Use the existing paper benchmark inputs:

- script reference:
  - [run_msk_30polyko_benchmark.sh](/mnt/pikachu/STAR-suite/scripts/paper/run_msk_30polyko_benchmark.sh)
- fastq root:
  - `/storage/MSK-perturb-comparison/msk30ko_full_3lib_20260304_095911/fastqs`
- gRNA feature ref:
  - `/mnt/pikachu/MSK-whitelists/ref_feature_geneBC.csv`
- LARRY feature ref:
  - `/mnt/pikachu/MSK-whitelists/ref_feature_larryBC.csv`
- genome:
  - `/storage/autoindex_110_44/bulk_index`

Again, remove `--readFilesCommand zcat` for the overnight run.

Suggested root:

- `/storage/solo_overnight_20260325/msk_perturb_no_bam/`

### B3. A375

Use the existing paper benchmark inputs:

- script reference:
  - [run_a375_benchmark.sh](/mnt/pikachu/STAR-suite/scripts/paper/run_a375_benchmark.sh)
- dataset root:
  - `/storage/A375`
- fastq root:
  - `/storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs`
- feature ref:
  - `/storage/A375/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv`
- whitelist:
  - `/storage/A375/3M-5pgex-jan-2023.txt`
- genome:
  - `/storage/autoindex_110_44/bulk_index`

Remove `--readFilesCommand zcat` for the overnight run.

Suggested root:

- `/storage/solo_overnight_20260325/a375_perturb_no_bam/`

### Metrics to collect for perturb track

- wall time
- max RSS
- speedup vs CR9
- filtered cell count
- mapping elapsed
- Solo timing lines
- phase breakdown:
  - genome load
  - feature assignment
  - mapping
  - Solo counting
  - PfMulti merge + calling
  - MEX writing
- `outs/filtered_feature_bc_matrix` parity vs CR reference where available
- `crispr_analysis` presence
- barcode Jaccard / cell Pearson / gene Pearson using the canonical parity script
- CRISPR match rate
- UMI Pearson where the parity tooling already emits it

## 100K Smoke-Fixture Family To Carry Forward

The benchmark plan should be backed by small stable fixtures so new changes can
be checked quickly before full overnight runs.

### Smoke fixture goals

Maintain 100K-class fixtures for:

1. no-BAM
2. sorted-BAM
3. Y-removal where biologically relevant

### Fixture families

#### UCSF perturb

Existing canonical 100K perturb fixture:

- `/storage/ucsf-2M/fixtures/ucsf2m_iPSC2_AALG2_100k_pfconfig`

Carry forward three test surfaces:

- `ucsf_perturb_100k_no_bam`
- `ucsf_perturb_100k_sorted_bam`
- `ucsf_perturb_100k_yremove`

Notes:

- `ucsf_perturb_100k_no_bam` can reuse the existing fixture root directly
- `ucsf_perturb_100k_yremove` should continue to use the UCSF GEX-only
  Y-removal smoke surface in
  [run_solo_yremove_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_solo_yremove_smoke.sh)
  and
  [run_ucsf_perturb_yremove_batch.sh](/mnt/pikachu/STAR-suite/scripts/run_ucsf_perturb_yremove_batch.sh)

#### A375 perturb

Create or preserve:

- `a375_perturb_100k_no_bam`
- `a375_perturb_100k_sorted_bam`

Suggested root:

- `/storage/A375/fixtures/a375_100k_20260325/`

#### MSK perturb

Create or preserve:

- `msk_perturb_100k_no_bam`
- `msk_perturb_100k_sorted_bam`

Suggested root:

- `/storage/MSK-perturb-comparison/fixtures/msk30ko_100k_20260325/`

Use the existing fixture makers as the starting point:

- [tests/multi_feature/create_fixture.sh](/mnt/pikachu/STAR-suite/tests/multi_feature/create_fixture.sh)
- [tests/multi_feature/create_fixture_downsampled.sh](/mnt/pikachu/STAR-suite/tests/multi_feature/create_fixture_downsampled.sh)
- [tests/downsample_fastq_gz.sh](/mnt/pikachu/STAR-suite/tests/downsample_fastq_gz.sh)

### What to store with each 100K fixture

Each carried fixture should include:

- source FASTQ provenance
- downsample seed
- exact read count target
- chemistry / whitelist choice
- whether BAM is expected
- whether Y-removal is expected
- a `README` or manifest with the canonical STAR command family

## Recommended Overnight Order

Run serially in this order:

1. UCSF GEX-only no-BAM — **everything except historical vanilla first:**
   - optimized STARsolo (internal gzip)
   - CR9 if rerunning rather than reusing reference
2. UCSF perturb no-BAM
3. MSK perturb no-BAM
4. A375 perturb no-BAM
5. UCSF GEX-only historical vanilla STARsolo (`7a7fb08` / `--historical-vanilla`)
   **last** — run **without** a short session/tool timeout (e.g. not 600s);
   allow multi-hour wall time so Solo can finish naturally.

If time remains after the primary matrix:

6. refresh or regenerate the 100K smoke fixtures
7. run the smoke fixtures on the optimized STAR build

## Script Follow-Up Needed

Several existing paper/benchmark scripts still hardcode:

```bash
--readFilesCommand zcat
```

Before using them for the overnight matrix, either:

- patch the scripts to omit `--readFilesCommand`, or
- create benchmark-only wrappers that remove the external decompressor path

Do not mix the old external-zcat timings with the new internal-gzip timings in
the same benchmark table.

## Expected Deliverables After Overnight Runs

- one benchmark summary table for UCSF GEX-only no-BAM
- one benchmark summary table for perturb-seq no-BAM across UCSF/MSK/A375
- updated `tests/ARTIFACTS.md` with all new roots
- if results are stable, promote the optimized build and command family into
  the benchmark scripts themselves
