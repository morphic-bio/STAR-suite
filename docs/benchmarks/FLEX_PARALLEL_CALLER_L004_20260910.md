# Parallel Flex sample calling: L004, 2026-09-10

The tag-aware integration previously called sample groups sequentially and
capped each group's EmptyDrops sampler at eight threads. The new scheduler
calls independent groups concurrently under one `--runThreadN` budget.
For L004, eight groups receive six workers each. Each sample's paired tags
remain in one model. Results retain whitelist order; joint observed-tag
occupancy runs only after every group completes.

## Results

Same L004 CBQ input: 1,823,648,323 read pairs, 48 threads, cloud instance
`i-06de289faa5d78117` (`m6id.12xlarge`), cold OS page cache. Both runs use the
compact half-khash cache, the complete 19,068-feature model, the 18,129-feature
filtered export, 100,000 simulations, BH FDR 0.01 and caller diagnostics.
No BAM, per-read decision sidecar, genome index loading or alignment.

| Metric | Saved serial control | Parallel groups |
|---|---:|---:|
| Total wall seconds | 522.28 | 428.95 |
| Cell calling/filtering seconds | 328 | 244 |
| Startup seconds | 18 | 18 |
| Read phase seconds | 110 | 101 |
| Entire post-read phase seconds | 384 | 299 |
| Peak RSS GiB | 75.11 | 75.73 |
| Called cells | 249,194 | 249,194 |

Cell calling/filtering is **1.34 times faster**. Total
runtime fell **17.9%** (1.22 times faster).
This is one execution per implementation, with the saved serial control from
[FLEX_HALF_KHASH_L004_20260910.md](FLEX_HALF_KHASH_L004_20260910.md).
Phase timestamps have one-second resolution; total wall includes cleanup.

All raw coordinates, all eight filtered matrices and all called barcode
identities match exactly. Every aggregate hash/read classification and the
per-sample stage summary also match. The comparison covers all
57 caller diagnostic files: OrdMag ranks, ambient profiles, candidate
p-values, decisions and bootstrap estimates are identical. Only the recorded
`mc_threads` field is excluded from diagnostic JSON equality (eight before,
six now). Joint occupancy still removes 3,348 calls and returns 249,194 cells.
Consequently, the previously measured CR concordance is unchanged.

## Per-group timing

These are overlapping durations, not additive wall times. Preparation copies
the group's sparse matrix and identities. Caller duration includes its
OrdMag/EmptyDrops calculations and diagnostic output. Worker budgets stay
with each group worker; this implementation does not redistribute a finished
group's budget to another group's running sampler.

| Sample | Workers | Prepare seconds | Caller seconds |
|---|---:|---:|---:|
| Glioblastoma_BC1-2 | 6 | 2.15 | 177.62 |
| Colorectal_BC3-4 | 6 | 6.47 | 225.87 |
| LungCancer2_BC5-6 | 6 | 7.36 | 232.87 |
| BreastCancer1_BC7-8 | 6 | 1.40 | 170.34 |
| LNReactive_BC9-10 | 6 | 2.72 | 187.71 |
| Kidney_BC11-12 | 6 | 2.46 | 195.83 |
| SkinMelanoma_BC13-14 | 6 | 3.68 | 209.53 |
| Endo_BC15-16 | 6 | 5.05 | 203.69 |

The gain is limited by work that already used the full bootstrap budget in
the serial caller. Each concurrent group now has six execution workers.
BreastCancer1 and Glioblastoma test no tail candidates, yet spend 170–178
seconds inside their callers. This points to the work before tail rescue as
a useful next profiling target; these timings do not isolate bootstrap,
ambient estimation and diagnostic output from one another.

## Reproducibility and validation

OrdMag historically selects bootstrap seeds and sample ranges from its thread
partition. Reducing that partition would change the estimator. The new
`maxConcurrentThreads` / `bootstrapWorkers` execution cap retains the previous
logical streams and chunk boundaries, and schedules them on fewer workers.
The `--soloCellFilterBootstrapThreads` setting therefore retains its previous
statistical meaning. EmptyDrops already seeds by simulation index and merges
integer tallies, so a changed MC worker count preserves its results.

A bounded task helper lets the calling thread participate, writes disjoint
result slots and joins every worker before propagating an exception. Worker
failure clears partial group results before joint occupancy can run.

- Clean Chromap-enabled STAR build based on `4c33c014` plus the previously
  validated half-khash work and this caller change.
- Five local executables passed: bootstrap parallelism, quality ties, the
  500-UMI floor, observed occupancy, and grouped caller integration.
- Bootstrap tests preserve exact estimates and membership for 48 configured
  streams, 1/6/13 execution workers, uneven bootstrap counts and custom seeds.
- Grouped fixture compares serial and parallel stage membership, p-values,
  calls, fused-tag behavior and empty-tag handling, and tests worker failure.
- Fresh full-STAR CBQ and BGZF fixtures match their saved controls exactly.
- One full L004 CBQ run proves count and caller diagnostic parity. Whole-L004
  BGZF was not repeated for this caller-only change.

The first fixture comparison expected MEX files in empty sample directories
created for gDNA metadata. The check was corrected to use actual matrices and
their barcode-file inventory, then continued from the completed CBQ fixture.
No STAR execution was repeated. The original path-error trace is retained.

## Source and artifacts

- Canonical implementation: `flex/source/libflex/FlexFilterTagAware.cpp`,
  `core/legacy/source/SoloFeature_flexfilter.cpp`, and shared `libscrna`.
- Local artifacts: `/mnt/pikachu/star_suite_paper/analysis/flex_caller_parallel_20260910/`.
- Cloud artifacts: `/scratch/flex_caller_parallel_20260910_v1/`.
- S3 bucket: `star-suite-320k-benchmark-alt-171440768238-us-west-2-20260904`,
  prefix `analysis-tools/flex_caller_parallel_20260910_v1/`.
- Fresh logs, diagnostics, source, tests, binary and validation hashes are
  archived there. Exact matrices reuse the existing durable archive under
  `analysis-tools/flex_half_production_20260910_v1/`.
- `BENCHMARK_COMPLETE.json` and `VALIDATION_COMPLETE.json` establish completion.
- Production binary SHA256: `dfa08ff2b26f1e4bb6a154527c5e878acda709ef0a7f743a10974c50710537ee`.

The validated binary is installed at `core/legacy/source/STAR`. Its predecessor
is saved as `STAR.before_parallel` in the local artifact directory. Source
changes remain uncommitted. The cloud instance remains running.
