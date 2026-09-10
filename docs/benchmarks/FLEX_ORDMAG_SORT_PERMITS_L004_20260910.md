# OrdMag sort reuse and thread permits: L004, 2026-09-10

The previous estimator copied and sorted the same bootstrap count vector for
every trial expected-cell value. It now sorts once per bootstrap replicate and
reuses the sorted vector for all trials. The trial grid, floating-point
rounding, inclusive cutoff, loss comparison and random draws are unchanged.
The public unsorted-input helper retains its original contract.

A shared permit pool also replaces fixed per-group worker shares. Every sample
coordinator reserves one worker; bootstrap and EmptyDrops tasks borrow spare
permits without blocking their coordinator. An active Monte Carlo sampler can
acquire returned permits between 64-iteration batches. Iteration seeds and
integer tallies preserve the exact p-values. Fused tags remain in one sample
model, results keep whitelist order, and joint occupancy follows every group.

## L004 benchmark

Same 1,823,648,323 L004 read pairs in CBQ, 48 threads, cold OS page cache,
`m6id.12xlarge` instance `i-06de289faa5d78117`. New arms run serially, once each,
with fresh outputs. The saved baseline is the preceding parallel-caller run.
All use the compact half-khash cache, 19,068 modeling features, 18,129 filtered
features, 100,000 simulations, BH FDR 0.01, and caller diagnostics. No BAM,
per-read decision sidecar, alignment or reference-index loading.

| Implementation | Caller seconds | Total seconds | Peak RSS GiB | Total CPU seconds |
|---|---:|---:|---:|---:|
| Previous fixed group shares | 244 | 428.95 | 75.73 | 14603.61 |
| Sort once, fixed group shares | 28 | 211.72 | 77.05 | 5944.11 |
| Sort once + shared permits | 18 | 201.90 | 76.02 | 6053.83 |

Sorting once reduces caller time **244 -> 28 seconds** (8.71 times faster).
Permits reduce the remaining caller time **28 -> 18 seconds** (35.7% less),
saving another 9.82 seconds overall.
Combined total runtime is **2.12 times faster**, or
52.9% less wall time than the previous 428.95-second run.
Compared with the earlier serial-group caller's 522.28 seconds, the combined
implementation takes 201.90 seconds. The sort change provides the largest gain.

Both new runs spend 18 seconds in startup and 100 seconds in read processing.
The full post-read phase falls from 299 to 84 to 73 seconds. Phase timestamps
have one-second resolution; total wall includes cleanup. These are single
executions, not repeated-trial estimates. Permit CPU time is slightly higher
than fixed-share sort-only CPU time despite lower wall time.

The permit pool records capacity 48, peak reserved 48, 832 acquisitions, and
all 48 permits returned. The six samplers launch 43, 44, 45, 46, 47 and 48
workers over their respective executions. Those are per-sampler totals, not
simultaneous allocations: the shared pool bounds the combined budget at 48.
A regression test also demonstrates a running executor growing from one
worker to four after another simulated group returns its permits.

## Exact regression results

Both new implementations match the saved control's complete raw matrix and
all eight filtered matrices, including axes and every coordinate. All
249,194 cell identities, aggregate read classifications, and per-sample stage
summaries are identical. Each arm compares all 57 caller diagnostic files:
OrdMag ranks, recovered-cell estimates, ambient profiles, candidate p-values
and decisions match. Only the `mc_threads` JSON field is excluded because the
execution limit changes. Occupancy still removes 3,348 calls. Existing CR
concordance is therefore unchanged; CR was not rerun or inspected for this task.

| Sample | Sort-only preparation s | Sort-only caller s | Permit preparation s | Permit caller s |
|---|---:|---:|---:|---:|
| Glioblastoma_BC1-2 | 2.13 | 1.13 | 2.13 | 0.95 |
| Colorectal_BC3-4 | 5.55 | 13.23 | 5.82 | 5.83 |
| LungCancer2_BC5-6 | 6.17 | 18.15 | 6.63 | 7.24 |
| BreastCancer1_BC7-8 | 1.43 | 0.79 | 1.39 | 0.62 |
| LNReactive_BC9-10 | 2.70 | 4.98 | 2.71 | 2.01 |
| Kidney_BC11-12 | 2.44 | 2.37 | 2.44 | 1.64 |
| SkinMelanoma_BC13-14 | 3.70 | 4.07 | 3.68 | 2.79 |
| Endo_BC15-16 | 5.05 | 9.22 | 5.29 | 3.70 |

Durations overlap and must not be added to obtain wall time. Total group wall
is 24.37 seconds with fixed shares and 13.93 seconds with permits; the surrounding
caller phase also includes common preparation and joint occupancy.

## Implementation and tests

- `OrdMagStage.cpp` uses one private sorted-input helper for all estimator trials.
- `ThreadPermits.h` provides a call-local budget and RAII permit/context guards.
  Every coordinator has a reserved worker; borrowing helpers never blocks
  progress. `ParallelTasks.h` joins all workers before returning or rethrowing.
- `EmptyDropsCRSampler.cpp` processes independent iteration batches with
  disjoint worker-local integer tallies. Worker counts cannot change draws.
- `FlexFilterTagAware.cpp` schedules groups using shared permits by default.
  Library comparisons can retain fixed shares through
  `FlexFilter::Config::useThreadPermits=false`. No new STAR CLI flag is needed.
- Clean Chromap-enabled builds for sort-only and final implementations, based
  on `4c33c014` plus the previously validated uncommitted hash/caller changes.
- Five sort-only test executables and seven final test executables pass.
  Coverage includes 90 exact comparisons against the frozen pre-change STAR
  estimator, uneven bootstrap partitions/custom seeds, both bootstrap stages,
  live permit borrowing, budget bounds, exception cleanup, nontrivial exact
  Monte Carlo tallies, queued groups, fused tags, empty tags, quality ties,
  the 500-UMI floor and joint occupancy.
- Each new binary passes fresh full-STAR CBQ and BGZF fixtures before its full
  L004 CBQ run. Whole-L004 BGZF was not repeated for this caller-only change.
  No application execution was repeated.

## Artifacts and installation

- Report: `docs/benchmarks/FLEX_ORDMAG_SORT_PERMITS_L004_20260910.md`.
- Local: `/mnt/pikachu/star_suite_paper/analysis/flex_ordmag_sort_permits_20260910/`.
- Cloud: `/scratch/flex_ordmag_sort_permits_20260910_v1/`.
- S3 bucket `star-suite-320k-benchmark-alt-171440768238-us-west-2-20260904`,
  prefix `analysis-tools/flex_ordmag_sort_permits_20260910_v1/`.
- Each arm has `BENCHMARK_COMPLETE.json` and `VALIDATION_COMPLETE.json`.
  New logs, diagnostics, tests, source manifests and binaries are archived.
  Large identical matrices reuse the durable `flex_half_production_20260910_v1`
  archive, with fresh hashes proving equivalence.
- Sort-only binary SHA256: `985580edc670c2599079aaec237102daf673090d0bca08d7577c786aa08cd3c2`.
- Installed final binary SHA256: `255821386fb760c2a4480186a7afa450e9d787a8bd22b5f5d98cb1246c027e6b`.
- The prior executable is saved locally as `STAR.before_sort_permits`.
  Code remains uncommitted, and the cloud instance remains running.

Preceding reports: [parallel groups](FLEX_PARALLEL_CALLER_L004_20260910.md),
[half-khash integration](FLEX_HALF_KHASH_L004_20260910.md).
