# Runbook: Streaming CB bucketing with spill — Flex post-map tail parallelization

Date: 2026-09-02. Branch: cut `core-bucket-tail` from the head of `core-bgzf-ingest`
(this work builds directly on the BGZF ingest module and its benchmark conventions).

## Why (measured, not speculative)

The BGZF/CBQ ingest work removed the front-end bottleneck; Amdahl now points at the
post-map tail. From the committed benchmark `jax_bgzf_20260901T221610Z` (full JAX,
2.011B pairs, 32 threads, 7:38.71 total):

| Phase | Wall | Threading today |
| --- | ---: | --- |
| Streaming (ingest + assignment + solo feed) | ~247 s | permit-scheduled, readerQ=0, healthy |
| `sumThreads` region (gather per-thread buffers) | ~56 s | none |
| `countCBgeneUMI` (group records by CB) | 90.8 s | none |
| direct hash-collapse core | 6.1 s | OpenMP (`collapseThreads`) — keep as-is |
| MEX write | 5 s | none |
| flexfilter (per-tag cell calling) | 28 s | none |
| collapse wrapper remainder (readInfo/stats) | ~13 s | none |

The tail (`processRecords`, 208.3 s) is ~45% of wall time and ~95% serial.
Threading audit: no threading constructs in `SoloFeature_sumThreads.cpp`,
`SoloFeature_countCBgeneUMI.cpp`, `SoloFeature_flexfilter.cpp`,
`SoloFeature_collapseUMIall.cpp`. OpenMP exists only in
`SoloFeature_collapseUMI_fromBridgeHash.cpp` (`#pragma omp parallel
num_threads(collapseThreads)`, `schedule(dynamic, 16)`) — use it as the house
pattern.

Separately, the 320k public run (`scffpe320k_bgzf_20260902T071608Z`) was
OOM-killed at 118.5 GB anon RSS with 4.34B of 7.3B pairs ingested (59%;
readerQ=0 throughout — the reader was never the problem). Memory model validated
on both datasets: **RSS ≈ ~30 bytes × kept reads** (JAX 1.68B kept ≈ 50 GB
observed; 320k ~6.4B kept → ~190 GB needed on a 126 GB box). Fixed structures
are small (hash cache 0.2 GB). The per-read accumulation is the whole problem.

One design solves both: bucket records by cell barcode as they are produced.

## Design

### 1. Streaming CB-hash bucketing (one primitive, two backends)

During the streaming phase, route every kept record into one of P buckets by CB
index (high bits of the CB index over the 737,280-entry whitelist; P
power-of-two, default 256). Producers append to per-worker, per-bucket segments;
a sealed segment is handed to the bucket store under an atomic claim — mirror
the `BgzfRangeReader::claim_work` frontier pattern (claim under lock, do the
heavy work outside it).

Packed record target ≤ 16 bytes: CB index (20 bits), UMI (24 bits for 12 bp),
gene/probe index (~18 bits), tag + status flags, remainder reserved. Byte-layout
documented in a header so the spill file is a defined format.

Backends behind one interface:
- **RAM mode**: sealed segments stay in memory (default when the dataset fits).
- **Spill mode**: sealed segments append to one file per bucket under a scratch
  directory. Sequential writes, overlapped with ingest; nothing is read back
  until the tail. No indexes needed — buckets are self-contained by
  construction.

Flags:
- `--soloBucketMode ram|spill|auto` (default `auto`)
- `--soloBucketMemGB <n>` — budget for `auto`: run in RAM until the packed
  accumulation crosses the budget, then flush all buckets to disk and continue
  in spill mode (one-way transition; must be output-equivalent to spill-from-start)
- `--soloBucketSpillDir <path>` (default: `outFileTmp`)
- `--soloBucketCount <P>` (default 256)

### 2. Bucket-parallel tail

Replace the monolithic gather + global CB grouping with: tail workers atomically
claim buckets, load (RAM segment list or sequential file read), group/sort
within the bucket, run the existing per-CB collapse logic, and emit
per-bucket partial matrices. CB ranges are disjoint across buckets, so the final
MEX is a concatenation in CB order — no merge step exists anywhere.
Overlap read-back with compute (claim bucket k+1's file read while processing k).

Consequences (do not implement as separate features — they fall out):
- `sumThreads` is eliminated for this path (nothing to gather).
- `countCBgeneUMI`'s global pass becomes a small per-bucket group.
- The OpenMP collapse core runs per-bucket unchanged.

### 3. flexfilter parallel across tags

Per-tag cell calling is independent — `omp parallel for` over tags with the
house pattern. Matters at 16 tags / 320k cells far more than on JAX's 4 tags.

### 4. Count-only bookkeeping audit

In `--flexNoAlign 1` count-only mode, audit the collapse wrapper's per-read
bookkeeping (`recordReadInfo` writes, stats passes): keep exactly what the
requested outputs need (tag tables when requested), skip the rest behind the
existing output-mode checks. No behavior change for runs that request
read-level outputs.

### Scope guards

- Scope to the fused Flex inline-hash path (where this tail runs). Classic
  STARsolo paths and all non-Flex behavior untouched; new behavior behind the
  flags above, default `auto` may be enabled for Flex only.
- Byte-identical outputs: matrices, per-sample outputs, stats, and the pipeline
  counters must equal the current path on identical inputs, across thread
  counts, and across `ram`/`spill`/`auto` (including the mid-run transition).
  UMI collapse is CB-local and sorts within CB, so equality is achievable —
  gates below enforce it.
- No htslib changes, no new dependencies, vendored code only.
- Permits unchanged: the tail is a single workload; a plain OpenMP pool is
  correct there. Do not thread the tail through the permit controller.

## Phases and expected-red TDD gates

Same discipline as the BGZF runbook: per-phase manifest at
`tests/bucket/PHASES.tsv`, tests SKIP until their enabling phase, gate runner
`tests/run_bucket_tests.sh`, `tests/bucket/IMPLEMENTED_PHASE` advanced only when
the phase's tests pass.

| Test | Enabling phase | Description |
| --- | --- | --- |
| B1 | 1 | Packed-record round-trip; bucket partition equals independent reference partition; atomic segment claims correct under 1/8/32 producer threads |
| B2 | 1 | Spill segment file round-trip byte equality with RAM segments |
| B3 | 2 | Bucket-parallel tail: matrices + per-sample outputs byte-identical to current path on gold fixture and JAX 800k downsample |
| B4 | 2 | Determinism across 1/8/32 tail threads and across bucket counts (64/256/1024) |
| B5 | 3 | Spill end-to-end equality, including the auto mid-run RAM→spill transition vs spill-from-start |
| B6 | 3 | flexfilter across-tag parallel equality (per-tag outputs identical to serial) |
| B7-off | 0 | Existing gzip/BGZF/CBQ regressions with bucketing disabled |
| B7-auto | 3 | Existing regressions with `auto` default active |

Phase 0: fixtures + an independent Python reference (synthetic multi-tag FASTQ
with known per-CB/per-UMI content; reference counter). Phase 4: benchmarks.

## Phase 4 benchmarks and acceptance criteria

Benchmark conventions exactly as `docs/benchmarks/bgzf_ingest_20260901/`: cold
page cache, serialized runs on a quiet machine, `/usr/bin/time -v`, command
lines + logs + preflight committed.

1. **JAX full triple** (gzip / BGZF / CBQ, 32 threads):
   - `processRecords` tail **< 45 s** (from 208.3 s)
   - streaming throughput regression **< 3%**
   - pipeline counters identical to the committed baselines
     (BGZF: 2,011,130,186 / 1,681,459,858 keep / 16,111,757 deny / 313,558,571 miss)
2. **320k scFFPE BGZF run** (`run_320k_benchmark.sh bgzf`, spill mode):
   - completes on this machine; peak RSS **< 110 GB**
   - per-sample outputs for all 16 tags produced
   - record wall/RSS artifacts in the benchmark directory
3. RAM-vs-spill overhead on JAX recorded (expected: spill visible cost ≈ 0 —
   writes overlap ingest, read-back overlaps tail compute).

## Anchors

- Tail call site: `core/legacy/source/Solo.cpp` → `SoloFeature::processRecords`
  (`SoloFeature_processRecords.cpp`)
- Serial stages: `SoloFeature_sumThreads.cpp`, `SoloFeature_countCBgeneUMI.cpp`,
  `SoloFeature_flexfilter.cpp`
- OpenMP exemplar: `SoloFeature_collapseUMI_fromBridgeHash.cpp`
- Atomic-claim exemplar: `core/legacy/source/input/BgzfRangeReader.cpp::claim_work`
- Benchmark wrappers: `docs/benchmarks/bgzf_ingest_20260901/run_jax_benchmark.sh`,
  `run_320k_benchmark.sh`
- OOM evidence: `docs/benchmarks/bgzf_ingest_20260901/scffpe320k_bgzf_20260902T071608Z.*`
