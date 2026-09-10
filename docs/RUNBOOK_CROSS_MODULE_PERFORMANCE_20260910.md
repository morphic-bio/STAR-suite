# Cross-module performance: process_features, BGZF, cell calling, bulk and SLAM

Date: 2026-09-10

Status: **ACTIVE; user authorized implementation on 2026-09-10 ("OK - let's start"). Phase 1 PF BGZF implementation and full feature-library comparison are complete; phase 2 is next.**

The user requested an ordered implementation plan for the performance
opportunities found after FLEX 1.9.0, including BGZF and an explicit priority
for `process_features`. This document is the execution record for that next
piece of work. Creating it does not launch datasets, provision cloud resources,
or publish another release.

## Outcomes and priority

1. Deliver parallel BGZF FASTQ ingestion to **process_features first**. Measure
   feature-library ingestion/assignment separately from the whole perturb run.
2. Reuse that input infrastructure for bulk, scRNA-seq, perturb GEX and SLAM.
3. Remove duplicated caller work, use bounded parallelism, and reduce matrix
   copies in scRNA/perturb; extend scheduling to multiplexed OCM samples.
4. Cache invariant SLAM probabilities, reuse fitted results, and fit genes in
   parallel without changing the model or output order.
5. Integrate the existing bulk VB component execution option and measure its
   runtime and memory effect with the current convergence settings.

Preserve feature assignment, barcode correction, UMI resolution, cell calls,
SLAM statistical models and bulk convergence rules. Do not use a statistical
change to produce a performance result. Report ingestion, model computation,
output and end-to-end time separately where they explain the result.

## Baseline and boundaries

| Item | Starting point |
| --- | --- |
| Canonical checkout | `/mnt/pikachu/STAR-suite`, keep on `master` |
| Current master at planning | `9e146eef20080108cc075a977938757f90221dbc` |
| Immutable v1.9.0 source tag | `fc59b1b20b9fd34c1176acebd7e1688d4a82c1ec` |
| Suggested implementation branch | `perf/cross-module-bgzf-callers` in an isolated checkout |
| Proposed local evidence root | `/mnt/pikachu/star_suite_paper/analysis/cross_module_performance_20260910/` |
| Previous benchmark instance | `i-06de289faa5d78117` is **terminated**; do not submit jobs to it |
| Retained root snapshot | `snap-0428a020b3a0fbbd4`; SSD data is in S3, not this snapshot |
| Durable FLEX evidence | See [completed release/retirement runbook](RUNBOOK_V1_9_0_FULL320K_RELEASE_RETIREMENT_20260910.md) |

Refresh branch and workspace state before implementation. Preserve unrelated
edits, currently in `docs/HANDOFF_FLEX_PERF_20260903.md` and the two modified
JAX comparison sources. Keep large inputs, binaries and outputs untracked.
Use explicit staged paths, no destructive reset/cleanup, and `--no-ff` merges
for shared-source changes. Do not move v1.9.0 or infer a new release version.

Follow `AGENTS.md` and the host's `AGENTS.local.md`. Cell Ranger remains a
clean-room comparator: results, logs, public documentation, MEX/BAM/HDF5 are
allowed; its code, installed modules, source archives and source-derived notes
are excluded. The performance comparator is the accepted STAR baseline;
running Cell Ranger again is unnecessary for exact-output changes.

The FLEX no-reference benchmark setting does **not** apply to bulk/scRNA/SLAM
alignment jobs: those need their matching reference. Isolated matrix callers,
SLAM histogram fitting and standalone VB finalization should not load a STAR
alignment index merely to reach their computation stage.

## What already exists, and what remains

| Finding | Current state | Work here |
| --- | --- | --- |
| OrdMag sort once per bootstrap draw | Shared `libscrna`, already shipped; bootstrap scRNA/perturb GEX reaches it | Preserve it; do not count it as a new port |
| EmptyDrops MC executor | Shared implementation supports multiple workers | scRNA and OCM set `mc_threads=0`, which means one worker; wire explicit budgets |
| FLEX sample permits | Group scheduling creates a bounded call-local pool | Extend the pattern to OCM; coordinate with existing core/PF scheduling |
| scRNA diagnostic OrdMag | Adapter recomputes bootstrap after the actual call | Return/reuse the original stage result |
| scRNA sparse matrix | STAR copies to separate arrays; C API reconstructs interleaved counts | Introduce validated non-owning views and remove redundant conversions |
| BGZF decoder | Existing block/range reader, CRC checks, leases and permit hooks | Reuse; do not create another decoder in PF |
| PF record APIs | Existing queued and worker-owned direct record-view entry points, used by CBQ | Add BGZF adapter and native dispatch to these paths |
| Core BGZF mapping | Restricted to paired FLEX; special BAM bridge excludes SLAM, batch and other modes | Generalize by capability, with tested module integration |
| SLAM fitting | Serial gene loops; invariant PMFs recomputed every iteration | Cache probabilities, reuse results, parallelize independent genes |
| Bulk VB | Already parallel; optional component execution exists in standalone finalizer | Expose/integrate the option in normal STAR and validate it |
| FLEX half-probe cache | FLEX-specific matching representation | No general cache port to alignment-based modules |

Measured context: FLEX caller-only sorting/permits reduced L004 calling from
244 to 28 to 18 seconds, with peak RSS approximately 75.7–76.0 GiB. A separate
half-probe cache change reduced peak RSS by 18.25 GiB. Neither establishes an
EmptyDrops memory saving for other modules. The bulk component option's older
results were hardware-dependent (roughly 3–8% on EPYC and 20% locally), and are
background evidence, not predictions for the new builds.

## Ordered implementation

Complete each phase's correctness gate before using its output as the next
baseline. Use separate commits for independently measurable changes. The
default order is below; change it only by recording the reason in the ledger.

| Phase | Deliverable | Dependency / acceptance |
| --- | --- | --- |
| 0 | Pin source, input, reference, output contracts and prior evidence | All execution follows this inventory |
| 1 | **PF BGZF adapter and input dispatch** | Exact feature counts and read identities; bounded queues and permits |
| 2 | PF direct worker integration and copy reduction | Phase 1; prove no additional FASTQ/queue bottleneck or count change |
| 3 | General core BGZF support across modules | Phases 1–2 reader contract; enable simple modes before complex ones |
| 4 | scRNA/perturb OrdMag reuse, then bounded MC threading | Two commits; exact cell identities, thresholds and MC tallies |
| 5 | scRNA/perturb matrix-view memory reduction | Phase 4; exact counts with smaller caller allocation/peak memory |
| 6 | OCM sample scheduling and shared caller permits | Phases 4–5; exact per-sample calls and total worker bound |
| 7 | SLAM probability caching, then fitted-result reuse | Two commits; same fits and output fields |
| 8 | Parallel SLAM gene fitting | Phase 7; bounded workers, same gene/output ordering |
| 9 | Bulk TranscriptVB component execution | Existing engine option; unchanged model/convergence, validated numeric differences |
| 10 | Combined regression, representative benchmarks, documentation/integration | All enabled phases accepted; retain individual attribution |

### Phase 0 — inventory and comparison contract

- Save Git status, exact source revision, compiler/build flags and effective
  generated defaults. Clean-build before testing changed compiled code.
- Locate existing completed controls by binary/input/config identity. Record
  what actually ran; old handoff statements are not completion evidence.
- Choose one local PF feature library with substantial input (CRISPR first,
  plus a small ADT/HTO fixture), one ordinary scRNA and one perturb GEX count
  matrix, a multi-sample OCM fixture, the public PE bulk fixture, and SLAM
  blank/treatment SE and PE fixtures. Record exact paths, hashes, read counts,
  references and output surfaces before launching any arm. Larger local
  performance inputs must be identified here rather than invented later.
- Prefer saved canonical counts/histograms/equivalence classes for isolated
  algorithm comparisons. Cache large MEX inputs once as a versioned binary
  sparse representation when useful; preserve axes, integer widths, feature
  inclusion, barcode namespace and original MEX content hashes. Conversion is
  outside timed caller work and is not permission to regenerate counts.
- Record known order-sensitive behavior: library auto-detection, online model
  updates, feature rescue modes, read sampling, output taps and per-file trims.
- Define numerical acceptance before timing. Reuse meaningful existing bands
  for order-sensitive VB behavior; do not replace count parity with correlation.

### Phase 1 — process_features BGZF first

Source entry points:

- `core/legacy/source/input/BgzfBlockReader.{h,cpp}` and
  `BgzfRangeReader.{h,cpp}`: detection, bounded inflate work and leased buffers.
- `core/legacy/source/PfMultiAssign.cpp`: feature-library dispatch and existing
  CBQ adapters; preserve its barcode/feature read-role mapping.
- `core/features/process_features/include/pf_api.h` and `src/pf_api.c`:
  `pf_read_record_view`, `pf_process_record_views()` and direct worker APIs.
- `core/features/process_features/src/assignBarcodes.c`: existing FASTQ
  producer/consumer behavior and feature counting kernel.

Implement in small steps:

1. Define an input-neutral BGZF-to-PF adapter using the existing decoder and
   PF record views. Keep C/C++ ownership explicit; build standalone PF without
   introducing a dependency on the full STAR executable or a genome index.
2. Pair streams by logical record ordinal and the established name policy.
   R1/R2/R3 block offsets and compressed lengths are independent; never equate
   byte ranges between mates. Support the PF barcode/feature/multiple-feature
   stream layouts actually selected by the library configuration.
3. Feed bounded batches through the existing PF stream API as a correctness
   integration. Retain backing buffers until all consumers finish with them;
   only claim copy elimination where lifetimes prove it.
4. Add consistent native selection and effective-mode logs in both embedded
   `pf-multi` and standalone PF entry points. Audit whether existing input
   options can express the choice; add a small explicit option only if needed.
5. Auto-detect genuine BGZF from header metadata. Preserve explicit user
   decompression commands. For mixed lanes, select a supported per-lane path
   or record a whole-library fallback. Forced BGZF mode must explain and fail
   on unsupported input; it must not silently become ordinary gzip.

Reuse existing `BgzfWorkPermitHooks` and PF/core permit interfaces. Preserve
PF auto-sized consumers (`--crAssignConsumerThreads -1`), search-thread policy
and dynamic scheduling; do not introduce another hard-coded small consumer
cap. With GEX and PF running together, account for inflate, assignment and
mapping work under the defined total compute budget. Record coordinator
accounting and reserve progress capacity so a producer/consumer cannot hold
every permit while waiting on its counterpart. Existing CBQ/FLEX behavior
must remain covered during this integration.

Gate: same ordered record identity/sequence/quality digest, same assigned and
unassigned totals, feature/barcode/UMI decisions, deduplicated sparse counts,
axes and relevant diagnostics as the existing reader. Include unequal mate
block layouts, records spanning blocks, multiple lanes, optional third stream,
CRC failure, truncated records, mate/name mismatch, long names/reads, worker
exceptions, backpressure and early cancellation. Reject unsupported capacity
limits before counting, rather than truncating fields.

### Phase 2 — remove PF handoff overhead

Connect leased BGZF batches to worker-owned PF processing, reusing the direct
API already used by CBQ. Prefer fixing a common record-view bottleneck so CBQ
also benefits. Preserve one feature-matching implementation; avoid copying its
logic into a BGZF-only consumer. Explicitly validate asynchronous ownership
and any temporary sequence/quality materialization inside `pf_api.c`.

Keep order-sensitive feature rescue policies on their existing supported path
until an equivalent dispatch is demonstrated. Use bounded memory and stable
merging of per-worker counts/statistics. Do not require an upfront whole-file
scan or an external BGZF index merely to begin reading.

Measure these distinct surfaces once per selected build/configuration:

- Decode/pair throughput and peak queue occupancy.
- Feature ingestion plus assignment/counting throughput and peak RSS.
- Whole perturb run with GEX and PF concurrent under its usual budget.

Compare the **same BGZF files** via the existing gzip path and native BGZF to
isolate reader benefit. A paired CBQ control can contextualize throughput.
Plain-gzip transcoding/reblocking is preparatory work and must be reported
separately. Do not claim PF is the largest winner until the measurements show it.

### Phase 3 — general core BGZF integration

Reuse the shared reader and pairing contract. Replace the FLEX-only activation
decision in `Parameters_openReadsFiles.cpp`/`mapThreadsSpawn.cpp` with explicit
capability checks; audit the bridges in `ReadAlignChunk_processChunks.cpp`.
Removing a guard alone is not implementation of the excluded mode.

Enable and validate in this order:

1. Paired bulk and scRNA/perturb GEX, ordinary single-pass alignment; test
   no-BAM counting and the output modes selected for that fixture.
2. Single-end bulk/scRNA and multi-lane inputs; handle legitimate barcode-read
   layouts rather than assuming every assay is a two-read FLEX library.
3. SLAM fixed-trim SE/PE, then automatic-trim detection/reopen/replay, per-file
   processing, file skipping and batch transitions.
4. Order-sensitive production surfaces: paired input-order output, Y/noY taps,
   transcriptome BAM, two-pass/SJ filtering and batch boundaries. Each mode
   needs explicit coverage before it is enabled. Unsupported combinations
   keep a logged safe fallback in auto mode and an actionable forced-mode error.

Use stable read/lane ordinals and bounded ordered delivery where downstream
behavior depends on input order. Preserve names, qualities, clipping, read
limits and restart behavior. Trace queue wait, inflate CPU and mapper CPU to
ensure faster decompression does not oversubscribe mapping or inflate memory.
Test FLEX regression as well as the new consumers.

### Phase 4 — scRNA and perturb caller CPU work

First commit: expose the actual OrdMag result/ambient membership required by
diagnostics through an internal result/trace interface. Remove the second
bootstrap call in `SoloFeature_emptyDrops_libscrna.cpp`. Preserve the existing
external C API contract or extend it compatibly; retain diagnostics from the
same calculation that determined calls.

Second commit: pass an explicit run-budget-derived MC worker limit through
`scrna_api.cpp` to `EmptyDropsCRSampler.cpp`. Current adapter zero means one
worker, not auto-detection. Keep logical bootstrap streams/seeds separate from
execution-worker allocation: changing a permit count must not change the
bootstrap sample. Preserve simulation-index seeds and integer MC reduction.
Single-sample calls can use the available budget after alignment; permits are
most useful when independent sample work overlaps.

Validate original and changed caller results on identical counts. Require
exact selected barcodes, OrdMag estimates/ranks/ties, ambient membership and
probabilities, p-values, BH decisions and diagnostics except explicit timing/
worker metadata. Include low counts, zero candidates, ties and uneven worker
budgets. Preserve scRNA/perturb parameters; do not import the FLEX 500-UMI
policy into another assay. Guide assignment and guide-calling models stay
outside this algorithm change.

### Phase 5 — scRNA and perturb caller memory

Introduce validated sparse count views that can describe STAR's existing
strided counts as well as split gene/count arrays. Avoid the adapter's copied
arrays followed by `scrna_api.cpp` reconstructing interleaved counts. Preserve
an owning fallback for callers that require it and public API compatibility.
Carry barcode/feature identity and MT rank metadata alongside the view.

Make index width, strides, lifetime and mutability explicit. Reject overflow
instead of truncating large matrices. Keep original counts available through
all likelihood and diagnostic computations. Measure live matrix bytes and
caller peak RSS as well as whole-process peak RSS; an earlier mapping peak
can hide a successful caller memory reduction.

Gate: exact phase-4 outputs across both owning and borrowed representations,
including empty/sparse/dense cells, gene strides, lifetime/error paths and
large-offset boundary checks. No unnecessary reread of huge MEX files.

### Phase 6 — OCM sample concurrency

Refactor the serial native sample loop in `OcmMultiMaterialize.cpp` into
bounded independent calls. Retain fused tag-to-sample grouping and final
sample order. Each active group reserves progress capacity and borrows spare
caller workers; pool capacity never becomes the sum of per-sample maxima.

Budget memory as well as CPU. Do not load every sample MEX simultaneously;
cap in-flight sample preparations based on their count footprint. Preserve
the established OCM materialization/namespace rules and avoid falling back to
standalone `soloCellFiltering` as a substitute OCM pipeline.

Gate: serial versus scheduled per-sample matrices, called barcodes, ranks,
ambient profiles, p-values and routing exact; test unequal sample sizes,
empty groups, more samples than workers and exception cleanup. Record permits
returned and peak in-flight matrix bytes. Audit other PF calls into the shared
C API individually: feature-only callers may use fixed expected-cell counts
and do not automatically receive an OrdMag-bootstrap speedup.

### Phase 7 — SLAM repeated arithmetic and result reuse

First commit: in `SlamSolver.cpp`, precompute old/new log probabilities for
each occupied `(n, conversions)` histogram entry at the fixed error/conversion
rates. Use the same floating-point expressions and histogram traversal order;
only the mixture weight changes during fitting. Cover boundary probabilities,
empty histograms, non-convergence and non-finite values using the established
behavior. Apply the same invariant-cache principle to `slam_vb_overdisp.cpp`
without changing dispersion, priors or enabling that model by default.

Second commit: calculate fitted results once and pass them to the standard
and optional GRAND-SLAM-format writers in `SlamQuant.cpp`. Reuse results only
for identical histograms, rates and model parameters. Debug histograms and
alternative-rate reports must not accidentally reuse incompatible fits.
Apply the common helper in `slam_requant` where appropriate. Keep output
format, precision and gene order stable.

Gate: exact fits and serialized outputs against the frozen calculation for
unchanged operation order; compare convergence and likelihood as well as NTR.
Test blank/treatment, both solvers, SE/PE-derived histograms, optional outputs
and repeated writing without recomputation. Time solver work separately from
alignment and output. Do not reuse stale caches after a parameter change.

### Phase 8 — SLAM gene parallelism

Run independent gene fits as bounded tasks into disjoint result slots, then
write in original gene order. Immutable fit parameters can be shared; keep
mutable caches worker-owned or prepare immutable tables before scheduling.
Use the available post-mapping budget and propagate exceptions after workers
are joined. A single fit's iteration and summation order must not depend on
which worker receives it.

Gate: same phase-7 outputs with one and multiple execution workers, including
skewed histogram sizes and empty genes; same read counts, conversions, coverage,
NTR, likelihood, cB counts and QC/trim decisions in the integrated SLAM fixture.
Report added cache/result memory; do not imply a FLEX-scale SLAM memory saving.

### Phase 9 — bulk TranscriptVB component execution

The implementation is already in `libem/vb_engine.cpp`; the standalone
`transcriptvb_finalize --component-parallel` sets `EMParams::component_partition`.
Normal STAR currently leaves it false. Add an explicit integration option,
validate it, and promote a default only after representative acceptance.

Give each connected component's transcripts one writer. This eliminates the
`threads × transcripts` accumulation buffer and its per-iteration reduction.
Preserve global convergence, priors, initialization, effective lengths and GC
update behavior. Keep `per_component_convergence` off: it changes how far
components converge and is a separate accuracy/speed tradeoff.

Use identical saved equivalence classes to isolate the engine, followed by
integrated bulk validation. Report actual allocation removal, engine and total
time, iteration count, transcript/gene NumReads and TPM, effective lengths and
tximport outputs. Accumulation order can change floating-point results; report
max absolute/relative errors, affected low-count transcripts and decision/
convergence differences. Existing parity thresholds must be pinned in phase 0,
not loosened after a failure. Older byte-identical examples are not a universal
guarantee. Include few/one-component cases, threads exceeding components and
an empty/no-support input.

Standalone no-index finalization already exists; reuse it for engine timing.
Keep static-range/distributed VB, different convergence rules and parallel GC
reduction changes as later work; their isolated historical speedups do not
establish end-to-end equivalence for this plan.

## Validation and benchmark protocol

### Required surfaces

| Path | Exact checks / numerical checks |
| --- | --- |
| Shared input | Logical record coverage exactly once, mate/name/quality identity, lanes, limits, stream digest, error behavior |
| PF | Feature/barcode/UMI assignment counters, corrected identities, raw and filtered matrices, per-library summaries |
| scRNA / perturb GEX | Raw count coordinates/axes, called-cell identities and all caller stages; preserve guide results in integrated perturb |
| OCM | Per-sample namespaces, grouping, raw/filtered counts, callsets and caller diagnostics |
| SLAM | Histograms, counts/conversions/coverage, fit outputs, trims/QC, standard and requested alternative output formats |
| Bulk | Read/alignment/count behavior, transcript/gene quantification, effective lengths and tximport; documented numerical bands for VB |
| Existing FLEX / CBQ | Relevant small reader/caller fixtures continue to match saved accepted outputs |

Use existing harnesses selectively: `tests/run_bgzf_ingest_tests.sh`,
`tests/run_cbq_pf_adapter_smoke.sh`, `tests/compare_pf_assign_outputs.py`,
`tests/run_pf_dynamic_permit_100k_smoke.sh`, `tests/emptydrops/` tests,
`tests/test_ocm_mex_materializer_tiny.sh`, `tests/run_ocm_materializer_memory_smoke.sh`,
`tests/run_slam_solver_test.sh`, `tests/run_slam_unit_tests.sh`,
`tests/run_slam_fixture_pe_parity.sh`, `tests/run_public_bulk_pe_smoke.sh` and
`tests/run_transcriptvb_chr22_parity_smoke.sh`. Extend these where they lack the
new behavior; an existing harness does not automatically test the new path.

Inspect scripts before launch for nested executions, default output deletion,
missing fixtures and skip-as-success behavior. For example, the old
`run_emptydrops_parity.sh` compares different caller backends and permits skips;
it is not the exact old/new shared-caller acceptance gate. Build a same-count,
same-model comparison instead. Do not run aggregating suites that silently
repeat an already executed fixture.

### Execution discipline

- Follow the host's identical-execution rule: saved results are reused by
  exact binary/input/argument/environment identity. A source fix and clean
  rebuild allow testing the changed build; a new output directory or timestamp
  alone does not authorize repeating the same execution. This runbook does
  not itself authorize repeat runs. Any required identical repeat follows the
  explicit user-authorization rule in `AGENTS.local.md`.
- Serialize benchmarks on a host. Do not time alongside other benchmarks,
  builds, compression, archive uploads or checksum sweeps. Use an isolated
  checkout and fresh output directories; preserve failed attempts.
- Start with bounded fixtures and isolated kernels, then one representative
  larger arm for each accepted change. Avoid a full experiment grid after
  every commit. Reuse checks whose executed binary and affected behavior have
  not changed; broaden validation only for new failures or unresolved risks.
- Prefer current local resources and saved inputs. No cloud reprovisioning,
  full 320K rerun, new release or infrastructure retirement is queued by this
  planning task. Record any later execution scope before launch.
- Log effective native/fallback path, decoder/worker budgets, queue occupancy,
  permit waits, phase wall/CPU times, overall wall/CPU time, peak RSS and input
  counts. Capture memory at caller preparation and fitting, not just process
  peak. Record cache policy, CPU/storage type and thread count.
- Keep conversions and reference preparation outside timed runs and report
  their cost separately. Match reference, input, clipping, statistical options
  and output settings between arms. BAM/sidecars are included only when the
  comparison or mode contract requires them, equally on both arms.
- A successful wrapper completion record and validated final outputs define
  completion. Progress logs or process disappearance alone do not.
- Treat single runs as single measurements, without confidence intervals.
  Use phase attribution/Amdahl's law to explain small total gains. Label
  ordinary gzip, BGZF and CBQ separately; file suffixes are not format proof.

### Evidence layout

Under the proposed local evidence root, retain:

```text
00_inventory/{source.json,inputs.json,references.json,acceptance.json}
execution_registry.jsonl
phase_N/commit.json
phase_N/arm_name/{argv.json,environment.json,stdout.log,stderr.log}
phase_N/arm_name/{BENCHMARK_COMPLETE.json,VALIDATION_COMPLETE.json,metrics.json}
phase_N/RESULTS.md
FINAL_REPORT.md
```

Each registry entry identifies binary SHA-256, source/effective-build inputs,
data/config hashes, environment, output path, execution status and whether the
result was reused. Preserve source-to-binary provenance, numerical comparisons
and failed attempts. Record artifact locations in `tests/ARTIFACTS.md` when
the corresponding phase creates them. Archive new evidence under its own
prefix when storage is selected; do not overwrite the completed FLEX archive.

## Phase 10 — combined acceptance and integration

After individual gates pass, run the selected combined regression once against
the accepted build and compare against the pinned controls. Report one table
per module with input format, exact workload, total and phase time, memory,
thread budget, correctness result and source identity. Separate measured
benefits from unmeasured opportunities; do not sum overlapping phase durations
or transfer FLEX's percentage gains to other modules.

Review and commit each accepted phase with its evidence summary. Keep shared
source history with `--no-ff` integration, complete repository-required CI for
the affected paths, and preserve unrelated workspace work. Coordinate a final
merge/push when that integration is in the execution scope; this runbook does
not specify a new release tag. Unsupported modes and unsuccessful changes
must be listed explicitly rather than marked implemented.

## Execution ledger

| Phase | Status | Commit / evidence / next action |
| --- | --- | --- |
| Runbook and source audit | **COMPLETE** | Current source inspected; PF given first BGZF priority |
| 0 — pin evidence and acceptance | COMPLETE FOR PF; LATER MODULES PENDING | Frozen source, hashed small/full PF inputs and references; pin later-module fixtures before their phases |
| 1 — PF BGZF input | COMPLETE | Full 9,748,584-pair exact parity; 28.2722 → 25.5696 s, RSS +37.8 MiB; see phase-1 handoff |
| 2 — PF direct worker/copy work | NOT STARTED | Depends on phase 1 exact parity |
| 3 — general core BGZF modes | NOT STARTED | Bulk/scRNA, then SE/SLAM and complex modes |
| 4 — scRNA CPU work | NOT STARTED | Separate OrdMag-reuse and MC-threading commits |
| 5 — scRNA matrix views | NOT STARTED | Owning/view parity and memory measurements |
| 6 — OCM sample permits | NOT STARTED | CPU and memory bounded |
| 7 — SLAM cache/result reuse | NOT STARTED | Two separately measured commits |
| 8 — SLAM gene parallelism | NOT STARTED | Same fits and output order |
| 9 — bulk component execution | NOT STARTED | Global convergence and numeric acceptance |
| 10 — combined acceptance/integration | NOT STARTED | Per-module final report and scoped integration |

## Supporting records

- [PF CBQ direct-range design and existing implementation](RUNBOOK_PF_CBQ_RANGE_COUNTING_20260531.md)
- [FLEX BGZF ingest contract](RUNBOOK_BGZF_PARALLEL_INGEST_20260901.md)
- [Shared OrdMag/permit measurements](benchmarks/FLEX_ORDMAG_SORT_PERMITS_L004_20260910.md)
- [Full FLEX 1.9.0 results](benchmarks/FULL320K_V190_20260910.md)
- [Standalone VB finalizer and component experiments](STANDALONE_VBEM_FINALIZER_20260807.md)
- [Public bulk fixture](PUBLIC_BULK_PE_FIXTURE_AND_REGRESSION_TESTS.md)
- [SLAM PE production modes](RUNBOOK_SLAM_PE_PRODUCTION.md)
- [SLAM bounded smoke](RUNBOOK_SLAM_PE_100K_SMOKE.md)
- [Perturb paper methodology](PAPER_BENCHMARK_METHODOLOGY.md)

## Execution update — 2026-09-10

Implementation is isolated at
`/mnt/pikachu/star_suite_paper/analysis/cross_module_performance_20260910/worktree`,
branch `perf/cross-module-bgzf-callers`, based on the pinned master above.
Phase 0 inventory is staged by implementation phase: PF references and both
small/full feature libraries are selected; later modules still require their
own fixture/reference manifests before execution. No new cloud resources are
involved.

Phase 1 deliberately attaches the adapter to the common PF FASTQ producer
ring, shared by the standalone CLI and embedded file APIs. It preserves the
existing feature kernel and leaves batch/direct-worker handoff removal to
phase 2. The C adapter returns borrowed record spans, including names for mate
validation, with leases retained until PF has copied them.

The 200,000-pair A375 check has exact keyed integer counts, barcode/feature
sets and assignment diagnostics: 8,423 barcodes, 11 features, 8,581 occupied
entries, 146,979 UMIs. Whole-process times were 2.0695 s (frozen baseline) and
2.0689 s (native BGZF); peak RSS was 210,688 KiB for both. This is **no measured
speed or memory improvement**. These are PF-only timings, excluding GEX,
reference alignment, cell calling and BGZF conversion. Both consumed the same
BGZF files under an eight-CPU affinity; inflaters and assignment shared six
permits in the native arm, leaving capacity for the two producer threads.

The frozen standalone CLI failed before processing reads with explicit lane
lists (`-a`, `--barcode_fastqs`, `--forward_fastqs`), so that failed arm supplies
no performance result. The successful comparison uses `pf_process_fastqs()`.
The standalone directory-input gates passed. Record the explicit-list CLI
failure separately; do not confuse it with a native BGZF regression.

See [PF BGZF phase-1 handoff](HANDOFF_PF_BGZF_PHASE1_20260910.md) and the evidence
root above for commands, hashes, completed outputs and remaining work. No
identical benchmark or fixture execution was repeated; two corrections to
comparison assertions reused existing output logs.

Full PF follow-up is complete: **9,748,584 paired feature reads**, exact
3,218,159 UMIs, 110,772 barcodes and 124,871 occupied matrix entries, plus
identical assignment diagnostics. Whole-process wall time was 28.2722 s for
the frozen baseline and 25.5696 s for native BGZF (9.6% lower in this single
comparison). RSS increased from 1,554,892 to 1,593,580 KiB (+37.8 MiB).
The reported read-processing phase was 22.25 → 21.26 s; whole-process gains
also include later PF work and should not be assigned entirely to decoding.
This is the **full feature library**, not a whole perturb GEX benchmark.
Phase 2 remains the next implementation step; do not rerun the accepted
baseline merely to establish it again.
