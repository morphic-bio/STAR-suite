# Runbook: hierarchical permits for A375 — completion between workloads, throughput within each workload

Date: 2026-09-11. Status: **implementation in progress; not released.** Base: branch `fix-larry-feature-regression-20260911`
(worktree `/mnt/pikachu/STAR-suite-larry-regression-20260911`), which carries the feature-learning revert
(`PfMultiProcess.cpp:2429`, learning always 100,000 reads). Nothing here is committed or released.

Owner's design clarification, 2026-09-11: MAP and FEATURE are independent workloads to finish. BGZF is a
component of each workload. The outer allocator minimizes the final completion tail; the inner allocator
maximizes its parent's processing rate within that parent's budget. Implement and validate this hierarchy
on A375 first. LARRY's separate cache investigation is not a prerequisite for this work.
The owner also confirmed that the finish-together heuristic already exists. Reuse that outer completion
balancing machinery; the new scheduling work is BGZF allocation within each parent's budget.

Owner's final metric ordering: establish correctness first, then minimize **total end-to-end runtime**,
then minimize the **tail after the first workload finishes**, then consider **average CPU utilization**.
Each workload includes its gather/merge steps. The last input read or the last worker processing a read
is not its finish boundary. Report processing, gather and full-workload timestamps separately; audit
the completion signals and ETA observations against those boundaries before accepting controller results.
For the same total duration, bringing the shorter workload's finish closer to the longer one's reduces
the tail. A reduction in decoder wait time alone does not outrank these objectives.

### Next performance audit, after correctness

The owner's Opus review identifies the exit gap (~27 s), consecutive gather work (~17 s), and
low-CPU startup/genome loading (~58 s) as the remaining A375 opportunities. It reports that permits
already give GEX approximately 98% of CPU during mapping. Retain these as priorities for measurement
before another scheduler change. They are not three proven, independently recoverable time savings.

The saved hierarchy trace supports near-saturation in **steady GEX-only mapping**: about 31.3 of
32 permits in use, or 98% of the permit budget. This is a permit-occupancy measurement; retain a
separate CPU measurement rather than silently equating an occupied permit with on-core CPU time.
Whole-job process CPU in the three fixed-four Step 4 controls averages 14.45–14.47 cores (~45.2%
of 32), including initialization, output and teardown. These denominators describe different phases.

Audit in the review's order:

1. **Exit gap:** measure from the successful-finish message to process exit and attribute object
   destruction/freeing and any remaining work. The approximately 27 s gap is visible in saved full
   wrapper timings; the wrapper's terminal exit remains the authoritative end-to-end boundary.
   In `a375_outer_fixed4`, sampled STAR CPU grows by about 26.8 CPU seconds during the final
   ~26.8 wall seconds while RSS declines. This supports a serial cleanup hypothesis, but does not
   identify the responsible destructor or justify bypassing cleanup.
2. **Gather work:** identify the two gathers behind the reported ~17 s, timestamp their dependencies,
   and determine which work can overlap. PF's fixed-four phase trace measures 8.60 s in thread-hash
   merge/cleanup and 2.92 s in UMI deduplication. They are included in its 24.89 s feature API phase.
   This feature work overlaps GEX mapping, so those durations cannot simply be subtracted from total
   wall. The saved log separately places Solo counting at 16:53:51–16:53:56 (~5 s), followed by
   pf-multi merge/finalize and guide calling at 16:53:56–16:54:09 (~13 s). That ~18 s consecutive
   post-map interval is consistent with the review's ~17 s. The filtered merge and guide calls consume
   GEX counts/cell calls, so establish dependencies before proposing overlap. Include both intervals
   when defining full output completion.
3. **Startup and genome loading:** attribute the reported ~58 s across initialization, whitelist/
   feature preparation, input counting and the index read. In `a375_outer_fixed4`, job start is
   16:51:28, genome loading starts 16:52:14, finishes 16:52:21, and mapping starts 16:52:22.
   That is ~54 s before mapping, with ~7 s inside the explicit genome-load interval. The broader
   startup period, rather than an assumed 58 s index read, is the useful investigation scope.

### Exit-table isolation (follow-up, 2026-09-11)

After the success message, STAR removes temporary files, closes input streams/joins readers, cleans up
its parameter registry and log streams, then destroys local/member objects as `main` returns. Mapping,
Solo/feature output and the RAchunk read state have already finished and been released in these A375 runs.
The main Parameters object still owns `CbCorrector`. At the time of the isolation probe its exact, H1
and ambiguous tables used `std::unordered_map` containers. Other barcode lookup tables already used
khash; this separate path had not received that conversion.

An isolated probe linked to the same clean `CbCorrector.o`, with the actual A375 whitelist, measured:

| Table/component | Entries | Construction/destruction wall |
| --- | ---: | ---: |
| Construct correction tables | 3,686,400 whitelist barcodes | 43.91 s |
| Destroy ambiguity map and its vectors | 42,938,880 keys | 4.41 s |
| Destroy H1 map | 88,270,080 keys | 7.47 s |
| Destroy exact map | 3,686,400 keys | 0.056 s |
| Destroy both whitelist string copies | 7,372,800 strings | 0.061 s |

The table uses the arm that creates/joins one thread before destruction, disabling glibc's single-thread
allocator shortcut. A fully single-threaded arm takes 11.43 s total destruction; the threaded-history arm
takes 12.00 s. Thus allocator single-thread mode alone does not explain the whole STAR exit gap. These
are isolated component measurements, not a breakdown of the actual 27–28 s full-process exit. Remaining
full-exit attribution requires phase/stack evidence.

Construction at ~44 s explains most of the measured pre-genome startup. The owner subsequently authorized
the khash conversion and contiguous ambiguity array together. Both are now implemented: exact/H1/range
tables use khash, and two construction passes populate one candidate array in the original order.
ASan/UBSan and preserved-implementation tests pass, including packed CBQ/FASTQ ordering, N expansion,
candidate limits and shared readers. No scheduler or correction policy was changed in this follow-up.

The new component probe measures 26.44 s construction, 0.163 s destruction and 3.31 GiB peak RSS,
versus 43.91 s, 12.00 s and 7.61 GiB. Full A375 confirms the benefit:

| Matched configuration | Total before → after | Peak RSS before → after | Exit gap before → after |
| --- | ---: | ---: | ---: |
| Fixed-four decoding | 192.10 → 145.47 s | 44.20 → 39.82 GiB | ~28.39 → ~1.10 s |
| Hierarchical balanced, four workers | 190.60 → 146.21 s | 44.20 → 39.80 GiB | ~28.45 → ~0.61 s |

Each new arm matches its corrected control exactly across all six MEX matrices, all three guide CSVs
and the biological feature totals. The two new arms also match each other. The changes recover almost
the full exit gap and 18–19 s of startup; mapping remains around 89–91 s. Feature API time increases by
about 2 s while overlapping GEX. These measurements support the total-runtime improvement, not a new
scheduler claim. The remaining exit is roughly a second, including log rounding and wrapper polling.
Full phase/CPU results and provenance:
[A375 storage benchmark](benchmarks/A375_CB_CORRECTOR_FLAT_HASH_20260911.md).

Artifacts: `/home/lhhung/pf_larry_regression_20260911/hierarchy/exit_cleanup_probe/` and its
`threaded_allocator/` arm retain diagnostic source, copied object/header/source, hashes, build commands,
execution manifests, timings and `summary.json`. No STAR production source was modified for these probes.
The subsequent production change and validations are under `hierarchy/cbcorrector_flat_hash/` and
`hierarchy/a375_flat_hash_{fixed4,balanced4}/` beneath the same artifact parent.

The current MAP permit-completion signal follows mapping-worker joins, before Solo gather/counting;
FEATURE's signal follows its assignment/gather/dedup/MEX phase. Audit full-phase finish separately
from releasing idle permits. Preserve these distinctions in tail calculations and ETA validation.

## 1. What the A375 tests established

Full A375 in-process (GEX 47,095,183 pairs + CRISPR 9,748,584 pairs, 32 threads). The native-BGZF/floor and
zcat experiments used the frozen v1.9.2 binary, with learning reported as
`STAR_PF_FEATURE_BOOTSTRAP_READS=100000`. The exact execution environment was not captured. Saved bootstrap
messages establish that learning was active, but do not independently establish that variable's value.
Do not assume `STAR_SOLO_NONFLEX_HASH_BRIDGE=1` was set on every historical run: two runs failed afterward
in GEX counting. The seven archived configurations in the reference manifest each report 9,274,336 feature
counts (March and the earlier patched worktree run: 9,274,337), with 349,645 unmatched reads.

| Configuration | Feature step (API timer) | Feature wait / work | Mapping | Total wall |
| --- | ---: | ---: | ---: | ---: |
| GEX native BGZF, permits on (P01 settings) | 64.2 s | 1,531 s / 40 s | 89 s | 173.69 s, **exit 1; excluded** |
| same, permits off | 48.6 s | - | - | 173.12 s, **exit 1; excluded** |
| GEX zcat, feature BGZF off, permits on | 20.8 s | 34 s / 26 s | - | 196 s |
| GEX zcat, feature BGZF off, permits off | 35.9 s | - | 96 s | 197 s |
| GEX 100,000 pairs, permits on | 10.6 s | 0.3 s / 22.5 s | - | - |
| GEX native BGZF, fixed feature floor 8 + FIFO | 23.4 s | 230 s / 39 s | 88 s | 188 s |
| GEX native BGZF, fixed feature floor 16 + FIFO | 27.7 s | 235 s / 39 s | 89 s | 190 s |
| Standalone process_features, same settings | 11.3-15.8 s | - | - | - |
| March 26 in-process (`a70a039-dirty`, GEX zcat) | 11.9 s | 0 blocked acquires | 86 s | 233 s |

Use the **working March baseline** to explain changes: feature API **11.8753 s**, mapping **86 s**.
The owner clarified that **end-to-end completion time is the objective**: a longer feature phase is
acceptable if allocating more capacity to GEX improves overall completion. The 11.9 s feature time is
a diagnostic reference, not an independent acceptance gate. Floor-8 is also an interim reference.
March's `BENCHMARK_SUMMARY.txt` records 233 s end to end. Its `Log.final.out` records mapping start
13:45:45; `Log.out` records feature finish 13:45:49 and mapping finish 13:47:11. Thus only 4 s of mapping
overlapped the feature-assignment phase, followed by an 82 s GEX-only tail. These are historical timing
comparisons across builds and input transports; require matched controls for causal attribution.

Interpretation and scope:
- The A375 results justify tackling BGZF permit sharing now. The feature step is much faster with a trivial
  GEX library and with the alternative input configuration; the floor-8/FIFO combination substantially
  improves the native-BGZF feature step without an observed mapping penalty in these rows.
- STAR's native BGZF decoders take MAP permits (`Parameters_openReadsFiles.cpp:455`); PF decoders and PF
  consumers take FEATURE permits (`PfMultiAssign.cpp:23,449-451`,
  `core/features/process_features/src/pf_bgzf_input.cpp:47-61`). They currently
  have no separate child-level admission policy within those domains.
- The table motivates the hierarchy. It does not separately quantify FIFO versus floors, or prove the
  cause of the floor-16 plateau. Further FIFO-only experiments are not a prerequisite for the owner's
  chosen hierarchy work: keep FIFO and parent settings fixed when testing the child allocator.
- Feature completion precedes mapping on A375. Its feature-step improvement is useful evidence for the
  inner allocator even when end-to-end wall time is unchanged. Do not promise an A375 total-wall speedup.
- The two `a375_inproc_bootstrap_permits_{on,off}` runs completed feature assignment but `time.txt` records
  exit status 1, their stdout lacks a successful finish, and neither has `Log.final.out`. The permits-on
  log reports `collapseUMIall_fromHash called but inlineHash_ is null`. Their feature-phase measurements
  remain diagnostic evidence; their total durations cannot be compared to a successful full pipeline or
  used to establish a 174–192-second normal timing-noise range.
- The zcat rows move GEX decompression CPU outside the permit pool, and also disable PF native BGZF. They
  are transport/accounting references, not full-budget performance targets or a pure one-variable control.
- The 100K-GEX run's 10.6416-second feature phase is an approximate low-contention reference for what the
  inner allocator might approach. Its much smaller GEX workload prevents treating it as a full-GEX target.
- Exact commands, artifact hashes, startup messages and completion evidence for seven configurations are
  now collected in [the reference manifest](ADAPTIVE_PERMIT_A375_REFERENCE_RUNS_20260911.json). All seven
  report the same assigned-read count despite different scheduling/reader settings. This supports strict
  same-build scheduler parity; it does not prove the cause of the one-read worktree difference.

### Recovered floor-8 reference

- Run and command source: `/storage/pf_regression_20260911/a375_inproc_featfloor8_fifo/Log.out`.
  This is the P01 command plus `--dynamicThreadFeatureFloor 8 --dynamicThreadFifoWaiters 1`, with both
  `--readFilesBgzfMode auto` and `--crAssignBgzfMode auto`. PF controller mode was the default `off`; ATAC
  controller was the default `0`. Copy the full captured argv from the manifest when building new arms.
- Logged executable: `/storage/paper_bench_v192_20260911/stage/STAR`, release `68eb92ff3cc107179d32ed52d95229b121debae9`.
  Its SHA256 measured during this review is `def647c0b027e7d60c661a2b001001e2a3a2516610ce10619669abef60495013`;
  this is a present-file measurement, not a separately captured historical execution hash.
- Feature API time: **23.3765 s**; total wall: **187.85 s**. Successful finish is recorded in `stdout.txt`
  and `Log.final.out` exists. The saved timing file has no explicit successful exit-code field.
- GEX startup records **16 inflater workers per mate, 32 total**, with `compute_budget=32`.
  PF startup in `time.txt` records **32 inflater workers, 2 producers, 31 consumers, 1 search thread**.
  The API record contains `bgzfReaderThreads=32`; distinguish a configured count from actual startup.
  `time.txt` also contains bootstrap finalization/replay messages, so learning is evidenced beyond counts.
- This run used the older binary plus reported environment override. It is not the hierarchy-off arm for
  a new worktree build. The pre-hierarchy worktree binary currently hashes to
  `c0e399de851638f6f3eeab6f176fe95b393dee19a0503636aa16d651a722d941`; implementation will change that hash.

## 2. What already exists

- **Permit pool**: one pool of `runThreadN` permits with domains MAP=0, FEATURE=1, ATAC=2
  (`ThreadControl.h:18`); borrowable per-domain floors (`mapPermitConfigureDomainFloors`, `ThreadControl.h:126`),
  FIFO waiters (`:143`), domain completion (`mapPermitMarkDomainComplete`, `:146`), per-domain wait/work/unit
  counters (`:241-259`). Floors default to 0, which leaves the legacy notify_one path that
  `parametersDefault` itself describes as a wakeup-fairness pathology under contention.
- **Work-balancing controller**: `SaturationPermitController.h`. Probes each domain's sustained occupancy
  (largest work estimate first), keeps it as a borrowable floor, and when capacity is short moves one floor
  permit per tick from the domain with the smallest ETA to the one with the largest (10% hysteresis), never
  above a domain's measured saturation. This is the "floors move with remaining work" behaviour.
- **But it only runs for multiome**: configured in `mapThreadsSpawn.cpp:946` only when
  `P.chromapAtac.enabled == 1 && P.dynamicThreadAtacController == 2`; its periodic observe/apply loop is the
  sampler thread in `star_chromap_orchestration.cpp:798-824` (observe ~1307, apply floors ~1339).
  GEX + feature runs (A375, MSK) never get it.
- **Work estimates**: controller takes `--dynamicThread{Map,Feature,Atac}WorkEstimate`. P01 recorded
  `feature.workUnits=9,806,352`, but these existing counters combine processed reads with inflated BGZF
  blocks (`BgzfRangeReader.cpp:382-386`). They cannot be used directly as completed-read counters.
  STAR already estimates reads per library for GEX and each feature
  library (`PfMultiProcess.cpp`, `estimateReadsAnchored`, ~882-944, logged as `pf-dynamic-controller:
  estimator[map|feature]`) but does not pass them to the controller.
- **Tests**: `tests/test_thread_control_permit_instrumentation.cpp`, `tests/dynamic_threads/`,
  `tests/run_dynamic_threads_tiny_fixture.sh`, `tests/run_dynamic_threads_variable_sequences_smoke.sh`,
  `tests/run_pf_dynamic_permit_100k_smoke.sh`. Exit invariant: `STAR.cpp:3580-3605`.

## 3. Two levels, two objectives

```text
Global compute budget
  MAP workload budget                 FEATURE workload budget
    BGZF decode + mapping               BGZF decode + feature assignment
```

**Outer objective: minimize the last workload's completion time.** For a workload i, let R_i be its remaining
read pairs and q_i(B_i) its observed completed-read rate at parent budget B_i. Its current drain-time estimate
is R_i / q_i(B_i), in seconds. Allocate the global budget to reduce the maximum drain time. Read counts alone
are not comparable computational costs, and equal completion times are not mandatory when a workload is
serial, input-limited or already saturated. BGZF has no separate outer ETA, floor or completion target.
The practical heuristic is to bring MAP and FEATURE's predicted finish times together. This is existing
policy (`SaturationPermitController::observe` and the PF ETA controller), not a new estimator to invent.

**Inner objective: maximize q_i(B_i).** Within a fixed parent budget, allocate capacity between decoding and
processing to maximize mapped/assigned read pairs per wall second. Decompressed bytes, decoder utilization
and queue occupancy are diagnostic signals, not the success metric. Finishing decompression early has no
independent value. The decoder is valuable only insofar as it feeds useful processing.

The two controllers must not respond to the same transient at the same time. Let the inner split settle
after a parent-budget change before the outer controller interprets the resulting processing rate. Do not
freeze a saturation estimate measured during bootstrap, reference preparation or a library transition.

### Budget and floor contract

- One live decoder consumes one permit from its parent, counted once against the global pool. It must not
  hold a parent permit while waiting to acquire another global permit.
- Parent floors protect the whole workload: decoding plus processing. A child does not gain an extra
  top-level reservation or more top-level priority by spawning additional inflater threads.
- Distinguish a parent floor from its current usable budget. A floor is a borrowable minimum; the outer
  allocator may lend idle capacity to a parent. Children can use that granted capacity, but cannot bypass
  the outer allocator to borrow independently. Sum of active parent permits never exceeds the global pool.
- Return loans at permit-release boundaries when the owning workload needs them; do not preempt a running
  block or strand a decoder needed to produce the next input buffer.
- Actual workload completion releases the parent's allocation. An estimate reaching zero, an empty queue,
  decoder EOF, or a pause between feature libraries is not completion. FEATURE publishes assignment-phase
  completion at `PfMultiProcess.cpp:3548`; MAP publishes mapping-worker completion at
  `mapThreadsSpawn.cpp:1042`. These release unused processing capacity. They do not both mark full
  workload completion: MAP still has Solo gather/counting afterward. Keep a separate full-workload
  boundary for the owner's tail metric; do not label worker completion as final cell/count output.
- Thread supply and permit admission are separate. Start the configured complement of decoder and processing
  workers when each reader/pipeline opens; keep excess workers blocked on demand/admission. Initial scope
  does not add threads dynamically. The admitted concurrency is bounded by both the parent's budget and
  live workers of that kind. Log configured and actually started counts, per mate/lane and in total, so a
  thread-supply ceiling cannot be mistaken for a controller failure or a measured CPU saturation point.

### Inner control signals and admission

- Measure completed read pairs, processing input-starvation time, decoder/processor permit occupancy and
  wait time, ordered input ready depth, in-flight decode work, and time blocked by full output buffers.
- A persistently low ready queue with processors waiting for input is evidence that decode capacity may
  need to increase. A full queue with ready processing work is evidence that decode capacity can decrease.
  Confirm a change using the parent's completed-read rate at the same budget.
- Distinguish missing input/I/O delay from insufficient inflate capacity; more decoder permits cannot repair
  a missing file, serialization bottleneck or slow storage by themselves.
- Account for ordered delivery: later ready blocks do not feed a consumer waiting for an earlier block.
  Preserve progress for the required block and both mates. Release permits before waiting for buffer space.
- Use bounded queues with low/high watermarks and hysteresis, plus throughput feedback. A permanent
  one-decoder cap or unconditional processor priority is not the design: either can make decoding limit the
  parent's throughput. Decoder concurrency may grow within the budget when it improves production.
- Use short inner observations and slower outer adjustments. Log the selected intervals and settling rule.
  The existing outer controller probes one workload at a time, assigning it most of the budget for two valid
  windows. Two workloads at ten seconds per window need at least forty seconds, already longer than the
  floor-8 A375 feature phase. Budget outer probing as well as inner settling for a short run; choose shorter
  outer windows or a bounded initial probe seeded from existing measurements, and record probe time spent.

## 4. Implementation order and gates

Make one scheduler change at a time. Clean-build before diagnosing a regression, instrument first, and use
timings as validation. Keep the learning size and assignment/calling rules fixed during scheduler comparisons. The singleton
correctness investigation below exposed separate buffer and learning-sample defects; validate their fixes
independently, then use the same corrected build in each scheduler arm.

**Step 1 — Introduce the parent/child accounting contract and diagnostics.** Separate processed read pairs
from inflated blocks/bytes, and split decode versus processing permit observations under MAP and FEATURE.
Record parent budgets/floors/loans and queue demand. Preserve existing admission while checking accounting.
Gate: balanced acquisitions/releases, exact global budget accounting, unchanged biological outputs. Cover
plain gzip and BGZF; gzip has no native inflate child. Record reader mode and the decoder/producer/consumer
threads actually started, including PF's existing startup log, with the corresponding live admission limits.
Record CPU used by any zcat/decompression subprocess separately from STAR and permit-accounted CPU; capture
it through the benchmark wrapper's child-process accounting even if the subprocess finishes before STAR.
Do not silently equate zcat CPU with zero decoding work. These measurements are inputs and evidence for the
hierarchy, not a standalone telemetry cleanup.

**Step 2 — Implement inner allocation at controlled parent settings.** Opt-in hierarchy first (proposed
`--dynamicThreadBgzfHierarchy 1`). Hold parent settings and FIFO identical between off/on controls. Implement
decode/processing admission under the parent's usable budget with bounded prefetch and throughput feedback.
Validate with synthetic decode-limited and processing-limited pipelines, ordered-mate delivery, a one-permit
budget, shrinking budgets, completion and cancellation. Demonstrate that capacity moves toward the limiting
stage, rather than merely reducing decoder thread count or its reported wait.

**Step 3 — Validate the inner hierarchy on full A375.** Use native BGZF for both arms and the same build,
inputs, learning policy, FIFO and parent settings for off/on controls. Use the recovered floor-8 argv above.
For the inner-allocation experiment, make the parent allocation explicit with `--dynamicThreadMapFloor 24`
and `--dynamicThreadFeatureFloor 8` in **both** arms. The historical floor-8 run left MAP's floor at zero;
it remains a reference. Holding 24/8 constant in the new matched pair isolates the child policy from a
change in parent floors and gives the inner controller a defined initial workload allocation.
Run one **new hierarchy-off floor-8 control on the clean implementation build**, then one hierarchy-on arm
on that identical binary; only the hierarchy flag changes. Neither historical v1.9.2 nor the older
`c0e399de...` worktree binary can substitute for the same-build off control. Recompute and save the actual
binary hash used for both arms. This is a new configuration/build comparison, not an identical rerun of
the archived floor-8 execution.

Assess processing rate, feature-step time, mapping time, total wall, memory and queue/permit traces. The
recorded 23.4-second feature step is a reference, not an exact cross-session cutoff. Require strict equality
of feature/GEX matrices and cell/guide calls between these same-build arms, zero leaked permits and no
deadlock. Investigate a greater-than-5% mapping or total-wall regression against the matched control.
The 5% rule is an investigation trigger, not a statistically established noise bound or an automatic
rejection. Inspect traces, phase boundaries and completed work before attributing a difference to the
scheduler. Do not schedule a duplicate control or confirmation run; an absolutely identical execution
requires the owner's explicit repeat authorization. A changed implementation may receive a new validation.
Owner's later acceptance clarification: a traced difference caused solely by an equally supported
ambiguous tie is acceptable. Establish identical input evidence and identify the tie before applying
this exception; unexplained singleton differences, lost reads, overlap or incorrect merging remain
failures to investigate. Preserve the strict comparison report even if an explained tie is accepted.
Use phase time and input starvation to locate contention, without requiring a shorter feature phase.
End-to-end time determines whether the allocation is useful. LARRY testing does not gate this A375 work.

Owner's added control: before judging adaptive allocation, use **four actual BGZF workers per parent**
with `--bgzfReaderThreads 4 --dynamicThreadBgzfHierarchy 0`, retaining MAP24/FEATURE8, FIFO and total32.
This parameter gives GEX two workers per mate (four total) and PF four across its lanes. Confirm startup
counts, rather than treating a permit cap as an actual thread count. Compare this fixed-worker control
with the adaptive hierarchy on the same frozen binary and biological inputs. A slower FEATURE phase
is acceptable if total wall improves. Keep source provenance pinned to the binary's build, even while
later controller implementation is in progress in the worktree.

Owner's subsequent clarification: four workers may already be near the optimum for A375. Remaining
competitive with that choice is a useful result; an adaptive policy need not beat a well-chosen fixed
allocation on every workload. The earlier FLEX permit result motivates the inner allocator, but does
not establish an A375 speedup. Record observed differences without calling a single comparison a
statistical equivalence test. In particular, the matched demand-correction build completed in
189.36 s with fixed four versus 190.40 s with adaptive admission and 32 available decoder workers,
with identical matrices and guide calls. This is about 0.5% slower for adaptive, not evidence that
adaptation is harmful or that four is optimal for other workloads.

**Step 4 — Connect the hierarchy to existing outer completion balancing.** After the inner allocator works,
retain the existing finish-together heuristic and adapt its integration where needed behind the proposed
`--dynamicThreadBalance 1`. This step is integration and accounting, not replacement of the outer policy.
Give it explicit active workloads: the current controller always probes ATAC, so MAP+FEATURE requires a real
active-domain mask. The current periodic sampler is owned by `StarChromapAtacAsyncRun::Impl`; create/start
the shared sampler on the GEX+FEATURE path as well, with explicit start, stop and join ownership.
Supply estimated total GEX pairs and the sum over feature libraries, with separate completed-pair counters
and explicit completion signals. The sampler (`star_chromap_orchestration.cpp:1225-1245`) currently computes
ETA as remaining work divided by an EWMA of observed work units per wall second, then passes that ETA to
the controller's `eta()` accessor. **The current formula mixes units:** total estimates are reads, while
both MAP and FEATURE's native-BGZF `workUnitsTotal` and rate include reads plus inflated blocks. Correct the
progress counter and rate together. This can distort remaining work and ETA; equality can report zero
prematurely. When the mixed counter exceeds the estimate, `updateEstimate` at lines 1177-1190 grows the
estimate, so an early zero ETA is not the inevitable outcome. Do not infer completion from either behavior.
A reads-per-permit-second rate describes efficiency and must incorporate allocation before becoming a
drain time in seconds. No replacement of the existing finish-together objective is needed.

Keep a single owner for parent budgets/floors. Reject either new hierarchy/balance mode during parameter
parsing unless `--dynamicThreadPfControllerMode off`; the PF controller changes the global permit target
at `PfMultiProcess.cpp:3104` and resets it at lines 3140/3151. Test the conflict rejection before launching
workers. For the initial GEX+FEATURE mode also require the ATAC controller off; preserve its legacy mode
when the new flags are off. Existing finish-together policy may be reused by the shared controller, but
the old PF retune thread must not run alongside it.
Control must function independently of whether verbose telemetry is enabled. The inner allocator reacts to
the new budget first; the outer allocator evaluates settled parent throughput. Preserve multiome behavior
when the new mode is off. Cover MAP-only, FEATURE-only, MAP+FEATURE, delayed work, sequential libraries and
completion before adding ATAC to the new mode.

**Step 5 — Validate the complete hierarchy on A375.** Compare outer balancing off/on with the inner policy
held fixed. Primary outcome is total wall and the remaining completion tail; feature time and per-workload
rates explain the allocation. Do not require the outer policy to reproduce a fixed feature-stage time if
it improves the overall finish. Require unchanged outputs and permit invariants as in Step 3.
Use three arms on the same new build, with four decoder workers per parent throughout: hierarchy off
and balance off; hierarchy on and balance off; hierarchy on and balance on. This retains the owner's
fixed-four control and separates inner admission from outer balancing.

## 5. Then LARRY

1. Build a MSK in-process downsample: GEX (plain gzip, zcat reader as in the ES wrapper) cut to a few million
   pairs, PolyIII and LARRY cut to ~2 million reads each, full feature references. Keep the learning revert
   fixed; distinguish inner hierarchy and outer balancing in comparisons. Compare feature step and total
   wall with standalone process_features on the same LARRY subset. A375 establishes the BGZF hierarchy;
   LARRY separately exercises the large feature-reference path and plain-gzip inputs. GEX zcat still consumes
   CPU outside the pool: retain separate subprocess CPU accounting and report the complete process-tree CPU
   budget. Do not present that arm as validation of a fully permit-accounted decoder hierarchy.
2. Full P02 (MSK ES), **learning/search-regression acceptance**: LARRY assignment within 30% of 1,442.9 s,
   PolyIII within 30% of 146.4 s, counts and calls matching April or every difference explained. These bounds
   primarily assess recovery from the learning/search defect. Attribute permit improvements through the
   matched scheduler controls and traces, not these historical bounds. Paper timings wait for a released build.

## 6. A375 test command

Use the recovered floor-8 command from `/storage/pf_regression_20260911/a375_inproc_featfloor8_fifo/Log.out`
(`##### Command Line`, also captured as argv in the reference manifest), with a
fresh `--outFileNamePrefix`, a copy of `multi_config.csv` for `--pfMultiConfig`, the worktree `STAR` binary
(learning revert built in, so no bootstrap environment variable), and:

```
export STAR_SOLO_NONFLEX_HASH_BRIDGE=1      # the A375 wrapper sets this; without it GEX counting fails
... --readFilesBgzfMode auto --crAssignBgzfMode auto \
    --dynamicThreadInterface 1 --dynamicThreadConstMapPermits 32 --dynamicThreadTelemetry 1 \
    --dynamicThreadFifoWaiters 1 --dynamicThreadMapFloor 24 --dynamicThreadFeatureFloor 8 \
    --dynamicThreadPfControllerMode off --dynamicThreadAtacController 0 \
    --dynamicThreadBgzfHierarchy 0
```

The hierarchy flag above is implemented. The command describes the new off control;
Step 3 changes only that flag to `1` for the on arm. Both use the same clean implementation build and built-in
learning revert. Capture the relevant environment explicitly; omit/unset `STAR_PF_FEATURE_BOOTSTRAP_READS`
in both arms so an inherited override cannot change their learning policy. Step 5 keeps the
hierarchy enabled and additionally varies the implemented `--dynamicThreadBalance` flag. Do not silently turn
off native BGZF or give its decoders unaccounted CPU capacity to obtain an improvement.

Before every new run, save a dev provenance manifest with actual binary SHA256, source revision/dirty patch,
complete argv and relevant environment, input/reference identities, host/thread budget, output path and
reason for the run. Save exit status and completion evidence afterward, including failed/interrupted runs.
Register artifact locations in `tests/ARTIFACTS.md`; retain a durable manifest for later central ingestion.

Record for every run: `stats.processing_time_sec` and `dynamicPermitDelta.{map,feature}.{waitNs,workNs}` from
`cr_assign/.../assignBarcodes.api_run.txt`; `started/finished mapping` and `started/finished pf-multi feature
assignment` (STAR stdout and `Log.out`); total wall (`/usr/bin/time`); `stats.txt` feature counts and unmatched;
GEX cells; parent and child decision logs; completed-pair rates, decode blocks/bytes, queue depth and input
starvation; parent budgets/floors/loans; per-child permit accounting; peak RSS; output comparison; exit status.
Include actual thread counts, reader modes, subprocess and whole-process-tree CPU, outer probing time and
inner settling windows. The historical manifest includes the floor-8, 100K-GEX and zcat reference bounds,
and labels failed native-BGZF totals as ineligible for full-run comparison.
Use matched controls for decisions and retain section 1 as historical context.

## 7. Constraints

- One run at a time; fresh output directories; clean builds for any regression diagnosis.
- Do not change the default permit behaviour until the owner accepts the A375 and LARRY results; keep the
  balancer opt-in until then.
- Clean room: never inspect Cell Ranger or cyto source.
- No commit, merge, push or release without the owner.
- Separate defects stay separate: the learning revert (counts, LARRY search cost) and the missing guide-call
  rows for EmptyDrops-rescued cells (`docs/HANDOFF_CRISPR_CALLS_MISSING_TAIL_CELLS_20260911.md` in the main
  checkout) are not part of this work.


## 8. Correctness findings during scheduler validation

The first Step 4 fixed-four controls completed in 189.36 s (hierarchy off), 191.59 s (inner only),
and 190.86 s (inner plus outer). Those small differences do not establish an adaptive speedup. GEX,
filtered feature matrices and guide-call tables agree, but raw feature matrices differ by one/two
singleton UMIs. They failed the strict parity gate and are retained as diagnostic evidence.

Two separate defects were isolated in process_features, before merge/gather:

1. **Zero-N fallback buffer corruption.** With `max_feature_n=0`, the old expansion-size expression
   shifts by a negative exponent and can allocate a zero-length buffer. Actual A375 reads showed their
   input sequences modified while recording a match position. The expression predates the scheduler;
   restoring the broad learning fallback exposed it. Give the zero-N case one expansion slot and check
   the N limit before writing the indices array. The actual two-read regression fails on the prior clean
   library and passes on the corrected clean library; ASan/UBSan and existing anchor tests pass.
2. **Lane-availability-dependent learning cutoff.** The existing learning phase has one consumer but
   takes whichever input lane is ready. This changes which reads get the broad learning search versus
   the later learned-position window. The two singleton reads match the short seven-base IL1B feature
   at positions 54/70; the learned position is 79 in both source-order controls, with a unique maximum.
   This is not an equally supported ambiguous tie. In a controlled 200,000-pair, two-lane fixture,
   reversing a two-second producer startup delay changes 215 raw matrix entries on the zero-N-fixed
   library. Waiting for strict round-robin lane order during learning gives identical histograms and
   matrices under both delays (147,638 UMIs each). After learning, normal ready-lane consumption resumes.
   The learning size remains 100,000; search/acceptance rules are unchanged.

Artifact root: `/home/lhhung/pf_larry_regression_20260911/hierarchy/`.
`umi_trace/` contains producer/consumer and pre/post-merge ledgers; `umi_trace_modes/` isolates input
ordering and records the failing UBSan case. `zero_n_fix/` contains the clean zero-N-only build and
regressions; `learning_order_fix/` contains both delayed-producer controls, source identities, commands,
input hashes and exact matrix comparisons. Diagnostics are distinct from timing measurements.
Full A375 validation of the two correctness fixes **passes** on clean binary
`ce9d8142f6e8f0f52d4755b669265166f4a08e673c65650be362e8ffe2577a3c`.
The three fresh scheduler configurations match exactly across **six matrices and three guide-call CSVs**,
including the previously differing raw features. All feature totals also agree: 9,275,612 assigned counts,
3,221,846 UMIs, 109,746 raw barcodes and 348,305 unmatched reads. Permit exit invariants pass.

| Corrected build, four decoder workers per parent | Total wall | Mapping | Feature API | Average CPU cores | Peak RSS |
| --- | ---: | ---: | ---: | ---: | ---: |
| Fixed four, hierarchy off | 192.10 s | 91 s | 24.48 s | 14.33 | 44.20 GiB |
| Inner hierarchy | 190.62 s | 89 s | 23.66 s | 14.41 | 44.19 GiB |
| Inner plus outer | 190.60 s | 89 s | 24.69 s | 14.43 | 44.20 GiB |

Adaptive settings are ~0.8% faster in these single matched observations. They remain competitive with
fixed four; this does not establish a repeatable speedup. Total wall includes all initialization, gathers,
output and exit. Tail from feature assignment/gather completion to Solo completion is 66/70/67 seconds,
respectively; the subsequent shared merge/guide phase is another 12–13 seconds. These full-phase timings
supersede read-worker-only tail measurements for the owner's metric.

Compared with the older Step 4 build, the correctness fixes change 421 raw feature entries, net +36 UMIs:
IL1B_sg2_HEK +269, RAB1A-2_MS -230, Non_Target-1_MS -3. Filtered feature matrices change at 329 entries
(net -28 UMIs). Both GEX matrices and all GEX cell identities are identical. Guide-call identities and
thresholds are unchanged; 178 per-cell rows have updated UMI counts. Do not describe this cross-build
change as a single raw singleton or unchanged feature matrices. The exact cross-build report is
`a375_correctness_cross_build_impact.json`; same-build scheduler equality is recorded separately in
`a375_corrected_{inner4,balanced4}_comparison.json` and `a375_corrected_summary.json`.

One outer-estimator limitation remains for the next phase audit: after FEATURE reads finish, its read
rate falls to zero during gather. The approximate-total/read-rate ETA then spikes (69.29 s in the balanced
trace), prompting a floor change from MAP28/FEATURE4 to MAP27/FEATURE5. Unused permits remain borrowable,
so this is not evidence of idle MAP capacity or a measured runtime regression. It does show why a
read-drain estimate is not a full-workload ETA. Account for the gather phase explicitly before claiming
that the outer estimator predicts complete workload finish times.

[Correctness handoff](HANDOFF_A375_HIERARCHY_CORRECTNESS_20260911.md) records the source changes and evidence.
No further performance change, commit or release has been made. LARRY remains the subsequent validation
stage; this A375 correctness result does not complete the full runbook.
