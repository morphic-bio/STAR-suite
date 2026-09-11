# A375 hierarchy and singleton correctness — 2026-09-11

Development worktree: `/mnt/pikachu/STAR-suite-larry-regression-20260911`, branch
`fix-larry-feature-regression-20260911`, base `68eb92ff3cc107179d32ed52d95229b121debae9`.
Changes are committed on the development branch and remain unreleased. No Cell Ranger/cyto source was inspected.
The controlling plan is [the reviewed runbook](RUNBOOK_ADAPTIVE_PERMIT_BALANCING_A375_20260911.md).

## Intended behavior and implementation

MAP and FEATURE are the outer workloads. The existing finish-together controller transfers parent
floors; each parent's BGZF admission optimizes completed parent read pairs per second. A decoder uses
one permit from its parent/global budget, not an additional independent allocation. Typed accounting
separates read pairs from inflated blocks. The inner policy tests one permit at a time, waits for
settling, measures parent throughput and rolls back harmful changes. Empty input with pending decode
work protects a borrowable decoder share; full inputs allow processors to use the capacity.

`--dynamicThreadBgzfHierarchy 1` enables this opt-in admission. `--dynamicThreadBalance 1` also enables
the existing outer policy on the MAP+FEATURE path. The old PF/ATAC controllers must be off; invalid
combinations are rejected before workers start. The new sampler runs independently of verbose
telemetry, excludes inactive ATAC from the outer policy and uses read-pair progress counters. Other
legacy controller modes remain covered by their existing tests.

PF phase diagnostics separate read processing, hash merge/cleanup, pending barcode rescue, UMI
collapse, filtering, output and cleanup. Their CPU field measures the coordinator thread only.
The owner's priority is correctness, then whole-process elapsed time, full-workload tail including
gathers, then average CPU. Do not optimize feature duration independently of total completion.

## Why the raw singleton changed

The original raw singleton was `IL1B_sg2_HEK / GGGTCAGCACCTCACT`, UMI `GACTAACCCAAG`.
A later arm also exposed `CGAAGGGGTAGGACTT`, UMI `CACCAGTCGGCA`. The feature is a seven-base
sequence (`TGAACCA`), unlike most of the full guide sequences. Relevant read pairs appear once at the
producer and once at the consumer in both native BGZF and gzip ledgers. The count differences already
exist before thread-hash merge, and deduplication has a unique single-read winner. They are not an
acceptable equally supported tie and there is no observed target-read overlap/loss at BGZF boundaries.

Two independent defects were isolated:

1. With zero permitted Ns, expansion buffer sizes use a negative shift exponent and can be zero.
   The actual fallback modifies input sequence bytes when storing match position. This code predates
   the scheduler (allocation/index expressions trace to January); enabling the old broad learning
   fallback exposes it. The fix provides one slot for the original sequence and checks N bounds before
   storing indices. The actual-read regression fails on the previous clean library, and passes with
   the fix, including ASan/UBSan. The correct match positions are 54 and 70.
2. Serialized feature learning still consumes whichever lane is ready. This changes which reads get
   the broad learning search versus the later learned-position window. Both controlled source-order
   arms learn a uniquely supported position 79 for the short feature, yet accept different singletons.
   A two-lane fixture holds the 200,000 input pairs constant and delays either producer by two seconds.
   The zero-N-only library differs at 215 raw entries. Deterministic round-robin consumption during
   learning yields identical histograms and all raw entries (147,638 UMIs) under both delays. After
   learning, normal opportunistic lane consumption resumes. Learning remains 100,000 reads; search and
   acceptance rules are unchanged.

## Evidence and validation

Artifacts: `/home/lhhung/pf_larry_regression_20260911/hierarchy/`.

- `step1/`, `step2/`, `step2_demand/`, `step4/`: frozen clean builds and controller/input/CLI tests.
  Tests include one-permit real BGZF input, ordered mates, cancellation, CRC/name/short-mate errors,
  demand reservation, legacy saturation policy, and actual outer transfers with telemetry disabled.
- `a375_inner_*`, `a375_demand_*`, `a375_fixed4`, `a375_outer_*`: all earlier controls and failed raw
  parity comparisons. Keep these; do not overwrite or relabel cross-build results as matched controls.
- `umi_trace/`, `umi_trace_modes/`: read, source-order and pre/post-merge ledgers; failing UBSan case.
- `zero_n_fix/`: clean zero-N-only binary, preserved library/source and before/after regression evidence.
- `learning_order_fix/`: opposite producer-delay controls, exact input hashes, source patches, build
  commands, binary hashes and exact matrix comparisons.
- Final clean development binary: `learning_order_fix/build/STAR`, SHA256
  `ce9d8142f6e8f0f52d4755b669265166f4a08e673c65650be362e8ffe2577a3c`.
  Clean PF+core build and zero-N/anchor tests pass. Full corrected scheduler controls completed under
  `a375_corrected_{fixed4,inner4,balanced4}` in 192.10/190.62/190.60 seconds. All six matrices, three
  guide-call CSVs and feature totals match exactly. All runs exited successfully without permit leaks.
- `a375_corrected_summary.json`, `a375_corrected_{inner4,balanced4}_comparison.json`: accepted same-build
  comparisons, complete runtime/CPU/RSS/phase data. Peak RSS is ~44.2 GiB in all three arms.
- `a375_correctness_cross_build_impact.json`: compared with the older Step 4 build, 421 raw feature
  entries change (net +36 UMIs), 329 filtered feature entries change (net -28). GEX matrices/cells,
  guide identities and thresholds stay identical; 178 per-cell guide rows have different UMI counts.
  Keep this consequence of correctness fixes separate from exact scheduler parity.

## Performance audit after correctness

The runbook records the Opus priorities: approximately 27 seconds from the successful-finish message
to actual exit, consecutive post-map gather/finalization, and low-CPU startup. Steady GEX-only mapping
already occupies about 31.3/32 permits. This is permit occupancy, not a measured 98% whole-run CPU.

The saved Step 4 fixed-four log has about 54 seconds before mapping, including only 7 seconds explicitly
labeled genome loading. Solo counting takes about 5 seconds, then pf-multi merge/finalize plus guide
calling takes about 13 seconds. PF's separate 8.6-second thread merge and 2.9-second dedup overlap MAP.
The exit gap consumes roughly one core while RSS falls: cleanup is a hypothesis requiring attribution.
Do not sum these reported durations as proven recoverable savings or remove cleanup without analysis.

Subsequent authorized storage work resolved most of this startup/exit cost: CbCorrector now uses khash
tables and a contiguous ambiguity-candidate array. Full A375 fixed-four wall falls 192.10 → 145.47 s,
balanced wall 190.60 → 146.21 s, with all six matrices and three guide tables exactly preserved.
Peak RSS falls ~44.2 → ~39.8 GiB and the exit gap falls ~28.4 → ~1 s. See the
[storage benchmark](benchmarks/A375_CB_CORRECTOR_FLAT_HASH_20260911.md) for the new binary, component
measurements and phase metrics. These later results supersede the unresolved exit hypothesis above;
the gather/ETA work below remains separate.

MAP's current permit-completion signal is after mapping-worker join, before Solo gather. FEATURE's
signal includes its assignment/gather/dedup/output, before the shared GEX+feature finalization. Releasing
idle processing permits and reporting a complete workload are separate events. Full output timings must
include Solo and the common downstream finalization; read-drain ETA alone does not predict gather cost.
The corrected balanced trace demonstrates this limitation: at 69.29 seconds, FEATURE read-rate decay
in gather inflates its estimated ETA and transfers one floor permit back to FEATURE. Idle floors are
borrowable, so this is not proof of reduced MAP CPU or a measured total-wall penalty. Treat explicit
phase-aware ETA as outstanding work in the phase audit, not an already solved completion estimator.

LARRY subset/full P02 validation remains pending in the runbook. Its plain-gzip/zcat input does not
exercise native BGZF admission, and zcat CPU is outside the permit pool. Missing tail-cell guide rows
and the non-Flex EmptyDrops minimum are separate defects, unchanged here. No commit/release or cloud
lifecycle work was performed.
