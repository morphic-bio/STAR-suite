# Cross-module performance implementation and regression results

Runbook: [ordered implementation](../RUNBOOK_CROSS_MODULE_PERFORMANCE_20260910.md).
Source baseline: `9e146eef20080108cc075a977938757f90221dbc`; implementation branch:
`perf/cross-module-bgzf-callers`. Evidence root:
`/mnt/pikachu/star_suite_paper/analysis/cross_module_performance_20260910`.
Final acceptance passed, with the TranscriptVB native-input restriction below.
The source is integrated into master with a no-fast-forward merge. The local host is an Intel i9-13900KF
with 128 GiB RAM; the recorded affinity and worker budgets identify each arm.

## Measured changes

All entries are single serialized measurements. Conversion is outside timing.
Whole-process memory can be dominated by input construction or the reference.
No cloud provisioning, full 320K rerun, or new release is part of this work.

| Workload | Frozen control | Changed implementation | Correctness / interpretation |
| --- | --- | --- | --- |
| Full A375 PF, 9,748,584 pairs, same BGZF inputs | 28.2722 s, 1,554,892 KiB RSS | Queued BGZF 25.5696 s; direct BGZF **18.0532 s**, 1,576,536 KiB RSS | Exact 3,218,159 UMIs, 110,772 CBs, 11 features, 124,871 nonzeros and assignment statistics; direct is 36.1% faster than control |
| Saved detected A375 GeneFull counts, 290,020 CBs, 7,533,915 nonzeros, 100K simulations | MC1: 24.0404 s, 380,400 KiB RSS | Borrowed strided view, MC8: **3.68962 s**, **296,876 KiB**; split view MC3: 10.5522 s, 296,152 KiB | Exact 1,188 cells, 233 tail candidates, likelihoods, p-values and decisions; logical bootstrap streams stay at four |
| SLAM solver, 1,500 deterministic histograms, both models | 0.377186 s | **0.093858 s** | Exact serialized fits, likelihoods, iterations and convergence; microbenchmark, not an end-to-end speed claim |
| Full saved public bulk EC evidence, 50,917,353 pairs | VB engine 69.61 s; 14,464,320-byte accumulation buffer | Component engine **68.52 s**; buffer eliminated | Byte-identical transcript/gene quantification, 1,357 iterations, GC update at 11; engine gain 1.6% |

The bulk finalizer's overall times were 95.6926 and 70.3726 seconds. Initial
transcriptome loading took 20.97 and 1.06 seconds respectively, and output also
varied. Those startup differences are excluded from the executor speed claim.
The full PF benchmark uses the full feature library; it is not a full GEX run.
The caller harness retains both split and strided source representations, so
its RSS comparison isolates the C API copy removal and does not measure the
additional removed STAR adapter copies.

## Integrated SLAM

Same real public FASTQs reblocked as BGZF, matched GRCh38-2024-A index, eight
CPU affinity slots, no BAM/sidecar, and identical model/output options. Blank
and treatment inputs each contain 100,000 reads. The PE case is a diagnostic
paired clone of treatment input, not an independent biological paired library.

| Input / model | Control whole wall / RSS KiB | Changed whole wall / RSS KiB | Parity |
| --- | --- | --- | --- |
| Blank SE, fixed trims | 19.759 s / 30,028,288 | 21.079 s / 30,042,144 | Eight output files and QC JSON exact |
| Treatment SE, automatic trims | 20.203 s / 30,348,336 | 19.020 s / 30,357,588 | Eight output files and QC JSON exact; both select 5′=10, 3′=15 |
| PE clone, overdispersed model | 35.919 s / 30,348,700 | 29.597 s / 30,363,504 | Eight output files and QC JSON exact |

GRAND-SLAM sample-prefix labels are normalized for comparison. The biological
blank/treatment fixtures establish output correctness; small whole-run timing
differences are not a general SLAM speedup. Per-gene fit results add one small
fixed-size result per reference gene and are invalidated on histogram, rate
or model changes. Standard and GRAND-SLAM writers reuse one fit pass.

## Input and caller compatibility

Native BGZF feeds ordered raw windows through STAR's established input parser.
It preserves FILE markers, names/qualities, read limits, lane transitions,
two-pass processing, transcriptome BAM, Y/noY taps, paired input-order output,
and scRNA two-/three-read layouts. PF uses independent mate block boundaries
and leased batches, with the established matching kernel. Chemistry detection,
three-stream input and order-sensitive rescue retain the queued PF path.

Inflation and mapping/assignment share CPU permits. The ordinary core path
also enforces this budget with dynamic retuning disabled. FLEX retains its
existing fused reader ownership. Decoder corruption, mate mismatch, truncated
records, capacity overflow, cancellation and worker failures have explicit
regression coverage. Failed input must not produce a final MEX.

OCM schedules sample preparations within a configurable estimated matrix
budget (`--ocmCellCallMaxMemory`, default 1 GiB). A sample larger than that
budget runs alone. This is an admission estimate, not an RSS limit. Sample
routing and output order remain stable, and empty groups produce empty calls.
The initial native tiny fixture matched all 60 output files at one/four workers.

## Final integrated and scheduling gates

| Workload | Control whole wall / RSS KiB | Final whole wall / RSS KiB | Acceptance |
| --- | --- | --- | --- |
| A375 perturb, 200K GEX pairs plus 200K guide pairs, eight CPUs | 34.870 s / 31,033,860 | 21.863 s / 31,058,404 | Six keyed matrices and 19 caller/guide diagnostic files exact; 269 GEX cells, four internal PF cells, 114 guide singlets |
| Bulk, 1M pairs, deterministic one-worker online learning, tximport output | 105.534 s / 31,136,140 | 93.153 s / 31,135,948 | All 226,005 transcript rows and both 38,606-gene outputs exact; 437 iterations in both arms |
| SLAM blank+treatment, per-file automatic trimming and reopen/skip | 40.203 s / 30,505,784 | 36.246 s / 30,523,572 | Eight scientific output/QC files exact |
| Larger OCM routing fixture, one versus three workers, 1,000 simulations | 32.651 s / 231,596 | 28.754 s / 263,416 | All 60 outputs exact; peak workers 1/3, all permits returned |

Whole-run timings include reference loading or MEX compression and are not
isolated estimates of executor gains. In particular, the accepted bulk input
uses the established gzip reader in both arms; the component option is the
changed executor. The eight-worker saved EC benchmark supplies its isolated
engine measurement. The first, rejected integrated bulk control took 1,095 s
because of a cold index load, so its wall time must never be used as the
baseline for a claimed native BGZF speedup.

The larger OCM fixture preserves every count in the detected A375 matrix and
adds synthetic tag suffixes to create groups of 116,008 / 87,006 / 87,006 / 0
columns, plus a 203,014-column union. Both schedules called 464 / 342 / 352 / 0
and 826 cells. The 100,000,000-byte admission budget allowed the oversized
110,539,072-byte union estimate to run alone, as specified. Peak RSS increased
with concurrency. This tests routing and scheduling on real count distributions;
it is not a biological OCM benchmark or a production 100K-simulation speed claim.

Final small gates cover static decoder/mapping budgets, core and SLAM batches,
fused FLEX ownership/count parity, the tag-aware caller, empty/one-component VB,
SLAM requant one/four workers, PF malformed direct batches and CBQ assignment.
Ordinary scRNA two-/three-read adapter parity uses the pinned synthetic fixture;
full caller arithmetic is tested on saved real GEX counts. The updated PF
adapter's final internal four-cell callset and outputs match its frozen control.

## Findings that changed implementation

1. **PF pre-MEX sparse extent:** the adapter populated sparse arrays but left
   `sparse_nnz` at zero. The new bounds checks correctly rejected this, and an
   integrated run silently skipped its nonfatal internal PF EmptyDrops step.
   The adapter now derives the occupied extent from offsets/counts without
   narrowing. The failed run is retained and cannot support a speed claim.
2. **Online bulk model sensitivity:** an eight-worker native-reader comparison
   had exact alignment summary counts but different upstream EC/fragment-length
   models (67,524 versus 67,516 ECs; FLD means 280.858 versus 280.729), before
   component VB ran. This failed the pinned numerical and iteration gates
   (438 versus 536 iterations). The same saved ECs produced byte-identical
   outputs. Native BGZF is therefore **disabled for TranscriptVB online model
   learning**: auto logs a fallback to the established reader; forced range
   rejects it. General bulk alignment/GeneCounts BGZF support remains available.
   `--quantVBComponentParallel 1` remains opt-in and keeps global convergence.
3. **Error/progress audit:** direct PF workers now report failures without a
   shared-buffer race, release permits on validation errors, check capacity
   before inspecting the last byte, and preserve probe-only early stopping.
   EOF batch recycling explicitly wakes waiting producers. Default quality
   synthesis excludes a sequence newline, including at the buffer boundary.
4. **Fixture limitations:** error-free toy SLAM reads cannot estimate trim
   variance; real public reads supplied the successful auto-trim gate. An old
   cached FLEX command included a legacy expected-cell override rejected by
   today's tag-aware caller; the corrected fixture removes that override.

No Cell Ranger code was inspected. This work changes execution and storage,
not the cell-calling estimator, feature inclusion, rescue thresholds, priors,
MC seeds, logical bootstrap streams, or biological concordance targets.

## Evidence and reproducibility

The evidence root retains frozen executables, build logs, input hashes, exact
commands, successful and failed wrapper completion records, and comparison
results. Main directories: `00_inventory`, `01_pf_bgzf_tests`,
`05_a375_pf_full`, `phase3_tests`, `phase456_tests`, `kernel_gates`, `vb_saved`,
`integrated_inputs`, `integrated_runs`, `final_gates`, `followup_integrated`,
`ocm_large`, `boundary_gate`, and `boundary_cbq`. The final registry is `execution_registry.jsonl`; acceptance is recorded in
`FINAL_ACCEPTANCE.json`, with `final_source_manifest.json` recording source and
frozen executable hashes. Existing controls are reused; changing only an output path
never authorizes a repeat. Repository test sources accompany the implementation.
