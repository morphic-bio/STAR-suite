# Feature search policy revert: full A375 validation

Date: 2026-09-11. Status: narrow source fix built and tested locally; uncommitted,
unreleased. A375 assignment recovery verified. A375 runtime regression remains
open. Full LARRY performance has not been validated with this patch.

For the current task, constraints, and ordered next steps, start with
[the next-agent handoff](HANDOFF_LARRY_A375_REGRESSION_NEXT_AGENT_20260911.md).

## Change and diagnosis

The only production source change is in `core/legacy/source/PfMultiProcess.cpp`:
restore `assignOpts.featureModeBootstrapReads = 100000`, including when an explicit
feature offset is supplied. This restores the historical feature-position learning
policy. It is unrelated to cell-calling bootstrap or UMI tie breaking.

Commit `864b04f2bcbcbd48dbddfae70f8255a6d54941ca` (2026-05-31,
"Add PF CBQ direct range mode") disabled learning when
`featureConstantOffset >= 0`. Both MSK and A375 use offset 0. This change predates
the July checkout used to bound the original handoff's regression search; that
checkout is not evidence of the exact source executed in April.

Without learning, the required-anchor path loops over features and searches their
anchors. LARRY has 245,979 features sharing one suffix anchor, so this repeats the
same search extensively. A clean API diagnostic found all 245,979 exact feature
keys in the hash. With BGZF disabled, the first 20,000 LARRY read pairs still
generated 244,880,980 resolver calls. This demonstrates repeated resolution in
the required-anchor path and shows that disabling BGZF does not remove it.
Reference self-hits verify population of the exact hash, but do not establish
real-read H0/H1 hit rates or exclude another large-cache defect. Full LARRY
validation remains outstanding.

The fix preserves existing environment overrides and the explicit split-read
`searchMode == "free"` override. No search optimization, reader implementation,
or permit scheduling change is included. An earlier draft anchor-grouping
optimization was discarded before building this binary.

## Full A375 test

The actual STAR GEX + CRISPR workflow completed with exit 0 and its wrapper-written
`BENCHMARK_SUMMARY.txt`. It processed all 9,748,584 feature read pairs, using 32
threads, shared dynamic permits, BGZF auto, feature offset 0, and no BAM. Inputs,
reference and options match the saved 1.9.2 P01 wrapper. This is diagnostic
validation of an unreleased patch, not a replacement paper benchmark row.

| Measurement | Saved March 26 (`a70a039-dirty`) | Saved v1.9.2 P01 | v1.9.2 + policy revert |
| --- | ---: | ---: | ---: |
| Assigned feature reads | 9,274,337 | 9,063,551 | 9,274,337 |
| Deduplicated feature UMIs | 3,221,810 | 3,191,995 | 3,221,811 |
| Feature API time, seconds | 11.8753 | 58.4219 | 59.3113 |
| Whole wrapper wall time, seconds | — | 192.34 | 195.24 |
| Peak RSS, GiB | — | 44.704 | 44.706 |
| GEX cells | — | 1,170 | 1,170 |

The feature matrix was compared by feature name and barcode, independent of row
and column order. Relative to v1.9.2, 3,032 entries changed, with a net gain of
29,816 UMIs and absolute difference of 29,918 UMIs. Assigned reads increased by
210,786.

Relative to March, exactly one entry differs: `IL1B_sg2_HEK` at barcode
`GGGTCAGCACCTCACT` has 1 UMI now versus 0 historically. That barcode is not in the
current filtered GEX cell set. The cause of this residual UMI is not established;
do not describe the raw feature matrices as identical.

Raw and filtered GEX `barcodes.tsv`, `features.tsv`, and `matrix.mtx` are all
byte-identical to saved v1.9.2 P01.

The log confirms feature-position learning finalized with 2 base offsets and 6
search positions. BGZF now uses the existing ring handoff because learning is
order-sensitive. The API time remains about 59 seconds: the A375 slowdown is
not resolved by this change. Permit wait counters alone do not establish a
separate scheduling defect. No scheduling change was tested.

## Provenance and artifacts

- Worktree: `/mnt/pikachu/STAR-suite-larry-regression-20260911`, branch
  `fix-larry-feature-regression-20260911`.
- Base revision: `68eb92ff3cc107179d32ed52d95229b121debae9` (v1.9.2).
- Binary: `core/legacy/source/STAR`, reporting that revision with `-dirty`.
- Binary SHA256: `c0e399de851638f6f3eeab6f176fe95b393dee19a0503636aa16d651a722d941`.
- Both core and process_features were cleaned before a full Chromap-enabled
  `make -C core/legacy/source -j8 STAR` build.
- Build logs: `/home/lhhung/pf_larry_regression_20260911/revert_{clean,build}.log`.
- New run: `/home/lhhung/pf_larry_regression_20260911/a375_revert_star/`.
  Contains `source.patch`, `binary.sha256`, `run_a375.sh`, `stdout.log`,
  `time.txt`, `completion.json`, `comparison.json` and `out/`.
- Comparator: `/home/lhhung/pf_larry_regression_20260911/compare_a375_revert.py`.
- Wrapper source:
  `/mnt/pikachu/star_suite_paper/analysis/paper_benchmark_refresh_20260910/wrappers/run_a375_benchmark.v192.CORRECTED.sh`.
- Saved v1.9.2 comparator: `/storage/paper_bench_v192_20260911/out/P01/`;
  wall/RSS evidence: `/storage/paper_bench_v192_20260911/logs/P01.time.txt`.
- Historical comparator: `/storage/A375/paper_bench_20260326_134444/`.
- Feature matrices under each output root:
  `cr_assign/CRISPR_Guide_Capture/A375_CRISPR_Guide_Capture_1/crispr/`.

The earlier 20K LARRY API run with only 1,000 learning reads is a mechanism
diagnostic, not a production-default or historical timing validation. The
production patch uses 100,000 learning reads. The original full LARRY timing and
matrix acceptance gates remain outstanding. Keep subsequent changes isolated
and test one defect at a time, as requested by the owner.

## Follow-up: instrumented process_features comparison

The owner requested direct cache-resolution instrumentation before using timing
as a final check, and requested earlier diffs because the July v1.4.3 tag is not
sufficient to bound the regression.

Clean, isolated process_features builds were made from `74b6e05^` (`f2ac1e2`,
before the May 29 changes), `v1.4.3` (`e4292f1`), and `68eb92f`. Each processed the
same first 200,000 A375 lane-1 read pairs, with one consumer and the gzip reader.
Each revision was tested once with learning disabled and once with the production
100,000-read learning threshold. Counters only were added; these runs are not
performance measurements. All six executions completed successfully.

**For each setting, all recorded counters match exactly across the three
revisions.** Selected counts:

| Counter | Learning disabled | Learning at 100,000 reads |
| --- | ---: | ---: |
| Feature-resolution calls, including replay | 200,000 | 203,665 |
| Successful feature resolutions, including replay | 183,876 | 191,469 |
| Per-read offset-cache hits | 0 | 1,537,396 |
| Offset-cache fills | 0 | 801,908 |
| Successful resolver returns rejected by the outer anchor feature-index check | 207,881 | 0 |
| Full feature-hash lookup calls | 0 | 0 |

The zero full-hash count is expected for this fixture: A375 has 11 guides, so
`simple_search` uses direct comparisons (threshold: fewer than 150 features).
Two guides differ from the common length, disabling Hamming prehash and the
uniform-length tiered path. This test does **not** rule out a defect in the large
LARRY hash. It establishes that the same A375 reads take the same resolution path
at these revisions when the configuration is held fixed.

The earlier diff audit found:

- Exact lookup, `simple_search`, `simple_hash_search`, `simple_hamming_search`,
  and `simpleCorrectFeature` bodies are unchanged from `5a20ef2` (March 17).
- March 18 changes include hot-hash handling, a fallback for Hamming distances
  greater than 2, bootstrap replay, and namespace fixes. The >2 fallback does not
  apply to the benchmark's Hamming ceiling of 1.
- `74b6e05` (May 29) synchronizes/serializes feature-position learning and fixes a
  match-position output. The single-consumer trace agrees before and after it;
  this is not a concurrency validation.
- `864b04f` (May 31) disables learning for explicit offsets in the STAR caller.
  This bypasses the offset cache and, for uniform-length LARRY, prevents reaching
  the learned-offset tiered path. It remains the specific policy change reverted.
- CATATAC's split-read `free` mode explicitly disables required-anchor matching
  and uses a different resolution branch. The user's report of unaffected
  multiomics is a useful control but has not yet been mapped to a particular run
  configuration. Small feature sets can also mask costs that grow with feature
  count; do not infer that all assays exercise the same cache.

Artifacts: `/home/lhhung/pf_larry_regression_20260911/cache_trace_a375/`, including
clean source snapshots, build logs, `input_manifest.json`, counter files,
`source_comparison.json`, and `earlier_source_comparison.json`. Instrumentation
setup: `/home/lhhung/pf_larry_regression_20260911/setup_cache_trace.py`.
No diagnostic counters were added to the production worktree.

## Follow-up: uninstrumented v1.4.3 feature-only timing

At the owner's subsequent request, a clean, uninstrumented v1.4.3
process_features library was run on **all 9,748,584 A375 feature read pairs**,
with 31 consumers, no GEX, no shared permits, no EmptyDrops, and no feature QC
plots. Learning was 0, matching the explicit-offset policy in the v1.4.3 tag and
the earlier current-version feature-only control. The v1.4.3 reader uses gzip.

- Read processing: **6.74 seconds**.
- Feature API processing: **9.663660049 seconds**.
- Whole process including initialization: **11.30 seconds**.
- Exit status: **0**.
- Assigned reads: **9,063,551**; deduplicated UMIs: **3,191,995**.
- Raw feature matrix is semantically identical to saved v1.9.2 P01.

The previous current-version feature-only control also completed quickly
(7.4115 seconds API time, BGZF auto). Reader modes differ, so do not present
these as a controlled version speedup ratio. Neither feature-only test reproduces
the approximately 59-second integrated A375 phase. The cause of that integrated
delay remains unproven; wait counters alone are insufficient.

Artifacts: `/home/lhhung/pf_larry_regression_20260911/a375_v143_feature_only/`:
`provenance.json`, `binary.sha256`, `command.json`, clean build logs, `stdout.log`,
`time.txt`, `completion.json`, `comparison.json`, and `out/crispr/`.

## Follow-up: provenance audit

The owner expects original paper runs to have used v1.4.3 and directed inspection
of the provenance folder. The available records contain conflicting version
metadata; retain the distinction between the intended baseline and observed
run provenance.

- `/mnt/pikachu/STAR-suite-provenance/records/paper-2026/perturbseq-full-depth.json`
  labels the software `1.7.0-candidate`, revision `c12ad552`, recorded August 8.
  Its A375/MSK elapsed times are 241/1,656 seconds and it links
  `STAR-suite/comparisons/paper_benchmarks_20260318`.
- That archive's `compiled_stats.tsv` points to
  `/mnt/pikachu/paper_bench_rerun_20260318_065211/{a375,msk_30polyko}/`.
  Both actual `Log.out` headers record **`ebc13622459ecc5c5d79cafceb887e4f3c3c30df`**,
  branch `upstream-benchmarking`, compiled March 18, with an empty dirty-file list.
  These are direct run records, not an inference from directory names.
- The historical comparators used by the regression handoff are different runs:
  March 26 A375 logs record `a70a03929b14331f841988df4ccb87ca26f9a770-dirty`;
  April 30 MSK ES logs record `cc84382afe51a5e06e76dd4073c23a3d5872e54f-dirty`.
  The dirty-file lists in those headers do not include process_features.
- The v1.4.3 release tag resolves to `e4292f1`, July 1. The new isolated timing
  test is explicitly pinned to that tag. The above historical logs do not prove
  those archived rows were executed with that release.

No public provenance record was edited. The aggregate software label needs an
explicit correction/reconciliation with the execution logs; it must not silently
replace the recorded execution revision.

### Correction: the March ledger is not the current run inventory

The owner correctly pointed out that tests continued after March. The aggregate
public record above must not be used to infer which later runs supplied the
current manuscript. `STAR-suite-provenance/README.md` explicitly describes a
curated public evidence ledger, not a project run registry.

Verified later records include:

- April 30 MSK ES full GEX + PolyIII + LARRY: 1,811 seconds / 30.2 minutes,
  33,226 cells. This is the current manuscript's MSK Perturb row and is linked
  by `PAPER_BENCHMARK_REFRESH_RUNBOOK_20260910.md`, rather than the March 18
  1,656-second MSK result in the public aggregate record.
- June 30 MSK ES GEX-only:
  `/storage/MSK-perturb-comparison/msk_gexonly_ES_modern_v1.4.2_20260630_122301/`.
  Its `Log.out` records v1.4.2, clean commit
  `996eea0ce97f127a9f25e459821e7eb76f26a942`; its wrapper summary records 1,894
  seconds / 31.57 minutes and 33,346 cells. This matches the current manuscript's
  31.6-minute scRNA row.
- July 25 A375 CR-compatible regression test:
  `/storage/A375/test_cr_compat_crispr_20260725_010102/`. Its log records v1.5.0,
  `26e6d05a48537761be330fd926b4bef9b68b3e8e-dirty`, 16 threads and BAM output.
  It is evidence of later testing, not an equivalent no-BAM paper timing.

The inspected A375/MSK benchmark wrappers save provenance with their outputs;
the September driver additionally saves binary hashes, arguments, and preflight
records. They do not automatically register these runs in either central
provenance repository. Thus absence from the curated ledger does not establish
absence of later tests. A375's latest equivalent full-workflow baseline still
needs to be identified from the run records; do not assume the March aggregate
is its latest validation.

The owner requested a follow-up fix for this provenance-registration gap on
2026-09-11. It is recorded under Priority 0 in [docs/todos](todos), with automatic
registration, durable retry, explicit receipts, historical backfill, and exact
paper-row/run linkage as acceptance requirements. Implementation is deferred
while the process_features regression remains the active task.
The owner clarified that development and production records may use separate
`dev/` and `production/` directories, but every run in both must record provenance.
The TODO specifies a common schema and searchable index across those directories.
