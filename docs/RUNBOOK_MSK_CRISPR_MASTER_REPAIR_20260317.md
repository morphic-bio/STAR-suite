# Runbook: MSK CRISPR Master Repair

**Date**: 2026-03-17  
**Status**: Proposed plan for review  
**Primary source**: `docs/HANDOFF_MSK_CRISPR_REGRESSION_DEBUG_20260317.md`

## Goal

Repair the MSK 30polyKO mixed-chemistry CRISPR path on the current `master`
checkout by restoring the missing `CR-Larry-perturb` behavior that is still
required for correctness, while avoiding a destructive overwrite of newer local
work already present in this tree.

## What We Know

- The paper benchmark regression is real and is not explained by GEX/index/poly-G
  differences. The bad run used the `master` binary and lost major chunks of the
  CR-Larry feature-assignment path.
- The current checkout is **not** plain `master`. It already contains local WIP
  in the same files, plus post-squash fixes on `master`.
- Because the worktree is dirty and the overlap is non-trivial, the repair
  should be a **manual backport/file-by-file overlay**, not a blind merge or
  cherry-pick of the 15 missing commits.
- GEX expression looks unaffected. The immediate problem is feature assignment
  and downstream CRISPR calling on the MSK mixed TRU/NXT/TRU run.
- The MSK 100k guide fixture does **not** currently point to a wrapper-level
  namespace flip:
  - guide library uses explicit `star_chemistry=nxt`
  - whitelist is 1-column `3M-february-2018_TRU.txt`
  - known-good and bad runs both record
    `effective_read_namespace=NXT`,
    `assignment_whitelist_namespace=TRU`,
    `translate_nxt_for_assign=yes`
- The first hard divergence on that fixture is inside feature matching:
  - reference PolyIII run: `Total feature counts 94394`, `Total barcodes 41418`,
    `Resolve_calls_total 0`
  - bad/current PolyIII run: `Total feature counts 31976`, `Total barcodes 14318`,
    `Resolve_calls_total 200782`
  This pushes the main P0 target into `assignBarcodes.c` guide matching rather
  than `PfMultiProcess.cpp` namespace plumbing.
- A simplified UCSF-style guide run on the same MSK 100k fixture recovered the
  catastrophic top-line regression without touching the matcher:
  - run shape: GEX + PolyIII only, no LARRY, explicit NXT `crWhitelist`,
    explicit CR assign knobs
  - PolyIII result: `Total feature counts 95733`, `Total barcodes 42084`,
    `Resolve_calls_total 200782`
  - final CRISPR calls rebounded from the bad run's near-zero signal
    (`1 feature=39`, `>1 features=3`) to a plausible distribution
    (`1 feature=6466`, `>1 features=223`), though still below the reference
    (`1 feature=8898`, `>1 features=459`)
- That means the current priority should shift from "matcher-only regression" to
  "multi-feature whitelist/control-path regression that can trigger the bad
  14k-barcode collapse on MSK".
- Small mixed-fixture tracing confirmed one real downstream namespace bug in the
  multi-feature path:
  - raw PolyIII guide barcodes overlapped filtered GEX barcodes only after the
    2-base NXT->TRU transform (`292` matches on the smoke fixture)
  - deferred `cr_assign/.../filtered` generation was matching feature
    `barcodes.tsv` against already-TRU-normalized GEX barcodes without first
    normalizing the feature-side namespace
  - this collapsed the final per-library PolyIII `filtered/` output from the
    expected `292` guide barcodes to `1`
  - patching `writeDeferredFilteredAssignOutput()` to normalize feature
    barcodes to TRU before subsetting restored the correct `292`-barcode smoke
    result
- Important: that deferred-filter bug was real, but it is **not sufficient** to
  explain the full MSK 100k regression by itself. On the same small fixture,
  the merged filtered feature matrix already carried the correct `292` CRISPR
  cells before this fix; the bug mainly affected the per-library `filtered/`
  artifacts and anything that depended on them directly.

## Current Tree vs `CR-Larry-perturb`

Already present or partially restored in the current checkout:

- Tiered hash search entry points exist in `assignBarcodes.c`.
- Bootstrap offset finalization exists.
- A brute-force bootstrap fallback already exists in `assignBarcodes.c`.
- `PfMultiProcess.cpp` already sets bootstrap/anchor options and has env
  overrides for bootstrap debugging.
- Async feature assignment launch is already present in `STAR.cpp` /
  `PfMultiProcess.cpp`.
- Consumer hot-path fixes on `master` are already newer than the missing branch
  in some places.

Still clearly missing relative to `CR-Larry-perturb`:

- Bootstrap replay buffer + some tiered-path bookkeeping/diagnostics in
  `core/features/process_features/src/assignBarcodes.c`
- `use_hot_hash` / `skip_heatmaps` config plumbing across `pf_api`,
  globals, and `PfMultiAssign`
- The later `PfMultiMerge` / `PfMultiProcess` MEX optimization stack
  (`pruneZeroCountFeatures`, vector remap, streaming gzip, in-memory GEX/MEX,
  Phase A/B split result object)

## Merge Strategy

- Do **not** run a raw `git merge` into the current dirty tree.
- Do **not** replace entire files with branch versions.
- Backport by topic, reconciling each change against:
  - current local WIP
  - newer `master` consumer-path fixes
  - current CLI/API options such as `probeOnly` and `skipQcOutputs`

## Proposed Changes

### P0: Must Restore for MSK Correctness

- [ ] Reproduce and isolate the bad MSK path as a multi-feature whitelist/control
  regression first.
  Concrete diff surface:
  - compare the bad mixed run against the simplified UCSF-style NXT run
  - trace what changes when `crWhitelist` is explicit NXT versus inherited TRU
    plus `translateNxt`
  - inspect how the multiple-feature path chooses filtered barcode inputs,
    namespace normalization, and whitelist mode when PolyIII is processed beside
    other libraries
  Validation surface:
  - PolyIII `stats.txt`
  - final `protospacer_calls_summary.csv`
  - raw/filtered guide barcode counts
  - `assignBarcodes.api_run.txt`

- [ ] Treat the MSK guide regression as an `assignBarcodes.c` matcher regression
  only after the simplified UCSF-style guide run path is accounted for.
  Concrete audit target:
  - compare the post-bootstrap fast path on current `master` against the
    known-good branch
  - current path uses `pf_search_hash_offsets()` /
    `pf_single_offset_hash_search()`
  - known-good branch used `pf_search_with_offsets()` /
    `pf_hamming_search_single()`
  - verify that the current helper still preserves the old fallback behavior
    when exact/hash probes miss, especially for PolyIII guide reads
  Validation surface:
  - PolyIII `stats.txt`
  - `Resolve_*` counters
  - raw/filtered guide barcode counts and MEX totals

- [ ] Backport the remaining bootstrap/tiered assignment logic from
  `CR-Larry-perturb` into
  `core/features/process_features/src/assignBarcodes.c`.
  Scope:
  - bootstrap replay buffer for reads that miss during bootstrap
  - `feature_mode_search_offsets_reset()`
  - guarded bootstrap histogram/finalize behavior from the branch
  - branch-side tiered-path safety helpers (`pf_bloom_maybe_*`,
    replay-drain path, tiered diagnostics where useful)
  Constraint:
  - preserve the newer consumer hot-path fixes already on `master`
  - preserve the current local brute-force bootstrap fallback unless testing
    shows the branch replay path supersedes it cleanly

- [ ] Backport the missing `process_features` config plumbing for
  `use_hot_hash` and `skip_heatmaps`.
  Files:
  - `core/features/process_features/include/globals.h`
  - `core/features/process_features/src/globals.c`
  - `core/features/process_features/include/pf_api.h`
  - `core/features/process_features/src/pf_api.c`
  - `core/features/process_features/src/main.c`
  - `core/features/process_features/src/mex_writer.c`
  - `core/features/process_features/src/assignBarcodes.c`
  - `core/legacy/source/PfMultiAssign.h`
  - `core/legacy/source/PfMultiAssign.cpp`
  Constraint:
  - keep current `probeOnly` and `skipQcOutputs` behavior on `master`
  - treat this as an overlay, not a rollback

- [ ] Reconcile `PfMultiAssign` option wiring with the current tree.
  Concrete changes:
  - wire `useHotHash` through to `pf_config_set_use_hot_hash()`
  - add `skipHeatmaps` and wire it through to `pf_config_set_skip_heatmaps()`
  - keep current `skipQcOutputs` support
  - keep current api-run summary fields and add branch fields only if they
    still match current permit-controller data

- [ ] Audit the mixed-chemistry namespace path for the MSK layout and patch only
  the deltas that are still missing on the current tree.
  Focus areas:
  - NXT PolyIII assignment against TRU GEX whitelist
  - `translateNxt` behavior through assign output and MEX merge
  - filtered-barcode normalization for mixed namespace runs
  - per-library `star_max_hamming`
  Rationale:
  - MSK is `GEX=TRU`, `gRNA=NXT`, `LARRY=TRU`
  - this is now back in scope because the simplified UCSF-style NXT run largely
    fixes the catastrophic MSK guide collapse
  Status:
  - per-library `star_whitelist` override is now implemented
  - deferred filtered assign output namespace normalization is now fixed for the
    smoke fixture
  - remaining work is the deeper guide-assignment regression on the full MSK
    dataset

### P1: Likely Useful but Not Required for First Correctness Pass

- [ ] Bring over the low-risk CR-compat behavior from `CR-Larry-perturb` in
  `core/legacy/source/PfMultiProcess.cpp`.
  Candidates:
  - `assignOpts.skipHeatmaps = true`
  - the `kPfConsumersPerProducer` producer/consumer apportioning logic
  - any branch-only logging that helps validate mixed-chem assignment
  Note:
  - this is more about runtime/QC noise than the primary MSK correctness bug

### P2: Defer Unless Post-Fix Runtime Is Still Bad

- [ ] Defer the large MEX optimization stack until after correctness is back.
  Files:
  - `core/legacy/source/PfMultiMerge.cpp`
  - `core/legacy/source/PfMultiMerge.h`
  - `core/legacy/source/PfMultiProcess.cpp`
  - `core/legacy/source/PfMultiProcess.h`
  - `core/legacy/source/STAR.cpp`
  - `core/legacy/source/SoloFeature.h`
  - `core/legacy/source/SoloFeature_outputResults.cpp`
  - `core/legacy/source/SoloFeature_processRecords.cpp`
  Features:
  - vector-remap merge
  - streaming gzip
  - pre-prune feature MEX
  - in-memory GEX/MEX caching
  - fuller Phase A/B split result object
  Rationale:
  - these are valuable, but they are not the first thing blocking MSK CRISPR
    correctness

## Confirmed Fixes

### Deferred Filtered Assign Output: Raw-vs-MEX Namespace Skew

Confirmed root cause:

- `assignBarcodes` raw outputs (`barcodes.txt`, `feature_per_cell.csv`) remain
  in native/raw barcode space.
- `PfMultiMexStub::copyBarcodesTsv()` rewrites `barcodes.tsv` into the
  whitelist-driven output namespace.
- `writeDeferredFilteredAssignOutput()` was subsetting through
  `PfMultiMerge::readMex(assignOut)`, which reads `barcodes.tsv`, then copying
  raw-side files separately.

That meant the deferred filtered writer could silently drop valid
`raw guide assignment ∩ Solo filtered GEX` cells whenever `barcodes.txt` and
`barcodes.tsv` diverged.

Applied fix in [PfMultiProcess.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/PfMultiProcess.cpp):

- deferred filtered assign output now subsets using raw `barcodes.txt`
- filtered `feature_per_cell.csv` is subset from the same raw barcode set
- filtered `barcodes.tsv` is regenerated from filtered `barcodes.txt` via
  `PfMultiMexStub::copyBarcodesTsv()`

Validation:

- Small newer-parameter fixture:
  `/tmp/msk_small_newpaperparams_fix_20260317`
  - filtered guide barcodes improved from `249` to `294`
  - filtered guide deduped counts improved from `299` to `355`
  - previously dropped cells such as `ACGTTCCGTAGCTGCC`,
    `AGACCCGGTCCCTGAG`, and `TACCTGCGTTCAGCTA` now survive the deferred
    filtered writer
  - all `raw guide assignment ∩ Solo filtered GEX` cells now survive
    (`missing_from_filtered = 0`)

- 100K newer-parameter fixture:
  `/tmp/msk_100k_ucsf_newpaperparams_fix_20260317`
  - no material change on this surface
  - filtered guide stats remain `6149` barcodes / `46826` deduped counts

Interpretation:

- this was a real late-stage namespace bug and is now fixed
- it explains the small-fixture filtered-guide dropouts
- it does not by itself explain the larger 100K assignment-level recovery gap

### 100K Guide Trace: NXT Barcode Mis-Rescued Against TRU Whitelist

Confirmed trace on the bad-path 100K `GEX + PolyIII` reproduction:
`/tmp/msk_bad_trace_ns_20260317`

Target barcode:

- raw guide read barcode in FASTQ: `AACGAAAAGAGCGTCG` (NXT)
- true translated whitelist barcode: `AACGAAATCAGCGTCG` (TRU)

Observed in the new targeted trace log:
`/tmp/msk_bad_trace_ns_20260317/namespace_trace.log`

- `NS_CHECK raw=AACGAAAAGAGCGTCG translated=AACGAAATCAGCGTCG translate_NXT=1 feat=18 match_pos=20 raw_hit=0 translated_hit=1`
- `NS_RESCUE ... n_matches=1`
- `cand[0]=AACGAAAAGAGCATCG mismatch_pos=12`
- `NS_DECISION ... decision=unique_rescue chosen=AACGAAAAGAGCATCG`

Interpretation:

- the read enters `checkAndCorrectBarcode()` in NXT space
- the code can see that the translated TRU form would hit the whitelist
- but the bad path does **not** translate before exact whitelist lookup
- instead it falls into the ordinary 1-mismatch rescue path and chooses the
  wrong TRU neighbor `AACGAAAAGAGCATCG`

Downstream consequence:

- the good/recovered run contains the expected cell barcode
  `AACGAAAAGAGCGTCG` with `(num_features=2, top_feature=18, total_deduped_umi=25)`
- the bad run does **not** contain that barcode
- instead the bad run contains the wrong translated output barcode
  `AACGAAATCAGCATCG` with `(num_features=2, top_feature=18, total_deduped_umi=24)`

Supporting outputs:

- good run:
  `/tmp/msk_100k_ucsf_guideflags_nolarry_20260317_1912/cr_assign/CRISPR_Guide_Capture/grna_de/PolyIII/feature_per_cell.csv`
- bad traced run:
  `/tmp/msk_bad_trace_ns_20260317/cr_assign/CRISPR_Guide_Capture/grna_de/PolyIII/feature_per_cell.csv`
- per-read corrected barcode trace:
  `/tmp/msk_bad_trace_ns_20260317/reads_trace.log`

What this means for the repair:

- the main 100K PolyIII failure is not just a late namespace leak
- it is also an early assignment-space error: NXT guide reads are being
  matched against a TRU whitelist via generic mismatch rescue instead of a
  deterministic NXT->TRU normalization step
- that gives silent misassignment, not just dropped barcodes

## Verification Gates

### Gate 0: Clean Rebuild First

Required before any regression validation:

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
```

### Gate 1: Small/Targeted Validation

- Run a mixed-chemistry smoke before the full paper benchmark.
  Suggested surfaces:
  - `tests/multi_feature/run_msk_multifeature_smoke.sh`
  - `tests/multi_feature/run_mixed_chemistry_filtered_smoke.sh`
- Validate using:
  - feature `stats.txt`
  - merged filtered barcodes overlap
  - final CRISPR/LARRY outputs
- Do **not** use only `resolve_resolved` counters as the success criterion.
  The bootstrap handoff already showed those counters are not the right parity
  surface.
- Do use `Resolve_calls_total` as a regression smoke metric for PolyIII.
  Current evidence:
  - known-good MSK 100k guide run: `0`
  - bad/current MSK 100k guide run: `200782`
  - simplified UCSF-style MSK guide run: `200782`
  This means `Resolve_calls_total` is useful as a smoke metric, but it is not
  sufficient by itself to explain the catastrophic MSK failure.

### Gate 2: Full MSK Benchmark

Run:

```bash
scripts/paper/run_msk_30polyko_benchmark.sh
```

Expected direction of change versus the bad `master` run:

- total feature counts move back toward the reference scale
  (`~100M`, not `~34M`)
- total deduped feature counts move back toward `~4.7M`
- gRNA-assigned cells recover dramatically
- zero-molecule filtered cells collapse from the catastrophic `28,550` level

### Gate 3: Regression Guard

After MSK repair:

- confirm GEX outputs are unchanged within expected tolerance
- confirm UCSF CRISPR path still behaves
- confirm mixed-chemistry filtered-barcode normalization still works

## Open Questions

- Is the replay buffer strictly required once the current local brute-force
  bootstrap fallback is present, or is it just the cleanest way to recover the
  original branch behavior?
- Do we want `useHotHash` enabled by default, opt-in only, or not at all for
  the first repair pass?
- Do we want to pull the `PfMultiMerge` optimization stack into the first repair,
  or keep the first pass tightly scoped to assignment correctness?

## Recommended Execution Order

1. Diff the bad MSK mixed run against the simplified UCSF-style NXT guide run.
2. Trace the multi-feature whitelist/namespace/filtering control path.
3. Patch the smallest control-path delta that reproduces the UCSF-style recovery
   under the normal MSK mixed run.
4. Only then revisit the PolyIII matcher path in `assignBarcodes.c` if the bad
   collapse still persists.
5. Clean rebuild.
6. Run targeted mixed-chem validation.
7. Run the full MSK paper benchmark.
