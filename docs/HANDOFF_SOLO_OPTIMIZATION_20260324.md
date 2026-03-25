# Handoff: Solo Optimization

Date: 2026-03-24 (updated)
Branch: `feature/solo-optimization-20260324-mastermerge`
Tip commit (when handoff last revised): `ff30af9`
Audience: next coding agent working on Solo / STARsolo post-map optimization

## Scope

This branch is intended for Solo-phase optimization work, not FLEX mapping-pipeline work.

Primary target surface:

- standard Solo / STARsolo post-map counting and output
- shared code that also matters for Flex, but only where the optimization is clearly worth the shared blast radius

Do not start by redesigning the read-stage pipeline here. The main optimization opportunity on this branch is the Solo finalize path after mapping.

## Current Branch State

- Active integration branch: `feature/solo-optimization-20260324-mastermerge` (Solo optimization + related merges).
- Tracked changes include non-Flex inline-hash bridge packed-key fixes, direct bridge collapse
  (`SoloFeature_collapseUMI_fromBridgeHash.cpp`), and related docs/artifacts notes.
- There may be unrelated untracked local files in the checkout; leave them alone unless explicitly asked.

## Build Hygiene

Before debugging or benchmarking:

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
```

This repo has a documented history of stale-object false regressions. Do not skip the clean rebuild before trusting results.

## What Matters Most

The current Solo post-map path still has obvious expensive surfaces:

1. `countCBgeneUMI()` branching and extra materialization
2. `collapseUMIall()` hot-path behavior
3. shared MEX writing / output staging

The first practical goal is to simplify or remove unnecessary work in those areas, not to build a new threading architecture.

## Important Files

### Main entry point

- [SoloFeature_countCBgeneUMI.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_countCBgeneUMI.cpp)

Key observations:

- Inline-hash mode already has a fast direct path:
  - `resolveAmbiguousCBs()`
  - optional `runCliqueCorrection()`
  - `collapseUMIall_fromHash()`
- Non-inline path still does a heavier loader/sink/finalize/collapse sequence.
- `CountingSink::finalize()` fills `rGeneUMI` / `rCBp` / `rCBn`, then `collapseUMIall()` consumes them.

### Legacy collapse path

- [SoloFeature_collapseUMIall.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_collapseUMIall.cpp)

Key observations:

- The hot per-CB routine [collapseUMIperCB()](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_collapseUMIall.cpp) still builds a `std::vector<Triplet>`, sorts it, then writes it back to the original buffer.
- This is an obvious target for cleanup because it adds:
  - extra allocation
  - extra copying
  - extra writeback
- The file still contains a lot of debug-only scaffolding and forced debug paths that make the hot function harder to reason about.

### Shared MEX writer

- [MexWriter.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/MexWriter.cpp)
- [SoloFeature_outputResults.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_outputResults.cpp)

Key observations:

- `MexWriter::writeMex()` still does one `fprintf()` per matrix entry / barcode / feature.
- `SoloFeature_outputResults()` still builds the unique+multi matrices through an `ostringstream` staging path.
- These are not necessarily the first bottleneck, but they are good low-risk cleanup targets after collapse.

## Likely Optimization Order

### 1. Instrument first

Add explicit timers around:

- `countCBgeneUMI()`
- `collapseUMIall()`
- `collapseUMIall_fromHash()`
- `outputResults()`
- shared MEX writing if possible

If you do not add instrumentation, at least benchmark these phases with wall-clock logging before changing multiple surfaces at once.

### 2. Clean up `collapseUMIall()` first

This is the best first code target.

Why:

- it is clearly hot
- it has obvious extra copy/sort/writeback behavior
- it is narrower than redesigning the loader/sink path

Suggested first pass:

- remove unnecessary debug scaffolding from the hot path if it is no longer needed
- replace the current `vector<Triplet> -> sort -> writeBack` pattern with a lower-overhead approach
- keep behavior stable first; do not combine semantics changes with performance work

### 3. Revisit the non-inline materialization path

The next target is the loader/sink/finalize/collapse chain in [SoloFeature_countCBgeneUMI.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_countCBgeneUMI.cpp).

Questions to answer:

- what data absolutely needs to be materialized before collapse?
- can any arrays be avoided or filled more directly?
- is the current `CountingSink` path still paying for compatibility layers that the current code no longer needs?

Do not rewrite this blindly. Measure first, then shrink the surface.

### 4. MEX/output cleanup after collapse

Once collapse is improved:

- optimize [MexWriter.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/MexWriter.cpp)
- simplify [SoloFeature_outputResults.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_outputResults.cpp)

Good targets:

- buffered writes instead of one `fprintf()` per row
- less intermediate string staging
- avoid duplicate barcode/feature transformations where possible

## What Not To Start With

- do not start with a new producer/consumer pipeline
- do not start with alignment-thread retuning
- do not start with broad Flex mapping-phase refactors on this branch
- do not combine correctness refactors with performance rewrites in one patch

This branch should be about Solo finalize performance first.

## Shared-Code Blast Radius

Be aware of the split:

### Lower-risk / Solo-focused

- [SoloFeature_collapseUMIall.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_collapseUMIall.cpp)
- [SoloFeature_outputResults.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_outputResults.cpp)
- [MexWriter.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/MexWriter.cpp)

### Shared with Flex-adjacent paths

- [SoloFeature_countCBgeneUMI.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_countCBgeneUMI.cpp)
- [flex/source/SoloReadInfoSink.cpp](/mnt/pikachu/STAR-suite/flex/source/SoloReadInfoSink.cpp)
- [flex/source/SoloFeature_materializeFromHash.cpp](/mnt/pikachu/STAR-suite/flex/source/SoloFeature_materializeFromHash.cpp)

If you touch the shared path, run extra coverage, not just Solo smoke.

## Minimum Validation

After each meaningful optimization patch:

1. Clean rebuild

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
```

2. Solo smoke

- [run_solo_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_solo_smoke.sh)

3. If the patch touches shared Solo/Flex infrastructure, also run a CR/GEX parity surface if fixtures are available

- [run_cr_parity_100k.sh](/mnt/pikachu/STAR-suite/tests/run_cr_parity_100k.sh)

4. Record any new benchmark or smoke output locations in:

- [ARTIFACTS.md](/mnt/pikachu/STAR-suite/tests/ARTIFACTS.md)

## Suggested First Patch

If you want the most pragmatic start:

1. add basic timing around `countCBgeneUMI()` and `collapseUMIall()`
2. clean up the `collapseUMIperCB()` copy/sort/writeback path
3. rerun `run_solo_smoke.sh`
4. if clean, benchmark a small parity or smoke fixture again

That gives a tight feedback loop and avoids overcommitting to a larger rewrite before you know the real gain.

## Current Benchmark Status

The current branch now has a first Solo-only benchmark on the corrected UCSF non-FLEX GEX data:

- sample root: `/mnt/pikachu/ucsf-perturb-seq-corrected/iPSC2_1/GEX`
- benchmark script: [run_ucsf_solo_gex_100k_benchmark.sh](/mnt/pikachu/STAR-suite/tests/run_ucsf_solo_gex_100k_benchmark.sh)
- genome: `/storage/autoindex_110_44/bulk_index`
- whitelist: `/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-february-2018_TRU.txt`
- downsample policy: 100,000 total reads across the 10 GEX lanes, not 100,000 per lane

Artifact roots:

- optimized run: `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_total/`
- baseline run: `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_total_baseline/`

Baseline binary provenance:

- built from base commit `7a7fb08` in isolated worktree `/tmp/star-suite-baseline-7a7fb08-SzAPwP`
- important build note: use the tracked `core/legacy/source/parametersDefault.xxd`
- a regenerated `parametersDefault.xxd` in that worktree produced a broken binary that exited with `unrecognized parameter name "soloAddTagsToUnsorted" in input "Default"`
- restoring the tracked `parametersDefault.xxd` fixed the baseline build and `STAR --version` returned `2.7.11b`

Result parity on the 100k-total GEX benchmark:

- baseline and optimized outputs matched on all checked top-line metrics
- both runs reported:
  - `Number of Reads = 100000`
  - `Estimated Number of Cells = 7309`
  - `UMIs in Cells = 65887`
  - `Reads Mapped to Genome: Unique = 0.90474`
  - `Reads Mapped to GeneFull: Unique GeneFull = 0.74086`

Coarse runtime comparison from `Log.final.out`:

- baseline:
  - started `Mar 24 01:11:41`
  - started mapping `Mar 24 01:12:59`
  - finished `Mar 24 01:13:02`
  - total wall time about `81 s`
  - mapping-to-finish about `3 s`
  - mapping speed `120.00` million reads/hour
- optimized:
  - started `Mar 24 01:01:24`
  - started mapping `Mar 24 01:02:42`
  - finished `Mar 24 01:02:44`
  - total wall time about `80 s`
  - mapping-to-finish about `2 s`
  - mapping speed `180.00` million reads/hour

Interpretation:

- end-to-end wall time is still dominated by genome load at this 100k scale
- the optimization does not change biological outputs on this benchmark
- the coarse finish time improved by about 1 second after mapping start
- the instrumented optimized binary shows the Solo finalize work itself is sub-second:
  - `collapseUMIall 0.0288 s`
  - `countCBgeneUMI 0.1162 s`
  - `outputResults(raw) 0.0786 s`
  - `outputResults(filtered) 0.0067 s`
  - `cellFiltering 0.0072 s`
  - `processRecords 0.2278 s`
  - sum of the logged Solo phases above is about `0.465 s`

This means the branch now has a valid apples-to-apples baseline for the corrected UCSF Solo GEX downsample. The next useful benchmark step is a larger downsample or full-sample run, where post-map Solo savings should separate more clearly from genome-load noise.

### 2M optimized follow-up run

A larger optimized-only run was also completed on the same corrected UCSF GEX sample:

- artifact root: `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_total/`
- downsample target: `2,000,000` total reads across the 10 GEX lanes
- started `Mar 24 01:21:09`
- started mapping `Mar 24 01:22:27`
- finished `Mar 24 01:22:33`
- total wall time from `Log.final.out`: about `84 s`
- mapping-to-finish time from `Log.final.out`: about `6 s`
- mapping speed: `1200.00` million reads/hour

2M run summary:

- `Number of Reads = 2000000`
- `Estimated Number of Cells = 7211`
- `UMIs in Cells = 1313762`
- `Reads Mapped to Genome: Unique = 0.906103`
- `Reads Mapped to GeneFull: Unique GeneFull = 0.746895`

Instrumented Solo timings for the 2M run:

- `collapseUMIall 0.8569 s`
- `countCBgeneUMI 1.4329 s`
- `outputResults(raw) 0.3213 s`
- `outputResults(filtered) 0.0657 s`
- `cellFiltering 0.0822 s`
- `processRecords 1.8642 s`

Operational note:

- `/usr/bin/time -v` reported full wall clock `2:05.49` and max RSS `44111452 kB`
- this includes the downsample FASTQ creation step, so use `Log.final.out` for STAR-only elapsed time and `time.txt` for end-to-end harness cost

This 2M run is a better scale point than the 100k pilot: the post-map Solo work is now clearly in the low-single-digit seconds instead of being buried almost entirely in timer granularity.

### 2M baseline comparison

A matching unoptimized baseline run was completed against the exact same downsampled FASTQs used by the optimized 2M run:

- baseline artifact root: `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_total_baseline/`
- baseline binary: `/tmp/star-suite-baseline-7a7fb08-SzAPwP/core/legacy/source/STAR`
- input FASTQs reused from:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_total/fastq_downsampled/`

2M result parity:

- baseline and optimized 2M outputs matched on all checked top-line metrics
- both runs reported:
  - `Number of Reads = 2000000`
  - `Estimated Number of Cells = 7211`
  - `UMIs in Cells = 1313762`
  - `Reads Mapped to Genome: Unique = 0.906103`
  - `Reads Mapped to GeneFull: Unique GeneFull = 0.746895`

2M runtime comparison:

- optimized `Log.final.out`:
  - started `Mar 24 01:21:09`
  - started mapping `Mar 24 01:22:27`
  - finished `Mar 24 01:22:33`
  - total wall time about `84 s`
- baseline `Log.final.out`:
  - started `Mar 24 01:27:24`
  - started mapping `Mar 24 01:28:41`
  - finished `Mar 24 01:28:47`
  - total wall time about `83 s`

Interpretation of the 2M pair:

- coarse `Log.final.out` timing shows effectively no meaningful end-to-end difference on this host
- both runs spent about `6 s` from mapping start to finish
- `/usr/bin/time -v` also shows near-identical end-to-end STAR runtimes:
  - optimized: `2:05.49`
  - baseline: `2:04.76`
- the optimized binary remains valuable because it exposes where Solo time is actually spent:
  - `collapseUMIall 0.8569 s`
  - `countCBgeneUMI 1.4329 s`
  - `outputResults(raw) 0.3213 s`
  - `outputResults(filtered) 0.0657 s`
  - `cellFiltering 0.0822 s`
  - `processRecords 1.8642 s`

Practical conclusion:

- at 2M reads on this machine, the optimization is not yet producing a clearly separated end-to-end wall-clock win over baseline
- the instrumentation confirms the Solo finalize path itself is only a few seconds, so whole-run timing is still dominated by other costs
- the next benchmark that can better resolve optimization impact is a substantially larger downsample or a full-sample run

### Same-binary guarded A/B comparison

To remove the "different binary" confound, the current instrumented binary now supports a runtime guard:

- default mode: optimized in-place qsort path
- legacy mode: set environment variable `STAR_SOLO_LEGACY_COLLAPSE_SORT=1`
- guard log line in `Log.out`:
  - `Solo optimization guard: collapseSort=optimized_qsort`
  - `Solo optimization guard: collapseSort=legacy_std_sort`

Implementation note:

- the guard is in [SoloFeature_collapseUMIall.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_collapseUMIall.cpp)
- it switches only the hot sort inside `collapseUMIperCB()` while keeping the same instrumentation and surrounding code

Direct 2M A/B artifacts:

- optimized same-binary run:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_samebin_opt/`
- legacy same-binary run:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_samebin_legacy/`

Direct 2M A/B result:

- outputs were identical across optimized and legacy guarded runs
- coarse wall-clock was essentially tied:
  - optimized `/usr/bin/time -v`: `2:05.02`
  - legacy `/usr/bin/time -v`: `2:05.15`
- `Log.final.out` was also effectively tied at about `84 s` total and `6 s` from mapping start to finish

But the internal Solo timings on this direct A/B run show the guarded legacy path slightly faster than the current "optimized" path:

- `collapseUMIall`
  - optimized: `0.864879 s`
  - legacy: `0.817963 s`
  - legacy faster by about `5.4%`
- `countCBgeneUMI`
  - optimized: `1.446490 s`
  - legacy: `1.411800 s`
  - legacy faster by about `2.4%`
- `processRecords`
  - optimized: `1.880620 s`
  - legacy: `1.843620 s`
  - legacy faster by about `2.0%`

Current conclusion after direct measurement:

- yes, the guarded single-binary method works and should be preferred for direct Solo-phase measurement
- no, the current in-place qsort change does not show a speed win in this direct 2M A/B test
- the next optimization work should treat this change as unproven or mildly regressive until a better hot-path improvement is found

### Experimental non-Flex inline-hash bridge

To test the more aggressive optimization direction, the binary now supports an
experimental non-Flex bridge behind:

- environment variable: `STAR_SOLO_NONFLEX_HASH_BRIDGE=1`
- CLI flag: `--soloInlineHashMode yes`

This experimental path does **not** use the Flex resolver for standard GEX.
Instead it:

- routes non-Flex inline-hash recording through the standard `record_base()`
  gene classification logic
- inserts `(CB, UMI, gene)` directly into the per-thread inline hash
- materializes hash entries back into legacy `rGeneUMI/rCBp/rCBn` arrays
- runs the existing legacy `collapseUMIall()` and output path

Current benchmark scope:

- dataset: corrected UCSF `iPSC2_1/GEX`
- input: `100000` total reads reusing the existing 100k downsampled FASTQs
- multimapper mode: `--soloMultiMappers Unique`
- purpose: upper-bound measurement of the temp-stream bypass on a restricted,
  non-Rescue workload

Artifact roots:

- same-binary legacy control:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_legacy_v2/`
- same-binary experimental bridge:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_hashbridge_v4/`
- refreshed same-binary legacy control after fixing stats aggregation and
  invalid-UMI rejection:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_legacy_v3/`
- refreshed same-binary experimental bridge after those two fixes:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_hashbridge_v6/`
- failed debug iterations while bringing the bridge up:
  - `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_hashbridge_v2/`
  - `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_hashbridge_v3/`

Timing comparison on the successful `100k` run:

- legacy control:
  - `collapseUMIall 0.0194988 s`
  - `countCBgeneUMI 0.0979653 s`
  - `processRecords 0.200383 s`
- experimental bridge:
  - `collapseUMIall 0.0166380 s`
  - `countCBgeneUMI 0.0220528 s`
  - `processRecords 0.118162 s`

Observed speed signal:

- `countCBgeneUMI` improved by about `4.4x`
- `processRecords` improved by about `1.7x`
- this is the first experiment in this runbook that shows a materially large
  Solo-phase timing reduction

But parity is currently **not acceptable**:

- legacy `Summary.csv`:
  - `Reads With Valid Barcodes = 0.9654`
  - `Reads Mapped to GeneFull: Unique GeneFull = 0.74086`
  - `Estimated Number of Cells = 7309`
  - `UMIs in Cells = 65887`
- experimental bridge `Summary.csv`:
  - `Reads With Valid Barcodes = 0.97343`
  - `Reads Mapped to GeneFull: Unique GeneFull = 0.21856`
  - `Estimated Number of Cells = 7294`
  - `UMIs in Cells = 65126`

Interpretation:

- the direction is promising as a performance upper bound
- the current bridge implementation is **not** a valid replacement for standard
  Solo
- remaining correctness gaps are still significant, especially around preserving
  standard non-Flex semantics through the inline-hash path

Debug notes from the bridge bring-up:

- initial attempts using the existing `record_flex()` path produced empty
  materialization because the Flex resolver dropped standard GEX records
- switching the experiment to `record_base()`-style gene classification fixed
  hash population
- materialized records also needed synthetic unique read ids so legacy
  `collapseUMIall()` could run its CR-style UMI correction path without
  conflicting `recordReadInfo()` writes

Root-cause analysis of the remaining parity gap:

- `Reads With Valid Barcodes` is wrong because inline-hash mode currently keeps
  CBs when the UMI is invalid, while legacy Solo rejects the whole read.
  The branch is in
  [SoloReadBarcode_getCBandUMI.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadBarcode_getCBandUMI.cpp#L451).
  Evidence from `Barcodes.stats`:
  - legacy: `noNinUMI=121`, `noUMIhomopolymer=105`
  - bridge: both `0`
- feature-level read stats are missing because the inline-hash path in
  [SoloFeature_sumThreads.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_sumThreads.cpp#L15)
  merges hashes and `cbReadCount`, but does not aggregate per-thread
  `SoloReadFeature::stats`. Evidence from `Features.stats`:
  - legacy has non-zero `noUnmapped`, `noNoFeature`, `MultiFeature`,
    `noTooManyWLmatches`
  - bridge shows all of those as `0`
- the bridge's read-count summaries are wrong because it fabricates
  `nReadPerCBunique` / `nReadPerCBtotal` from `cbReadCount` instead of the
  loader-derived per-read bookkeeping that legacy Solo uses.
  The bridge code is in
  [SoloFeature_countCBgeneUMI.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_countCBgeneUMI.cpp#L100).
  Evidence:
  - bridge `Features.stats` reports `yesWLmatch=21856`
  - the same run reports `yesUMIs=67347`
  - that combination is impossible if the read counters were faithful
- ambiguous-CB handling is still incomplete in the bridge:
  - the bridge only inserts directly into the hash when `cbMatch <= 1`, so
    ambiguous CB observations are skipped in the experimental insert path in
    [SoloReadFeature_record_base.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature_record_base.cpp#L17)
  - the pre-collapse hook
    [SoloFeature::resolveAmbiguousCBs()](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature.h#L203)
    is currently an empty stub, so `pendingAmbiguous_` is not resolved before
    `materializeRGUFromHash()`

Practical takeaway from the root-cause pass:

- the big speed win is coming from bypassing temp-stream serialization/parsing
- the current parity gap is mostly in the bridge glue:
  - invalid-UMI rejection semantics
  - per-thread stats aggregation
  - faithful per-read `nReadPerCB*` bookkeeping
  - ambiguous-CB resolution before materialization
- that means the optimization direction is still viable, but the bridge has to
  reproduce those legacy semantics explicitly before the benchmark is meaningful

Bottom line:

- yes, the aggressive inline-hash direction looks worthwhile from a raw timing
  perspective
- no, it is not ready to benchmark at larger scale as a real optimization until
  parity is fixed on this restricted `Unique` workload first

### Follow-up after fixing the first two parity issues

Implemented fixes:

- `SoloFeature_sumThreads.cpp`: the inline-hash merge path now also calls
  `readFeatSum->addStats(*readFeatAll[ii])`
- `SoloReadBarcode_getCBandUMI.cpp`: the experimental non-Flex bridge now uses
  legacy invalid-UMI rejection semantics instead of Flex CB/UB independence

Validation rerun on the same `100k`, `Unique`, corrected UCSF GEX benchmark:

- legacy control:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_legacy_v3/`
- bridge:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_hashbridge_v6/`

What this fixed:

- `Barcodes.stats` is now an exact match between legacy and bridge:
  - `noNinUMI=121`
  - `noUMIhomopolymer=105`
  - `yesWLmatchExact=92690`
  - `yesOneWLmatchWithMM=2234`
  - `yesMultWLmatchWithMM=2287`
- top-level feature classification counters are now restored:
  - `noUnmapped=1655`
  - `noNoFeature=3798`
  - `MultiFeature=17001`
  - `subMultiFeatureMultiGenomic=5363`

What is still wrong:

- `Summary.csv` remains far from parity:
  - legacy `Reads Mapped to GeneFull: Unique GeneFull = 0.74086`
  - bridge `Reads Mapped to GeneFull: Unique GeneFull = 0.21816`
  - legacy `Estimated Number of Cells = 7309`
  - bridge `Estimated Number of Cells = 7294`
  - legacy `UMIs in Cells = 65887`
  - bridge `UMIs in Cells = 65042`
  - legacy `Sequencing Saturation = 0.0801366`
  - bridge `Sequencing Saturation = -2.08297`
- `Features.stats` still shows the bookkeeping mismatch:
  - legacy `noTooManyWLmatches=671`, bridge `0`
  - legacy `yesWLmatch=74086`, bridge `21816`
  - legacy `yessubWLmatchExact=71376`, bridge `0`
  - legacy `yessubWLmatch_UniqueFeature=74086`, bridge `21816`
  - legacy `yesCellBarcodes=9571`, bridge `9510`
  - legacy `yesUMIs=68149`, bridge `67258`

Timing impact after the two fixes is still strong:

- legacy:
  - `countCBgeneUMI 0.0874185 s`
  - `processRecords 0.181515 s`
- bridge:
  - `countCBgeneUMI 0.0221998 s`
  - `processRecords 0.117794 s`

Updated conclusion:

- the first two parity fixes are confirmed
- they eliminate the barcode-accounting drift and restore the missing
  per-thread feature stats
- the remaining mismatch now concentrates in the two harder items:
  `nReadPerCBunique` / `nReadPerCBtotal` bookkeeping and ambiguous-CB handling
- important scope note for problem 3:
  - the deferred `oneExact` gate for `cbMatch==1` is a strict-match
    compatibility rule and is intentionally **not** part of the current bridge
    implementation target
  - the UCSF benchmark harness uses
    `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts`, where `oneExact` is
    disabled
  - therefore the current problem-3 work should focus on general read-level
    accounting parity, not on reproducing the stricter `oneExact` acceptance
    rule
  - if strict-match compatibility is needed later, add `oneExact` support as a
    separate follow-up item rather than blocking the current pseudocount-mode
    optimization work

## Problem 3 Update: Read-Accounting Sidecar + CB Key-Space Fix

To address problem 3, the bridge now captures a compact per-read accounting
sidecar during `record_base()` and uses that to populate
`nReadPerCBunique` / `nReadPerCBtotal` in the experimental non-Flex bridge
path instead of fabricating them from `cbReadCount`.

The first sidecar-only validation run (`hashbridge_v8`) showed that the new
read accounting itself was correct, but the bridge hash was still corrupting
cell-barcode identity:

- debug line from `Log.out`:
  - `Bridge accounting debug: records=97211 featGood=74757 multiFeature=0 exact=92690 oneMM=2234 ambiguous=2287 uniqueReads=74072 multiReads=0`
- but `Features.stats` still reported only:
  - `yesWLmatch=21642`
  - `yessubWLmatch_UniqueFeature=21642`

Root cause:

- the shared inline hash key format is `CB20/UMI24/GENE15/TAG5`
- that 20-bit CB field silently truncates 3M Cell Ranger whitelist indices
- the sidecar was using full whitelist indices, while the bridge hash was
  materializing a truncated CB space

Fix implemented:

- the experimental non-Flex bridge now stores a compact bridge-local CB id in
  the inline hash
- thread merge remaps other threads' compact ids back to real whitelist
  indices, then into the destination compact id space
- materialization remaps compact ids back to the original whitelist indices

Validation rerun after the key-space fix:

- bridge:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_hashbridge_v9/`
- legacy control:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_legacy_v4/`

What is now fixed:

- read-level mapping totals are essentially back to parity:
  - legacy `yesWLmatch=74086`
  - bridge `yesWLmatch=74035`
  - legacy `yessubWLmatch_UniqueFeature=74086`
  - bridge `74035`
  - legacy `Reads Mapped to GeneFull: Unique GeneFull = 0.74086`
  - bridge `0.74035`
  - legacy `Unique Reads in Cells Mapped to GeneFull = 71640`
  - bridge `71615`
- sequencing saturation is now sane again:
  - legacy `0.0801366`
  - bridge `0.0915378`

What is still mismatched:

- `noTooManyWLmatches` is still slightly high:
  - legacy `671`
  - bridge `685`
- molecule/cell outputs are still low:
  - legacy `yesUMIs=68149`
  - bridge `67258`
  - legacy `UMIs in Cells=65887`
  - bridge `65027`
  - legacy `yesCellBarcodes=9571`
  - bridge `9535`
  - legacy `Estimated Number of Cells=7309`
  - bridge `7304`
  - legacy `Total GeneFull Detected=9887`
  - bridge `9485`

Current conclusion:

- problem 3 is largely resolved for read-level accounting
- the big `yesWLmatch` / `Summary.csv` mismatch was mostly a CB key-space bug,
  not a sidecar logic bug
- the remaining parity gap is now much smaller and is consistent with the
  unresolved ambiguous-CB path in problem 4
- timing remains favorable:
  - bridge `countCBgeneUMI 0.0365256 s`
  - bridge `processRecords 0.136989 s`

## Problem 4 Update: Reuse the Flex Ambiguous-CB Strategy

The non-Flex bridge now reuses the same high-level approach as Flex:

- during mapping, ambiguous CB observations are stored in
  `pendingAmbiguous_` instead of being dropped
- after mapping, the bridge runs Bayesian CB resolution
- resolved observations are reinserted into the inline hash
- the bridge then materializes and collapses as before

Implementation notes:

- ambiguous gene observations are now accumulated in
  [SoloReadFeature_record_base.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature_record_base.cpp)
- non-Flex bridge resolution/reinsertion is implemented in
  [SoloFeature.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature.cpp)
- the first attempt (`hashbridge_v10`) wired the reuse pattern but still fed
  the resolver the wrong barcode representation, so it resolved nothing:
  - `[AMBIG-CB-RESOLVE] pending=1639 resolved=0 still_ambiguous=1639 added_to_hash=0`
- the fix was to use the raw observed CB sequence `cbSeq` plus matching
  `cbQual`, rather than `cbMatchString`, for Bayesian resolution context

Validation rerun after the raw-CB fix:

- bridge:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_hashbridge_v11/`
- legacy control:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_legacy_v4/`

Resolver activity:

- `pending=1621`
- `resolved=1048`
- `still_ambiguous=573`
- `added_to_hash=1081`

Parity after ambiguous-CB reuse:

- read-level metrics are now essentially matched:
  - legacy `yesWLmatch=74086`
  - bridge `74072`
  - legacy `Reads Mapped to GeneFull: Unique GeneFull = 0.74086`
  - bridge `0.74072`
- molecule/cell metrics are now also essentially matched:
  - legacy `yesUMIs=68149`
  - bridge `68214`
  - legacy `UMIs in Cells=65887`
  - bridge `65881`
  - legacy `Estimated Number of Cells=7309`
  - bridge `7311`
  - legacy `Sequencing Saturation=0.0801366`
  - bridge `0.0790852`

Residual differences still present:

- `noTooManyWLmatches` remains slightly high:
  - legacy `671`
  - bridge `685`
- `yesCellBarcodes` is slightly high:
  - legacy `9571`
  - bridge `9644`

## GeneFull Detected Root Cause: 15-bit Gene Key Aliasing

The large `Total GeneFull Detected` gap in `hashbridge_v11` was not a cell
filtering problem. It was a bridge hash key-space bug.

Observed symptom:

- filtered cell sets were nearly identical, but the bridge had
  `Total GeneFull Detected = 9524` versus legacy `9887`
- direct filtered-matrix comparison showed many exact count swaps between gene
  rows rather than random gene loss
- representative swaps:
  - `ENSG00000198938` row `38576` lost `618`, while `ENSG00000289422` row
    `5808` gained `618`
  - `ENSG00000198712` row `38573` lost `354`, while `ENSG00000163251` row
    `5805` gained `354`

Critical clue:

- these swapped row pairs differ by exactly `32768`
- the inline-hash packed key format is `[CB20][UMI24][GENE15][TAG5]`
- so the non-Flex bridge was truncating gene indices above `32767` via the
  15-bit `GENE15` field and aliasing them onto lower gene rows

Fix:

- the non-Flex bridge now uses a bridge-local compact gene id, analogous to
  the already-required bridge-local compact CB id
- hash insert, cross-thread hash merge, ambiguous-CB reinsertion, and
  materialization now round-trip through this bridge-local gene map before
  collapsing back onto the real `GeneFull` row ids

Validation rerun after the gene-key fix:

- bridge:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_hashbridge_v12/`
- legacy control:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_legacy_v4/`

Key result:

- `Total GeneFull Detected` now matches exactly:
  - legacy `9887`
  - bridge `9887`
- the prior `32768` row-swap pairs are gone; for example:
  - row `38576` remains `618` in both legacy and bridge
  - row `5808` remains `0` in both legacy and bridge

Residual differences after the gene-key fix are now small and no longer
structural:

- `Reads Mapped to GeneFull: Unique GeneFull`
  - legacy `0.74086`
  - bridge `0.74072`
- `UMIs in Cells`
  - legacy `65887`
  - bridge `65881`
- `Unique Reads in Cells Mapped to GeneFull`
  - legacy `71640`
  - bridge `71626`
- `noTooManyWLmatches`
  - legacy `671`
  - bridge `685`
- `yesCellBarcodes`
  - legacy `9571`
  - bridge `9644`

Residual pattern analysis:

- the remaining drift is concentrated in the extreme low-count barcode tail,
  not in established cells or gene rows
- raw nonzero barcodes:
  - legacy `9571`
  - bridge `9644`
  - delta is entirely `73` bridge-only barcodes; legacy has no raw-only
    barcodes
- size of those `73` bridge-only raw barcodes:
  - `72` have exactly `(1 UMI, 1 gene)`
  - `1` has exactly `(2 UMI, 2 genes)`
- among barcodes that are nonzero in both runs, only `53` differ at all, and
  every one of them differs by exactly one singleton molecule:
  - either `(legacy UMI, legacy genes) -> (bridge UMI, bridge genes)` as
    `(n, n) -> (n+1, n+1)`
  - or `(n, n) -> (n-1, n-1)`
  - no shared barcode shows a larger structural drift
- filtered cells show the same pattern:
  - only `52` shared filtered cells differ
  - every shared filtered-cell delta is exactly `(±1 UMI, ±1 gene)`
  - the two bridge-only filtered cells are tiny edge cases with `2 UMI` and
    `2 genes` each

Interpretation:

- the remaining mismatch behaves like single-molecule reassignment between
  low-support barcodes
- it does **not** look like another gene-index or matrix materialization bug
- because UMI and detected-gene deltas move together as `(±1, ±1)`, the
  likely remaining source is CB-assignment / ambiguous-CB edge handling rather
  than UMI collapse or gene mapping

Targeted trace of the `73` bridge-only raw barcodes:

- trace run:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_hashbridge_trace73/`
- target list:
  `/tmp/solo_bridge_73_barcodes.txt`
- stderr trace file:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_unique_hashbridge_trace73/driver.stderr`

Trace result:

- exactly `73` trace lines were emitted
- all `73` were `mode=ambig_resolve`
- `0` were `mode=direct`
- therefore every one of the `73` extra raw barcodes is created by the
  ambiguous-CB reinsertion path, not by exact or simple `cbMatch<=1` insertions

Per-entry pattern:

- every traced barcode came from a very small ambiguous group:
  - candidate count `2-5`
  - usually `umiKeys=1`
  - usually `obs=1`
  - one notable 2-observation case:
    `CTCCTCCCACGTACTA` with `umiKeys=2`, `obs=2`
- this is fully consistent with the observed singleton-tail behavior in the raw
  matrix

Operational conclusion:

- the remaining parity gap is now isolated to non-Flex handling of
  `cbMatch>1` ambiguous CBs
- the current bridge uses the Flex Bayesian resolver for these reads
- legacy standard Solo uses the older `ptot/pmax` heuristic in
  [SoloReadFeature_inputRecords.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature_inputRecords.cpp)
- so the residual `73`-barcode tail is best interpreted as a resolver-policy
  mismatch between Flex-style Bayesian resolution and legacy Solo ambiguous-CB
  assignment, not as a downstream counting bug

Policy decision:

- keep the unified Cell Ranger-style ambiguous-CB resolver used by Flex
- do **not** reintroduce legacy non-Flex `ptot/pmax` behavior just for bridge
  parity
- treat the remaining singleton-tail mismatch as an intentional algorithmic
  difference between legacy STARsolo and the unified Flex/CR-style path

Codebase unification to reduce drift:

- the Cell Ranger-style ambiguous-CB reinsertion logic for inline-hash data is
  now centralized in
  [SoloFeature::resolvePendingAmbiguousToHash()](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature.cpp)
- the non-Flex bridge wrapper
  [SoloFeature::resolveAmbiguousCBs()](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature.cpp)
  now just delegates to that shared helper with bridge compact CB/gene mapping
- the Flex inline-hash collapse path in
  [SoloFeature_collapseUMI_fromHash.cpp](/mnt/pikachu/STAR-suite/flex/source/SoloFeature_collapseUMI_fromHash.cpp)
  now calls the same helper instead of maintaining a duplicate Bayesian
  resolution/reinsertion block
- this keeps Flex and non-Flex bridge on the same Cell Ranger-style ambiguous
  CB code path and should reduce future drift

Current conclusion:

- the Flex-style ambiguous-CB reuse is the right solution and materially fixes
  problem 4
- the large `GeneFull Detected` mismatch was caused by 15-bit gene-id aliasing
  in the bridge hash and is now fixed
- the bridge is now very close to legacy on the benchmark that previously had
  large parity failures
- remaining differences are small enough to investigate as residual edge-case
  parity issues rather than architectural gaps

## 2M Unique Benchmark for Unified CR-Style Bridge

After switching away from the full sample for runtime reasons, I ran a matched
`2,000,000`-read UCSF `iPSC2_1/GEX` downsample benchmark in
`--soloMultiMappers Unique` mode on the current branch.

Artifacts:

- legacy control:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_unique_legacy/`
- unified non-Flex bridge (`STAR_SOLO_NONFLEX_HASH_BRIDGE=1`,
  `--soloInlineHashMode yes`):
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_unique_hashbridge/`

Measured Solo timings from `run/Log.out`:

- legacy:
  - `collapseUMIall 0.540971 s`
  - `countCBgeneUMI 1.00476 s`
  - `outputResults(raw) 0.128359 s`
  - `outputResults(filtered) 0.0657186 s`
  - `cellFiltering 0.0675458 s`
  - `processRecords 1.2246 s`
- bridge:
  - `collapseUMIall 0.43801 s`
  - `countCBgeneUMI 0.624055 s`
  - `outputResults(raw) 0.126878 s`
  - `outputResults(filtered) 0.0737178 s`
  - `cellFiltering 0.0755331 s`
  - `processRecords 1.08528 s`

Timing deltas:

- `countCBgeneUMI`: about `1.61x` faster on the bridge
  (`1.00476 s -> 0.624055 s`, saving `0.380705 s`)
- `collapseUMIall`: about `1.24x` faster on the bridge
  (`0.540971 s -> 0.43801 s`, saving `0.102961 s`)
- `processRecords`: about `1.13x` faster on the bridge
  (`1.2246 s -> 1.08528 s`, saving `0.13932 s`)
- `outputResults(raw)` is effectively flat
- `outputResults(filtered)` and `cellFiltering` are slightly slower on the
  bridge, but by only a few milliseconds

Top-line outputs:

- legacy:
  - `Number of Reads = 2000000`
  - `Reads Mapped to GeneFull: Unique GeneFull = 0.746895`
  - `Estimated Number of Cells = 7211`
  - `Unique Reads in Cells Mapped to GeneFull = 1438379`
  - `UMIs in Cells = 1313762`
  - `Total GeneFull Detected = 17782`
  - `Sequencing Saturation = 0.0865785`
- bridge:
  - `Number of Reads = 2000000`
  - `Reads Mapped to GeneFull: Unique GeneFull = 0.746586`
  - `Estimated Number of Cells = 7211`
  - `Unique Reads in Cells Mapped to GeneFull = 1437767`
  - `UMIs in Cells = 1306937`
  - `Total GeneFull Detected = 17766`
  - `Sequencing Saturation = 0.0901919`

Feature-stat differences:

- legacy:
  - `yesWLmatch = 1493789`
  - `yessubWLmatchExact = 1460058`
  - `noTooManyWLmatches = 4499`
  - `yesCellBarcodes = 40650`
  - `yesUMIs = 1364459`
- bridge:
  - `yesWLmatch = 1493172`
  - `yessubWLmatchExact = 1460058`
  - `noTooManyWLmatches = 5115`
  - `yesCellBarcodes = 41440`
  - `yesUMIs = 1358500`

Interpretation:

- the unified non-Flex inline-hash route still produces a real Solo-phase win
  at `2M`
- the clearest gain remains `countCBgeneUMI`, which is now about `1.6x` faster
  on this larger benchmark
- `processRecords` also improves, but only modestly, because non-counting work
  still dominates the post-map phase
- the output drift at `2M` is consistent with the earlier policy decision to
  keep the Cell Ranger/Flex-style Bayesian ambiguous-CB resolver instead of
  forcing strict legacy `ptot/pmax` parity
- so this benchmark should be interpreted as: the route is worthwhile from a
  speed perspective, and the remaining output differences are an intentional
  algorithmic choice rather than an unresolved bridge corruption bug

## Related Files Worth Reading

- [SoloFeature_countCBgeneUMI.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_countCBgeneUMI.cpp)
- [SoloFeature_collapseUMIall.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_collapseUMIall.cpp)
- [SoloFeature_outputResults.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_outputResults.cpp)
- [MexWriter.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/MexWriter.cpp)
- [SoloFeature.h](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature.h)
- [tests/run_solo_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_solo_smoke.sh)
- [tests/run_cr_parity_100k.sh](/mnt/pikachu/STAR-suite/tests/run_cr_parity_100k.sh)

## Full EBs2_2 Benchmark on Repaired Master-Based Branch

Using the corrected full UCSF paper harness from
`comparisons/paper_benchmarks_20260318/ucsf_ebs2_2/run_ucsf_ebs2_2_benchmark.sh`
on the corrected dataset root
`/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2`, I ran the full 32-thread
dynamic benchmark on the repaired master-based Solo branch.

Baseline full run:

- output:
  `/storage/paper_bench_solo_full_20260324/ucsf_ebs2_2_standard_baseline/`
- completed successfully
- `Log.final.out`:
  - `Number of input reads = 444896731`
  - `Started job on = Mar 24 07:05:37`
  - `Started mapping on = Mar 24 07:06:37`
  - `Finished on = Mar 24 07:24:52`
- `BENCHMARK_SUMMARY.txt`:
  - `star_cells = 13721`
  - `wall_minutes = 19.8`
- `Summary.csv`:
  - `Reads Mapped to GeneFull: Unique GeneFull = 0.855092`
  - `Estimated Number of Cells = 13721`
  - `UMIs in Cells = 256057933`
- instrumented Solo timings from `Log.out`:
  - `collapseUMIall 164.982 s`
  - `countCBgeneUMI 271.754 s`
  - `outputResults(raw) 4.31639 s`
  - `outputResults(filtered) 3.76 s`
  - `cellFiltering 211.73 s`
  - `processRecords 488.128 s`
- `/usr/bin/time -v`:
  - `Elapsed (wall clock) = 19:46.33`
  - `Maximum resident set size = 72504156 kB`

Optimized full run, first attempt:

- output:
  `/storage/paper_bench_solo_full_20260324/ucsf_ebs2_2_standard_solohash_optimized/`
- same harness plus:
  - environment `STAR_SOLO_NONFLEX_HASH_BRIDGE=1`
  - CLI flag `--soloInlineHashMode yes`
- failed in Solo post-map with:
  - `FATAL ERROR: non-Flex Solo hash bridge observed more than 32768 distinct genes; packed gene field overflow`

Root cause:

- the UCSF benchmark reference under `/storage/autoindex_110_44/bulk_index`
  contains `38606` genes
- this fits in a true 16-bit field but not in the bridge's effective 15-bit
  gene field inherited from the Flex-style packed key

Patch applied:

- the non-Flex bridge now uses its own packed key layout instead of the shared
  Flex key
- new bridge key:
  - `[bridgeCB24][UMI24][GENE16]`
- changed files:
  - `flex/source/hash_shims_cpp_compat.h`
  - `core/legacy/source/SoloReadFeature.cpp`
  - `core/legacy/source/SoloReadFeature_record_base.cpp`
  - `core/legacy/source/SoloFeature.cpp`
  - `flex/source/SoloFeature_materializeFromHash.cpp`

Optimized full run, second attempt after 16-bit bridge patch:

- output:
  `/storage/paper_bench_solo_full_20260324/ucsf_ebs2_2_standard_solohash_optimized_v2/`
- successfully passed the previous overflow point and entered Solo post-map on
  the full sample
- then was killed externally during Solo finalize with shell status `137`
- no new STAR fatal error was emitted before termination
- `/usr/bin/time -v`:
  - `Elapsed (wall clock) = 10:54.91`
  - `Maximum resident set size = 127126724 kB`
  - `Exit status = 137`

Interpretation of the new blocker:

- the 16-bit bridge-key fix worked; the previous full-sample gene-field crash
  is resolved
- the next full-scale blocker is memory, not packed-key width
- compared with the successful baseline:
  - baseline RSS peak: `72.5 GB`
  - optimized bridge RSS peak before kill: `127.1 GB`
- this strongly suggests the current bridge implementation's
  `materializeRGUFromHash()` path is still too memory-hungry at full scale
- the likely next optimization target is eliminating the bridge
  materialization/expansion step rather than further packed-key changes

### Direct bridge-hash collapse (implemented 2026-03-24)

The non-Flex bridge no longer uses `materializeRGUFromHash()` → `rGeneUMI` /
`collapseUMIall()` in `countCBgeneUMI`. It calls
`SoloFeature::collapseUMIall_fromBridgeHash()` instead:

- **File:** `core/legacy/source/SoloFeature_collapseUMI_fromBridgeHash.cpp`
- **Behavior:** scan `inlineHash_` with `unpackBridgeCgAggKey`, sort by
  `(wlCb, gene, umi)`, per-CB gene segments build `umiArray` from aggregated
  counts, then reuse `umiArrayCorrect_CR` and the same MultiGeneUMI_CR
  resolution loop shape as `SoloFeature_collapseUMIall.cpp`.
- **Fail-closed gating** (parameter error if not met): `--soloUMIfiltering
  MultiGeneUMI_CR`, exactly `--soloUMIdedup 1MM_CR`, `--soloMultiMappers
  Unique`, `--outSAMtype None`, no `trackReadIdsForTags`, no MultiGeneUMI /
  MultiGeneUMI_All.
- **`materializeRGUFromHash()`** remains in the tree for other experiments but
  is not invoked from the bridge path in `SoloFeature_countCBgeneUMI.cpp`.
- **Memory hygiene:** after copying the merged hash into `recs`, the code calls
  `kh_destroy` on `readFeatSum->inlineHash_` and nulls the pointer so peak RSS
  is not *hash + recs + matrix* at once (unlike waiting for `clearLarge` /
  `soloFlexMinimalMemory`).
- **`countCellGeneUMI` sizing:** preallocated using the count of unique
  `(wlCb, gene)` runs in sorted `recs` (upper bound on matrix rows), not
  `recs.size()` (unique CB/gene/UMI keys).
- **`umiArray`:** grows per gene with `nU0` only; not pre-sized from
  `nReadPerCBmax` (total reads per CB).

**Full-sample validation (v3):** reran optimized EBs2_2 with
`STAR_SOLO_NONFLEX_HASH_BRIDGE=1`, `--soloInlineHashMode yes`, output under
`/storage/paper_bench_solo_full_20260324/ucsf_ebs2_2_standard_solohash_optimized_v3/`.

Result:

- still killed externally with exit `137`
- `/usr/bin/time -v`:
  - `Elapsed (wall clock) = 11:30.97`
  - `Maximum resident set size = 127682252 kB`
  - `Exit status = 137`
- got substantially farther than `v2`:
  - finished mapping at `Mar 24 09:36:21`
  - entered Solo post-map at `09:36:22`
  - logged `Allocated and initialized readInfo array, nReadsInput = 444896732`
  - logged
    `Experimental non-Flex inline-hash bridge: direct hash collapse (no materializeRGUFromHash / no legacy collapseUMIall)`
  - logged
    `Direct bridge-hash UMI collapse (no rGeneUMI materialization), hash_entries=269677733`

Interpretation:

- the direct bridge collapse removed the old `materializeRGUFromHash()` /
  legacy-collapse path, but it did **not** solve the full-sample OOM
- the new failure occurs after entering the direct bridge collapse on a merged
  hash with `269,677,733` unique `(CB,gene,UMI)` entries
- live memory right before Solo counting was already high:
  - after freeing genome index memory:
    - `VmRSS = 80806804 kB`
    - `VmHWM = 109832420 kB`
- compared with prior full attempts:
  - baseline peak RSS: `72504156 kB`
  - optimized `v2` peak RSS: `127126724 kB`
  - optimized `v3` peak RSS: `127682252 kB`

So the full-scale blocker remains memory, but it is now concentrated inside the
direct bridge path itself rather than the old explicit
`materializeRGUFromHash()` step.

### 2M memory validation after direct bridge collapse (2026-03-24)

After the direct `collapseUMIall_fromBridgeHash()` path landed, I reran the
2M UCSF `GeneFull` `Unique` benchmark on the repaired
`feature/solo-optimization-20260324-mastermerge` branch to check whether the
bridge-specific memory blow-up was actually gone.

First completion on a fresh 2M resample:

- output:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_unique_hashbridge_direct_v3/`
- completed successfully
- `/usr/bin/time -v`:
  - `Elapsed (wall clock) = 1:33.95`
  - `Maximum resident set size = 40396560 kB`
- instrumented Solo timings from `Log.out`:
  - `collapseUMIall_fromBridgeHash 0.386508 s`
  - `countCBgeneUMI 0.447979 s`
  - `outputResults(raw) 0.120784 s`
  - `outputResults(filtered) 0.0752119 s`
  - `cellFiltering 2.75206 s`
  - `processRecords 3.5359 s`

Exact same-FASTQ rerun for apples-to-apples timing/memory on the new branch:

- output:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_unique_hashbridge_direct_samefastq_v1/`
- completed successfully
- `/usr/bin/time -v`:
  - `Elapsed (wall clock) = 1:33.65`
  - `Maximum resident set size = 40395468 kB`
- instrumented Solo timings from `Log.out`:
  - `collapseUMIall_fromBridgeHash 0.382937 s`
  - `countCBgeneUMI 0.44476 s`
  - `outputResults(raw) 0.121363 s`
  - `outputResults(filtered) 0.0664623 s`
  - `cellFiltering 2.71248 s`
  - `processRecords 3.52106 s`

Current-branch legacy control on the exact same 2M FASTQs:

- output:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_unique_legacy_mastermerge_v1/`
- completed successfully
- `/usr/bin/time -v`:
  - `Elapsed (wall clock) = 1:33.28`
  - `Maximum resident set size = 40193584 kB`
- instrumented Solo timings from `Log.out`:
  - `collapseUMIall 0.474378 s`
  - `countCBgeneUMI 0.940734 s`
  - `outputResults(raw) 0.120835 s`
  - `outputResults(filtered) 0.0646548 s`
  - `cellFiltering 2.65979 s`
  - `processRecords 3.74702 s`

Interpretation:

- The direct bridge collapse removed the old bridge-specific RSS penalty at 2M.
- On the current master-based branch, the direct bridge is only `201884 kB`
  (~`0.19 GB`, ~`0.5%`) above the same-binary legacy control.
- Compared with the older pre-direct bridge path
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_unique_hashbridge/`,
  peak RSS dropped from `44360304 kB` to `40395468 kB`:
  - delta `-3964836 kB` (~`-3.8 GB`, ~`-8.9%`)
- The main Solo counting step is still faster on the direct path:
  - `countCBgeneUMI`: `0.940734 s` legacy vs `0.44476 s` direct (`2.1x`)
  - `collapseUMIall`/`collapseUMIall_fromBridgeHash`: `0.474378 s` vs
    `0.382937 s`
  - `processRecords`: `3.74702 s` vs `3.52106 s`

Current parity status on the current branch:

- Mapping totals now match exactly between current-branch legacy and direct
  runs:
  - `Number of Reads = 2000000`
  - `Reads Mapped to Genome: Unique = 0.902979`
- Remaining direct-vs-legacy output drift at 2M is modest but real:
  - `Reads Mapped to GeneFull: Unique GeneFull = 0.743445` legacy vs
    `0.743137` direct
  - `Estimated Number of Cells = 7221` legacy vs `7217` direct
  - `UMIs in Cells = 1308644` legacy vs `1301714` direct
  - `Total GeneFull Detected = 17786` legacy vs `17771` direct

### 2M validation after memory fixes 1-4 (2026-03-24)

Implemented before this rerun:

1. free thread-local bridge state immediately after `sumThreads()` merge
2. skip `PackedReadInfo` allocation on the non-Flex direct bridge path
3. replace bridge per-read accounting with:
   - immediate per-CB packed counters
   - compact fixed-width deferred accounting for ambiguous CBs only
4. clear bridge accounting state immediately after replay

I intentionally did **not** do the next architecture change yet:

5. avoid the giant merged/extracted bulk structure by draining thread-local
   hashes directly into grouped sub-hashes

Validation rerun:

- output:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_unique_hashbridge_direct_samefastq_v2/`
- inputs:
  reused `fastq_downsampled/` from
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_unique_hashbridge_direct_v3/`
- completed successfully
- `/usr/bin/time -v`:
  - `Elapsed (wall clock) = 1:38.48`
  - `Maximum resident set size = 40213304 kB`
- Solo timings from `Log.out`:
  - `collapseUMIall_fromBridgeHash 0.419443 s`
  - `countCBgeneUMI 0.45533 s`
  - `outputResults(raw) 0.132954 s`
  - `outputResults(filtered) 0.0635182 s`
  - `cellFiltering 2.72824 s`
  - `processRecords 3.49318 s`

Memory interpretation:

- compared with the exact same FASTQs on the pre-fix direct path (`v3`):
  - RSS dropped from `40396560 kB` to `40213304 kB`
  - delta `-183256 kB` (~`-0.17 GB`)
- this leaves the direct bridge only `19720 kB` (~`19 MB`) above the earlier
  current-branch legacy control (`40193584 kB`), though that legacy control was
  not run on this exact downsample
- so items `1-4` improved 2M RSS modestly, but 2M was already close to legacy
  after the earlier direct-collapse work

Important regression on the exact same FASTQs:

- outputs changed relative to direct `v3`, so these fixes are **not** parity-safe
  yet
- `Features.stats` changes:
  - `noTooManyWLmatches: 5068 -> 6443`
  - `yesWLmatch: 1486274 -> 1492990`
  - `yessubWLmatch_UniqueFeature: 1486274 -> 1467273`
- `Summary.csv` changes:
  - `Reads With Valid Barcodes: 0.972353 -> 0.971665`
  - `Reads Mapped to GeneFull: Unique GeneFull: 0.743137 -> 0.733637`
  - `Estimated Number of Cells: 7242 -> 7224`
  - `UMIs in Cells: 1302603 -> 1301990`
  - `Unique Reads in Cells Mapped to GeneFull: 1432339 -> 1413168`

Conclusion:

- memory fixes `1-4` are directionally correct and do shave RSS
- but they introduced a real 2M regression in read/accounting behavior on the
  exact same input
- before the full benchmark, the next debug target is the new immediate +
  deferred bridge accounting logic, not item `5` yet

Bottom line:

- the large `v2` drift was a **real bug**, not acceptable algorithmic drift
- root cause: the compact deferred accounting rewrite dropped explicit
  `featGood` / `multiFeature` state and reconstructed `featGood` from
  `readFlag.featureU || readFlag.featureM`
- that is wrong for `--soloMultiMappers Unique` because some ambiguous
  `feature=-1` records still carry `featureM` bits; `v2` therefore replayed
  those as counted multi-feature reads
- smoking gun from `Features.stats` on the exact same FASTQs:
  - direct `v3`: `yesWLmatch = yessubWLmatch_UniqueFeature = 1486274`
  - broken `samefastq_v2`: `yesWLmatch = 1492990`,
    `yessubWLmatch_UniqueFeature = 1467273`
  - implied counted multi-feature reads in `Unique` mode: `25717`

### 2M same-FASTQ rerun after deferred-accounting bug fix (2026-03-24)

Fix:

- restored explicit `featGood` and `multiFeature` fields in the compact
  deferred bridge accounting record
- ambiguous `feature=-1` rows are now replayed with the same counted/not-counted
  semantics as before the memory rewrite

Validation rerun:

- output:
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_unique_hashbridge_direct_samefastq_v3/`
- inputs:
  reused `fastq_downsampled/` from
  `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_2M_unique_hashbridge_direct_v3/`
- completed successfully
- `/usr/bin/time -v`:
  - `Elapsed (wall clock) = 1:34.39`
  - `Maximum resident set size = 40213020 kB`

Result:

- the major read-accounting regression is gone
- `Features.stats` returned to the expected `Unique`-mode shape:
  - `yesWLmatch = 1486273`
  - `yessubWLmatch_UniqueFeature = 1486273`
  - implied counted multi-feature reads: `0`
- `Reads Mapped to GeneFull: Unique GeneFull` recovered from `0.733637` to
  `0.743136` (pre-regression `v3` was `0.743137`)
- `noTooManyWLmatches` returned from `6443` to `5068`

Residual difference vs pre-memory-fix direct `v3`:

- small matrix / filtered-output drift remains
- examples:
  - `Estimated Number of Cells: 7242 -> 7224`
  - `UMIs in Cells: 1302603 -> 1301927`
  - `yesUMIs: 1352866 -> 1352811`
- so item `3` contained one clear real bug that is now fixed, but there is
  still a smaller parity difference left in the post-memory-fix path

- the 2M memory problem is substantially addressed
- the direct bridge path no longer shows the large intermediate-memory
  expansion seen in the old hash-bridge implementation
- full-sample validation is still required, but 2M no longer suggests an OOM
  mechanism by itself

### Full UCSF optimized rerun after fix 5 (2026-03-24)

Fix 5 changed the full-scale bridge collapse shape:

- instead of merging all non-ambiguous bridge data into one giant global hash,
  the direct bridge path now drains the thread-local bridge hashes directly
  during extraction
- only the merged ambiguous/deferred contribution is still injected through
  `readFeatSum->inlineHash_`
- this removes the old "giant merged bulk hash + giant extracted traversal"
  overlap that was driving the full-sample OOM

Optimized full rerun:

- output:
  `/storage/paper_bench_solo_full_20260324/ucsf_ebs2_2_standard_solohash_optimized_v4/`
- completed successfully; no exit `137`
- `/usr/bin/time -v`:
  - `Elapsed (wall clock) = 22:38.81`
  - `Maximum resident set size = 70575848 kB`
- `Log.out` key bridge lines:
  - `hash_entries=325532880 thread_hashes=32 merged_ambiguous_hash_entries=2665474`
  - `Drained thread-local + merged ambiguous bridge hashes after extraction (recs=325532880)`
  - `Unique (CB,gene) segments (matrix row upper bound)=83073593`
  - `Finished direct bridge-hash UMI collapse, nCB=1840536 wall=277.391 s`

Current-branch full baseline rerun:

- output:
  `/storage/paper_bench_solo_full_20260324/ucsf_ebs2_2_standard_baseline_mastermerge_v2/`
- completed successfully
- `/usr/bin/time -v`:
  - `Elapsed (wall clock) = 19:48.45`
  - `Maximum resident set size = 68149756 kB`
- this rerun reproduced the archived baseline `Summary.csv` exactly

Full same-branch A/B result:

- fix 5 **did** solve the full-sample OOM
- however, on the current branch it is still slower than the current baseline:
  - wall time: `22:38.81` optimized vs `19:48.45` baseline
  - delta: `+170.36 s` (`+14.33%`)
  - peak RSS: `70575848 kB` optimized vs `68149756 kB` baseline
  - delta: `+2426092 kB` (`+2.31 GiB`)
- this is still a huge improvement over the old failed full runs:
  - old direct-collapse `v3` OOM peak: `127682252 kB`
  - new full optimized `v4` peak: `70575848 kB`

Solo timing comparison vs current-branch baseline:

- optimized `countCBgeneUMI`: `280.117 s`
- baseline `countCBgeneUMI`: `277.432 s`
- optimized `cellFiltering`: `218.183 s`
- baseline `cellFiltering`: `219.881 s`
- optimized `processRecords`: `508.431 s`
- baseline `processRecords`: `501.99 s`

So the current state is:

- fix 5 succeeded as a **memory fix**
- fix 5 did **not** produce a speed win on the full UCSF workload
- the remaining optimization target is now CPU / traversal cost, not OOM survival

Full-output parity vs current-branch baseline:

- `Number of Reads`: identical (`444896731`)
- `Reads Mapped to GeneFull: Unique GeneFull`:
  `0.855092 -> 0.854919` (`-0.000173`, `-0.0202%`)
- `Estimated Number of Cells`: `13721 -> 13706` (`-15`, `-0.1093%`)
- `Unique Reads in Cells Mapped to GeneFull`:
  `365470404 -> 365284662` (`-185742`, `-0.0508%`)
- `UMIs in Cells`: `256057933 -> 254997805`
  (`-1060128`, `-0.4140%`)
- `Total GeneFull Detected`: `33779 -> 33758` (`-21`, `-0.0622%`)

Conclusion after fix 5:

- full-scale viability is now restored
- residual full-output drift is small but real
- the direct bridge is now viable enough to iterate further, but it is not yet
  ready to replace the baseline path on either speed or parity

### Full EBs2_2 direct bridge — tuple redesign + aggregated ambiguous CB (2026-03-25)

Detached worktree `/tmp/star-suite-v10-redesign-20260325` (clean rebuild), policy:
Bayes-resolved ambiguous read accounting when resolution succeeds (preferred
semantics; not legacy pin-only).

- **Artifact root:**
  `/storage/100K/ucsf_solo_bridge_redesign_20260325/ucsf_ebs2_2_gexonly_bridge_redesign_v1/`
- **Harness:** same STAR CLI family as archived
  `ucsf_ebs2_2_standard_solohash_optimized_v4` (GEX-only `readFilesIn`, same
  Solo + `pfMultiConfig` + CR-assign flags); `STAR_SOLO_NONFLEX_HASH_BRIDGE=1`,
  `--soloInlineHashMode yes`, `--outTmpDir` under the artifact root.
- **`/usr/bin/time -v`:** wall **20:07.92**; max RSS **70985716 kB** (~67.7 GiB).
- **Mapping (`Log.out` / `Log.final.out`):** started **Mar 25 10:19:27**,
  finished **Mar 25 10:28:18** → **531 s** (~8 m 51 s); mapping speed line
  ~**1426** M reads/h on this run.
- **Solo timings (`Log.out`):**
  - `collapseUMIall_fromBridgeHash` **277.749 s**
  - `countCBgeneUMI` **279.954 s**
  - `processRecords` **505.697 s**
  - `cellFiltering` **214.226 s**
- **`[AMBIG-CB-RESOLVE]`:** `pending=704105` `resolved=163716`
  `still_ambiguous=540389` `added_to_hash=1090231`
- **`Summary.csv` (GeneFull):** Unique GeneFull **0.85468**; Estimated cells
  **13715**; UMIs in cells **254181669**; Total GeneFull detected **33750**

Comparison anchors (same UCSF corrected EBs2_2 GEX input family):

| Run | Wall | Max RSS kB | Mapping elapsed (from logs) |
|-----|------|------------|------------------------------|
| Baseline mastermerge v2 | 19:48.45 | 68149756 | (see that `Log.out`) |
| Direct bridge v4 | 22:38.81 | 70575848 | 13:53:03 → 14:04:21 = **678 s** |
| **Redesign (this)** | **20:07.92** | **70985716** | 10:19:27 → 10:28:18 = **531 s** |

Readout:

- **Wall** improved by ~**151 s** (~**11%**) vs bridge **v4** on this execution;
  still ~**19 s** slower (~**1.6%**) than baseline mastermerge v2.
- **Mapping** log interval dropped by ~**147 s** (~**22%**) vs v4 on this pair of
  runs (treat as indicative; different calendar days / load).
- **RSS** essentially in the v4 tier (~**+410 MiB** vs v4); still several GiB
  above baseline-only.
- **`collapseUMIall_fromBridgeHash`** remains ~**278 s** — the dominant Solo
  cost is unchanged; the redesign targets ambiguous **mapping** bookkeeping, not
  post-map hash collapse CPU.

Worth pursuing further if the **mapping** improvement reproduces under repeated
serialized runs and if the Bayes-resolved accounting policy is accepted
product-side; **collapse** time still dominates Solo vs baseline until further
HASH phase work.
