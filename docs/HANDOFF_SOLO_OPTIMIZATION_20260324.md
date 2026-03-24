# Handoff: Solo Optimization

Date: 2026-03-24
Branch: `feature/solo-optimization-20260324`
Base commit: `7a7fb08`
Audience: next coding agent working on Solo / STARsolo post-map optimization

## Scope

This branch is intended for Solo-phase optimization work, not FLEX mapping-pipeline work.

Primary target surface:

- standard Solo / STARsolo post-map counting and output
- shared code that also matters for Flex, but only where the optimization is clearly worth the shared blast radius

Do not start by redesigning the read-stage pipeline here. The main optimization opportunity on this branch is the Solo finalize path after mapping.

## Current Branch State

- Branch created from `perturb` as `feature/solo-optimization-20260324`
- No tracked code changes yet on this branch from me
- There are unrelated untracked local files in the checkout; leave them alone unless explicitly asked

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

## Related Files Worth Reading

- [SoloFeature_countCBgeneUMI.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_countCBgeneUMI.cpp)
- [SoloFeature_collapseUMIall.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_collapseUMIall.cpp)
- [SoloFeature_outputResults.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_outputResults.cpp)
- [MexWriter.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/MexWriter.cpp)
- [SoloFeature.h](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature.h)
- [tests/run_solo_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_solo_smoke.sh)
- [tests/run_cr_parity_100k.sh](/mnt/pikachu/STAR-suite/tests/run_cr_parity_100k.sh)
