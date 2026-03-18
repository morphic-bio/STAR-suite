# Handoff: Flex Tiered Hash Pre-Filter Post-Master Merge

**Date**: 2026-03-18  
**Worktree**: `/mnt/pikachu/STAR-suite-flex-h01-cache`  
**Branch**: `feature/flex-hash-prefilter-postmaster-20260318`  
**Base merge commit**: `98a8026` (`master` merged into prior Flex hash branch)  
**Status**: merged onto current `master`, not yet rebuilt or revalidated after merge

## What Has Been Done

The Flex H0/H1 pre-filter prototype is implemented in the feature worktree.

Main code paths:

- runtime cache loader and classifier:
  - [FlexHashScreen.cpp](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/FlexHashScreen.cpp)
  - [FlexHashScreen.h](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/FlexHashScreen.h)
- runtime hook before alignment:
  - [ReadAlign_oneRead.cpp](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/ReadAlign_oneRead.cpp)
- parameter plumbing / auto-discovery / disable flag:
  - [ParametersSolo.cpp](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/ParametersSolo.cpp)
  - [ParametersSolo.h](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/ParametersSolo.h)
  - [parametersDefault](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/parametersDefault)
- Flex-side keep/deny record handling:
  - [SoloReadFeature_record_flex.cpp](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/flex/SoloReadFeature_record_flex.cpp)
- cache builder / pilot harness:
  - [flex_h01_pilot.py](/mnt/pikachu/STAR-suite-flex-h01-cache/scripts/flex_h01_pilot.py)
  - [run_flex_h01_pilot.sh](/mnt/pikachu/STAR-suite-flex-h01-cache/scripts/run_flex_h01_pilot.sh)
  - [run_flex_hash_screen_internal_100k.sh](/mnt/pikachu/STAR-suite-flex-h01-cache/tests/run_flex_hash_screen_internal_100k.sh)

Runtime behavior implemented on this branch:

- hash screen is enabled by default for Flex when a cache is available
- environment or discovered cache path is accepted
- `--no-hash-screen yes` forces the legacy align-everything path
- read screening happens before normal alignment
- `KEEP` can short-circuit into the Flex counting path
- `DENY` currently corresponds to the `NEG_PROBE_AMBIG` class

## What Was Validated Before the Master Merge

The important validation was done before merging current `master` into this
branch, and those results are still the reference point.

Primary summary doc:

- [FLEX_H01_DENY_CLASS_ANALYSIS_20260315.md](/mnt/pikachu/STAR-suite-flex-h01-cache/docs/FLEX_H01_DENY_CLASS_ANALYSIS_20260315.md)

Key milestones already achieved:

1. 100-probe hybrid H0/H1 pilot:
   - exact raw-MEX parity on the 100K fixture for the selected probe set
2. full probe-set cache build:
   - `8,070,044` records
   - binary cache artifact:
     `/storage/downsampled_100K/SC2300771/results/flex_h01_full_cache_20260315_153914/reclassified/sequence_cache.bin`
3. deny-class analysis:
   - only `NEG_PROBE_AMBIG` looked safe as a hard deny
   - `NEG_NO_GENE` and `NEG_NO_CANDIDATES` were not safe terminal classes
4. internal Flex hash-screen integration:
   - after fixing the sample-whitelist surface, the internal 100K run was
     byte-identical to legacy on keyed raw MEX

Most important exact-parity result:

- artifact root:
  `/storage/downsampled_100K/SC2300771/results/flex_hash_screen_internal_full16_20260315_203745/`
- same-binary keyed raw-MEX compare:
  - `619448` exact entries
  - `0` mismatches
- timing:
  - hash-on: `50.48s`
  - legacy: `59.67s`

That is the strongest proof point for the approach so far.

## What Happened On 2026-03-18

Current `master` was merged into the previous Flex branch:

- old branch head: `feature/flex-h01-cache-instrumentation`
- merge commit: `98a8026`

During the merge:

- working changes were stashed
- `master` merged cleanly
- stash reapply produced conflicts in:
  - [tests/ARTIFACTS.md](/mnt/pikachu/STAR-suite-flex-h01-cache/tests/ARTIFACTS.md)
  - [parametersDefault.xxd](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/parametersDefault.xxd)
- both conflicts were resolved:
  - `tests/ARTIFACTS.md` kept both the UCSF bootstrap notes from `master` and
    the Flex hash-screen artifact entries from this branch
  - `parametersDefault.xxd` was regenerated from merged
    [parametersDefault](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/parametersDefault)

After that merge, this new branch was created:

- `feature/flex-hash-prefilter-postmaster-20260318`

## Current Branch State

This branch is intentionally still dirty. The merge and stash reapply restored
the Flex work, but the post-merge branch has **not** yet been rebuilt or
revalidated.

Currently modified tracked files include:

- [Makefile](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/Makefile)
- [Parameters.cpp](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/Parameters.cpp)
- [ParametersSolo.cpp](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/ParametersSolo.cpp)
- [ParametersSolo.h](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/ParametersSolo.h)
- [ReadAlign.h](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/ReadAlign.h)
- [ReadAlign_oneRead.cpp](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/ReadAlign_oneRead.cpp)
- [ReadAlign_outputAlignments.cpp](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/ReadAlign_outputAlignments.cpp)
- [SoloReadFeature.h](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/SoloReadFeature.h)
- [SoloReadFeature_record_shared.h](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/SoloReadFeature_record_shared.h)
- [Stats.cpp](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/Stats.cpp)
- [Stats.h](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/Stats.h)
- [SoloReadFeature_record_flex.cpp](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/flex/SoloReadFeature_record_flex.cpp)
- [parametersDefault](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/parametersDefault)
- [parametersDefault.xxd](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/parametersDefault.xxd)
- [ARTIFACTS.md](/mnt/pikachu/STAR-suite-flex-h01-cache/tests/ARTIFACTS.md)

Relevant untracked files restored from the Flex work:

- [FlexHashScreen.cpp](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/FlexHashScreen.cpp)
- [FlexHashScreen.h](/mnt/pikachu/STAR-suite-flex-h01-cache/core/legacy/source/FlexHashScreen.h)
- [FLEX_H01_DENY_CLASS_ANALYSIS_20260315.md](/mnt/pikachu/STAR-suite-flex-h01-cache/docs/FLEX_H01_DENY_CLASS_ANALYSIS_20260315.md)
- [run_flex_hash_screen_internal_100k.sh](/mnt/pikachu/STAR-suite-flex-h01-cache/tests/run_flex_hash_screen_internal_100k.sh)

There is also a safety backup stash still present:

- `stash@{0}` = `pre-merge-flex-h01-20260318`

## What Needs To Be Done Next

### 1. Rebuild on the merged branch

Mandatory first step:

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
```

Do not trust any regression signal until this merged branch builds cleanly.

### 2. Run the 100K internal parity smoke again

Use:

- [run_flex_hash_screen_internal_100k.sh](/mnt/pikachu/STAR-suite-flex-h01-cache/tests/run_flex_hash_screen_internal_100k.sh)

Primary check:

- confirm the exact-parity full-16-whitelist result still holds after the
  `master` merge

Success criterion:

- keyed raw-MEX exact match against same-binary legacy control

### 3. Verify parameter and startup behavior

Recheck:

- cache auto-discovery
- `--no-hash-screen yes`
- missing-cache fallback to legacy path
- logging in `Log.out`

This matters because `master` brought substantial parameter and `pf` changes.

### 4. Clean up branch contents

Decide what to do with:

- untracked test outputs like `-/BC001/`
- `scripts/__pycache__/`
- any restored helper docs or scripts that should be committed vs ignored

### 5. Move beyond 100K

If the 100K rerun still passes:

1. run the merged branch on the 2M fixture
2. capture keyed raw/filtered MEX parity against legacy
3. compare elapsed time and alignment-call reduction

### 6. Prepare for upstreaming

Before bringing this back toward the main tree:

- collapse the branch into the minimal required code/files
- keep the offline cache builder external
- keep the runtime path small and auditable
- preserve `--no-hash-screen yes` as the legacy escape hatch

## Practical Recommendation

Treat the branch as having one strong proof point and one major missing step:

- strong proof point:
  - exact internal 100K parity already achieved before the merge
- missing step:
  - no clean rebuild + rerun has been done since merging `master`

So the immediate next agent should not redesign anything. The correct next step
is simply:

1. rebuild
2. rerun the 100K internal parity harness
3. only then decide whether any merge fallout exists
