# Handoff: Bootstrap + Tiered Hash Activation for Feature Assignment

**Date**: 2026-03-17  
**Branch**: `CR-Larry-perturb` (ahead of origin by 4 commits)  
**Last commit**: `2f0bd9c` — "pf: restore tiered hash search in bootstrap path, add safe merge policy"  
**Uncommitted changes**: Yes — bootstrap wiring + brute-force fallback (see below)

## Context

Feature assignment performance is dominated by the slow brute-force
`checkAndCorrectFeature` path, which scans every possible offset × every
feature for every read. Two fast paths exist but were not activated from STAR:

1. **Bootstrap-learned tiered hash** — first 100K reads use brute force to
   learn per-feature offsets, then all subsequent reads use O(n_offsets) hash
   probes via `pf_search_hash_offsets`.
2. **Integer-key prehash tables** — `feature_hamming_le1_hash`/`le2_hash` use
   `uint64_t`/`seq128_t` keys instead of heap-allocated strings.

Both are implemented and compiled into the binary. The integer-key conversion
works correctly. The problem is the bootstrap was never activated from STAR's
`PfMultiProcess.cpp` — it was added in commit `b0db247` but lost in the
destructive squash merge `b26a66e`.

## Update: 2026-03-17 Verification

The apparent `848 resolved` regression is **not** a read-assignment regression.

On a clean rebuild of the current tree, the UCSF 10K/20K fixture run at:

- `/tmp/ucsf_bootstrap_regression_current_20260317_080129`

produced:

- `Resolve outcomes total resolved=848 no_hit=34770 too_many_n=195 ...`
- but `guides/stats.txt` shows `Total_unmatched_reads 5648`

That means the actual assigned-read count is approximately:

- `20000 - 5648 = 14352`

which is already in the expected `~14K/20K` range for the old non-bootstrap
surface.

So the key conclusion is:

- `Resolve calls total` / `Resolve outcomes total resolved` are **not** the
  correct validation surface for this path.
- Those counters are dominated by `simpleCorrectFeature` /
  `checkAndCorrectFeature` call sites and do **not** represent final assigned
  reads after the exact/anchor fast paths.
- The real parity surface for this fixture is `guides/stats.txt` and, for
  stronger validation, raw/filtered MEX outputs.

An exploratory restoration of the `b0db247` replay-style bootstrap flow
did **not** improve the UCSF fixture on the correct surface:

- `/tmp/ucsf_bootstrap_regression_patched_20260317_080624`
- `guides/stats.txt`: `Total_unmatched_reads 5692`

So the immediate action item is **not** “restore replay because of the 848
number.” The next action item is to validate bootstrap/tiered-hash behavior on
assignment outputs, not resolve-call counters.

## Update: 2026-03-17 Direct Hash Validation

The remaining concern was whether the gRNA hash lookup itself was broken or
whether reads were simply reaching it with the wrong offset surface.

To test that directly, a small debug override was added in
`PfMultiProcess.cpp` so STAR can be run with environment overrides for:

- `STAR_PF_USE_FEATURE_ANCHOR_SEARCH`
- `STAR_PF_REQUIRE_FEATURE_ANCHOR_MATCH`
- `STAR_PF_FEATURE_BOOTSTRAP_READS`

This makes it possible to force the small UCSF fixture through the exact
non-bootstrap tiered-hash path without changing CLI semantics.

### Forced offset 0, bootstrap off, anchor off

Run:

- `/tmp/ucsf_hash_offset0_20260317`

Command surface:

- `STAR_PF_FEATURE_BOOTSTRAP_READS=0`
- `STAR_PF_REQUIRE_FEATURE_ANCHOR_MATCH=0`
- `STAR_PF_USE_FEATURE_ANCHOR_SEARCH=0`
- `--crAssignFeatureOffset 0`

Observed in `run.log`:

- three immediate `[DBG_MISS] off=0 ... sub=AAGCAGTGGTATCAACGCAG...`
- checkpoint trace:
  - `tiered_hit=0`
  - `tiered_miss=8192`

Observed in `guides/stats.txt`:

- `Total_unmatched_reads 5999`
- `Total deduped feature counts 13151`

Interpretation:

- the hash code is running
- it is probing offset `0`
- it is missing because the read starts with the TSO/preamble, not the guide

This argues **against** a broken integer encoding or hash-table lookup.

### What offset does the working path actually use?

From the working fallback surface in:

- `/storage/test_10k_trace/run/cr_assign/CRISPR_Guide_Capture/EBs2_2_AALG2_CRISPR_Guide_Capture_0/guides/feature_sequences.txt`

weighted match-position histogram:

- `27`: `45`
- `28`: `394`
- `29`: `205`
- `30`: `2874`
- `31`: `9860`
- `32`: `817`
- `33`: `12`
- `34`: `3`

So the true assignment surface for this library is concentrated at:

- `30-32`, dominated by `31`

not `0`.

### Forced offset 31, bootstrap off, anchor off

Run:

- `/tmp/ucsf_hash_offset31_20260317`

Command surface:

- same env overrides as above
- `--crAssignFeatureOffset 31`

Observed in `run.log`:

- immediate `[DBG_HIT] off=30 feat=... d=0 sub=GGAGGCTCCGAACTCGCCGG...`
- checkpoint trace:
  - `tiered_hit=5484`
  - `tiered_miss=2708`

Observed in `guides/stats.txt`:

- `Total_unmatched_reads 5658`
- `Total deduped feature counts 13472`

This is very close to the working baseline (`5648` unmatched, `13481`
deduped counts).

### Direct conclusion

For this small UCSF fixture:

- the tiered hash implementation itself is **working**
- the integer-key encoding is **not** the primary bug
- the reads fall through when the hash path is driven at the wrong offset
  surface
- the effective gRNA offset for this library is around `31`, with `30/+/-1`
  also contributing

So the real problem on small fixtures is:

- bootstrap never finalizes at `100000` reads
- default small-run routing does not promote the learned/observed `30-32`
  offset surface into the tiered-hash fast path soon enough

not a broken hash table.

## What Was Restored (Uncommitted)

### 1. C++ wiring (restored from `b0db247`)

Three files need the bootstrap plumbed through STAR → pf_api:

**`PfMultiAssign.h`** — Added fields to `AssignOptions`:
```cpp
bool useFeatureAnchorSearch = false;
bool requireFeatureAnchorMatch = false;
int featureModeBootstrapReads = 0;
bool useHotHash = false;
```

**`PfMultiProcess.cpp`** — Set defaults:
```cpp
assignOpts.useFeatureAnchorSearch = true;
assignOpts.requireFeatureAnchorMatch = true;
assignOpts.featureModeBootstrapReads = 100000;
```

**`PfMultiAssign.cpp`** — `applyAssignOptions()` wires to pf_api:
```cpp
pf_config_set_use_feature_anchor_search(cfg, 1);
pf_config_set_require_feature_anchor_match(cfg, 1);
pf_config_set_feature_mode_bootstrap_reads(cfg, options.featureModeBootstrapReads);
```

The pf_api setters already exist (`pf_api.c` lines 401-411).

### 2. Bootstrap brute-force fallback (new, at line ~4626 of `assignBarcodes.c`)

The bootstrap phase used only anchor-guided matching (`simpleCorrectFeature`
at anchor-derived positions). For UCSF CRISPR data, many reads have a 30bp
TSO preamble (`AAGCAGTGGTATCAACGCAG...`) before the guide, so anchor
matching fails to locate guides (anchors not at expected position relative to
offset 0).

Added: when the anchor passes find no match during bootstrap, fall back to
`checkAndCorrectFeature` (full brute-force scan). If it finds a match, record
the offset in the histogram via `feature_mode_record()`:

```c
if (!bestFeature && feature_mode_bootstrap_done == 0) {
    int bf_feat = checkAndCorrectFeature(sequence, features, maxHammingDistance,
                                          nThreads, &bf_hamming, bf_match_seq,
                                          max_feature_n, &bf_ambiguous, &bf_match_pos, stats);
    if (!bf_ambiguous && bf_feat > 0 && bf_hamming <= maxHammingDistance) {
        bestFeature = (uint32_t)bf_feat;
        feature_mode_record(bf_feat, (int)bf_match_pos);
    }
}
```

## Current Status / Open Problem

The brute-force fallback runs but produces fewer matches than expected:

| Path                          | Resolved | Total calls | Rate (reads/sec) |
|-------------------------------|----------|-------------|-------------------|
| Non-bootstrap (old)           | 14,001   | 20,000      | 13,900            |
| Bootstrap anchor-only         | 706      | 30,023      | 327,800           |
| Bootstrap + brute-force (new) | 848      | 35,813      | 37,700            |

Expected: the brute-force fallback should recover similar match counts to the
non-bootstrap path (~14K/20K). The 848 matches suggest something is
suppressing the brute-force from firing or succeeding.

Update: this interpretation was too strict. The `848` figure does **not**
track final assigned reads. The current clean run already lands at
`14352/20000` assigned reads on `guides/stats.txt`, so the assignment-level
regression signal is not present on this fixture.

### Likely root causes to investigate

1. **`Resolve calls total 35,813` for 20,000 reads** — the count exceeds the
   number of reads. This comes from `simpleCorrectFeature` calls in the
   anchor passes (each calls `pf_stats_record_resolve`). The brute-force
   fallback's `checkAndCorrectFeature` adds more. But the resolve outcomes
   show `no_hit=34,770` — most calls produce no hit. Need to verify the
   brute-force fallback is actually being reached for the reads that the
   anchor path misses.

2. **Guard condition `!bestFeature`** — the brute-force only fires when
   `bestFeature == 0`. But the anchor passes might set `bestFeature` to a
   non-zero value with a hamming distance > `maxHammingDistance` (a "soft
   failure"). Check whether `bestFeature` is being set to a feature with
   `bestHamming > maxHammingDistance`, which would prevent the fallback from
   running but still be treated as no-match later.

3. **`feature_mode_bootstrap_done == 0` guard** — with only 20K reads and
   `bootstrap_reads=100000`, bootstrap never finishes. All reads stay in the
   bootstrap phase. This is expected for the 10K test fixture but not an issue
   for full runs (millions of reads).

4. **Double-counting in resolve stats** — the `simpleCorrectFeature` calls
   within the anchor passes also record resolve stats. This inflates the
   resolve total and may confuse the match rate calculation. The 848 resolved
   may be only from the brute-force path, while the anchor path contributes
   separately. Need to separate the two sources.

5. **`checkAndCorrectFeature` behavior with `maxN=0`** — when
   `--crAssignFeatureN 0`, `checkAndCorrectFeature` calls
   `checkSequenceAndCorrectForN` which rejects reads with ANY N in the guide
   region. The full read includes Ns at various positions. In the non-bootstrap
   path, the code extracts the guide subsequence before calling
   `checkAndCorrectFeature`. In the bootstrap fallback, we pass the full
   `sequence` pointer — `checkAndCorrectFeature` handles the full read and
   internally scans at all offsets, so this should be OK. But verify.

## Debugging approach for next agent

1. **Add a counter** for how many times the brute-force fallback fires vs
   how many reads enter the bootstrap phase. Print at end:
   ```c
   static unsigned long long bf_attempted = 0, bf_hit = 0;
   ```
   This tells you whether the guard conditions are suppressing the fallback.

2. **Compare head-to-head**: run the same 10K fixture through the
   non-bootstrap path (set `featureModeBootstrapReads = 0` and
   `requireFeatureAnchorMatch = false`) to confirm the brute-force still
   produces 14K matches.

3. **Once the bootstrap phase is producing matches**: create a larger test
   fixture (e.g., 200K reads) so the bootstrap can finalize and switch to the
   tiered hash. Verify the tiered hash then gets hits.

4. **Full E2E benchmark**: run the corrected `EBs2_2` dataset
   (`/storage/EBs2_2_output_nvme/`) with bootstrap enabled. Target: 22-25 min
   total wall time (down from 61 min without tiered hash).

## Test infrastructure

- **10K fixture**: `/storage/test_10k_trace/` — 10K reads per lane×2 lanes
  from `EBs2_2_AALG2` guide CRISPR Guide Capture library.
- **Config**: `/storage/test_10k_trace/pf_multi_config.csv`
- **Feature ref**: `/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv`
  (548 features, 20bp, pattern `(BC)GTTTNAGAGCTAAGC`).
- **Full dataset**: `/storage/EBs2_2_output_nvme/` (EBs2_2 on NVMe).
- **Run script for full benchmark**: `scripts/run_ucsf_perturb_yremove_batch.sh`
  (but parameters may need updating — check the `--crAssignFeatureOffset`).

## Smoke test command (10K)

```bash
rm -rf /storage/test_10k_trace/run /storage/test_10k_trace/tmp
mkdir -p /storage/test_10k_trace/run

/mnt/pikachu/STAR-suite/core/legacy/source/STAR \
  --runThreadN 4 --genomeDir /storage/autoindex_110_44/bulk_index \
  --readFilesIn /dev/null /dev/null --readFilesCommand zcat \
  --outFileNamePrefix /storage/test_10k_trace/run/ \
  --outTmpDir /storage/test_10k_trace/tmp \
  --outSAMtype None --soloType None \
  --pfMultiConfig /storage/test_10k_trace/pf_multi_config.csv \
  --crChemistry NXT --crOutputChemistry TRU \
  --crWhitelist /home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt \
  --crMinUmi 3 --crAssignMaxHamming 1 --crAssignFeatureOffset -1 \
  --crAssignLimitSearch -1 --crAssignMinCounts 0 \
  --crAssignMaxBarcodeMismatches 5 --crAssignFeatureN 0 --crAssignBarcodeN 1 \
  --crAssignConsumerThreads 2 --crAssignSearchThreads 1 \
  --crAssignSkipQcOutputs 1 \
  --crFeatureRef /mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv \
  2>&1 | grep -E "TRACE|Finished|bootstrap|Finalized|offset|Resolve"
```

## Key code locations

| File | Lines | What |
|------|-------|------|
| `assignBarcodes.c` | 4312-4654 | Bootstrap block (bootstrap_reads > 0) |
| `assignBarcodes.c` | 4327-4361 | Post-bootstrap tiered hash (the fast path) |
| `assignBarcodes.c` | 4420-4624 | Bootstrap phase anchor matching |
| `assignBarcodes.c` | 4626-4641 | **NEW**: brute-force fallback during bootstrap |
| `assignBarcodes.c` | 87-96 | `feature_mode_record()` — histogram accumulation |
| `assignBarcodes.c` | 101-170 | `feature_mode_finalize()` — picks mode offsets, builds search_offsets |
| `assignBarcodes.c` | 3565-3660 | `pf_single_offset_hash_search()` — integer-key tiered d0→d1→d2 |
| `assignBarcodes.c` | 3661-3760 | `pf_search_hash_offsets()` — iterates learned offsets |
| `PfMultiProcess.cpp` | ~1847 | Bootstrap params set (restored) |
| `PfMultiAssign.cpp` | ~348-356 | `applyAssignOptions` wiring (restored) |
| `PfMultiAssign.h` | ~47-50 | AssignOptions fields (restored) |
| `pf_api.c` | 401-411 | Config setters (already existed) |
| `star_feature_call.cpp` | 177-179 | Reference: compat_perturb mode already sets bootstrap=100K |

## Debug trace code to remove before final commit

`assignBarcodes.c` contains several `[TRACE]` and `[DBG_HIT]`/`[DBG_MISS]`
fprintf statements added during investigation. These should be removed or
gated behind a `PF_DEBUG` env var before the final commit. Search for:
- `_trace_tiered_hit`, `_trace_tiered_miss`, `_trace_legacy`, etc.
- `[TRACE tid=`, `[DBG_HIT]`, `[DBG_MISS]`
- Static `_pf_dbg_hit`, `_pf_dbg_miss` counters in `pf_search_hash_offsets`

## Prior related handoff docs

- `docs/HANDOFF_INTEGER_KEY_PREHASH_REFACTOR_20260304.md` — integer-key conversion
- `docs/HANDOFF_PF_CONSUMER_OPTIMIZATIONS_20260317.md` — consumer threading fixes
- `docs/HANDOFF_NXT_AUTODETECT_FIX_AND_REGRESSION_20260316.md` — auto-detect bugs
- `AGENTS.md` section "Safe Merge Policy" — to prevent future squash-merge regressions

## Build

```bash
# Clean rebuild (mandatory after any source change)
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
```

Note: `make clean` in `core/legacy/source` does NOT clean
`libprocess_features.a`. If `assignBarcodes.c` changes are not picked up,
also run:
```bash
rm -f core/features/process_features/libprocess_features.a
```
