# Handoff: PF Cumulative Hash + Per-Library Lookup State (2026-03-02)

## Scope
This handoff covers the implementation of cumulative feature prehash lookups for `process_features`, with lookup state moved from globals to per-library (`feature_arrays`) ownership, plus regression runs on UCSF 2M and UCSF full dynamic E2E.

Primary goals implemented:
- Per-library lookup state for multi-library correctness.
- Cumulative hash semantics (`<=1`, `<=2`) so one lookup can return best distance and ambiguity status.
- Dynamic-ceiling pruning (`bestHamming`) in alternate-position searches.

## Code Changes
### Files modified
- `core/features/process_features/include/common.h`
- `core/features/process_features/include/globals.h`
- `core/features/process_features/src/io.c`
- `core/features/process_features/src/assignBarcodes.c`
- `core/features/process_features/src/barcode_match.c`
- `core/features/process_features/src/globals.c`

### Implementation summary
1. `feature_arrays` now owns lookup state:
- `feature_hamming_le1_hash` (cumulative `<=1`)
- `feature_hamming_le2_hash` (cumulative `<=2`)
- `feature_no_ambiguity_le1` (per-feature flags)
- `feature_no_ambiguity_le2` (per-feature flags)
- enable flags for each hash

2. Prehash payload now stores:
- feature id (1-based)
- best Hamming distance (`0/1/2`)
- ambiguity bit at the best distance

3. Prehash construction in `io.c`:
- exact entries included in cumulative tables
- `<=1` table includes exact + one-mismatch
- `<=2` table includes exact + one-mismatch + two-mismatch
- tie handling updates ambiguity bit and per-feature no-ambiguity flags

4. Lookup path in `simple_hamming_search`:
- exact check first
- single cumulative hash query based on current ceiling (`<=2` preferred when allowed)
- payload-driven return of feature and distance
- ambiguous payload returns no-call for that query
- fallback to full scan remains

5. Dynamic ceiling propagation:
- `simpleCorrectFeature` now passes current `bestHamming` as query ceiling when iterating N-alternates
- anchor/offset fallback loops in `process_feature_sequence` pass reduced ceiling as better matches are found

6. Memory lifecycle:
- per-library hashes/flexible arrays are freed in `free_feature_arrays`
- global cleanup only clears `feature_code_hash` (`clear_feature_lookup_hashes`)

## Build Status
- STAR rebuilt cleanly:
  - `make -C core/legacy/source clean`
  - `make -C core/legacy/source -j8 STAR`
- `process_features` tools built cleanly:
  - `make -C core/features/process_features -j8 tools`

Note:
- A parallel build race can transiently fail `assignBarcodes`/`call_features` link if `../libscrna/libscrna.a` is not yet built. Re-running tools after STAR/libscrna build resolves this.

## Validation Runs
## 1) UCSF 2M quick A/B regression (current code, same inputs)
Artifact root:
- `/tmp/ucsf2m_cumhash_ab_20260302_002130`

Run modes:
- `prehash_off`: `--feature_prehash_max_hamming 0`
- `prehash_on`: `--feature_prehash_max_hamming 2`
- both with `--feature_prehash_max_entries 50000000`
- input: `/storage/ucsf-2M/guides/iPSC2_1_AALG2`
- `--max_reads 1000000`

Results:
- `stats.txt`: identical
- `features.txt`: identical
- raw and filtered barcode sets: identical sets
- raw and filtered matrix content: canonically identical
  - order differs in text files (`barcodes.txt`, `matrix.mtx`, `feature_per_cell.csv`) but canonical matrix by barcode is equal

Timing:
- off: wall `0:35.76`, max RSS `329104 KB`
- on: wall `0:35.85`, max RSS `381236 KB`

## 2) UCSF full dynamic E2E (same policy)
New run:
- `/storage/ucsf-full/bench_20260218_dynamic_first/runs/star_full_dynamic_32x32_cumhash_perfeature_20260302_002345`

Reference run:
- `/storage/ucsf-full/bench_20260218_dynamic_first/runs/star_full_dynamic_32x32_prehash_20260301_222331`

Timing:
- new: wall `20:50.45`, max RSS `44021648 KB`
- ref: wall `20:53.24`, max RSS `44019664 KB`

Parity vs reference:
- `Log.final.out` key mapping metrics: identical
- `Solo.out/GeneFull/{Summary.csv,Features.stats,UMIperCellSorted.txt}`: identical
- CRISPR outputs in `outs/crispr_analysis/`: identical
- raw/filtered MEX outputs in `outs/*_feature_bc_matrix/`: identical

## 3) UCSF 2M STAR dynamic check against older fixture
New run:
- `/storage/ucsf-2M/star_runs/star_2m_dynamic_cumhash_perfeature_20260302_001811`

Reference fixture:
- `/storage/ucsf-2M/star_runs/star_2m_dynamic_regress_exactoffset_20260301_093812`

Observed:
- key mapping metrics match
- CRISPR and matrix files differ from this older fixture

Interpretation:
- this older fixture does not represent a strict parity anchor for current tree state; use the same-code A/B run above for algorithmic regression validation.

## Operational Notes For Next Agent
1. If extending optimization logic, keep cumulative-hash semantics (`<=h`) intact; distance-only tables are not enough for single-lookup correctness.
2. Keep lookup state in `feature_arrays`; do not revert to globals if multi-library support is required.
3. For parity checks on text MEX files, use canonical matrix comparison keyed by barcode/feature rather than file-order `cmp`.
4. If needed, add optional debug output of payload decode path in `simple_hamming_search` for targeted divergence investigations.

## Follow-up Policy Update (hash-first, 2026-03-02)
In response to policy refinement, lookup order in `simple_hamming_search` was changed to:
1. probe highest available cumulative hash first (for `maxH=1`, `<=1` hash; for `maxH=2`, `<=2` hash),
2. fast-return no-call on hash miss,
3. decode payload on hit and return directly when unique.

Additionally, alternate-path checks now reduce query ceiling more aggressively when a unique best hit is already known:
- unique best at distance 1 with no ambiguity metadata -> remaining checks forced to exact (`maxH=0`)
- unique best at distance 2 with no ambiguity metadata -> remaining checks forced to `<=1`

This behavior is implemented in:
- `core/features/process_features/src/assignBarcodes.c`

Quick regression smoke (UCSF guides, `max_reads=100000`) showed:
- canonical matrix parity (`off` vs `on`) preserved
- `stats.txt`/`features.txt` identical
- timing essentially unchanged (`4.26s` off vs `4.30s` on)

## Policy Corrections vs Legacy (N-alternative ambiguity)
This section documents deliberate correctness corrections in feature matching after `N` substitution.

Legacy behavior mixed "not found" and "ambiguous" in a single `0` state and used conservative tie handling that can under-count valid assignments.

Corrected policy:
1. **Strict improvement promotion only**
   - A `best_feature=0` state is promoted only when a candidate has strictly smaller Hamming distance (`new_h < best_h`).
   - No equal-distance promotion.
2. **Strict improvement resets ambiguity**
   - When a strictly better candidate appears, previous ambiguity state is cleared.
3. **Equal-distance ambiguity rule**
   - At equal best distance, ambiguity is set when competing hits map to different features.
   - Same-feature equal-distance support is not forced ambiguous by default.
4. **Ambiguous treated as not-found for promotion**
   - Ambiguous signals do not get promoted as a match by themselves.

Rationale:
- The legacy approach can reject valid assignments by treating all equal-distance multi-alt evidence as ambiguous, even when support collapses to the same feature.
- With `N` as an unknown base, the correct feature-level treatment is to preserve unique best-feature evidence and only mark ambiguity on true feature competition.

Implementation:
- `core/features/process_features/src/assignBarcodes.c`
  - shared resolver `resolve_feature_over_n_alts(...)`

## Current Working Tree (relevant files)
Modified and not committed:
- `core/features/process_features/include/common.h`
- `core/features/process_features/include/globals.h`
- `core/features/process_features/src/io.c`
- `core/features/process_features/src/assignBarcodes.c`
- `core/features/process_features/src/barcode_match.c`
- `core/features/process_features/src/globals.c`
