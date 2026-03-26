# Handoff: Solo Benchmark Status (2026-03-26)

## Scope

This handoff is the current source of truth for:

- UCSF GEX-only legacy vs current STAR timing
- UCSF GEX-only CR9 parity and historical-surface comparison
- UCSF / MSK perturb benchmark reruns with Solo bridge enabled
- current `zcat` vs native-gzip conclusion

It is intended to stop the next agent from re-debugging already-settled issues.

## Repo / Branch

- Repo: `/mnt/pikachu/STAR-suite`
- Branch: `feature/flex-optimization-using-solo-20260325`

## Settled Facts

### 1. The UCSF GEX-only wrapper is not the blocker

Script:

- `/mnt/pikachu/STAR-suite/scripts/paper/run_ucsf_gexonly_no_bam_benchmark.sh`

On this branch, it already has the intended two-mode split:

- `--historical-vanilla`
- `--modern-optimized`

Known-good invariants already verified:

- exact-one mode enforcement
- historical mode uses `--readFilesCommand zcat`
- modern mode does not use `--readFilesCommand zcat`
- wrapper removes `${OUTDIR}/tmp` rather than pre-creating it

If someone reports the old broken behavior, assume wrong checkout / stale tree
first.

### 2. The `7a7fb08` historical binary can parse and run

The historical UCSF GEX-only run completed successfully once built with the
documented hygiene:

- `make clean`
- regenerate `parametersDefault.xxd`
- build `core/features/libscrna`
- build `STAR`

Completed historical run:

- `/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_vanilla_7a7fb08_retry2/`

This run was not a STAR failure. Earlier confusion came from short timeouts and
watching partial logs.

### 3. Native gzip is currently slower than `zcat` on this host

This is now a practical conclusion, not a hypothesis.

For UCSF GEX-only optimized current STAR:

- native gzip:
  - `/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_optimized_current_retry2/`
  - wall `15.78 min`
- `zcat`:
  - `/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_optimized_current_zcat_20260326/`
  - wall `13.75 min`

Conclusion:

- keep `zcat` in the benchmark wrappers for now
- internal gzip can be tuned later as a separate task

### 4. The missing perturb win was caused by missing Solo bridge flags

The earlier perturb wrappers did not actually exercise the optimized Solo path.

They were patched to add:

- `export STAR_SOLO_NONFLEX_HASH_BRIDGE=1`
- `--soloInlineHashMode yes`

Patched scripts:

- `/mnt/pikachu/STAR-suite/scripts/paper/run_ucsf_ebs2_2_benchmark.sh`
- `/mnt/pikachu/STAR-suite/scripts/paper/run_msk_30polyko_benchmark.sh`
- `/mnt/pikachu/STAR-suite/scripts/paper/run_a375_benchmark.sh`

All three were also returned to `zcat` for now.

A375 was intentionally left with its existing chemistry-specific settings:

- no forced `--clip3pPolyG yes`
- keep `--soloStrand Unstranded`

## Current Valid Benchmark Numbers

### UCSF GEX-only

Legacy STAR (`7a7fb08`):

- run:
  `/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_vanilla_7a7fb08_retry2/`
- wall: `23.92 min`
- filtered cells: `13,872`

CellGenI-style historical STAR (`7a7fb08`, latest TRU whitelist, `GeneFull` only):

- run:
  `/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_cellgeni_historical_7a7fb08_truwhitelist_genefullonly_20260326/`
- wall: `26.83 min`
- filtered cells: `13,847`

Current STAR, optimized, `zcat`:

- run:
  `/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_optimized_current_zcat_20260326/`
- wall: `13.75 min`
- filtered cells: `13,723`
- `collapseUMIall_fromBridgeHash 9.55009 s`
- `countCBgeneUMI 9.75733 s`
- `processRecords 236.669 s`

Current STAR, optimized, native gzip:

- run:
  `/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_optimized_current_retry2/`
- wall: `15.78 min`
- filtered cells: `13,728`

Conclusion:

- current optimized + `zcat` is the best validated UCSF GEX-only surface
- current optimized + `zcat` is faster than legacy by `10m10s`
- current optimized + `zcat` is faster than the CellGenI-style historical rerun
  by `13m05s`

### UCSF perturb

Old baseline:

- run:
  `/storage/paper_bench_solo_full_20260324/ucsf_ebs2_2_standard_baseline/`
- wall: `19.8 min`
- filtered cells: `13,721`

New bridge + `zcat` rerun:

- run:
  `/storage/solo_overnight_20260326/solo_bridge_zcat_reruns_20260326/ucsf_perturb_bridge_zcat_v2/`
- wall: `16.4 min`
- filtered cells: `13,719`
- `collapseUMIall_fromBridgeHash 9.69803 s`
- `countCBgeneUMI 9.91244 s`
- `processRecords 244.088 s`

Conclusion:

- UCSF perturb improved by `3m21s`

### MSK perturb

Native-gzip comparison run:

- run:
  `/storage/solo_overnight_20260326/msk_perturb_no_bam/msk30ko_no_bam_nativegzip_20260326_081323/`
- wall: `28.6 min`
- filtered cells: `30,567`

New bridge + `zcat` rerun:

- run:
  `/storage/solo_overnight_20260326/solo_bridge_zcat_reruns_20260326/msk_perturb_bridge_zcat_v2/`
- wall: `25.0 min`
- filtered cells: `30,557`
- `collapseUMIall_fromBridgeHash 12.5934 s`
- `countCBgeneUMI 12.9358 s`
- `processRecords 271.427 s`

Conclusion:

- MSK improved by `3m31s`

### A375

- No new bridge + `zcat` rerun is completed yet.
- Leave the current README A375 numbers alone until that run exists.

## Current Valid Parity Numbers

### Valid now

The perturb parity table in `README.md` is still the current valid one:

- A375:
  - Jaccard `0.976`
  - gene Pearson `0.975`
  - cell Pearson `1.000`
  - CRISPR match `100%`
- UCSF perturb:
  - Jaccard `0.976`
  - gene Pearson `0.995`
  - cell Pearson `1.000`
  - CRISPR match `98.9%`
- MSK:
  - Jaccard `0.942`
  - gene Pearson `0.994`
  - cell Pearson `1.000`
  - CRISPR match `99.4%`

### UCSF GEX-only vs CR9

Correct CR9 reference:

- `/storage/cr9_ebs2_2_benchmark_20260318/cr9_ebs2_2`

| Arm | Wall | Filtered cells | Cell Jaccard vs CR9 | Barcode Pearson vs CR9 | Gene Pearson vs CR9 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Current optimized `zcat` | `13:45` | `13,723` | `0.976483` | `0.999946` | `0.994885` |
| Legacy local vanilla `7a7fb08` | `23:55` | `13,872` | `0.985771` | `0.998653` | `0.886463` |
| CellGenI-style historical `7a7fb08` | `26:50` | `13,847` | `0.989120` | `0.999949` | `0.963561` |
| CR9 reference | `58 min` | `13,760` | `1.000000` | `1.000000` | `1.000000` |

Interpretation:

- current optimized `zcat` is best on runtime and gene parity
- both historical surfaces have stronger filtered-cell overlap with CR9 than
  current, but materially worse gene-level agreement
- CellGenI-style historical is a clear improvement over the older local legacy
  surface on gene parity, but still behind current

### UCSF GEX-only filtered-cell overlap across all arms

| Pair | Cells A | Cells B | Common | Jaccard |
| --- | ---: | ---: | ---: | ---: |
| Current vs legacy local vanilla | `13,723` | `13,872` | `13,540` | `0.963358` |
| Current vs CellGenI historical | `13,723` | `13,847` | `13,549` | `0.966336` |
| Legacy local vanilla vs CellGenI historical | `13,872` | `13,847` | `13,809` | `0.992739` |
| Current vs CR9 | `13,723` | `13,760` | `13,578` | `0.976483` |
| Legacy local vanilla vs CR9 | `13,872` | `13,760` | `13,717` | `0.985771` |
| CellGenI historical vs CR9 | `13,847` | `13,760` | `13,728` | `0.989120` |

## GEX-only Parity Helper Status

Helper:

- `/mnt/pikachu/STAR-suite/scripts/run_ucsf_gexonly_gex_parity_vs_cr.sh`

Current behavior:

- `CR_RUN` is now required explicitly
- the helper no longer silently defaults to the stale UCSF CR path
- the validated CR9 comparison for UCSF GEX-only now points at the real
  `EBs2_2` CR9 run:
  `/storage/cr9_ebs2_2_benchmark_20260318/cr9_ebs2_2`

Related narrow handoff:

- `/mnt/pikachu/STAR-suite/docs/HANDOFF_UCSF_GEXONLY_PARITY_FIX_20260326.md`

## README Status

`README.md` has already been updated with:

- UCSF perturb `16.4 min`
- MSK perturb `25.0 min`
- UCSF GEX-only current vs legacy vs CellGenI timing table
- UCSF GEX-only CR9 parity table and pairwise filtered-cell overlap
- updated speedups
- updated phase breakdowns
- note that native gzip was slower on the comparison surface
- replacement of the misleading UCSF GEX-only parity note with valid CR9
  results

Do not overwrite those timing updates with older numbers.

## What The Next Agent Should Do

### Task 1: finish A375 bridge + zcat rerun

Run the patched A375 wrapper serially and record the artifact root.

Target script:

- `/mnt/pikachu/STAR-suite/scripts/paper/run_a375_benchmark.sh`

Expected behavior:

- use `zcat`
- use `STAR_SOLO_NONFLEX_HASH_BRIDGE=1`
- use `--soloInlineHashMode yes`
- keep A375’s existing strand / poly-G choices

After completion:

- update `tests/ARTIFACTS.md`
- update `README.md` if the timing materially changes

### Task 2: carry the GEX-only result set into the durable benchmark records

1. Add the finalized UCSF GEX-only parity artifact paths to
   `tests/ARTIFACTS.md`.
2. If a future rerun changes the best current surface, rerun all three UCSF
   GEX-only comparisons against the same CR9 reference:

   ```bash
   CR_RUN=/storage/cr9_ebs2_2_benchmark_20260318/cr9_ebs2_2 \
   STAR_RUN=/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_optimized_current_zcat_20260326 \
   /mnt/pikachu/STAR-suite/scripts/run_ucsf_gexonly_gex_parity_vs_cr.sh
   ```

3. Keep the helper on an explicit `CR_RUN`; do not restore the stale implicit
   default.

## Files To Read First

- `/mnt/pikachu/STAR-suite/README.md`
- `/mnt/pikachu/STAR-suite/tests/ARTIFACTS.md`
- `/mnt/pikachu/STAR-suite/docs/HANDOFF_UCSF_GEXONLY_PARITY_FIX_20260326.md`
- `/mnt/pikachu/STAR-suite/scripts/run_ucsf_gexonly_gex_parity_vs_cr.sh`
- `/mnt/pikachu/STAR-suite/scripts/paper/run_ucsf_ebs2_2_benchmark.sh`
- `/mnt/pikachu/STAR-suite/scripts/paper/run_msk_30polyko_benchmark.sh`
- `/mnt/pikachu/STAR-suite/scripts/paper/run_a375_benchmark.sh`

## One-Line Summary

The timings are now in good shape, the perturb bridge reruns are real, `zcat`
is the current benchmark default, and UCSF GEX-only now has a validated CR9
comparison across current, legacy, and CellGenI-style historical surfaces.
