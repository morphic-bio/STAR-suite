# Handoff: CRISPR guide calls missing for EmptyDrops-rescued cells (A375, STAR Suite 1.9.2)

Date: 2026-09-11. Status: **open, not debugged.** Separate from the LARRY feature-assignment slowdown
(`HANDOFF_LARRY_FEATURE_ASSIGNMENT_REGRESSION_20260911.md`).

## Summary

On A375 (paper arm P01), STAR Suite 1.9.2 writes a guide-call row only for OrdMag primary cells. The 35
cells added by EmptyDrops rescue are in the filtered GEX matrix but have **no row** in
`outs/crispr_analysis/protospacer_calls_per_cell.csv`. Six of them are cells Cell Ranger 9 calls a guide
for, so CRISPR call parity drops from 1,083/1,083 (v1.4.x, March) to 1,077/1,083. Their guide UMIs are
present and match Cell Ranger; only the call row is missing.

## Facts

| | v1.4.x (March, `/storage/A375/paper_bench_20260326_134444`) | 1.9.2 (`/storage/paper_bench_v192_20260911/out/P01`) |
| --- | --- | --- |
| GEX filtered cells | 1,188 | 1,170 = 1,135 OrdMag primaries + 35 EmptyDrops rescues |
| Rows in `protospacer_calls_per_cell.csv` | 1,188 (1,086 with a call, 102 None) | **1,135** (1,077 with a call, 58 None) |
| GEX cells with no call row | 0 | **35, exactly the 35 EmptyDrops rescues** |
| CR9 guide calls reproduced | 1,083 / 1,083 exact | 1,077 exact, **6 CR-called cells have no STAR row** |
| Guide UMI Pearson on common rows | 0.999994 | 0.999996 |

- Primary / rescue split is from `Solo.out/GeneFull/filtered/EmptyDrops/backend_debug_candidates.tsv`
  (`is_simple_cell`, `is_tail_candidate`). Primaries start at 3,662 UMIs; the 35 rescued cells have
  506-3,531 GEX UMIs.
- The six CR-called cells without a STAR row (CR call, CR UMIs / STAR guide UMIs in
  `outs/filtered_feature_bc_matrix`):
  `AAGCAAATCAGAGAAC` Non_Target-1_MS 73/73, `ATAGTGCGTACACTCC` RAB1A-2 217/244,
  `ATGGCTATCCGCATAA` Non_Target-1_MS 136/138, `CCTGTAGCAAAGTCTC` Non_Target-1_MS 32/32,
  `CTTGTGTTCATGCAGC` Non_Target-1_MS 96/97, `GGCCTGGGTGTGGTCA` Non_Target-1_MS 118/121.
  All six had the same call in the v1.4.x run.
- The guide library's own EmptyDrops set is not the cause: `cr_assign/.../crispr/filtered_barcodes.txt`
  has 742 cells on 1.9.2 and 743 on v1.4.x.
- Comparator: Cell Ranger 9.0.1 `multi`, `/storage/A375-CR-9.01/1k_CRISPR_5p_gemx_count_refmatch_2024a_fullraw`
  (count-layout view for the parity script: `/storage/paper_bench_v192_20260911/comparators/A375_CR9_multi_as_count`).
  Parity reports: `/mnt/pikachu/star_suite_paper/analysis/paper_benchmark_refresh_20260910/v192/P01/`.

## Not yet known

- Whether 1.9.1 already does this. P01 was never run on 1.9.1. 1.9.2 changed the non-Flex candidate floor
  (100 -> 500) and added a separate primary floor (`primary_umi_min`, 1 for non-Flex); 1.9.1 changed the
  shared caller to reuse the OrdMag result. Either could leave the guide caller reading the primary set
  (OrdMag `passingIndices`) instead of the final filtered set.

## Suggested approach

1. Run A375 once with the frozen v1.9.1 binary
   (`/mnt/pikachu/benchmark_build_v191_20260911/src/core/legacy/source/STAR`, sha256 `13a55ced...`) using the
   P01 command (`/storage/paper_bench_v192_20260911/out/P01/multi_config.csv`, wrapper
   `analysis/.../wrappers/run_a375_benchmark.v192.CORRECTED.sh` with `STAR_BIN` overridden). About 3 min.
   Count call rows against GEX cells.
2. Find which cell set the CR-compat CRISPR calling step takes as input, and make it the final filtered
   GEX callset (primaries plus rescues), as in v1.4.x.
3. Acceptance: one call row per filtered GEX cell; A375 CR9 call parity back to 1,083/1,083 exact; no
   change to GEX matrices or callsets.

## Constraints

- Do not modify the frozen build trees or `/storage/paper_bench_v192_20260911`; work in a new worktree.
- Commit messages carry no AI attribution; do not push or release without the owner.
- Paper benchmarks must use a released build; P01 is rerun after the fix is released.
