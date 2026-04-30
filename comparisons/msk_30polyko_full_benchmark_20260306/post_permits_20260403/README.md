# MSK 30polyKO DE — post-permits-fix STAR rerun (2026-04-03)

This is the canonical post-permits-fix DE benchmark, paired with the 2026-03-06
Cell Ranger 9 runs (CR is deterministic w.r.t. flags and data, so reusing the
same CR outputs is correct). It supersedes the pre-permits DE wall-time number
in the parent `../README.md` (42 min STAR / 4.0× speedup) — the permits fix
brought STAR DE wall time down from 42 min to **26.9 min**.

## Headline numbers (post-permits-fix)

| Metric | Value |
|---|---|
| STAR cells (Solo GeneFull filtered) | 33,095 |
| Cell Ranger cells (gRNA run, filtered) | 32,256 |
| Common cells | 32,248 |
| Barcode Jaccard | 0.9742 |
| Per-barcode UMI Pearson (Gene Expression) | 0.999903 |
| Per-feature Pearson (Gene Expression, all common) | 0.994554 |
| CRISPR set-equivalent calls | 98.04% (23,063 / 23,525) |
| CRISPR call UMI Pearson | 0.999708 |
| **STAR wall time** | **26.9 min** (1,612 s) |
| Cell Ranger total wall (gRNA 58 min + LARRY 110 min) | 168 min |
| **Speedup STAR vs CR total** | **6.24×** |

## Provenance

- STAR run: `/storage/MSK-perturb-comparison/paper_bench_emptydrops_guarded_redo_20260403_214718`
- CR gRNA run (paired):
  `/storage/MSK-perturb-comparison/cr_full_grna_30crispr_20260306_173247`
- CR LARRY run (paired):
  `/storage/MSK-perturb-comparison/canonical_tru_seq_cr_larry_20260306_052040`
- Same staged FASTQs as the 2026-03-06 canonical run (`msk30ko_full_3lib_20260304_095911/fastqs`)

## Files in this directory

| File | Description |
|---|---|
| `BENCHMARK_SUMMARY.txt` | STAR wall time / cells / config |
| `PARITY_vs_CR9.txt` | Per-feature-type parity (GEX, CRISPR, CRISPR call concordance) |
| `star_Log.final.out` | STAR alignment summary |
| `star_pf_multi_config.csv` | STAR pfMultiConfig (3-library) |
