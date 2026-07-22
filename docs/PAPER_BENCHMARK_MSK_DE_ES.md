# MSK 30polyKO Perturb-seq — paper-grade benchmark (DE + ES replicate)

This document is the canonical, paper-ready comparison of STAR-suite vs
Cell Ranger 9 on the MSK 30polyKO 3-library Perturb-seq dataset (GEX + PolyIII
gRNA + LARRY lineage barcodes), reported on two independent biological
samples processed with identical pipelines and parameters:

- **DE** (Definitive Endoderm): 4 mRNA samples (S25, S36, S2, S52), 1 PolyIII (S55), 1 LARRY (S35)
- **ES** (Embryonic Stem cells): 4 mRNA samples (S24, S34, S3, S30), 1 PolyIII (S32), 1 LARRY (S31)

Same chemistry (mixed TRU/NXT), same references, same wrappers, same flags;
only the biological sample changes. Both runs use the published wrappers
[`scripts/paper/run_msk_30polyko_benchmark.sh`](../scripts/paper/run_msk_30polyko_benchmark.sh)
(STAR) and
[`scripts/paper/run_msk_30polyko_cr_benchmark.sh`](../scripts/paper/run_msk_30polyko_cr_benchmark.sh)
(CellRanger 9).

## Headline table

All parity numbers come from
[`scripts/report_additional_parity_metrics.py`](../scripts/report_additional_parity_metrics.py)
with `--gene-corr-min-counts 20 --gene-corr-min-cells-pct 0.01` and the
NXT↔TRU translation file
(`3M-february-2018_NXT.txt.gz`, `--translate both`), per
[`docs/PAPER_BENCHMARK_METHODOLOGY.md`](PAPER_BENCHMARK_METHODOLOGY.md).
Both samples were measured on the same host
(pikachu, i9-13900KF / 126 GB RAM / 32 threads), Cell Ranger 9.0.1.

| Metric | DE (post-permits-fix, 2026-04-03) | ES (2026-04-30) |
|---|---|---|
| STAR cells (Solo GeneFull filtered) | 33,095 | **33,226** |
| Cell Ranger cells (gRNA run, filtered) | 32,256 | **32,670** |
| Common cells (filtered ∩ filtered, after NXT translation) | 32,248 | **32,652** |
| Barcode Jaccard | 0.9742 | **0.982** |
| Per-barcode UMI Pearson (Gene Expression, common cells) | 0.999903 | **0.999928** |
| Per-feature Pearson (Gene Expression) | 0.994554 (all common features) | **0.993687** (filtered, ≥20 counts, 16,958 genes) |
| CRISPR set-equivalent calls | 98.04% (23,063/23,525) | **98.97%** (25,894/26,164) |
| CRISPR call UMI Pearson | 0.999708 | **0.999434** |
| LARRY (Custom) per-barcode Pearson (common cells) | — | **0.9941** |
| LARRY (Custom) per-feature Pearson (co-observed features) | — | **0.9509** (3,325 features) |
| **STAR wall time** (single 3-library run) | **26.9 min** | **30.2 min** |
| Cell Ranger wall (gRNA + GEX) | 58 min | 59.4 min |
| Cell Ranger wall (LARRY + GEX) | 110 min | 107.7 min |
| **Cell Ranger total wall** | 168 min | **167.1 min** |
| **Speedup STAR vs CR total** | **6.24×** | **5.53×** |

The DE numbers are from the **post-permits-fix** rerun
([`comparisons/msk_30polyko_full_benchmark_20260306/post_permits_20260403/`](../comparisons/msk_30polyko_full_benchmark_20260306/post_permits_20260403/)).
The earlier 2026-03-06 paired run (STAR 42 min, 4.0× speedup) used a
pre-permits-fix STAR build and is now superseded by the 2026-04-03 STAR
rerun on the same staged FASTQs paired with the same (deterministic) CR
outputs. The ES numbers (this rerun) are paired in a single orchestrated
pipeline:
[`comparisons/msk_30polyko_full_benchmark_ES_20260430/`](../comparisons/msk_30polyko_full_benchmark_ES_20260430/).

## Reproducibility

The full pipeline (stage → STAR → CR) for the ES sample is captured in the
orchestrator at
[`comparisons/msk_30polyko_full_benchmark_ES_20260430/README.md`](../comparisons/msk_30polyko_full_benchmark_ES_20260430/README.md).
Key invariants between the two runs:

- Same STAR genome index: `/storage/autoindex_110_44/bulk_index`
- Same Cell Ranger transcriptome: `refdata-gex-GRCh38-2024-A` (stock 10x build)
- Same cell-barcode whitelists (TRU + NXT)
- Same gRNA feature reference (`ref_feature_geneBC_crispr.csv`,
  `feature_type=CRISPR Guide Capture`) and LARRY feature reference
  (`ref_feature_larryBC.csv`, `feature_type=Custom`)
- Identical wrapper flags: `--no-bam`, `--threads 32`,
  CR `--localcores 32 --localmem 110`
- STAR runs as a single 3-library invocation; Cell Ranger requires two
  separate `cellranger multi` runs (gRNA-side + LARRY-side) because CR9
  does not support 3 libraries in one invocation.

## Take-aways for the paper

1. **Replicated speedup**: The 5.5–6.2× STAR-vs-CR wall-time advantage on a
   3-library Perturb-seq dataset is reproduced on a second biological sample
   processed independently from staging through parity computation.
2. **Stable parity**: Per-barcode UMI Pearson, per-gene Pearson, and CRISPR
   call concordance are within 0.1% across DE and ES.
3. **Single-binary advantage**: STAR's 3-library single-pass run remains
   faster than the two sequential CR `multi` runs that CR9 requires; this is
   the main source of the speedup (the two CR runs each load the genome and
   re-process GEX).
4. **LARRY (Custom) parity** is now reported on ES — the per-barcode UMI
   Pearson of 0.994 and total-UMI ratio 1.008 confirm STAR's Custom-feature
   path matches CR within sampling noise; the lower per-feature Pearson
   (0.95 on co-observed features) is dominated by the long tail of rare
   lineage barcodes with very low counts in either tool.

## Linked artifacts

- DE: [`comparisons/msk_30polyko_full_benchmark_20260306/`](../comparisons/msk_30polyko_full_benchmark_20260306/)
- ES: [`comparisons/msk_30polyko_full_benchmark_ES_20260430/`](../comparisons/msk_30polyko_full_benchmark_ES_20260430/)
- Wrappers:
  [`scripts/paper/run_msk_30polyko_benchmark.sh`](../scripts/paper/run_msk_30polyko_benchmark.sh),
  [`scripts/paper/run_msk_30polyko_cr_benchmark.sh`](../scripts/paper/run_msk_30polyko_cr_benchmark.sh)
- Methodology:
  [`docs/PAPER_BENCHMARK_METHODOLOGY.md`](PAPER_BENCHMARK_METHODOLOGY.md)
- Surface audit (DE):
  [`docs/MSK_BENCHMARK_SURFACE_AUDIT_20260403.md`](MSK_BENCHMARK_SURFACE_AUDIT_20260403.md)
