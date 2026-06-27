# Paper Benchmarks: STAR-suite Multi-Feature Performance (2026-03-18)

**Date**: 2026-03-18
**Branch**: `multi-feature` (post MSK CRISPR master repair + PfMultiMerge streaming optimization)
**Host**: pikachu (i9-13900KF, 126 GB RAM, 32 threads)
**STAR version**: 2.7.11b (compiled 2026-03-18T06:31:31+00:00)

## Datasets

| Dataset | Libraries | Chemistry | Reads | Expected cells |
|---|---|---|---|---|
| A375 1k CRISPR 5' GemX | GEX + CRISPR (2) | TRU (single-column WL) | 47M | ~1,200 |
| UCSF EBs2_2 Perturb-seq | GEX + CRISPRa guides (2) | NXT→TRU (2-column WL) | 445M | ~14,000 |
| MSK 30polyKO | GEX + gRNA + LARRY (3; 245,979 LARRY barcodes) | Mixed TRU/NXT (per-library WL) | 669M | ~30,000 |

## Wall-Time Summary

| Dataset | Threads | BAM | Wall time | Cells | writeCombinedMex (raw/filt) |
|---|---|---|---|---|---|
| A375 | 32 | none | **4.0 min** (241s) | 1,187 | 2.0s / 1.5s |
| UCSF EBs2_2 | 32 | none | **19.0 min** (1141s) | 13,721 | 19.5s / 15.0s |
| MSK 30polyKO | 32 | none | **27.6 min** (1656s) | 30,567 | 35.1s / 24.4s |

## Phase Breakdown

### A375 (47M reads, 2 libraries, 38K + 11 features)

| Phase | Start | End | Duration |
|---|---|---|---|
| Genome load | 06:52:18 | 06:53:02 | 44s |
| Feature assignment | 06:53:02 | 06:53:34 | 32s |
| Mapping | 06:53:10 | 06:54:31 | 81s |
| Solo counting | 06:54:32 | 06:55:41 | 69s |
| PfMulti merge + CRISPR calling | 06:55:41 | 06:55:51 | 10s |

### UCSF EBs2_2 (445M reads, 2 libraries, 38K + 548 features)

| Phase | Start | End | Duration |
|---|---|---|---|
| Genome load | 07:00:28 | 07:01:16 | 48s |
| Feature assignment | 07:01:16 | 07:05:33 | 4m 17s |
| Mapping | 07:01:16 | 07:09:47 | 8m 31s |
| Solo counting | 07:09:48 | 07:17:34 | 7m 46s |
| PfMulti merge + CRISPR calling | 07:17:34 | 07:18:54 | 1m 20s |

### MSK 30polyKO (669M reads, 3 libraries, ~38K genes + 30 gRNA features + 245,979 LARRY barcodes)

| Phase | Start | End | Duration |
|---|---|---|---|
| Genome load | 07:29:39 | 07:30:27 | 48s |
| Feature assignment | 07:30:27 | 07:50:17 | 19m 50s |
| Mapping | 07:30:27 | 07:45:09 | 14m 42s |
| Solo counting | 07:45:10 | 07:54:39 | 9m 29s |
| PfMulti merge + CRISPR calling | 07:54:39 | 07:56:38 | 1m 59s |

Note: Feature assignment and mapping run concurrently via `dynamicThreadInterface`.

## Parity vs CellRanger 9

All parity metrics computed with `scripts/report_additional_parity_metrics.py`
using `--gene-corr-min-counts 20 --gene-corr-min-cells-pct 0.01` per
`docs/PAPER_BENCHMARK_METHODOLOGY.md`. CR9 references use
`refdata-gex-GRCh38-2024-A` (gencode v44, mkref 8.0.0).

| Dataset | Cells (STAR / CR) | Jaccard | Gene Pearson | Cell Pearson | CRISPR match | UMI Pearson |
|---|---|---|---|---|---|---|
| A375 | 1,187 / 1,162 | 0.976 | 0.975 (15,673 genes) | 0.9995 (1,160 BCs) | 100% (1,083/1,083) | 1.000 |
| UCSF EBs2_2 | 13,721 / 13,760 | 0.976 | 0.995 (18,061 genes) | 1.000 (13,571 BCs) | 98.9% (11,902/12,038) | 0.999 |
| MSK 30polyKO | 30,567 / 32,256 | 0.942 | 0.994 (17,460 genes) | 1.000 (30,481 BCs) | 99.4% (23,210/23,341) | 1.000 |

### CR9 Reference Runs

| Dataset | CR9 run ID | Reference | CR9 wall time |
|---|---|---|---|
| A375 | `1k_CRISPR_5p_gemx_count_refmatch_2024a_fullraw` | `refdata-gex-GRCh38-2024-A` | ~15 min |
| UCSF EBs2_2 | `cr9_ebs2_2` (run 2026-03-18) | `refdata-gex-GRCh38-2024-A` | 58 min |
| MSK 30polyKO | `cr9_starindex_grna` (GEX+gRNA only) | `refdata-gex-GRCh38-autoindex11044` | ~58+110 min |

Full parity output files: `{a375,ucsf_ebs2_2,msk_30polyko}/parity_vs_cr9.txt`.

## Regression Check vs Previous Benchmarks (2026-03-17)

### A375: No regression

| Metric | Old (Mar 17) | New (Mar 18) |
|---|---|---|
| Wall time | 245s | 241s (-2%) |
| Cells | 1,187 | 1,187 |
| Filtered MEX delta (old vs new) | — | 1 count / 24.4M total |
| CRISPR calls (old vs new) | — | 1187/1188 lines identical (1 off by 1 UMI) |
| GEX vs CR ref (filtered, common BCs) | delta 19,284 | delta 19,284 (identical) |
| CRISPR vs CR ref (filtered, common BCs) | delta 48,574 | delta 48,575 (off by 1) |

### MSK 30polyKO: Major CRISPR improvement (namespace fix)

| Metric | Old (Mar 17) | New (Mar 18) |
|---|---|---|
| Wall time | 2516s (with BAM) | 1656s (no BAM) |
| Cells | 30,520 | 30,567 (+47, ED variance) |
| CRISPR: cells with 0 molecules | **28,550 (93.5%)** | **233 (0.8%)** |
| CRISPR: cells with 1+ call | **691** | **23,546** |
| Filtered MEX barcodes (old vs new) | 30,520 | 30,567 |
| Common barcodes | — | 30,499 |
| Total feature count delta (old vs new) | — | +5.0M (from recovered CRISPR) |

The old run suffered from the NXT namespace bug: gRNA barcodes were matched in
the wrong namespace, causing 93.5% of cells to show zero CRISPR guide molecules.
The new run, with per-library whitelist support and deterministic namespace
normalization, correctly assigns guides to 23,546 cells.

### UCSF EBs2_2: First run on this branch (no prior baseline)

Establishes the baseline for 2-library NXT perturb-seq with CRISPRa v2 guides.

## Solo GEX Statistics

### A375

| Metric | Value |
|---|---|
| Reads | 47,095,182 |
| Valid barcodes | 89.9% |
| Sequencing saturation | 21.9% |
| Uniquely mapped | 74.1% |
| Cells | 1,187 |
| Median UMI/cell | 17,885 |
| Median genes/cell | 5,562 |

### UCSF EBs2_2

| Metric | Value |
|---|---|
| Reads | 444,896,731 |
| Valid barcodes | 97.6% |
| Sequencing saturation | 29.8% |
| Uniquely mapped | 96.6% |
| Cells | 13,721 |
| Median UMI/cell | 15,431 |
| Median genes/cell | 5,223 |

### MSK 30polyKO

| Metric | Value |
|---|---|
| Reads | 668,705,043 |
| Valid barcodes | 95.8% |
| Sequencing saturation | 31.8% |
| Uniquely mapped | 93.3% |
| Cells | 30,567 |
| Median UMI/cell | 8,408 |
| Median genes/cell | 3,933 |

## CRISPR GMM Calling Summary

### A375 (11 guides, minUMI=10)

| Metric | Value |
|---|---|
| Total cells | 1,187 |
| Cells with 0 molecules | 16 |
| Cells with no call | 85 |
| Cells with 1 feature | 1,051 |
| Cells with >1 features | 35 |

### UCSF EBs2_2 (548 CRISPRa guides, minUMI=3)

| Metric | Value |
|---|---|
| Total cells | 13,721 |
| Cells with 0 molecules | 435 |
| Cells with no call | 1,191 |
| Cells with 1 feature | 1,183 |
| Cells with >1 features | 10,912 |

### MSK 30polyKO (30 gRNA guides, minUMI=2)

| Metric | Value |
|---|---|
| Total cells | 30,567 |
| Cells with 0 molecules | 233 |
| Cells with no call | 6,788 |
| Cells with 1 feature | 20,972 |
| Cells with >1 features | 2,574 |

## Key Code Changes Since Previous Benchmarks

1. **Namespace/whitelist correctness** (primary): Per-library `star_whitelist`
   in pfMultiConfig, deterministic NXT→TRU normalization in assignBarcodes and
   PfMultiMerge. Fixes the MSK gRNA zero-molecule bug.

2. **PfMultiMerge streaming optimization** (secondary): Direct gzip writes via
   `gzwrite` replacing the write-plaintext-then-recompress cycle. Vector-based
   O(1) barcode remapping replacing `std::map`. `unordered_map` for barcode
   lookups. Zero-count feature pruning.

3. **assignBarcodes fastHamming hardening**: `pf_hamming_search_fasthamming`
   brute-force fallback only activates for `maxHammingDistance > 2`, preventing
   performance regression on standard prehash tiers.

## PE Bulk Benchmark (Integrated STAR-suite vs External Stepwise Pipeline)

Benchmark script: `scripts/paper/run_pe_bulk_feature_benchmark.sh`
Dataset: JAX PE (21033-09-01-13-01\_S1\_L007), full sample on `/storage`, 32 threads.
STAR index: `/storage/autoindex_110_44/bulk_index`

Important: the archived March 2026 tables below were collected in parity-QC
mode. The integrated arm emitted `Aligned.toTranscriptome.out.bam` and ran
integrated Salmon QC, so the archived "Total" values include a parity artifact
that normal STAR-suite production does not require. Current production-mode
benchmarking uses internal TranscriptVB only on the integrated arm; use
`--parity-qc` only when Salmon parity artifacts are intentionally needed.

### With Y-removal (2026-03-10)

| Step | Integrated | External |
|---|---|---|
| STAR (trim + align + Y-split + TranscriptVB) | 29.20s | — |
| Salmon QC | 31.51s | 31.03s |
| Decompress | — | 13.44s |
| Trimvalidate | — | 16.89s |
| STAR align | — | 49.94s |
| remove\_y\_reads | — | 14.03s |
| **Total** | **60.71s** | **125.33s** |
| **Speedup** | **2.1x** | — |

### Without Y-removal (2026-03-18)

| Step | Integrated | External |
|---|---|---|
| STAR (trim + align + TranscriptVB) | 5.69s | — |
| Salmon QC | 31.00s | 30.49s |
| Decompress | — | 12.79s |
| Trimvalidate | — | 14.36s |
| STAR align | — | 29.65s |
| **Total** | **36.69s** | **87.29s** |
| **Speedup** | **2.4x** | — |

Note: The low 5.69s integrated STAR time in the no-Y run reflects a warm page cache
(genome loaded by the preceding downsampled stage). The Y-removal run's 29.20s
is a cold-cache measurement and is more representative for single-invocation use.

### Quantification Parity (no Y-removal, storage)

| Comparison | Transcript Pearson | Gene Pearson |
|---|---|---|
| TranscriptVB vs integrated Salmon | 0.995 | 0.997 |
| Integrated Salmon vs external Salmon | 1.000 | 0.997 |
| TranscriptVB vs external Salmon | 0.995 | 0.997 |

Artifacts: `/tmp/pe_bulk_feature_benchmark_no_yremove_20260318_144657/`

### Rerun Recipe Notes (2026-06-27)

The paper recipe now treats STAR-suite internal TranscriptVB as the production
quantifier. The integrated arm does not emit TranscriptomeSAM or run integrated
Salmon QC by default. Those parity artifacts are enabled only with
`--parity-qc`; this wrapper does not use a FIFO/streaming Salmon contract. The
integrated arm passes `--transcriptomeFasta` alongside `--quantVBgcBias 1` so
the GC-bias and effective-length model use the same transcriptome FASTA as
Salmon.

For headline speedups, use `scripts/paper/run_pe_bulk_feature_benchmark.sh` and
keep both arms on the wrapper's matched output mode. Use `--integrated-only`
when refreshing only the STAR-suite production timing and the external control
does not need to be rerun.

A corrected full PPARG no-Y STAR-suite-only production run on `/storage`
measured integrated trim+align+sorted BAM+internal TranscriptVB at `8:54.52`.
It did not emit integrated TranscriptomeSAM or run integrated Salmon QC. Run
root: `/storage/JAX_PE/results/pparg_prod_benchmark_no_y_20260627_172349/`.

An earlier PPARG no-Y sanity rerun with a lean unsorted-BAM serial comparator
measured STAR-suite integrated at `7:18.06`, Salmon QC at `1:00.08`, and the
lean serial chain at `9:50.17` (`1.35x`). That diagnostic confirms the current
parity path but is not an apples-to-apples replacement for the production
benchmark because the integrated arm still emitted the parity transcriptome BAM.

## File Inventory

```
paper_benchmarks_20260318/
├── README.md                          (this file)
├── compiled_stats.tsv                 (machine-readable summary)
├── scripts/                           (analysis & comparison tools)
│   ├── compare_feature_mex.py         (MEX-level parity comparison)
│   ├── compute_parity_metrics.py      (Jaccard/Pearson/Spearman/CRISPR)
│   ├── report_additional_parity_metrics.py  (canonical parity script)
│   ├── run_pe_bulk_feature_benchmark.sh     (PE bulk feature benchmark)
│   └── gather_pe_bulk_external_tools.sh     (external tool runner)
├── a375/
│   ├── BENCHMARK_SUMMARY.txt
│   ├── Log.final.out
│   ├── pf_multi_config.csv
│   ├── RUN_COMMAND.sh                 (exact STAR invocation)
│   ├── run_a375_benchmark.sh          (benchmark script snapshot)
│   ├── star_solo_summary.csv
│   ├── protospacer_calls_summary.csv
│   ├── protospacer_umi_thresholds.csv
│   ├── parity_vs_cr9.txt             (canonical parity metrics)
│   ├── phase_timings.txt
│   └── dynamic_thread_telemetry.txt
├── ucsf_ebs2_2/
│   ├── BENCHMARK_SUMMARY.txt
│   ├── Log.final.out
│   ├── pf_multi_config.csv
│   ├── RUN_COMMAND.sh
│   ├── run_ucsf_ebs2_2_benchmark.sh
│   ├── star_solo_summary.csv
│   ├── protospacer_calls_summary.csv
│   ├── protospacer_umi_thresholds.csv
│   ├── parity_vs_cr9.txt             (canonical parity metrics)
│   ├── phase_timings.txt
│   └── dynamic_thread_telemetry.txt
├── msk_30polyko/
│   ├── BENCHMARK_SUMMARY.txt
│   ├── Log.final.out
│   ├── pf_multi_config.csv
│   ├── RUN_COMMAND.sh
│   ├── run_msk_30polyko_benchmark.sh
│   ├── star_solo_summary.csv
│   ├── protospacer_calls_summary.csv
│   ├── protospacer_umi_thresholds.csv
│   ├── parity_vs_cr9.txt             (canonical parity metrics)
│   ├── phase_timings.txt
│   └── dynamic_thread_telemetry.txt
└── pe_bulk/
    ├── BENCHMARK_SUMMARY_yremove.txt  (with Y-removal, 2026-03-10)
    ├── BENCHMARK_SUMMARY_no_yremove.txt (without Y-removal, 2026-03-18)
    └── comparison_metrics_no_yremove.tsv
```

## Full Run Outputs

Complete STAR outputs (Solo MEX, outs/filtered_feature_bc_matrix, cr_assign,
crispr_analysis, logs) are archived at:

```
/mnt/pikachu/paper_bench_rerun_20260318_065211/
├── a375/           (313 MB)
├── ucsf_ebs2_2_standard/  (2.8 GB)
└── msk_30polyko/   (4.9 GB)
```

## Reproducibility

Each subdirectory contains a `RUN_COMMAND.sh` with the exact STAR invocation.
The benchmark scripts are also included for the full wrapper (FASTQ discovery,
multi-config generation, output validation).

Build: `make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR`
Branch: `multi-feature` at commit used for this run.
