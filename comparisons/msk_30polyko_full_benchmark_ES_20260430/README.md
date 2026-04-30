# MSK 30polyKO Full Benchmark — ES sample (STAR vs Cell Ranger 9)

**Date**: 2026-04-30
**Dataset**: MSK 30-guide polyKO Perturb-seq, ES sample (Embryonic Stem cells)
**Source FASTQs**: `/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_ES_DE_XM` — `*_ES_*_R[12]_001.fastq.gz`
**Host**: pikachu (i9-13900KF, 126 GB RAM, 32 threads)

This is the second full-dataset benchmark in the MSK 30polyKO series; the
companion DE-sample benchmark is in
[`msk_30polyko_full_benchmark_20260306`](../msk_30polyko_full_benchmark_20260306/).
Both runs use the same chemistry, references, and tool versions; the only
change is the sample (DE → ES).

## Summary

All parity numbers below come from `scripts/report_additional_parity_metrics.py`
with `--gene-corr-min-counts 20 --gene-corr-min-cells-pct 0.01` and the
NXT↔TRU translation file
(`3M-february-2018_NXT.txt.gz`, `--translate both`), matching the DE-precedent
methodology. CRISPR concordance excludes STAR rows with `feature_call=None`.

| Metric | DE (precedent) | ES (this run) |
|---|---|---|
| STAR cells (Solo GeneFull filtered) | 30,497 | **33,226** |
| Cell Ranger cells (gRNA run, filtered) | 32,256 | **32,670** |
| Common cells (filtered ∩ filtered, after NXT translation) | 30,417 | **32,652** |
| Barcode Jaccard (STAR ∩ CR / STAR ∪ CR) | 0.94 | **0.982** |
| Per-barcode UMI Pearson (filtered_vs_filtered) | 0.9999 | **0.999928** |
| Per-gene Pearson (filtered, min_counts=20, min_cells_pct=0.01) | 0.993 | **0.993687** |
| CRISPR set-equivalent calls (eval cells in both, % match) | 98.5% (22,200/22,531) | **98.97%** (25,894/26,164) |
| CRISPR call UMI Pearson | 0.999 | **0.999434** |
| LARRY (Custom) per-barcode Pearson (32,637 common cells) | — | **0.9941** |
| LARRY (Custom) per-feature Pearson (3,325 co-observed features) | — | **0.9509** |
| **STAR wall time (single 3-lib run)** | 26.9 min | **30.2 min** |
| CR wall time (gRNA + GEX) | 58 min | **59.4 min** |
| CR wall time (LARRY + GEX) | 110 min | **107.7 min** |
| **CR total wall time** | 168 min | **167.1 min** |
| **Speedup STAR vs CR total** | 4.0× | **5.5×** |

## Run inputs

| Item | Path |
|---|---|
| GEX FASTQs | `/storage/MSK-perturb-comparison/msk30ko_full_3lib_ES_20260430_054814/fastqs/mRNA` (4 samples: S24, S34, S3, S30) |
| gRNA (PolyIII) FASTQs | `.../fastqs/PolyIII` (S32) |
| LARRY FASTQs | `.../fastqs/LARRY` (S31) |
| STAR genome | `/storage/autoindex_110_44/bulk_index` |
| CR transcriptome | `/storage/A375-CR-9.01/refdata-gex-GRCh38-2024-A` (stock 10x GRCh38-2024-A) |
| Cell barcode whitelist (TRU) | `/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt` |
| gRNA whitelist (NXT) | `/storage/scRNAseq_output/whitelists/3M-february-2018_NXT.txt` |
| gRNA feature ref (CR) | `/mnt/pikachu/MSK-whitelists/ref_feature_geneBC_crispr.csv` (`feature_type=CRISPR Guide Capture`) |
| LARRY feature ref | `/mnt/pikachu/MSK-whitelists/ref_feature_larryBC.csv` (`feature_type=Custom`) |
| Cell Ranger | `/home/lhhung/cellranger-9.0.1/bin/cellranger` |
| STAR | `core/legacy/source/STAR` (master @ commit captured in `RUN_COMMAND.sh`) |

## STAR run

- Wall time: **30.2 min** (1,811 s), 32 threads, `--outSAMtype None`
- Input reads: 626,412,313
- Uniquely mapped: 93.16%
- Multi-mapped: 3.56% (54.8M reads enter CR-compat rescue, 22.0M rescued)
- Unmapped (too short): 1.69%
- Estimated cells: 33,226 (Solo GeneFull EmptyDrops_CR)
- Median UMI / cell: 7,850; Median GeneFull / cell: 3,735; Mean reads / cell: 11,906
- Single-pass 3-library mode: GEX + PolyIII (CRISPR Guide Capture) + LARRY (Custom)
- Wrapper: [`scripts/paper/run_msk_30polyko_benchmark.sh`](../../scripts/paper/run_msk_30polyko_benchmark.sh)

See `star_Log.final.out`, `star_solo_summary.csv`, `star_pf_multi_config.csv`.

## Cell Ranger 9 runs

CR9 cannot combine 3 libraries in a single `multi` invocation, so two
sequential `cellranger multi` runs were executed against the same staged
fastqs:

| Run | Wall | Output |
|---|---|---|
| `cr9_gex_grna` (GEX + CRISPR Guide Capture) | 59.4 min (3,563 s) | `cr_gex_grna/cr9_gex_grna/outs/` |
| `cr9_gex_larry` (GEX + Custom)              | 107.7 min (6,460 s) | `cr_gex_larry/cr9_gex_larry/outs/` |
| **Total**                                   | **167.1 min**       | |

CR command (from `cr_grna_time.log`):
```
cellranger multi --id cr9_gex_grna --csv multi.csv --localcores 32 --localmem 110
```

Both `multi.csv` configs use `create-bam,false` and
`no-secondary-analysis,true`. The wrapper that generated them is
[`scripts/paper/run_msk_30polyko_cr_benchmark.sh`](../../scripts/paper/run_msk_30polyko_cr_benchmark.sh).

See `cr_grna_multi_config.csv`, `cr_larry_multi_config.csv`,
`cr_grna_time.log`, `cr_larry_time.log`.

## Parity detail

### GEX (Gene Expression)

Authoritative numbers from `report_additional_parity_metrics.py` (NXT↔TRU
translation, `--translate both`). Both CR runs measure GEX on the same library;
the gRNA-side run is the canonical comparison.

- **Common barcodes (filtered_vs_filtered, after translation)**: 32,652 (33,226 STAR / 32,670 CR)
- **Per-barcode UMI Pearson**: **0.999928** (Spearman 0.999930)
- **Per-gene Pearson** (filtered, min_counts=20, ≥1% common cells; 16,958 filtered genes): **0.993687** (Spearman 0.991760)
- **Per-gene Pearson on all 38,606 common genes (no filter)**: **0.993769**
- CR filtered total UMIs: 262,847,658; STAR: 265,865,263 (ratio 1.011)
- Raw matrix per-barcode Pearson across 2,021,731 common raw barcodes: **0.999984**

Secondary script (`compute_parity_metrics.py`, no min-count filter) reports
**Barcode Jaccard 0.982** for the filtered ∩ filtered set.

### CRISPR Guide Capture (gRNA / PolyIII)

- CR rows: 26,185; STAR non-None rows: 27,138; common evaluated rows: **26,164**
- **Set-equivalent calls (excluding STAR None)**: **25,894 / 26,164 = 98.97%**
- **CRISPR call UMI Pearson**: **0.999434** (Spearman 0.997804)
- 29 features common (CR has 29 / STAR 30 — STAR includes one additional non-targeting feature)
- Per-feature UMI Pearson on 29 common: **0.9999** (from `compute_parity_metrics.py`)
- Total assigned UMIs: CR 12,774,427; STAR 13,402,328 (ΔSTAR = +627,901)

### LARRY (Custom)

Compared from the LARRY-side CR run (`cr9_gex_larry`):

- 245,979 features in both refs; 3,325 features observed with non-zero counts in both
- **Per-barcode UMI Pearson**: **0.9941** (32,637 common cells)
- **Per-feature Pearson on co-observed features**: **0.9509** (Spearman 0.9585)
- CR total LARRY UMIs: 10,506,672; STAR: 10,588,558 (ratio 1.008)

The lower per-feature Pearson on LARRY is expected: LARRY is a 245k-feature
sparse barcode space; most features have very few counts in either tool, so
the correlation is dominated by sampling noise on rare lineage barcodes.
Per-barcode and total-UMI parity are excellent.

## Reproducing this benchmark

End-to-end orchestration is `scripts/paper/run_msk_30polyko_benchmark.sh` +
`scripts/paper/run_msk_30polyko_cr_benchmark.sh`, with shared `MSK_FASTQ_ROOT`:

```bash
STAGING=/storage/MSK-perturb-comparison/msk30ko_full_3lib_ES_<ts>
# (1) Stage fastqs into ${STAGING}/fastqs/{mRNA,PolyIII,LARRY} with the
#     <library>_S<n>_L<lane>_R<read>_001.fastq.gz naming used by both tools.

# (2) STAR
MSK_FASTQ_ROOT=${STAGING}/fastqs \
MSK_OUTDIR=/storage/MSK-perturb-comparison/paper_bench_ES_<ts> \
  scripts/paper/run_msk_30polyko_benchmark.sh --no-bam --threads 32

# (3) CellRanger 9 (two sequential `multi` invocations)
MSK_FASTQ_ROOT=${STAGING}/fastqs \
MSK_OUTDIR=/storage/MSK-perturb-comparison/cr_paper_bench_ES_<ts> \
  scripts/paper/run_msk_30polyko_cr_benchmark.sh --no-bam --threads 32 --localmem 110

# (4a) Authoritative parity (matches PAPER_BENCHMARK_METHODOLOGY.md). Requires
#      a CR-standard `outs/{raw,filtered}_feature_bc_matrix` layout; CR9 `multi`
#      hides those under per_sample_outs/, so symlink the wrapper layout first:
WRAP=/tmp/cr_es_grna_wrapper && mkdir -p "$WRAP/outs"
ln -s <CR_GRNA>/outs/multi/count/raw_feature_bc_matrix             "$WRAP/outs/raw_feature_bc_matrix"
ln -s <CR_GRNA>/outs/per_sample_outs/cr9_gex_grna/count/sample_filtered_feature_bc_matrix \
      "$WRAP/outs/filtered_feature_bc_matrix"
ln -s <CR_GRNA>/outs/per_sample_outs/cr9_gex_grna/count/crispr_analysis "$WRAP/outs/crispr_analysis"

python3 scripts/report_additional_parity_metrics.py \
  --cr-run "$WRAP" --star-run <STAR> \
  --barcode-translation /home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt.gz \
  --translation-direction left-to-right --translate both \
  --gene-corr-min-counts 20 --gene-corr-min-cells-pct 0.01

# (4b) Secondary breakdown (per-feature-type Jaccard, LARRY)
python3 scripts/paper/compute_parity_metrics.py \
  --cr-mex   <CR_GRNA>/outs/per_sample_outs/cr9_gex_grna/count/sample_filtered_feature_bc_matrix \
  --star-mex <STAR>/outs/filtered_feature_bc_matrix \
  --cr-crispr   <CR_GRNA>/outs/per_sample_outs/cr9_gex_grna/count/crispr_analysis/protospacer_calls_per_cell.csv \
  --star-crispr <STAR>/outs/crispr_analysis/protospacer_calls_per_cell.csv \
  --label MSK_30polyKO_ES_gRNA --tsv-out parity_grna.tsv
```

## Files in this directory

| File | Description |
|---|---|
| `parity_metrics.txt` | **Authoritative** parity (NXT translation, gene-corr filter; matches paper methodology) |
| `parity_grna.txt` / `parity_grna.tsv` | Per-feature-type breakdown (Jaccard, GEX, CRISPR) — gRNA run |
| `parity_larry.txt` / `parity_larry.tsv` | Per-feature-type breakdown — LARRY run (Custom features) |
| `star_Log.final.out` | STAR alignment summary |
| `star_solo_summary.csv` | STAR Solo GeneFull summary |
| `star_pf_multi_config.csv` | STAR pfMultiConfig (3-library) |
| `cr_grna_multi_config.csv` | CR `multi` config for GEX + gRNA |
| `cr_larry_multi_config.csv` | CR `multi` config for GEX + LARRY |
| `cr_grna_time.log` / `cr_larry_time.log` | CR `/usr/bin/time -v` logs |

## Data provenance

- STAR run: `/storage/MSK-perturb-comparison/paper_bench_ES_20260430_054814`
- CR runs:  `/storage/MSK-perturb-comparison/cr_paper_bench_ES_20260430_054814/{cr_gex_grna,cr_gex_larry}`
- Staging:  `/storage/MSK-perturb-comparison/msk30ko_full_3lib_ES_20260430_054814/fastqs`
- Pipeline logs: `.../msk30ko_full_3lib_ES_20260430_054814/_logs/{orchestrator,star,cr}.log`
