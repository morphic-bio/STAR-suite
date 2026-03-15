# MSK 30polyKO Full Benchmark: STAR vs Cell Ranger 9

**Date**: 2026-03-06  
**Dataset**: MSK 30-guide polyKO Perturb-seq (NXT chemistry, 3 libraries: GEX + gRNA + LARRY lineage)  
**Host**: pikachu (32 threads, 110 GB RAM)

## Summary

| Metric | Value |
|---|---|
| STAR cells | 30,497 |
| Cell Ranger cells | 32,256 |
| Barcode Jaccard | 0.94 |
| Gene Pearson (filtered, 17,448 genes) | 0.993 |
| Cell Pearson (per-barcode UMI) | 1.000 |
| CRISPR set-equivalent calls | 98.5% (22,200/22,531) |
| CRISPR UMI Pearson | 0.999 |
| STAR wall time | 42 min (single run, all 3 libraries) |
| CR wall time (gRNA + GEX) | 58 min |
| CR wall time (LARRY + GEX) | 110 min |
| CR total wall time | 168 min |
| Speedup | 4.0× |

## Configuration

### Shared

- Reference: `/storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-crstar/star`
- GEX reads: 668,705,043 (28 FASTQ files across 4 samples, 8 lanes)
- Feature reference: 30 CRISPR guides (8 bp)
- Chemistry: NXT (barcode translation via `3M-february-2018_NXT.txt.gz`)
- Threads: 32

### Cell Ranger 9.0.1

- Config: `cr_multi_config.csv`
- `localcores,32`, `localmem,110`

### STAR 2.7.11b (STAR-suite)

- Config: `star_pf_multi_config.csv`
- Key flags:
  - `--soloFeatures Gene GeneFull`
  - `--soloCrGexFeature GeneFull`
  - `--soloCrMultimapRescue yes`
  - `--soloMultiMappers Rescue`
  - `--crChemistry auto` (NXT inferred from whitelist filename)
  - `--crMinUmi 2`
  - `--crAssignMaxHamming 5`
  - `--dynamicThreadInterface 1` (concurrent alignment + feature calling)
  - `--outSAMtype BAM Unsorted`

## GEX Parity Detail

From `report_additional_parity_metrics.py` (with NXT↔TRU barcode translation):

- **Per-barcode UMI Pearson** (filtered vs filtered, 30,417 common barcodes): **0.9999**
- **Per-gene total Pearson** (17,448 filtered genes, min 20 counts, min 1% cells): **0.993**
- UMI totals on common barcodes: STAR 266,445,024 vs CR 264,048,892 (ratio 1.009)

## CRISPR Parity Detail

- CR evaluated rows: 22,715; STAR evaluated rows: 30,497 (23,541 non-None)
- Common evaluated rows (excl STAR None): 22,531
- Set-equivalent calls: 22,200 / 22,531 = **98.5%**
- Singleton exact match rate: **99.99%**
- UMI sum Pearson on common: **0.999**
- Per-feature UMI Pearson (filtered): **0.9999**
- Mean cell-wise call Jaccard: **0.997**

## STAR Alignment Statistics

- Input reads: 668,705,043
- Uniquely mapped: 92.11%
- Multi-mapped: 3.84%
- Unmapped (too short): 2.95%
- Median UMI per cell: 8,315
- Median genes per cell: 3,911
- Estimated cells: 30,497 (EmptyDrops_CR: 29,287 simple + 1,210 rescued)

## Files in This Directory

| File | Description |
|---|---|
| `parity_metrics.txt` | Full output from `report_additional_parity_metrics.py` |
| `crispr_call_concordance.json` | CRISPR call-level concordance summary |
| `crispr_guide_parity.json` | Per-guide UMI count parity |
| `cr_multi_config.csv` | Cell Ranger multi config |
| `star_pf_multi_config.csv` | STAR pfMultiConfig |
| `cr_time.log` | Cell Ranger `/usr/bin/time` output |
| `star_Log.final.out` | STAR alignment summary |
| `star_solo_summary.csv` | STAR Solo GeneFull summary |

## Data Provenance

- STAR run: `/storage/MSK-perturb-comparison/canonical_tru_seq_20260306_052040/star_3lib`
- CR run: `/storage/MSK-perturb-comparison/cr_full_grna_withcalls_20260306_095454`
- Paper artifacts: `/storage/MSK-perturb-comparison/paper_artifacts/msk_grna_parity_30guide_20260306`
