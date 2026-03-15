# A375 Perturb-seq Full Benchmark: STAR vs Cell Ranger 9

**Date**: 2026-03-15  
**Dataset**: 10x 1k CRISPR 5' GEM-X (A375), full FASTQs  
**Host**: pikachu (32 threads, 128 GB RAM)

## Summary

| Metric | Value |
|---|---|
| STAR cells | 1,191 |
| Cell Ranger cells | 1,162 |
| Barcode Jaccard | 0.98 |
| Gene Pearson (filtered) | 0.975 |
| Cell Pearson (per-barcode UMI) | 1.000 |
| CRISPR exact match | 100% (1,083/1,083) |
| CRISPR UMI Pearson | 1.000 |
| STAR wall time | 6 min 48 s |
| Cell Ranger wall time | 14 min 36 s |
| Speedup | 2.1× |

## Configuration

### Shared

- Reference: `/storage/A375-CR-9.01/refdata-gex-GRCh38-2024-A` (GRCh38-2024-A, STAR 2.7.1a index)
- GEX FASTQs: `/storage/A375-CR-9.01/full_fastqs/gex/` (47,095,182 reads, 2 lanes)
- CRISPR FASTQs: `/storage/A375-CR-9.01/full_fastqs/crispr/` (9,748,584 reads, 2 lanes)
- Feature reference: `/storage/A375/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv` (11 guides)
- Whitelist: `/storage/A375/3M-5pgex-jan-2023.txt` (3.7M 5' GEX barcodes, TRU namespace)
- CRISPR min-UMI: 10
- BAM output: disabled
- Threads: 32

### Cell Ranger 9.0.1

- Config: `cr_multi_config.csv`
- `create-bam,false`
- `localcores,32`
- `localmem,128`

### STAR 2.7.11b (STAR-suite)

- Config: `pf_multi_config.csv`
- Key flags:
  - `--soloStrand Unstranded` (required for this dataset)
  - `--soloFeatures GeneFull`
  - `--soloCrGexFeature genefull`
  - `--soloCrMultimapRescue yes`
  - `--clipAdapterType CellRanger4 --clip3pPolyG yes`
  - `--soloCellFilter EmptyDrops_CR`
  - `--crChemistry TRU`
  - `--crMinUmi 10`
  - `--dynamicThreadInterface 1` (permits: alignment + feature calling concurrent)
  - `--dynamicThreadTelemetry 1`
  - `--outSAMtype None`

## GEX Parity Detail

From `report_additional_parity_metrics.py` (no barcode translation needed for TRU):

- **Per-barcode UMI Pearson** (filtered vs filtered, 1,162 common barcodes): **0.9995**
- **Per-gene total Pearson** (15,677 filtered genes, min 20 counts, min 1% cells): **0.975**
- UMI totals on common barcodes: STAR 21,376,229 vs CR 21,361,149 (ratio 1.0007)

## CRISPR Parity Detail

- CR evaluated rows: 1,083; STAR evaluated rows: 1,191 (1,086 non-None)
- Common evaluated rows: 1,083
- Set-equivalent calls: 1,083 / 1,083 = **100%**
- UMI sum Pearson on common: **0.9999**
- Total assigned UMIs: STAR 3,049,220 vs CR 3,001,352 (delta +47,868)

## STAR Alignment Statistics

- Input reads: 47,095,182 (GEX only; CRISPR routed through pf-multi)
- Uniquely mapped: 71.80%
- Multi-mapped: 10.71%
- Unmapped (too short): 16.57%
- Reads in cells mapped to GeneFull: 33,399,180
- Median UMI per cell: 18,213
- Median genes per cell: 5,360
- Estimated cells: 1,191 (EmptyDrops_CR: 1,137 simple + 54 rescued)

## Timing Breakdown (STAR)

| Phase | Time |
|---|---|
| Genome loading | ~1 min |
| Mapping | 2 min 28 s |
| Solo counting + EmptyDrops | 1 min 21 s |
| Feature assignment (CRISPR) | ~40 s |
| GMM calling | < 1 s |
| **Total wall time** | **6 min 48 s** |

## Files in This Directory

| File | Description |
|---|---|
| `parity_metrics.txt` | Full output from `report_additional_parity_metrics.py` |
| `cr_multi_config.csv` | Cell Ranger multi config used |
| `pf_multi_config.csv` | STAR pfMultiConfig used |
| `cr_run.log` | Cell Ranger run log |
| `star_run.log` | STAR run log (stdout/stderr) |
| `star_Log.final.out` | STAR alignment summary |
| `star_solo_summary.csv` | STAR Solo GeneFull summary |

## Reproduction

```bash
# Cell Ranger
cellranger multi --id cr_a375_32t_nobam --csv cr_multi_config.csv

# STAR (after building STAR-suite)
STAR --runThreadN 32 \
  --genomeDir /storage/A375-CR-9.01/refdata-gex-GRCh38-2024-A/star \
  --readFilesIn <GEX_R2>,<CRISPR_R2> <GEX_R1>,<CRISPR_R1> \
  --readFilesCommand zcat \
  --outSAMtype None \
  --clipAdapterType CellRanger4 --clip3pPolyG yes \
  --alignEndsType Local --chimSegmentMin 1000000 \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloCBwhitelist /storage/A375/3M-5pgex-jan-2023.txt \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
  --soloMultiMappers Unique --soloCellFilter EmptyDrops_CR \
  --soloCbUbRequireTogether no --soloStrand Unstranded \
  --soloFeatures GeneFull --soloCrGexFeature genefull \
  --soloCrMultimapRescue yes \
  --pfMultiConfig pf_multi_config.csv \
  --crWhitelist /storage/A375/3M-5pgex-jan-2023.txt \
  --crChemistry TRU --crMinUmi 10 \
  --crFeatureRef /storage/A375/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv \
  --dynamicThreadInterface 1

# Parity check
python3 scripts/report_additional_parity_metrics.py \
  --cr-run <CR_RUN> --star-run <STAR_RUN> \
  --gene-corr-min-counts 20 --gene-corr-min-cells-pct 0.01
```

## Notes

- STAR routes CRISPR FASTQ reads through the pf-multi assignment pipeline,
  not through genome alignment. STAR's "Number of input reads" (47M) reflects
  GEX reads only; CR reports 47M GEX + 9.75M CRISPR = 56.8M total.
- `--soloStrand Unstranded` is required for this A375 5' GEX dataset.
  Using `Forward` drops GeneFull mapping to ~11% (strandedness mismatch).
- The `--crChemistry TRU` flag resolves the whitelist namespace for
  feature barcode assignment. Prior to commit db68dd2, whitelists without
  "NXT" or "TRU" in the filename required a symlink workaround.
