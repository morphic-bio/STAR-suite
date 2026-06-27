# Paper Benchmark Methodology

This document records the canonical methodology for all perturb-seq parity
benchmarks reported in the paper. It exists because a parameter inconsistency
between the A375 and UCSF/MSK paper scripts caused a 0.975 → 0.943 Gene
Pearson regression that took significant effort to diagnose.

## Root Cause of the A375 Regression

The A375 paper script (`scripts/paper/run_a375_benchmark.sh`) was missing
`--soloCrMultimapRescue yes`, which the UCSF and MSK scripts both had. This
single omission caused:

- 1.7M reads shifted from unique to multi-mapped (no rescue to recover them)
- 12% fewer UMIs on common barcodes (18.7M vs 21.4M)
- Gene Pearson dropped from 0.975 to 0.943

A secondary factor was a different genome index (autoindex `bulk_index` with
`sjdbOverhang 100` and `cellrangerLegacyGtfFilter No`) vs the original CellRanger 9
pre-built reference (`sjdbOverhang 90`, legacy GTF filter). The original CR9
reference at `/storage/A375-CR-9.01/` was deleted, and the `crstar` symlink index
at `/storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-crstar/star/` has
15 broken symlinks pointing to that deleted path.

### Timeline

| Date | Event |
|---|---|
| 2026-03-15 | Original A375 benchmark: STAR vs local CR9 run, Gene Pearson **0.975** |
| 2026-03-17 | Paper script created, inadvertently dropped `--soloCrMultimapRescue yes` |
| 2026-03-17 | Re-run produced Gene Pearson **0.943** (also used autoindex instead of CR9 ref) |
| 2026-03-17 | Investigation traced root cause to missing rescue flag + different index |

### Verification

Running the old STAR build (commit c30a106, same era as the 0.975 benchmark)
against the CR9 reference with the official parity script also yields ~0.943
Gene Pearson — confirming the script/index are the variables, not the STAR build.

## Canonical STAR Parameters (All Perturb-seq Benchmarks)

Every paper benchmark script MUST include the following CR-compat parameters.
Dataset-specific values are noted.

### Shared across all datasets

```
--clipAdapterType CellRanger4
--alignEndsType Local
--chimSegmentMin 1000000
--soloType CB_UMI_Simple
--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12
--soloBarcodeReadLength 0
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
--soloUMIfiltering MultiGeneUMI_CR
--soloUMIdedup 1MM_CR
--soloMultiMappers Unique
--soloCellFilter EmptyDrops_CR
--soloCbUbRequireTogether no
--soloFeatures GeneFull
--soloCrGexFeature genefull
--soloCrMultimapRescue yes          ← CRITICAL, do not omit
--crAssignMaxHamming 1
--crAssignFeatureOffset <dataset>   ← 0 for A375, -1 for UCSF
--crAssignLimitSearch -1
--crAssignMinCounts 0
--crAssignMaxBarcodeMismatches 5
--crAssignFeatureN 0
--crAssignBarcodeN 1
--crAssignConsumerThreads -1
--crAssignSearchThreads 1
--crAssignSkipQcOutputs 1
--dynamicThreadInterface 1
--dynamicThreadConstMapPermits <threads>
--dynamicThreadTelemetry 1
```

### Dataset-specific parameters

| Parameter | A375 | UCSF EBs2_2 | MSK 30polyKO |
|---|---|---|---|
| `--soloStrand` | `Unstranded` | `Forward` | `Forward` |
| `--clip3pPolyG` | omit (not NovaSeq) | `yes` | `yes` |
| `--crChemistry` | `TRU` | `auto` | `auto` |
| `--crOutputChemistry` | n/a (TRU native) | `TRU` | n/a |
| `--crWhitelist` | single-col TRU | 2-col NXT | per-library |
| `--crMinUmi` | `10` | `3` | `2` |
| `--crAssignFeatureOffset` | `0` | `-1` | `0` |

### Why `--soloCrMultimapRescue yes` matters

CellRanger 9 uses a multimap rescue algorithm that reassigns multi-mapped reads
to unique genes when one mapping is exonic and the others are intergenic or
intronic. Without this flag, STAR leaves those reads as multi-mapped and they are
excluded from Solo counting (which defaults to `--soloMultiMappers Unique`).

For A375, rescue recovered ~1.77M reads in the original benchmark:
- Gene-vs-non-gene fast path: 817K
- Exonic winner (Phase 3): 949K
- Intronic fallback (Phase 3): 3K

## Canonical Parity Script

All paper parity numbers MUST use the same script and parameters:

```bash
python3 scripts/report_additional_parity_metrics.py \
  --cr-run <CR_DIR> \
  --star-run <STAR_DIR> \
  --gene-corr-min-counts 20 \
  --gene-corr-min-cells-pct 0.01 \
  --translate <none|both>
```

- `--translate none` for TRU-native datasets (A375) where no barcode translation
  is needed.
- `--translate both` for NXT datasets (UCSF EBs2_2) that use 2-column whitelists.
- `--gene-corr-min-counts 20` and `--gene-corr-min-cells-pct 0.01` are the
  standard gene-level filtering thresholds.

### Gene filtering logic

The script uses a strict AND condition for gene filtering
(`summarize_gene_corr`, lines 398–403):

```
left_sum >= min_counts
AND right_sum >= min_counts
AND left_cell_n >= min_cells_abs
AND right_cell_n >= min_cells_abs
```

where `min_cells_abs = ceil(len(common_barcodes) * min_cells_pct)`.

Do NOT use ad-hoc scripts with OR-based filtering — those produce inflated
Gene Pearson values by including low-expression genes that are noise-dominated.

### CellRanger reference outputs

For parity, the CR reference must be a **CellRanger 9** run, not the
10x-published reference downloads (which may be CellRanger 8.0.0 or earlier).

For A375, the 10x-published `sample_filtered_feature_bc_matrix` was generated by
CellRanger 8.0.0 (different from CR9). Using it as the reference produces
different parity numbers (0.943 Gene Pearson vs 0.975 against actual CR9 output).
Always verify the `software_version` in the MTX header.

## Genome Index

The original A375 benchmark used the CellRanger 9 pre-built reference at
`/storage/A375-CR-9.01/refdata-gex-GRCh38-2024-A/star/` (now deleted).

Current benchmarks use `/storage/autoindex_110_44/bulk_index` (sjdbOverhang 100,
gencode v44, `cellrangerLegacyGtfFilter No`). This index produces a ~3%
lower unique mapping rate (70.3% vs 73.6%) because it includes more gene models,
shifting reads from unique to multi-mapped. With `--soloCrMultimapRescue yes`,
most of these are recovered.

**Action needed**: The CR9 reference needs to be re-downloaded or the autoindex
needs to be validated as producing equivalent parity when rescue is enabled.

## PE Bulk Benchmark

The PE bulk benchmark (`scripts/paper/run_pe_bulk_feature_benchmark.sh`) compares
integrated STAR-suite (trimming + alignment + Y-removal + internal TranscriptVB)
against an external stepwise pipeline (trimvalidate + STAR + remove_y_reads +
Salmon).

### Quantification contract

- The production STAR-suite arm is internal TranscriptVB, enabled with
  `--quantMode TranscriptVB --quantVBgcBias 1`.
- The integrated recipe must pass the same transcriptome FASTA used by Salmon via
  `--transcriptomeFasta`, so GC-bias and effective-length calculation use the
  pinned reference sidecar/input.
- External Salmon is part of the external stepwise production baseline, because
  that arm needs a separate quantifier after STAR emits `Aligned.toTranscriptome.out.bam`.
- Integrated TranscriptomeSAM emission and integrated Salmon QC are opt-in
  parity artifacts enabled with `--parity-qc`. They are not part of the
  STAR-suite production timing.
- Internal TranscriptVB collects transcript evidence during the STAR run, then
  merges/finalizes ECs, GC/effective-length state, and VB/EM convergence after
  alignment EOF. Do not describe this as external Salmon running concurrently.
- Headline speedups must be measured with this wrapper and the same output mode
  on both arms. The integrated total is the timed STAR-suite production command;
  the external total is decompress + trimvalidate + STAR TranscriptomeSAM +
  optional remove_y_reads + Salmon. Lean sanity checks using direct Trim Galore
  or unsorted BAM are valid diagnostics, but their ratios are not
  apples-to-apples with the paper benchmark.

### Modes

- `--yremove` (default): Both arms perform Y-chromosome removal. Integrated uses
  `--emitNoYBAM yes --emitYNoYFastq yes`; external uses `remove_y_reads` on
  trimmed FASTQs.
- `--no-yremove`: Y-chromosome removal disabled on both arms. Useful for
  datasets/analyses where Y-removal is not relevant.
- `--integrated-only`: Run only the STAR-suite production arm. Use this for
  STAR-suite timing refreshes when the external control does not need to be
  rerun. In default mode this is equivalent to `--skip-external --skip-compare`
  and does not enable parity artifacts.

### Parity comparison

When `--parity-qc` is enabled, transcript and gene-level Pearson/Spearman
correlations are computed by `tests/transcriptvb/compare_salmon_star.py`. The
benchmark produces a `comparison_metrics.tsv` with three comparisons per stage:
1. Integrated TranscriptVB vs Integrated Salmon (internal consistency)
2. Integrated Salmon vs External Salmon (pipeline effect on same tool)
3. Integrated TranscriptVB vs External Salmon (the headline comparison)

Do not run parity-QC mode for headline wall-time claims. It deliberately emits
and rereads an integrated transcriptome BAM that normal STAR-suite production
does not require.

### 2026-06-27 PPARG timing reruns

The corrected production-mode STAR-suite-only PPARG no-Y rerun on `/storage`
measured integrated trim+align+sorted BAM+internal TranscriptVB at `8:54.52`.
It did not emit integrated TranscriptomeSAM and did not run integrated Salmon
QC, matching the normal STAR-suite production arm. Run root:
`/storage/JAX_PE/results/pparg_prod_benchmark_no_y_20260627_172349/`.

An earlier PPARG no-Y sanity rerun before the production-mode wrapper correction
measured STAR-suite integrated trim+align+TranscriptVB plus TranscriptomeSAM
emission at `7:18.06`; Salmon QC on the integrated transcriptome BAM added
`1:00.08`.
A lean serial comparator using Trim Galore, upstream STAR, unsorted transcriptome
BAM, and Salmon took `9:50.17`, giving `1.35x` for that diagnostic setup. This
run confirmed high current Salmon parity (NumReads Pearson `0.999979` on all
transcripts; `0.999980` at sum>=1000) but should not replace the archived paper
speedup because its serial comparator is lighter than the paper wrapper and that
earlier integrated arm was still parity-artifact enabled.

## Paper Scripts

| Dataset | Script | Status |
|---|---|---|
| A375 | `scripts/paper/run_a375_benchmark.sh` | Fixed (rescue added) |
| UCSF EBs2_2 | `scripts/paper/run_ucsf_ebs2_2_benchmark.sh` | OK |
| MSK 30polyKO | `scripts/paper/run_msk_30polyko_benchmark.sh` | OK (but crstar index broken) |
| PE Bulk | `scripts/paper/run_pe_bulk_feature_benchmark.sh` | OK (supports --yremove/--no-yremove) |

## Broken Index: crstar

`/storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-crstar/star/` has
15 broken symlinks pointing to the deleted `/storage/A375-CR-9.01/` directory.
The MSK paper script defaults to this path and will fail until repaired.
