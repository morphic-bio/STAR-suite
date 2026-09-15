# GEX Compatibility Notes (STAR vs Cell Ranger)

## Scope

These notes capture the GEX parity settings and the main reasons Pearson/Spearman can drop when comparing STARsolo outputs to Cell Ranger filtered MEX outputs.

## Reference/Input Parity Scope

For CellRanger-style reference preparation (`--autoIndex Yes --cellrangerStyleIndex Yes`), parity guarantees apply to formatted input files:

- `${genomeDir}/cellranger_ref/genome.fa`
- `${genomeDir}/cellranger_ref/genes.gtf`

For `2024-A`, updated filtering mode (`--cellrangerLegacyGtfFilter No`, or `Auto` with `2024-A`) includes 47 `protein_coding_LoF` genes and excludes 33 chrY PAR genes versus legacy mode.

Index binaries (`Genome`, `SA`, `SAindex`) are not guaranteed byte-identical across STAR builds/versions, even with identical FASTA/GTF inputs.

## CR-Parity Parameter Set (Recommended)

Use these settings when the goal is to match Cell Ranger GEX counts as closely as possible:

- `--soloFeatures GeneFull`  
  Cell Ranger 7.0 and later count intronic reads by default (`include-introns`
  true, per 10x's documentation), so `GeneFull` is the matching surface. Use
  `Gene` (exons only) only against a Cell Ranger run made with introns excluded
  or an older Cell Ranger release.
- `--soloStrand` matching the library  
  `Forward` for 3' libraries; `Reverse` for 10x 5' R2-only libraries (Cell
  Ranger chemistry `SC5P-R2*`, e.g. A375), because read 2 is antisense to the
  transcript. `Unstranded` drops reads where opposite-strand genes overlap and
  counts antisense reads.
- `--soloCrMultimapRescue yes`  
  Cell Ranger-matched multimapper rescue (default `no`).
- `--clip3pPolyG yes` on NovaSeq/NextSeq data  
  Poly-G tails otherwise inflate specific genes (default `auto` trims only with
  `--clipAdapterType CellRanger4`).
- `--soloMultiMappers Unique`  
  CR does not use EM/Uniform multimapper redistribution for the standard GEX matrix.
- `--soloCellFilter EmptyDrops_CR`  
  Run EmptyDrops on GEX only; apply the filtered barcodes to features later.
  - **Backend note**: EmptyDrops_CR now uses libscrna by default.
  - Use `--soloEmptyDropsLegacy yes` to force the legacy STAR EmptyDrops_CR path.
- `--soloUMIdedup 1MM_CR` and `--soloUMIfiltering MultiGeneUMI_CR`  
  Keep CR-like UMI handling.
- `--soloCbUbRequireTogether no`  
  Used in recent CR-like parity runs; default is `yes` now, so set explicitly if you want parity with those runs.
- `--soloCrGexFeature genefull`  
  Ensures CR-compat merge uses the GeneFull MEX (errors if GeneFull is missing).

The full option set used for the STAR Suite 1.9.4 manuscript benchmarks is in
[PAPER_BENCHMARK_METHODOLOGY.md](PAPER_BENCHMARK_METHODOLOGY.md) Section 1.6.

Optional stricter multimapper handling:
- `--outFilterMultimapNmax 1` (if you want to drop all multimappers at alignment time).

### Example STARsolo invocation (key flags only)

```bash
STAR \
  --soloFeatures GeneFull \
  --soloStrand Forward \
  --soloMultiMappers Unique \
  --soloCrMultimapRescue yes \
  --soloCellFilter EmptyDrops_CR \
  --soloUMIdedup 1MM_CR \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloCbUbRequireTogether no \
  --soloCrGexFeature genefull
```

## Comparison Script Settings

The gene-level concordance reported for STAR Suite 1.9.4 is Spearman over every
gene in both annotations (zero-count genes included) together with Pearson on
raw per-gene totals, over the cells both tools called: `spearman_all_genes` and
`pearson_all_genes` from `scripts/report_additional_parity_metrics.py`. Its
`--gene-corr-min-counts` / `--gene-corr-min-cells-pct` thresholds only affect
the `*_filtered_genes` fields (and the `pearson`/`spearman` aliases), which are
not the reported metric. See
[PAPER_BENCHMARK_METHODOLOGY.md](PAPER_BENCHMARK_METHODOLOGY.md) Section 1.5.

The older comparison helper `tests/compare_a375_star_mex.py` uses filtering thresholds that directly affect correlations:

- Default `--min-counts 20` and `--min-cells-pct 0.01`
- Earlier runs that reported higher Pearson (e.g., 0.95+ range) used **≥10 counts** in the filter step.

If you want to reproduce those older correlations, pass:

```bash
python3 tests/compare_a375_star_mex.py \
  <CR_MEX_DIR> <STAR_MEX_DIR> \
  --feature-types "Gene Expression" \
  --min-counts 10 \
  --min-cells-pct 0.01
```

## Common Sources of Pearson Drops

1. **Gene vs GeneFull mismatch**  
   - Cell Ranger 7.0 and later count introns by default, which matches `GeneFull`. Comparing `Gene` (exons only) against such a run lowers the correlations; `Gene` matches only a Cell Ranger run with introns excluded.
2. **Multimapper mode**  
   - `Unique` aligns best with CR. `EM` or `Uniform` inflate totals and reduce Pearson.
3. **Filtered vs raw MEX mismatch**  
   - CR’s filtered MEX has empty drops removed. STAR’s raw MEX is unfiltered unless `--soloCellFilter` is used.
4. **Barcode harmonization**  
   - Comparisons strip `-1` suffixes and only use common barcodes.
5. **Filter thresholds in the comparison script**  
   - `--min-counts` and `--min-cells-pct` change which genes are included in correlations.
6. **Dataset/version mismatch**  
   - 2024-A CR outputs vs `/storage/A375/outputs/unpacked` runs can differ; full-depth vs downsample also changes results.
7. **Strand mismatch**  
   - A 5' R2-only library run with `--soloStrand Unstranded` instead of `Reverse` loses reads at overlapping opposite-strand genes and counts antisense reads (A375 gene Spearman 0.952 with `Unstranded`, 0.988 with `Reverse`, against Cell Ranger 9.0.1).

## Notes on CR-Compat Merge (`--soloCrGexFeature`)

- `--soloCrGexFeature gene` forces the CR-compat merge to use the Gene MEX.
- If Gene output is missing, CR-compat merge errors; `auto` falls back to GeneFull.
- If you need both Gene and GeneFull outputs, run `--soloFeatures Gene GeneFull`.

## References

- `plans/a375_star_mex_comparison_summary.md` (Gene vs GeneFull correlation summary)
- `tests/a375_gex_cr_like_comparison_results.md` (latest full-depth CR-like run)
- `tests/compare_a375_star_mex.py` (comparison thresholds and filters)
