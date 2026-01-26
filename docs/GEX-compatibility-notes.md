# GEX Compatibility Notes (STAR vs Cell Ranger)

## Scope

These notes capture the GEX parity settings and the main reasons Pearson/Spearman can drop when comparing STARsolo outputs to Cell Ranger filtered MEX outputs.

## CR-Parity Parameter Set (Recommended)

Use these settings when the goal is to match Cell Ranger GEX counts as closely as possible:

- `--soloFeatures Gene`  
  CR GEX counts are exons-only. `GeneFull` adds introns and typically lowers Pearson vs CR.
- `--soloMultiMappers Unique`  
  CR does not use EM/Uniform multimapper redistribution for the standard GEX matrix.
- `--soloCellFilter EmptyDrops_CR`  
  Run EmptyDrops on GEX only; apply the filtered barcodes to features later.
- `--soloUMIdedup 1MM_CR` and `--soloUMIfiltering MultiGeneUMI_CR`  
  Keep CR-like UMI handling.
- `--soloCbUbRequireTogether no`  
  Used in recent CR-like parity runs; default is `yes` now, so set explicitly if you want parity with those runs.
- `--soloCrGexFeature gene`  
  Ensures CR-compat merge uses the Gene MEX (errors if Gene is missing).

Optional stricter multimapper handling:
- `--outFilterMultimapNmax 1` (if you want to drop all multimappers at alignment time).

### Example STARsolo invocation (key flags only)

```bash
STAR \
  --soloFeatures Gene \
  --soloMultiMappers Unique \
  --soloCellFilter EmptyDrops_CR \
  --soloUMIdedup 1MM_CR \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloCbUbRequireTogether no \
  --soloCrGexFeature gene
```

## Comparison Script Settings

The comparison helper `tests/compare_a375_star_mex.py` uses filtering thresholds that directly affect correlations:

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
   - CR uses exons-only (`Gene`). `GeneFull` adds introns and usually lowers Pearson.
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

## Notes on CR-Compat Merge (`--soloCrGexFeature`)

- `--soloCrGexFeature gene` forces the CR-compat merge to use the Gene MEX.
- If Gene output is missing, CR-compat merge errors; `auto` falls back to GeneFull.
- If you need both Gene and GeneFull outputs, run `--soloFeatures Gene GeneFull`.

## References

- `plans/a375_star_mex_comparison_summary.md` (Gene vs GeneFull correlation summary)
- `tests/a375_gex_cr_like_comparison_results.md` (latest full-depth CR-like run)
- `tests/compare_a375_star_mex.py` (comparison thresholds and filters)
