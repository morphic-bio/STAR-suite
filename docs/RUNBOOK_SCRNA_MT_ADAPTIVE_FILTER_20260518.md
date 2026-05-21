# Runbook: Adaptive Mitochondrial-% QC Guard for scRNA Downstream

Date: 2026-05-18

## Purpose

Replace the fixed `mt_pct <= 5%` cell-QC cutoff in the scRNA downstream
pipeline with a **combined strict-floor + MAD soft guard** that adapts to
samples whose baseline mitochondrial-% is biologically elevated.

This runbook covers two sets of STAR-suite helper scripts:

- **New implementation** — scripts the downstream pipeline calls going forward.
- **Conversion** — scripts that retrofit already-built h5ad releases.

## Background

The downstream filter is `n_genes ∈ [min_genes, effective_max_genes]` **and**
`mt_pct <= 5%`. The n_genes bounds are already per-sample adaptive
(`median + n_mad·MAD`, see `compute_adaptive_qc_threshold.py`); the mt% cutoff
was a single fixed 5%.

For the MSK 30KO ES sample the **median mt% of called singlets is ~5.07%** —
the median cell already sits at the cutoff, so the fixed 5% cap bisects a
healthy population. ES is undifferentiated/proliferating and naturally carries
a higher baseline mt%. Measured filter cause for ES: of ~18.7k cells dropped,
~93% fail `mt% > 5` while only ~10% fail the n_genes cap. The fixed mt% cutoff,
not the MAD-based gene cap, is the dominant — and overly aggressive — cut.

## Method

A cell passes the mt% guard when

```
mt_pct <= mt_threshold
mt_threshold = max(mt_floor, median(mt_pct) + n_mad · MAD(mt_pct))
```

- `mt_floor` = 5.0% — strict guard; any cell at/below it is always kept, so the
  rule never tightens below the historical cutoff.
- `median + n_mad·MAD` — per-sample soft guard; `n_mad` = 3, **raw** MAD
  (no 1.4826 scaling), matching the n_genes convention. Raw MAD × 3 is already
  on the strict side of the field convention (3 *scaled* MADs in scuttle/OSCA).
- Threshold and median/MAD are computed over **STAR-called singlets**
  (`obs['singlet']`), then applied to every barcode — same population basis as
  the n_genes adaptive filter.
- The strict floor is folded into the single `mt_threshold` value, so the
  combined guard is just `mt_pct <= mt_threshold`.

Worked examples (MSK 30KO, singlets):

| Sample | median mt% | raw MAD | median+3·MAD | mt_threshold | effect |
| --- | ---: | ---: | ---: | ---: | --- |
| ES   | 5.06 | 0.88 | 7.69 | **7.69** | MAD guard active; +~12.4k cells |
| PP1  | 1.47 | 0.52 | 3.02 | **5.00** | floored; clean sample unchanged |
| S6_1 | 4.42 | 1.86 | 10.01 | **10.01** | MAD guard active; sample flagged |

This is a **one-parameter model**: a single global `n_mad` shared across all
samples; `median`/`MAD` are robust summaries, not fitted quantities. No
per-sample tuning, no elbow/curve fitting — by design, to avoid overfitting.

## Scripts (in `scripts/`)

| File | Set | Role |
| --- | --- | --- |
| `scrna_mt_adaptive.py` | shared | Core: `compute_mt_threshold`, `mt_sample_flag`, `qc_filter_mask`, `apply_mt_adaptive_filter`, `merge_threshold_json`. Imported by both sets so they compute identically. |
| `apply_adaptive_mt_filter.py` | new implementation | Computes the threshold and recomputes `filter` / `singlet_filtered` on the unfiltered downstream h5ad; merges mt keys into the threshold JSON. |
| `generate_qc_histogram_mt_adaptive.py` | new impl + conversion | QC histogram with two mt rejection lines (floor + adaptive). Replaces `generate_qc_histogram.py`. |
| `convert_h5ad_mt_adaptive.py` | conversion | Retrofits a built sample dir: recomputes obs, rewrites `final_counts.h5ad`, rebuilds `filtered_counts.h5ad` + `default_singlet_filtered_counts.h5ad`, updates the JSON. |

## New keys in `adaptive_qc_threshold.json`

`mt_pct_median`, `mt_pct_mad`, `mt_pct_n_mad`, `mt_pct_floor`,
`mt_pct_raw_threshold`, `mt_pct_threshold`, `mt_pct_threshold_was_floored`,
`mt_pct_cells`, `mt_pct_source`, `min_genes`, `max_genes`,
`filter_cells_strict_mt5`, `filter_cells_mt_adaptive`,
`singlet_filtered_cells_strict_mt5`, `singlet_filtered_cells_mt_adaptive`,
and the sample-flag keys
`mt_pct_flag`, `mt_pct_flag_high_pct`, `mt_pct_flag_high_fraction`,
`mt_pct_flag_high_fraction_limit`. The existing n_genes keys are preserved.

## New / changed obs columns

- `filter` — recomputed with the adaptive mt guard.
- `singlet_filtered` — recomputed (`singlet & filter`).
- `filter_strict_mt5`, `singlet_filtered_strict_mt5` — the **legacy** strict-5%
  columns, preserved once for audit and rollback.
- `adata.uns['mt_adaptive_filter']` — the full threshold record.

## New-implementation procedure

The pipeline orchestrator is
`scripts/run_scrna_downstream_gene_full_velocyto.sh`. The orchestrator now runs
the adaptive mt guard immediately after `combineFilters.py`. `combineFilters.py`
itself is **not** modified — it still runs with the legacy 5% cutoff and its
`filter` column is overridden afterward.

1. Keep `compute_adaptive_qc_threshold.py` and `combineFilters.py` as they are.
2. **Insert** `apply_adaptive_mt_filter.py` after `combineFilters.py` and before
   `postprocess_downstream_filters.py`:

   ```bash
   python3 scripts/apply_adaptive_mt_filter.py \
     --input-h5ad     "${downstream_dir}/unfiltered_counts.h5ad" \
     --threshold-json "${downstream_dir}/adaptive_qc_threshold.json" \
     --mt-floor 5.0 --n-mad 3.0
   ```

3. `postprocess_downstream_filters.py` runs unchanged — it builds
   `filtered_counts.h5ad` / `default_singlet_filtered_counts.h5ad` from the
   recomputed `filter`.
4. **Swap** the QC-histogram call from `generate_qc_histogram.py` to:

   ```bash
   python3 scripts/generate_qc_histogram_mt_adaptive.py \
     --input-h5ad     "${downstream_dir}/unfiltered_counts.h5ad" \
     --output-dir     "${downstream_dir}" \
     --threshold-json "${downstream_dir}/adaptive_qc_threshold.json"
   ```

`apply_adaptive_mt_filter.py` requires `mt_pct`, `n_genes`, and `singlet` obs
columns (all produced by `combineFilters.py`) and reads the n_genes bounds
(`min_genes`, `effective_max_genes`) from the threshold JSON. The unchanged
`postprocess_downstream_filters.py` step also expects `non_empty` when it builds
the exported filtered views.

For GPU-side post-MEX execution, `scripts/run_remote_multiome_post_mex_rsync.sh`
and `scripts/run_remote_scrna_downstream_rsync.sh` stage the same adaptive-mt
helpers to the remote checkout before invoking the downstream wrapper.

## Conversion procedure (retrofit existing releases)

Conversion needs each sample's **`final_counts.h5ad`** (all barcodes) — the
delivery-tree `filtered_counts.h5ad` holds only survivors and cannot be
un-filtered. `final_counts.h5ad` lives in the production-sample archive, e.g.
`/mnt/pikachu/msk30ko-production-sample-archive/MSK-05-13-26-large-files/<sample>/downstream_genefull_velocyto_cellbender/`.

1. **Dry-run** to preview thresholds and cell deltas (writes nothing):

   ```bash
   python3 scripts/convert_h5ad_mt_adaptive.py --dry-run \
     --sample-dir /mnt/pikachu/msk30ko-production-sample-archive/MSK-05-13-26-large-files/{ES,DE,PP1,PP2,S5_1,S5_2,S6_1,S6_2}/downstream_genefull_velocyto_cellbender
   ```

2. **Convert.** In place (keeps `*.pre_mt_adaptive` backups) or to a new tree
   with `--output-dir`:

   ```bash
   python3 scripts/convert_h5ad_mt_adaptive.py \
     --sample-dir <dir1> <dir2> ... \
     --output-dir /mnt/pikachu/msk30ko-mt-adaptive-rebuild
   ```

   Per sample this rewrites `final_counts.h5ad`, rebuilds
   `filtered_counts.h5ad` + `default_singlet_filtered_counts.h5ad`, and updates
   `adaptive_qc_threshold.json`.

3. **Regenerate QC graphs** for each converted sample:

   ```bash
   python3 scripts/generate_qc_histogram_mt_adaptive.py \
     --input-h5ad     <out>/final_counts.h5ad \
     --output-dir     <out> \
     --threshold-json <out>/adaptive_qc_threshold.json
   ```

4. **Re-stage** the rebuilt `filtered_counts.h5ad`,
   `default_singlet_filtered_counts.h5ad`, `adaptive_qc_threshold.json`, and
   `gene_quantile_histogram.{png,html}` into the delivery tree
   (`/mnt/pikachu/msk30ko-h5ad-qc-delivery/MSK-05-13-26-large-files/<sample>/...`),
   then re-upload.

`DE_GemX` is not in the production-sample archive alongside the other eight —
it was processed in an earlier separate trial run and never re-run (the trial
output was good). Convert it from that run directory:

```
/storage/MSK-perturb-comparison/msk30ko_DE_GemX_velocyto_downstream_trial_20260513_190549/samples/DE_GemX/downstream_genefull_velocyto_cellbender/
```

That directory holds `final_counts.h5ad` and `adaptive_qc_threshold.json`, so
`convert_h5ad_mt_adaptive.py --sample-dir <that dir>` works directly. The
delivery-tree `DE_GemX` h5ads were copied from this run (byte-size identical).

## QC graph

`generate_qc_histogram_mt_adaptive.py` replaces the single 5% MT line with
**two** horizontal lines on the mitochondrial-% axis:

- **MT floor** (purple dashed) — the strict 5% guard.
- **MT adaptive** (dark-red dotted) — `median + n_mad·MAD`.

The effective cut is the higher of the two and is labelled `(effective)` plus
annotated in the on-plot caption.

`gene_quantile_histogram.{png,html}` is always emitted. In addition, when the
sample-level flag fires (`mt_pct_flag = true`), the script emits a second,
**conditional** artifact `mt_quantile_histogram.{png,html}` — a standalone mt%
distribution over singlets with the floor / adaptive lines and the high-mt
review band drawn on it. This file is **not** produced for unflagged samples,
so flagged samples have one extra file to stage/upload to the delivery tree.

## Sample-level flag

`mt_sample_flag` raises `mt_pct_flag = true` when more than 10% of singlets
exceed 20% mt%. This is a **human-review tripwire, not a per-cell cut** — it
catches degraded samples where the MAD soft guard could drift high (e.g. S6_1,
~18% of singlets above 20% mt → flagged). Flagged samples should be reviewed
against their `mt_quantile_histogram.{png,html}` (emitted only when flagged)
before release; the filter still runs.

## Validation

- Re-run the provider cell-call concordance after conversion; ES recovery
  should move up toward the provider call. Some overshoot is expected and
  acceptable (provider counts are exon-only; lower-complexity high-mt cells are
  intentionally still filtered by the n_genes guard).
- Check `mt_pct_threshold_was_floored`: `true` for clean samples (PP1/PP2),
  `false` for ES/S5_*/S6_*.
- Confirm `filter_cells_mt_adaptive >= filter_cells_strict_mt5` for every
  sample (the guard never removes cells the 5% cutoff kept).
- Review any sample with `mt_pct_flag = true`.

## Rollback

The legacy strict-5% masks survive every conversion as `filter_strict_mt5` /
`singlet_filtered_strict_mt5`, and the originals are kept as
`*.pre_mt_adaptive` unless `--no-backup` was passed. To revert, restore from
the backups or rebuild the views from `filter_strict_mt5`.
