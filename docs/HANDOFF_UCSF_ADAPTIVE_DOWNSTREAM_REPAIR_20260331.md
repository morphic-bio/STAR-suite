# UCSF Adaptive Downstream Repair 2026-03-31

This note records the downstream QC repair applied to the corrected UCSF
production rerun rooted at:

- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009`

## Problem

The corrected UCSF STAR rerun recovered the expected GEX and CRISPR surfaces,
but the downstream `h5ad` packaging drifted back to the old fixed Seurat-like
QC thresholds:

- `min_genes=200`
- `max_genes=2500`
- `mt_pct_cutoff=5`

This happened in the downstream layer, not in STAR or CellBender. For
`EBs1_1`, the fixed threshold path reduced:

- `14350` STAR-called cells
- `12500` singlets after `scDblFinder`
- to only `285` `singlet_filtered` cells

The direct evidence is in:

- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/samples/downstream_prep_logs_downstream_genefull_velocyto_cellbender/EBs1_1.log`

## Scripts Changed

The downstream defaults were switched back to adaptive QC so future samples do
not silently reuse the obsolete `max_genes=2500` path:

- `scripts/run_scrna_downstream_gene_full_velocyto.sh`
- `scripts/run_remote_cellbender_scan.sh`
- `scripts/run_ucsf_downstream_batch.sh`

Specifically:

- `ADAPTIVE_FILTER` now defaults to `1`
- the adaptive path computes `max_genes = median(n_genes in STAR singlets) + 3 * MAD`
- the wrapper still keeps `min_genes=200` and `mt_pct_cutoff=5` unless passed explicitly

## Repair Procedure

For already completed downstream outputs, the intended repair flow is:

1. Re-run downstream prep with adaptive QC against the existing STAR `run/`:
   - `scripts/run_scrna_downstream_gene_full_velocyto.sh --run-dir ... --output-dir ... --adaptive-filter`
2. Reuse the existing CellBender result:
   - `downstream_genefull_velocyto_cellbender/cellbender/cellbender_counts.h5`
3. Reapply the saved CellBender layer to rebuilt `counts.h5ad` and `unfiltered_counts.h5ad`
4. Propagate that layer into rebuilt `filtered_counts.h5ad` and `default_singlet_filtered_counts.h5ad`
5. Rebuild `final_counts.h5ad` from the repaired unfiltered object
6. Regenerate `final_counts.summary.txt`

This preserves CellBender provenance and avoids re-running denoising when only
the QC mask was wrong.

## Pilot Validation

`EBs1_1` was used as the repair pilot. The adaptive rerun produced:

- `star_called_cells = 14350`
- `singlet_cells = 12580`
- `effective_max_genes = 7529`
- `QC-only filtered cells = 11593`
- `default singlet-filtered cells = 10652`

These values are recorded in:

- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/samples/EBs1_1/downstream_genefull_velocyto_cellbender/adaptive_qc_threshold.json`
- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/samples/EBs1_1/downstream_genefull_velocyto_cellbender/adaptive_repair_20260331_085330.log`

The pre-repair fixed-threshold summaries were preserved as timestamped text
backups in the same output directory:

- `summary.fixed_threshold_*.txt`
- `final_counts.fixed_threshold_*.summary.txt`

## Operational Note

The remote watcher for this run is:

- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/REMOTE_SCAN_RELAUNCH_FIX_20260331_072949.log`

After the default change in `run_scrna_downstream_gene_full_velocyto.sh`, newly
prepared downstream samples will inherit adaptive QC even if the already-running
watcher process itself was launched without `--adaptive-filter`.
