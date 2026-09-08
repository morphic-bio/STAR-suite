# Streaming per-cell MEX statistics

`flex_mex_cell_stats` computes per-cell raw-count and detected-gene statistics
for one or more explicit cell cohorts while streaming a single MEX matrix. It
is intended for repeated tranche analysis where loading the whole raw matrix
would waste time and memory.

Cell identity is exact and tag-aware. The reader removes only a terminal `-1`
from barcode IDs; it never truncates a `CB16+TAG8` identity to `CB16`. This is
important when the same CB16 occurs under more than one tag.

Build and test:

```bash
make -C docs/benchmarks/jax_matrix_20260904/analysis flex_mex_cell_stats
make -C docs/benchmarks/jax_matrix_20260904/analysis test-flex-mex-cell-stats
```

Run on either a combined STAR raw MEX or a Cell Ranger per-sample raw MEX:

```bash
analysis/flex_mex_cell_stats \
  --mex-dir RUN/raw_feature_bc_matrix \
  --cohort shared=cohorts/shared.txt \
  --cohort star_only=cohorts/star_only.txt \
  --cohort cr_only=cohorts/cr_only.txt \
  --input-label star \
  --out-prefix results/star.raw
```

The MEX directory may contain plain or gzip-compressed `matrix.mtx`,
`barcodes.tsv`, and `features.tsv` files. Legacy `genes.tsv` is also accepted.
The program validates the matrix dimensions and declared number of entries,
but retains only one small statistics record per requested cell plus a column
lookup vector; it never retains sparse matrix entries.

Outputs:

- `.cells.tsv` has one row per cohort cell, in cohort-file order, with exact
  found/missing status, total counts, detected genes, and MatrixMarket `nnz`.
- `.summary.tsv` has cohort totals and deterministic nearest-rank quantiles
  (`min`, `q25`, median, `q75`, `p90`, `p95`, `p99`, `max`) over cells found on
  the MEX barcode axis. Missing cells are reported separately and excluded
  from quantiles.

Canonical 10x MEX files contain one positive coordinate per detected
gene/cell, so `detected_genes` and `nnz` agree. If a noncanonical matrix stores
explicit zero values, `nnz` includes those coordinates while
`detected_genes` does not.
