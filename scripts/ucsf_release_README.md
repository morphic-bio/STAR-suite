# UCSF Perturb-seq Release

**Date**: 2026-03-31
**Pipeline**: STAR-suite (GeneFull + Velocyto, Y-removal)
**Denoising**: CellBender 0.3.2 (GPU)
**Samples**: 16 (10 EBs, 6 iPSC)
**Feature type**: CRISPR Guide Capture (hCRISPRa-v2, AALG2 pattern)

## Quick start

For most analyses, use `filtered_counts.h5ad` per sample. It contains
singlet-filtered, QC-filtered cells with spliced/unspliced/ambiguous RNA
velocity layers and a CellBender-denoised count layer.

```python
import anndata
adata = anndata.read_h5ad("EBs1_1/filtered_counts.h5ad")
```

If you need to apply your own cell-calling or QC thresholds, use
`unfiltered_counts.h5ad` which contains all barcodes with boolean masks in
`.obs` to recover any subset.

## Directory layout

```
README.md                                    This file
samples.tsv                                  Per-sample manifest (see below)

<sample_id>/
    filtered_counts.h5ad                     Singlet + adaptive-QC filtered cells
    unfiltered_counts.h5ad                   All barcodes (cells + empty droplets)

    qc/
        summary.txt                          Detailed AnnData summary report
        adaptive_qc_threshold.json           Computed adaptive QC thresholds
        gene_quantile_histogram.png          Gene-count distribution QC plot
        gene_quantile_histogram.html         Interactive version (Plotly)

    features/
        crispr/
            protospacer_calls_per_cell.csv   Per-cell CRISPR guide assignments
            protospacer_calls_summary.csv    Guide-level assignment summary
            protospacer_umi_thresholds.csv   UMI calling thresholds (CSV)
            protospacer_umi_thresholds.json  UMI calling thresholds (JSON)
            filtered_feature_library.h5ad    CRISPR-only count matrix (filtered cells)
            raw_feature_library.h5ad         CRISPR-only count matrix (all barcodes)

    cellbender/
        cellbender_counts.h5                 CellBender denoised count matrix
        cellbender_counts.pdf                CellBender training/elbow QC plot
        cellbender_counts_report.html        CellBender detailed HTML report
        cellbender_counts_metrics.csv        CellBender summary metrics

    cellranger_compat/
        filtered_feature_bc_matrix/          CellRanger-style filtered GEX + CRISPR
        raw_feature_bc_matrix/               CellRanger-style raw GEX + CRISPR
        filtered_velocyto_feature_bc_matrix/ CellRanger-style filtered spliced/unspliced
        raw_velocyto_feature_bc_matrix/      CellRanger-style raw spliced/unspliced
```

## samples.tsv columns

| Column | Description |
|--------|-------------|
| `sample_id` | Sample identifier |
| `filtered_cells` | Cells in `filtered_counts.h5ad` (singlet + adaptive QC) |
| `singlet_cells` | Singlets after scDblFinder doublet removal |
| `star_called_cells` | Cells called by STAR (EmptyDrops_CR algorithm) |
| `effective_max_genes` | Adaptive max-genes threshold used for QC filtering |
| `median_n_genes` | Median genes per singlet (basis for adaptive threshold) |
| `mad_n_genes` | MAD of genes per singlet |
| `min_genes` | Minimum genes threshold (fixed at 200) |
| `n_mad` | Number of MADs above median for max-genes (fixed at 3.0) |

## h5ad file description

Both `filtered_counts.h5ad` and `unfiltered_counts.h5ad` share the same
schema. `filtered_counts.h5ad` is the subset of `unfiltered_counts.h5ad`
where `obs["singlet_filtered"] == True`.

### X (main count matrix)

CSR sparse matrix, float32. Total UMI counts per cell from STAR GeneFull
quantification (exonic + intronic reads mapped to gene bodies).

### Layers

| Layer | Type | Description |
|-------|------|-------------|
| `spliced` | int64 sparse | Exonic reads (mature mRNA) |
| `unspliced` | int64 sparse | Intronic reads (pre-mRNA) |
| `ambiguous` | int64 sparse | Reads that cannot be unambiguously assigned to spliced or unspliced |
| `denoised` | float64 sparse | CellBender ambient-RNA-removed counts (from `cellbender_counts.h5`) |

Note: `X = spliced + unspliced + ambiguous` for each cell.

### obs (cell-level annotations)

**Cell-calling and QC:**

| Column | Type | Description |
|--------|------|-------------|
| `is_cell` | bool | STAR cell call (EmptyDrops_CR algorithm) |
| `non_empty` | int64 | Non-empty barcode flag from STAR |
| `n_genes` | int64 | Number of genes detected (non-zero UMI count) |
| `total_counts` | float32 | Total UMI count |
| `mt_counts` | float32 | Mitochondrial gene UMI count |
| `mt_pct` | float32 | Mitochondrial fraction (%) |
| `filter` | bool | Passes QC thresholds (`min_genes <= n_genes <= effective_max_genes` and `mt_pct <= 5`) |

**Doublet detection (scDblFinder):**

| Column | Type | Description |
|--------|------|-------------|
| `doublet` | int64 | 1 = predicted doublet, 0 = singlet |
| `doublet_scores` | float64 | scDblFinder doublet probability score |
| `singlet` | bool | `is_cell and not doublet` |

**Combined filter:**

| Column | Type | Description |
|--------|------|-------------|
| `singlet_filtered` | bool | `singlet and filter` (the default usable-cell mask) |

**CRISPR guide assignments (library-specific, long-prefix columns):**

These columns use the prefix
`CRISPR_Guide_Capture_<sample>_AALG2_CRISPR_Guide_Capture_1__`:

| Suffix | Type | Description |
|--------|------|-------------|
| `num_features` | int64 | Number of distinct guides detected in this cell |
| `num_umis` | int64 | Total CRISPR guide UMIs |
| `feature_call` | category | Assigned guide(s), pipe-delimited if multiple |
| `is_featured` | bool | At least one guide assigned |
| `feature1_count` | int64 | UMI count of the top guide |
| `feature2_count` | int64 | UMI count of the second guide |
| `feature_call_category` | category | `single`, `multi`, `ambiguous`, or `none` |

**CRISPR guide assignments (generic aliases):**

Shortened aliases for the library-specific columns above, provided for
convenience when there is a single CRISPR feature library:

| Column | Type | Description |
|--------|------|-------------|
| `is_featured` | bool | At least one guide assigned |
| `feature_call` | category | Assigned guide(s) |
| `feature_call_num_features` | int64 | Same as library-specific `num_features` |
| `feature_call_num_umis` | int64 | Same as library-specific `num_umis` |
| `best_feature` | category | Top guide by UMI count |
| `feature1_count` | int64 | UMI count of the top guide |
| `feature2_count` | int64 | UMI count of the second guide |
| `feature_call_category` | category | `single`, `multi`, `ambiguous`, or `none` |

### var (gene-level annotations)

| Column | Type | Description |
|--------|------|-------------|
| `gene_symbols` | category | HGNC gene symbol |
| `feature_types` | category | Always `Gene Expression` |

Index (`var_names`) is the Ensembl gene ID.

### uns (unstructured metadata)

| Key | Description |
|-----|-------------|
| `feature_libraries` | List of integrated feature library IDs |
| `feature_library_generic_alias` | Mapping from generic alias columns to library-specific columns |
| `gene_expression_feature_kind` | Feature type identifier (`Gene Expression`) |
| `gene_expression_source` | Source matrix path |
| `velocyto_source` | Source Velocyto matrix path |

## Adaptive QC methodology

Each sample uses an adaptive maximum-genes threshold instead of a fixed
cutoff:

```
effective_max_genes = median(n_genes among singlets) + 3 * MAD(n_genes among singlets)
```

This prevents the high gene-complexity of these perturb-seq samples from
being erroneously filtered. The fixed `max_genes=2500` threshold used in
standard Seurat-like QC would discard the majority of real cells in this
dataset.

Thresholds per sample are recorded in `qc/adaptive_qc_threshold.json`.

## CellBender denoising

Ambient RNA removal was performed with CellBender 0.3.2 (GPU mode) using
the `unfiltered_counts.h5ad` as input. The denoised count matrix is
available both as:

- The `denoised` layer in the h5ad files (integrated by barcode matching)
- The standalone `cellbender/cellbender_counts.h5` file

## CellRanger-compatible outputs

The `cellranger_compat/` directory contains standard 10x-format sparse
matrices (MEX format: `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`)
compatible with CellRanger-based downstream tools (Seurat `Read10X`,
Scanpy `read_10x_mtx`, etc.).

- `filtered_*` directories contain STAR-called cells only
- `raw_*` directories contain all barcodes
- `*_velocyto_*` directories contain spliced/unspliced/ambiguous matrices

## Provenance

- **STAR binary**: `STAR.release` (GeneFull + Velocyto mode)
- **Genome index**: GRCh38 (GENCODE v44)
- **Y-removal**: Enabled (`--emitNoYBAM yes --emitYNoYFastq yes`)
- **Cell calling**: STAR EmptyDrops_CR algorithm
- **Doublet detection**: scDblFinder (R, via Docker)
- **CRISPR assignment**: STAR CellRanger-compatible protospacer calling
- **Whitelist**: 3M-february-2018 TRU (GEX/Solo), NXT translated (CR)
- **CellBender**: 0.3.2, GPU mode, remote offload to 10.159.4.53
- **Threads**: 24 STAR threads
- **Input root**: Corrected UCSF perturb-seq symlinked FASTQs
