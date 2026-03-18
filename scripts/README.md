# scripts/

Suite-level helper scripts for STAR-suite. These are standalone tools that
support preflight validation, QC, parity checking, downstream processing, and
fixture management. They are not part of the compiled STAR binary; they run
independently with Python 3 or Bash.

For module-specific scripts see `core/legacy/scripts/README.md` and
`flex/scripts/README.md`.

## Preflight and QC

| Script | Language | Purpose |
|--------|----------|---------|
| `preflight_library_pairing.py` | Python | FASTQ preflight: detects barcode chemistry (TRU/NXT), infers library type (GEX/GUIDE), and pairs libraries by top-barcode Jaccard overlap. Includes duplicate-file merging, confidence (S/N) tables, and post-hoc name-based sanity checks. Designed for Perturb-seq datasets where filenames may be mislabeled. See `docs/RUNBOOK_FASTQ_PREFLIGHT_LIBRARY_PAIRING_20260317.md` for the full spec and validated results. |
| `compute_adaptive_qc_threshold.py` | Python | Computes an adaptive `max_genes` QC threshold from STAR-called singlets using median + N*MAD. Used in downstream filtering pipelines to avoid hardcoded cutoffs. |
| `generate_qc_histogram.py` | Python | Builds QC histograms (genes, UMIs, mito fraction) from unfiltered downstream AnnData output for visual inspection. |
| `run_trim_qc_fastq.sh` | Bash | Wrapper that runs the `trim_qc_fastq` binary for per-sample FASTQ trimming and quality checks. |
| `run_trim_qc_fastq_batch.sh` | Bash | Batch runner: invokes `run_trim_qc_fastq.sh` on every sample directory under a root path. Logs to timestamped files. |

## Parity and Benchmarking

| Script | Language | Purpose |
|--------|----------|---------|
| `report_additional_parity_metrics.py` | Python | Computes parity metrics for STAR vs Cell Ranger: gene-level Pearson correlation, cell-level correlation, and CRISPR feature-call set-match rates. This is the primary numeric parity tool referenced in the paper benchmarks. |
| `compare_barcode_sets.py` | Python | Compares two barcode sets (e.g., STAR vs CR filtered barcodes) and reports overlap, Jaccard, and symmetric-difference statistics. Useful for quick cell-list parity checks. |
| `run_gex_feature_parity_checks.sh` | Bash | End-to-end GEX feature parity: runs STAR and Cell Ranger, then compares gene counts and correlations. |
| `run_ucsf_call_parity_recovery.sh` | Bash | UCSF CRISPR call-only parity recovery: compares STAR vs Cell Ranger GMM-based guide calling on the UCSF dataset. |
| `a375_feature_sweep.sh` | Bash | Runs `assignBarcodes` parameter sweeps on A375 CRISPR FASTQs with configurable barcode/UMI lengths and search settings. |
| `a375_make_downsample_tiers.sh` | Bash | Builds tiered downsample FASTQ fixtures for A375 GEX and CRISPR via symlinks and the process_features downsample script. |

## Downstream Processing

| Script | Language | Purpose |
|--------|----------|---------|
| `build_gene_full_velocyto_h5ad.py` | Python | Builds a `counts.h5ad` from GeneFull Cell-Ranger-style MEX output, packaging raw Velocyto spliced/unspliced/ambiguous layers alongside gene counts. |
| `prepare_velocyto_mex.py` | Python | Packages STARsolo Velocyto raw outputs into raw/filtered MEX directories with total and per-layer matrices. |
| `integrate_feature_library.py` | Python | Merges a feature library (e.g., CRISPR guides) into an existing downstream h5ad, adding the feature counts as a new layer or `obsm` entry. |
| `postprocess_downstream_filters.py` | Python | Rewrites downstream filtered outputs to produce both QC-only and default singlet-filtered versions of the AnnData. |
| `propagate_anndata_layer.py` | Python | Copies an AnnData layer from a full source object into a target subset by matching `obs_names`, preserving target cell order. |
| `filter_fastq_by_readnames.py` | Python | Filters paired FASTQ files (R1/R2) to keep only reads whose names appear in a supplied list. Useful for creating targeted test fixtures. |

## Flex

| Script | Language | Purpose |
|--------|----------|---------|
| `flex_h01_pilot.py` | Python | External 100-probe H0/H1 pilot for Flex: generates probe lists, probe FASTAs, and handles feature-matrix post-processing for the Flex validation pipeline. |
| `run_flex_cr_config.sh` | Bash | Runs STAR-Flex using a CellRanger-format config, with inputs rendered by `render_flex_inputs_from_cr_config.py`. |

## Fixtures and Utilities

| Script | Language | Purpose |
|--------|----------|---------|
| `download_public_bulk_fixture.sh` | Bash | Downloads a fixed public human bulk RNA-seq paired-end fixture from GEO/SRA (default: GSE88509 / GSM2344101). Used by bulk RNA-seq smoke tests. |
| `check_links.sh` | Bash | Validates a symlink manifest: checks that each link exists, is a symlink, and its target is present. |
| `sync_links.sh` | Bash | Creates symlinks from a tab-separated manifest (`link_path<TAB>target_path`). |

## Subdirectories

| Directory | Purpose |
|-----------|---------|
| `codespaces/` | GitHub Codespaces setup helpers (lightweight FASTQ copy, env config). |
| `docker/` | Docker image build and smoke-test runner scripts. See main `README.md` § Docker. |
| `lib/` | Shared shell library functions sourced by other scripts. |
| `paper/` | Scripts specific to paper figures and benchmark tables. |
| `release/` | Release packaging and distribution helpers. |
| `flex_compat/` | Flex compatibility / migration utilities. |
| `ucsf_parity/` | UCSF-specific parity check scripts. |
