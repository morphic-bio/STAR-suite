# scripts/

Suite-level helper scripts for STAR-suite. These are standalone tools that
support preflight validation, QC, parity checking, downstream processing, and
fixture management. They are not part of the compiled STAR binary; they run
independently with Python 3 or Bash.

For module-specific scripts see `core/legacy/scripts/README.md` and
`flex/scripts/README.md`.

## Production Recipe Boundary

Morphic/JAX/MSK production launchers and remote downstream/CellBender helpers
are moving to `/mnt/pikachu/morphic-recipes`. The STAR-suite copies of those
entrypoints are compatibility launchers that delegate through
`MORPHIC_RECIPES_ROOT` when present. New production recipe work should be made
in `morphic-recipes`; keep STAR-suite focused on core processing, shared
regression tools, and generic fixtures.

## Preflight and QC

| Script | Language | Purpose |
|--------|----------|---------|
| `preflight_library_pairing.py` | Python | FASTQ preflight: detects barcode chemistry (TRU/NXT), infers library type (GEX/GUIDE), and pairs libraries by top-barcode Jaccard overlap. Includes duplicate-file merging, confidence (S/N) tables, and post-hoc name-based sanity checks. Designed for Perturb-seq datasets where filenames may be mislabeled. See `docs/RUNBOOK_FASTQ_PREFLIGHT_LIBRARY_PAIRING_20260317.md` for the full spec and validated results. |
| `preflight_whitelist_family.py` | Python | FASTQ preflight: samples R1 cell barcodes and scores them against multiple whitelist families/namespaces, e.g. February-2018 TRU/NXT vs May-2023 GEM-X TRU/NXT. Emits summary/rates/JSON reports and exits non-zero on expected-whitelist mismatches. |
| `run_msk_30ko_fastq_preflight.sh` | Bash | MSK 30KO wrapper around `preflight_whitelist_family.py`; validates `docs/MSK_30KO_FASTQ_MANIFEST.tsv` against February-2018 and May-2023/GEM-X whitelist families before STAR production runs. |
| `run_msk_40ko_fastq_preflight.sh` | Bash | MSK 40KO wrapper around `preflight_whitelist_family.py`; validates `docs/MSK_40KO_FASTQ_MANIFEST.tsv` and confirms the February-2018 TRU/NXT chemistry assignment used by the 40KO production wrapper. |
| `compute_adaptive_qc_threshold.py` | Python | Computes an adaptive `max_genes` QC threshold from STAR-called singlets using median + N*MAD. Used in downstream filtering pipelines to avoid hardcoded n_genes cutoffs. |
| `scrna_mt_adaptive.py` | Python | Shared adaptive mitochondrial-percentage QC helpers: computes the MT strict-floor + MAD soft guard and rewrites downstream AnnData filter masks. |
| `apply_adaptive_mt_filter.py` | Python | Applies the adaptive MT percentage guard to `unfiltered_counts.h5ad` after `combineFilters.py`, preserving legacy strict-5% masks for audit. |
| `generate_qc_histogram_mt_adaptive.py` | Python | Builds QC histograms with both MT floor and adaptive MT threshold lines from `adaptive_qc_threshold.json`. |
| `convert_h5ad_mt_adaptive.py` | Python | Retrofitting tool for existing downstream h5ad releases: applies the adaptive MT guard, rebuilds filtered views, and updates threshold JSON. |
| `generate_qc_histogram.py` | Python | Legacy fixed-MT QC histogram generator for non-adaptive downstream runs. |
| `run_trim_qc_fastq.sh` | Bash | Wrapper that runs the `trim_qc_fastq` binary for per-sample FASTQ trimming and quality checks. |
| `run_trim_qc_fastq_batch.sh` | Bash | Batch runner: invokes `run_trim_qc_fastq.sh` on every sample directory under a root path. Logs to timestamped files. |

## Parity and Benchmarking

| Script | Language | Purpose |
|--------|----------|---------|
| `report_additional_parity_metrics.py` | Python | Computes parity metrics for STAR vs Cell Ranger: gene-level Pearson correlation, cell-level correlation, and CRISPR feature-call set-match rates. This is the primary numeric parity tool referenced in the paper benchmarks. |
| `compare_barcode_sets.py` | Python | Compares two barcode sets (e.g., STAR vs CR filtered barcodes) and reports overlap, Jaccard, and symmetric-difference statistics. Useful for quick cell-list parity checks. |
| `run_gex_feature_parity_checks.sh` | Bash | End-to-end GEX feature parity: runs STAR and Cell Ranger, then compares gene counts and correlations. |
| `run_ucsf_gexonly_gex_parity_vs_cr.sh` | Bash | GEX-only STAR `Solo.out/GeneFull` vs CR9 Gene Expression MEX; calls `report_additional_parity_metrics.py --skip-feature-call-parity` (no STAR `crispr_analysis`). |
| `run_ucsf_call_parity_recovery.sh` | Bash | UCSF CRISPR call-only parity recovery: compares STAR vs Cell Ranger GMM-based guide calling on the UCSF dataset. |
| `a375_feature_sweep.sh` | Bash | Runs `assignBarcodes` parameter sweeps on A375 CRISPR FASTQs with configurable barcode/UMI lengths and search settings. |
| `a375_make_downsample_tiers.sh` | Bash | Builds tiered downsample FASTQ fixtures for A375 GEX and CRISPR via symlinks and the process_features downsample script. |
| `run_jax_scrnaseq02_ocm_oracle_smoke.sh` | Bash | Harness for the JAX scRNAseq02 `25E32-L3` OCM oracle: stages a 2M read-pair downsample or full-depth FASTQ symlinks, renders an existing-binary STAR GeneFull+Velocyto/Y-removal smoke command with Velocyto low-memory spill defaults, and optionally runs the Cell Ranger multi layout validator. |
| `run_jax_scrnaseq02_ocm_composite_smoke.sh` | Bash | Comparator and timing harness for JAX scRNAseq02 OCM: defaults to native STAR `--ocmMultiBarcodeMode flex` (`CB16+OCM_TAG8` before correction/counting, with production `--soloFeatures GeneFull Velocyto` and native per-sample EmptyDrops after OCM split), native streaming `--ocmMultiBamSplit yes` with per-sample unsplit/Y/noY BAMs, Y/noY BAM/FASTQ split outputs (`STAR_YREMOVE=yes`), no barcode/UMI BAM tags or final CB/UB replay (`STAR_BAM_CBUB_TAGS=no`), no `GX/GN` BAM gene tags (`STAR_BAM_GXGN_TAGS=no`), and native OCM GeneFull/Velocyto materialization. Reuses completed Cell Ranger 9 multi artifacts by default for comparison; only runs Cell Ranger when `--run-cr` is explicit. Python materialization is retained for selectable STAR comparison features and fallback production routing. |
| `validate_jax_scrnaseq02_ocm_oracle.py` | Python | Validates STAR OCM multi-compatibility outputs against the `25E32-L3` Cell Ranger oracle bundle: `outs/multi`, `outs/per_sample_outs`, `cells_per_tag.json`, per-sample filtered barcode counts, oracle barcode recall, OCM tag-proportion warnings, and optional oracle-precision overlap. It does not assert GeneFull matrix-value parity because the oracle is exonic `Gene`. |
| `run_multiome_cell_call_from_arc.sh` | Bash | Thin driver for the experimental multiome cell-calling path: extracts an ARC barcode table, builds `gex_evidence.tsv` and `atac_evidence.tsv`, runs the multiome combiner, and optionally compares the output against ARC labels in one outdir. |
| `run_multiome_cell_call_external_gex_from_arc.sh` | Bash | Higher-level experimental runner: extracts a `Gene Expression` MEX from ARC `raw_feature_bc_matrix`, runs `scrna_simpleed`, and then invokes `run_multiome_cell_call_from_arc.sh` with the resulting GEX evidence surface. |
| `run_multiome_cell_call_harness_from_arc.sh` | Bash | Sequential harness around the ARC-driven multiome runner: executes a small set of named `scrna_simpleed` / combiner profiles against one completed ARC run and writes a single `harness_summary.tsv` for side-by-side comparison. |
| `extract_cr_feature_type_mex.py` | Python | Extracts a single `feature_type` surface such as `Gene Expression` from a Cell Ranger `raw_feature_bc_matrix` and rewrites it as a smaller MEX directory. |

## Downstream Processing

| Script | Language | Purpose |
|--------|----------|---------|
| `build_gene_full_velocyto_h5ad.py` | Python | Builds a `counts.h5ad` from GeneFull Cell-Ranger-style MEX output, packaging raw Velocyto spliced/unspliced/ambiguous layers alongside gene counts. |
| `package_star_genefull_mex.py` | Python | Packages STARsolo `Solo.out/GeneFull` raw/filtered MEX into CellRanger-style `outs/raw_feature_bc_matrix` and `outs/filtered_feature_bc_matrix` directories for downstream h5ad/MuData construction. |
| `prepare_velocyto_mex.py` | Python | Legacy repair/backfill packager for older STAR runs. Current STAR binaries write raw/filtered Velocyto `outs/` MEX internally. |
| `compare_velocyto_mex.py` | Python | Exact diff between two STAR runs: `--mode solo` (Velocyto `*.mtx` + mandatory `Solo.out/Gene/raw` axes), `--mode packaged`, `--mode all`, or `--mode genes` (`Solo.out/Gene` + `GeneFull` raw/filtered only; Phase 6). |
| `report_velocyto_sorted_replay_rss.py` | Python | Reads `Log.out` for VmRSS after sorted-replay materialization and integrated-hash spill lines; prints `MAX_VM_RSS_KB` plus `PER_LOG_MAX_VM_RSS_KB[<subdir>]` per input log. |
| `run_star_velocyto_canonical.sh` | Bash | Pinned UCSF perturb STAR for `Gene`/`GeneFull`/`Velocyto` parity; passes through `STAR_VELOCYTO_INTEGRATED_HASH` (Stage 2 CB-bucket path with deterministic replay). Refuses non-fresh `--out-prefix` unless `UCSF_VELOCYTO_REUSE_STAR_OUTDIR=1`. |
| `save_velocyto_baseline.sh` | Bash | Canonical STAR snapshot for `UCSF_VELOCYTO_BASELINE_OUTDIR` / `--mode all` parity; expects native Velocyto `outs/` MEX. |
| `run_ucsf_perturb_yremove_batch.sh` | Bash | Corrected-layout UCSF perturb production runner. Uses `STAR.release` by default, emits Y/noY BAM + FASTQ outputs, verifies native raw/filtered Velocyto MEX, optionally runs downstream h5ad generation + CellBender, and prepares per-sample Globus transfer manifests/helpers. |
| `run_ucsf_downstream_batch.sh` | Bash | Batch wrapper over `run_scrna_downstream_gene_full_velocyto.sh` for existing UCSF production roots. |
| `run_scrna_downstream_gene_full_velocyto.sh` | Bash | Downstream GeneFull + Velocyto packaging: builds `counts.h5ad`, applies adaptive n_genes + MT percentage QC, writes filtered/singlet-filtered h5ads, integrates feature libraries, and optionally adds a CellBender denoised layer plus `final_counts.h5ad`. |
| `run_msk_40ko_pipeline_from_manifest.py` | Python | Manifest-driven MSK 40KO production wrapper modeled on MSK 30KO: stages GEX/PolyIII/LARRY FASTQs, writes `pf_multi_config.csv`, runs CR-compatible STAR GeneFull+Velocyto with native Velocyto MEX output, optionally submits remote downstream/CellBender jobs, and transfers BAMs through Globus. See `docs/RUNBOOK_MSK_40KO_PRODUCTION.md`. |
| `build_multiome_mudata.py` | Python | Builds RNA+ATAC `.h5mu` outputs from STAR downstream RNA h5ad or GEX MEX plus ATAC peak MEX, preserving Velocyto layers when required. Supports sparse-smoke `--allow-empty-barcode-intersection` and records STAR/Chromap/Y-removal provenance in `mdata.uns["multiome"]`. |
| `run_star_multiome_lane_smoke.sh` | Bash | STAR/Chromap sample runner for 10x Multiome FASTQs. Runs GeneFull+Velocyto GEX with Y/noY outputs, native ATAC barcode parsing (`bc:8:23:-` by default), low-memory Chromap sidecar/BAM materialization, native ATAC peak MEX, optional remote post-MEX handoff, and `LOCAL_MEX_READY.txt`. |
| `run_remote_multiome_post_mex_rsync.sh` | Bash | Preferred multiome remote boundary: stages only RNA MEX plus local ATAC peak MEX/metrics to the GPU server, runs RNA downstream/CellBender and MuData remotely, syncs downstream h5ad and `.h5mu` outputs back, excludes remote `mudata_venv/`, and writes `REMOTE_POST_MEX_READY.txt`. |
| `run_multiome_mudata_smoke.sh` | Bash | Production-shaped multiome smoke wrapper: splits ARC MEX, regenerates missing GeneFull/Velocyto MEX packaging, runs remote GeneFull+Velocyto+CellBender downstream when configured, and writes ARC-only or STAR-RNA/ARC-ATAC MuData outputs. |
| `run_jax_multiome01_production.sh` | Bash | Manifest-driven JAX_Multiome01 production wrapper. Builds the sample manifest from the metadata workbook, runs each sample to the local STAR/Chromap MEX boundary, queues serialized remote post-MEX workers, and can trigger Globus large-file upload/cleanup per completed sample. |
| `upload_jax_multiome01_large_files_globus.sh` | Bash | Per-sample Globus uploader for JAX_Multiome01 large files. Uploads raw FASTQs, generated Y/noY FASTQs, and BAMs by default; waits for Globus task success before deleting generated local BAM/FASTQ files while preserving raw input FASTQs. |
| `normalize_multiome_atac_barcode_fastq.py` | Python | Legacy fallback for materializing ARC ATAC barcode reads before Chromap. Production uses native `--chromapAtacReadFormat` in STAR/Chromap instead. |
| `integrate_feature_library.py` | Python | Merges a feature library (e.g., CRISPR guides) into existing downstream h5ads, adding per-library obs annotations. Supports `--ambient-fdr-calls-csv` to import ambient-FDR guide call support columns such as `guide_fdr_num_umis`, `guide_fdr_min_called_umi`, and `guide_fdr_min_qvalue`. |
| `postprocess_downstream_filters.py` | Python | Rewrites downstream filtered outputs to produce both QC-only and default singlet-filtered versions of the AnnData. |
| `propagate_anndata_layer.py` | Python | Copies an AnnData layer from a full source object into a target subset by matching `obs_names`, preserving target cell order. |
| `filter_fastq_by_readnames.py` | Python | Filters paired FASTQ files (R1/R2) to keep only reads whose names appear in a supplied list. Useful for creating targeted test fixtures. |

## Flex

| Script | Language | Purpose |
|--------|----------|---------|
| `flex_h01_pilot.py` | Python | Flex H0/H1/H2 tooling: probe lists/FASTAs, MEX→sequence-cache post-processing, **H2** subcommands (`h2-make-synth-fastq`, `h2-build-cache-from-mex`, `h2-write-binary-cache`) for two-mismatch variants with `cache_class=3` KEEP rows. |
| `run_flex_h02_pilot.sh` | Bash | End-to-end H2 pilot: synthetic FASTQ from H0 seeds in an FH01SEQ1 cache → STAR-Flex (BAM+GX) → MEX-derived TSV → optional `h2_keep_only.bin`. See script header for sharding env vars. |
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
