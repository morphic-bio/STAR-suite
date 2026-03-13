# Follow-Up: scRNA Downstream MEX Reuse and Velocyto Constraints

Date: 2026-03-13
Source branch: `scrna-downstream`
Related handoff:
- `docs/HANDOFF_SCRNA_DOWNSTREAM_MEX_20260313.md`

## Objective

Record what exists today for downstream reuse of STAR-suite MEX outputs and
what blocks post hoc velocyto recovery from the current UCSF perturb run.

## scRNA-seq Repo Surfaces Reused

Repo inspected:

- `/mnt/pikachu/scRNA-seq`

Relevant downstream surfaces found there:

1. `workflows/scRNA_seq_features/widgets/scRNA_seq_features/read_counts/Dockerfiles/read_counts.py`
   - Detects Cell Ranger / STARsolo / kallisto / fry layouts.
   - For Cell Ranger-style outputs, reads:
     - `outs/filtered_feature_bc_matrix/`
     - `outs/raw_feature_bc_matrix/`
   - Marks filtered barcodes in `obs['is_cell']`.
   - For STARsolo layouts, expects `Solo.out/Velocyto/raw/{spliced,unspliced,ambiguous}.mtx`
     and attaches them as layers.

2. `workflows/scRNA_seq_features/widgets/scRNA_seq_features/read_counts/Dockerfiles/findValidCells.R`
   - Runs EmptyDropsCellRanger unless `SKIP_ED` is set.
   - Runs `scDblFinder` on the non-empty droplets.
   - Writes:
     - `non_empty_barcodes.txt`
     - `doublet_barcodes.txt`
     - `filtered_barcodes_with_scores.txt`

3. `workflows/scRNA_seq_features/widgets/scRNA_seq_features/read_counts/Dockerfiles/combineFilters.py`
   - Adds barcode-level annotations into `AnnData.obs`.
   - Computes:
     - `non_empty`
     - `doublet`
     - `doublet_scores`
     - `n_genes`
     - `mt_counts`
     - `total_counts`
     - `mt_pct`
     - `singlet`
     - `filter`
     - `singlet_filtered`
   - Writes:
     - `unfiltered_counts.h5ad`
     - `filtered_counts.h5ad`
   - With `--adaptive_filter`, also writes:
     - `gene_quantile_histogram.html`
     - `gene_quantile_histogram.png`

4. `scripts/run_flex_master.sh`
   - Contains an embedded Python MEX-to-`counts.h5ad` converter.
   - This path is convenient, but it is a minimal loader:
     - reads `matrix.mtx[.gz]`, `barcodes.tsv[.gz]`, `features.tsv[.gz]`
     - stores only the first column of `features.tsv` as `var_names`
   - It does not preserve the feature name / feature type columns from the
     3-column 10x `features.tsv.gz`.

## Implications For UCSF Perturb MEX

The UCSF perturb MEX is not gene-only.

From the original handoff:

- `filtered_feature_bc_matrix/features.tsv.gz` has three columns:
  - feature ID
  - feature name
  - feature type
- It contains both:
  - `Gene Expression`
  - `CRISPR Guide Capture`

This matters for downstream reuse:

- The minimal loader embedded in `run_flex_master.sh` loses feature-type
  metadata unless it is extended.
- The legacy `read_counts.py` path is a better starting point if we want to
  preserve the ordinary Cell Ranger-style import behavior.
- Any downstream workflow must not assume that the input MEX is gene-only.

## Best Current Reuse Path

For ordinary downstream QC/filtering from the finished UCSF outputs:

1. Start from one sample directory such as:
   - `/mnt/pikachu/ucsf-perturb-yremove_20260311_224631/samples/EBs1_1_AALG1/run`

2. Reuse the scRNA-seq workflow components in this order:
   - `read_counts.py` to create `counts.h5ad`
   - `findValidCells.R` to produce non-empty and doublet calls
   - `combineFilters.py --adaptive_filter` to create filtered `h5ad` outputs
     and QC plots

3. Join guide calls from:
   - `outs/crispr_analysis/protospacer_calls_per_cell.csv`

4. Keep sample identity in metadata before merging samples.

## Velocyto Findings For The UCSF Run

The current UCSF run was not executed with `Velocyto` as a STARsolo feature.

Observed in:

- `/mnt/pikachu/ucsf-perturb-yremove_20260311_224631/samples/EBs1_1_AALG1/run/Log.out`

Key command-line finding:

- `--soloFeatures GeneFull`

There is no local evidence of:

- `Solo.out/Velocyto/`
- `--soloFeatures ... Velocyto`

## BAM Availability Constraint

The current local UCSF sample directories contain `Aligned.out.bam`, but these
files are header-only stubs after cleanup, not full alignment payloads.

Observed across all 20 samples:

- `Aligned.out.bam` size: `4310` bytes

`samtools idxstats` on a sample reports only header references with zero mapped
records, consistent with payload cleanup.

This matches the original handoff note that BAM/YFASTQ payloads were
transferred/cleaned.

## Conclusion: Can We Recover Velocyto Without Re-Running STAR?

For the local artifacts currently present under:

- `/mnt/pikachu/ucsf-perturb-yremove_20260311_224631`

the answer is no.

Reason:

- the run only produced `GeneFull` locally
- no local `Solo.out/Velocyto` outputs exist
- the retained `Aligned.out.bam` files are header-only cleanup stubs, not full
  alignments

## What Would Make Post Hoc Velocyto Recovery Possible

Any one of the following would reopen the path:

1. Recover the full aligned BAM payloads that were transferred off this local
   run directory, then generate velocyto counts from those BAMs.

2. Re-run STAR/STARsolo with `Velocyto` enabled in `--soloFeatures`.

3. If there is another retained run directory elsewhere with full BAMs and the
   same reference/chemistry, use that as the input surface instead of the local
   cleaned UCSF directory.

## Suggested Next Branch Goal

The next feature branch should focus on one of these two paths:

- add a documented post hoc velocyto workflow from full retained BAMs
- add first-class STAR-suite support for keeping or emitting `Velocyto` outputs
  in the UCSF/CR-compat execution path

## Update On `velocyto-posthoc` Branch

The branch now includes a working path to package STARsolo `Velocyto` outputs
into Cell Ranger-style raw/filtered MEX directories.

Added:

- `scripts/prepare_velocyto_mex.py`
  - reads `Solo.out/Velocyto/raw/{spliced,unspliced,ambiguous}.mtx`
  - reuses `Solo.out/Gene/{raw,filtered}` barcode axes
  - writes:
    - `outs/raw_velocyto_feature_bc_matrix/`
    - `outs/filtered_velocyto_feature_bc_matrix/`
    - `outs/velocyto_feature_bc_matrix_manifest.json`

- `tests/run_perturb_velocyto_mex_smoke.sh`
  - runs the UCSF 100K perturb fixture with:
    - `--soloFeatures Gene GeneFull Velocyto`
    - `--soloCrGexFeature genefull`
  - validates that:
    - raw/filtered velocyto MEX directories are created
    - `matrix.mtx.gz == spliced + unspliced + ambiguous`
    - filtered barcodes are a subset of raw barcodes

## UCSF 100K Validation Result

Validated successfully on:

- `/storage/ucsf-2M/fixtures/ucsf2m_iPSC2_AALG2_100k_pfconfig`

Smoke output:

- `/mnt/pikachu/STAR-suite/tests/perturb_velocyto_mex_smoke_output_20260313_061919`

Observed summary:

- raw velocyto features: `38606`
- raw velocyto barcodes: `3686400`
- filtered velocyto barcodes: `3261`
- raw total nnz: `3145`
- filtered total nnz: `3145`

## UCSF 2M Validation Result

Validated successfully on the larger staged UCSF 2M set:

- GEX:
  - `/storage/ucsf-2M/GEX/iPSC2_1_AALG2`
- guides:
  - `/storage/ucsf-2M/guides/iPSC2_1_AALG2`
- pfMulti config:
  - `/storage/ucsf-2M/star_runs/star_baseline_iPSC2_1_AALG2_1M_20260217_065351/multi_config.csv`

Smoke output:

- `/mnt/pikachu/STAR-suite/tests/perturb_velocyto_mex_smoke_output_20260313_064248`

Observed summary:

- raw velocyto features: `38606`
- raw velocyto barcodes: `3686400`
- filtered velocyto barcodes: `6832`
- raw total nnz: `30436`
- filtered total nnz: `26987`

This rerun also passed with multithreading enabled:

- `PERTURB_VELOCYTO_MEX_THREADS=8`

## Root Cause Of Initial Failure

The first multithreaded UCSF 100K run failed with:

- `[ERROR] Conflicting CB for readId=0 existing=229421 new=3402852`

Root cause:

- `GeneFull` does not carry `readIndex` in this configuration
- the loader therefore generated synthetic read IDs starting at zero per thread
- `CountingSink` incorrectly treated those per-thread synthetic IDs as global
  keys and rejected legitimate records from another thread

Fix applied:

- `flex/source/SoloReadInfoSink.cpp`
  - only enforce the per-read CB conflict guard when a real `readIndex` is
    present
  - skip that guard for synthetic thread-local IDs

Result:

- the UCSF 100K perturb velocyto smoke now completes successfully with
  `PERTURB_VELOCYTO_MEX_THREADS=4`

## STAR-Cell Downstream Path

For downstream QC on the packaged GeneFull + velocyto output, do not run
external EmptyDrops on the 100K downsampled smoke. It is too sparse for that
path and STAR already provides the cell calls we want to honor.

Implementation added in this branch:

- `scripts/build_gene_full_velocyto_h5ad.py`
  - builds `counts.h5ad` with GeneFull counts and raw velocyto layers
- `scripts/run_star_cell_doublets.R`
  - writes `non_empty_barcodes.txt` from STAR-called cells (`obs.is_cell`)
  - runs `scDblFinder` on STAR-called cells only
  - writes `doublet_barcodes.txt` and `filtered_barcodes_with_scores.txt`
- `scripts/run_scrna_downstream_gene_full_velocyto.sh`
  - uses `outs/filtered_feature_bc_matrix` as the STAR-called cell set
  - runs the R/Python downstream steps inside `biodepot/scrna-matrices:latest`
  - sets writable `NUMBA` and matplotlib cache dirs inside the output tree

100K smoke rerun:

- command:
  - `scripts/run_scrna_downstream_gene_full_velocyto.sh --run-dir tests/perturb_velocyto_mex_smoke_output_20260313_061919/run --output-dir tests/perturb_velocyto_mex_smoke_output_20260313_061919/downstream_genefull_velocyto_100k_starcells --min-genes 1 --mt-pct-cutoff 100 --adaptive-filter`
- output:
  - `tests/perturb_velocyto_mex_smoke_output_20260313_061919/downstream_genefull_velocyto_100k_starcells`

Observed result on the 100K STAR-cell smoke:

- STAR-called cells: `4082`
- `scDblFinder` doublets: `248`
- singlets after doublet removal: `3834`
- cells passing QC filter: `3118`
- singlets passing QC filter: `3113`

Artifacts written:

- `counts.h5ad`
- `unfiltered_counts.h5ad`
- `filtered_counts.h5ad`
- `non_empty_barcodes.txt`
- `doublet_barcodes.txt`
- `filtered_barcodes_with_scores.txt`
- `gene_quantile_histogram.html`
- `gene_quantile_histogram.png`

Notes:

- The 100K fixture is extremely sparse; the adaptive smoke above used
  `--min-genes 1` and `--mt-pct-cutoff 100` only to exercise the path.
- With these data, adaptive `max_genes` collapsed to `1` because the singlet
  median and MAD were both `1` and `0`, respectively.

## CellBender Integration Method

The existing scRNA-seq CellBender path does not replace `X` directly. It does:

- `remove_noise.sh`
  - runs `cellbender remove-background`
  - default pipeline input is `filtered_counts.h5ad`
  - writes `cellbender/cellbender_counts.h5`
- `addCounts.py`
  - loads the original h5ad
  - loads CellBender output with `cellbender.remove_background.downstream.anndata_from_h5`
  - maps CellBender barcode rows back onto the original AnnData by barcode name
  - adds the denoised counts as a new sparse layer, default layer name:
    - `denoised`
  - writes:
    - `final_counts.h5ad`

This is the methodology used by:

- `/mnt/pikachu/scRNA-seq/scripts/run_flex_master.sh`
- `/mnt/pikachu/scRNA-seq/scripts/rerun_cellbender_sequential.py`

## UCSF 100K CellBender Check

I ran the existing CellBender container script (`remove_noise.sh` inside
`biodepot/cellbender:0.3.2`) against:

- `tests/perturb_velocyto_mex_smoke_output_20260313_061919/downstream_genefull_velocyto_100k_starcells/filtered_counts.h5ad`

Result:

- CellBender failed before producing `final_counts.h5ad`
- it warned that only `3117` barcodes were present and that the input appeared
  prefiltered
- failure terminated in CellBender prior estimation with:
  - `IndexError: index -100 is out of bounds for axis 0 with size 0`

Interpretation:

- the scRNA-seq integration method is clear and reusable
- but the 100K STAR-cell smoke is too sparse/prefiltered for CellBender itself
- for a real CellBender run we should feed it a true raw unfiltered barcode
  matrix/h5ad, then add the denoised layer back into the AnnData object

## Feature-Library Integration

The older scRNA-seq `gather_features.py` workflow is not directly compatible
with current CR-compat `cr_assign` outputs:

- the legacy script expects `features_matrix.mtx`
- current `cr_assign` libraries expose `matrix.mtx`
- the GEX-aligned CRISPR calls now come from:
  - `outs/crispr_analysis/protospacer_calls_per_cell.csv`
- that call file provides:
  - `num_features`
  - `feature_call`
  - `num_umis`
- it does **not** provide the top-two feature-count fields that the older
  `best_feature` / `feature1_count` / `feature2_count` logic used

Implementation added in this branch:

- `scripts/integrate_feature_library.py`
  - runs inside `biodepot/gather_features:latest`
  - discovers one feature library at a time from `run/cr_assign/*/*/<sample>`
  - writes per-library feature-matrix h5ads:
    - `feature_libraries/<library_id>/raw_feature_library.h5ad`
    - `feature_libraries/<library_id>/filtered_feature_library.h5ad`
  - annotates `counts.h5ad`, `unfiltered_counts.h5ad`, `filtered_counts.h5ad`
    (and `final_counts.h5ad` if present) with library-scoped obs columns
- `scripts/run_scrna_downstream_gene_full_velocyto.sh`
  - now auto-discovers feature libraries after QC / optional CellBender
  - runs the feature-library helper in the gather container
  - regenerates a single `summary.txt` after feature annotations are added

Current counts obs columns added for the UCSF CRISPR Guide Capture library:

- `CRISPR_Guide_Capture_iPSC2_1_AALG2_CRISPR_Guide_Capture_1__is_featured`
- `CRISPR_Guide_Capture_iPSC2_1_AALG2_CRISPR_Guide_Capture_1__feature_call`
- `CRISPR_Guide_Capture_iPSC2_1_AALG2_CRISPR_Guide_Capture_1__num_features`
- `CRISPR_Guide_Capture_iPSC2_1_AALG2_CRISPR_Guide_Capture_1__num_umis`

When exactly one compatible feature library is present, the wrapper also adds
generic aliases:

- `is_featured`
- `feature_call`
- `feature_call_num_features`
- `feature_call_num_umis`

Safety rule:

- if more than one `CRISPR_Guide_Capture` library is present but STAR only
  exposes the single global
  `outs/crispr_analysis/protospacer_calls_per_cell.csv`, the wrapper exits
  instead of guessing the library mapping

100K validation:

- command:
  - `scripts/run_scrna_downstream_gene_full_velocyto.sh --run-dir tests/perturb_velocyto_mex_smoke_output_20260313_061919/run --output-dir tests/perturb_velocyto_mex_smoke_output_20260313_061919/downstream_genefull_velocyto_100k_starcells_featurelibs_v2 --min-genes 1 --mt-pct-cutoff 100 --adaptive-filter`
- output:
  - `tests/perturb_velocyto_mex_smoke_output_20260313_061919/downstream_genefull_velocyto_100k_starcells_featurelibs_v2`

Observed result:

- `counts.h5ad`: `4082 x 38606`
- `unfiltered_counts.h5ad`: `4082 x 38606`
- `filtered_counts.h5ad`: `3118 x 38606`
- filtered cells with CRISPR calls:
  - `1337`
- feature-library raw h5ad:
  - `8368 x 548`
- feature-library filtered h5ad:
  - `2725 x 548`

Artifacts written:

- `feature_libraries/iPSC2_1_AALG2_CRISPR_Guide_Capture_1/raw_feature_library.h5ad`
- `feature_libraries/iPSC2_1_AALG2_CRISPR_Guide_Capture_1/filtered_feature_library.h5ad`
- `feature_libraries/iPSC2_1_AALG2_CRISPR_Guide_Capture_1/manifest.json`
