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
