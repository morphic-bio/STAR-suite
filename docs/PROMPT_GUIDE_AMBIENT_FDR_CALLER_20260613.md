# Prompt - Implement Ambient-FDR Guide Calling In STAR-Suite

Copy this into a new Codex agent.

```text
You are working in STAR-suite on the ambient-FDR guide caller for Perturb-seq.

Start by reading:

- AGENTS.md
- AGENTS.local.md, if present
- docs/HANDOFF_GUIDE_AMBIENT_FDR_CALLER_20260613.md
- docs/feature_barcodes.md, especially the CRISPR feature calling section
- docs/CRISPR_FEATURE_CALLING_IMPLEMENTATION_SUMMARY.md

Repository / branch setup:

- The default branch is master.
- Do not work directly in a dirty checkout if unrelated local changes are
  present.
- Prefer a clean worktree from master:

  cd /mnt/pikachu/STAR-suite
  git fetch origin
  git worktree add /mnt/pikachu/STAR-suite-guide-ambient-fdr \
    -b feature/guide-ambient-fdr origin/master
  cd /mnt/pikachu/STAR-suite-guide-ambient-fdr

- If the handoff docs are only present in the original dirty checkout, copy them
  into the worktree before implementing:

  cp /mnt/pikachu/STAR-suite/docs/HANDOFF_GUIDE_AMBIENT_FDR_CALLER_20260613.md docs/
  cp /mnt/pikachu/STAR-suite/docs/PROMPT_GUIDE_AMBIENT_FDR_CALLER_20260613.md docs/

Goal:

Implement an optional ambient-FDR guide caller in STAR-suite for Perturb-seq.
This complements, but must not replace, the existing Cell Ranger-compatible GMM
caller and current `outs/crispr_analysis/` outputs.

Core behavior:

- Use raw CRISPR guide counts from `outs/raw_feature_bc_matrix` to estimate the
  ambient guide background from non-cell barcodes.
- Use only the final filtered cell barcode set for q-values and calls.
- Do not compute or store q-values for rejected/non-cell barcodes.
- Do not emit dense cell-by-guide q-value matrices. Dense is out of scope.
- Emit sparse q-values only for observed guide counts in filtered cells.
- Treat missing cell-guide entries as non-calls; statistically their q-value is
  1 and they are not materialized.
- Apply BH using `n_filtered_cells * n_guides` as the number of tests, while
  sorting only observed-entry p-values.
- Keep memory O(nnz + n_cells + n_guides), never O(n_cells * n_guides).

Statistical contract:

For guide g and final called cell c:

ambient_rate[g] = raw guide UMIs in non-cell barcodes for g /
                  raw guide UMIs in non-cell barcodes over all guides

cell_depth[c] = total guide UMIs in c

null[c,g] = Poisson(ambient_rate[g] * cell_depth[c])

pvalue[c,g] = P(X >= observed_count[c,g] | null[c,g])

Then apply Benjamini-Hochberg over the full filtered-cell by guide universe.

Implementation targets:

- Put the shared statistical implementation in process_features, not duplicated
  in STAR core.
- Relevant files to inspect first:
  - core/features/process_features/src/call_features.c
  - core/features/process_features/include/call_features.h
  - core/features/process_features/src/call_features_main.c
  - core/legacy/source/star_feature_call.cpp
  - core/legacy/source/PfMultiProcess.cpp
  - core/legacy/source/CrMultiProcess.cpp
  - core/legacy/source/Parameters.cpp
  - core/legacy/source/Parameters.h
  - core/legacy/source/parametersDefault
- Current automatic GMM path:
  - STAR writes `outs/raw_feature_bc_matrix`.
  - STAR writes `outs/filtered_feature_bc_matrix`.
  - `runCrisprFeatureCalling(filteredOutDir, ...)` reads the filtered MEX.
  - It filters to `CRISPR Guide Capture`.
  - It calls `cf_process_mex_dir_gmm()`.

Suggested CLI:

Integrated STAR:

--crGuideCaller gmm|ambient-fdr|both|none
--crGuideFdr 0.01
--crGuideFdrMinUmi 1
--crGuideFdrEmitQvalues sparse|none

star_feature_call:

--guide-caller gmm|ambient-fdr|both|dominant
--guide-fdr 0.01
--guide-fdr-min-umi 1
--raw-mex-dir DIR
--filtered-barcodes FILE

Default behavior:

- Existing CR-compatible profiles must remain GMM by default.
- `--compat-perturb` should continue to imply GMM unless the user explicitly
  requests `ambient-fdr` or `both`.
- New perturb-QC/Multiomics/CAT-ATAC profiles may use `both`, but do not change
  existing parity defaults without a clear test.

Required outputs:

Keep existing files unchanged:

outs/crispr_analysis/protospacer_calls_per_cell.csv
outs/crispr_analysis/protospacer_calls_summary.csv
outs/crispr_analysis/protospacer_umi_thresholds.csv
outs/crispr_analysis/protospacer_umi_thresholds.json

Add ambient-FDR files under:

outs/crispr_analysis/ambient_fdr/

Required:

guide_fdr_calls_per_cell.csv
guide_fdr_summary.json
guide_ambient_rates.tsv
guide_qvalues.mtx
guide_qvalues_barcodes.tsv
guide_qvalues_features.tsv

Recommended per-cell call columns:

cell_barcode
num_features
feature_call
num_umis
min_qvalue
num_features_at_default_fdr
call_status
default_fdr
caller

Validation:

- Add a tiny synthetic test with a raw guide MEX containing called cells plus
  empty/non-cell barcodes.
- Validate ambient rates, BH monotonicity, sparse-only q-values, no q-values for
  rejected barcodes, and predictable calls as the FDR threshold changes.
- Verify `--crGuideCaller gmm` preserves current GMM outputs.
- Build at least the affected process_features / feature-barcodes tools.
- If touching STAR core parameters or PfMulti/CrMulti integration, do a clean
  rebuild before debugging any failure:

  make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR

Important guardrails:

- Do not materialize dense q-value matrices.
- Do not score all raw barcodes.
- Do not change guide counting or EmptyDrops behavior.
- Do not rename or replace existing CR-compatible output files.
- Do not commit generated binaries, test outputs, or large data.
- Do not revert unrelated local changes.

When finished:

- Summarize implementation and changed files.
- Report exact tests run.
- Note any residual limitations, especially numeric edge cases or missing
  dataset-level validation.
```
