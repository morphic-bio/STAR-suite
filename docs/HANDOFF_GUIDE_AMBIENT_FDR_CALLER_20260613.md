# HANDOFF - Ambient-FDR Guide Calling For STAR Perturb-Seq

Date: 2026-06-13

Update after QC default work: the implemented integrated STAR default is
`--crGuideCaller auto`, which preserves the root Cell Ranger-compatible GMM
outputs and writes the ambient-FDR sidecar under `outs/crispr_analysis/ambient_fdr/`
when raw and filtered guide MEX are available. Use explicit `--crGuideCaller gmm`
for strict GMM-only parity output.

## Agent Start And Branching

The primary checkout at `/mnt/pikachu/STAR-suite` may be dirty with unrelated
Larry/CBQ work. Do not overwrite, revert, or stage unrelated changes.

Recommended start for a new implementation agent:

```bash
cd /mnt/pikachu/STAR-suite
git fetch origin
git worktree add /mnt/pikachu/STAR-suite-guide-ambient-fdr \
  -b feature/guide-ambient-fdr origin/master
cd /mnt/pikachu/STAR-suite-guide-ambient-fdr
```

If this handoff was created in the dirty primary checkout and is not present in
the new worktree, copy it in before starting implementation:

```bash
cp /mnt/pikachu/STAR-suite/docs/HANDOFF_GUIDE_AMBIENT_FDR_CALLER_20260613.md \
  /mnt/pikachu/STAR-suite-guide-ambient-fdr/docs/
cp /mnt/pikachu/STAR-suite/docs/PROMPT_GUIDE_AMBIENT_FDR_CALLER_20260613.md \
  /mnt/pikachu/STAR-suite-guide-ambient-fdr/docs/
```

Read `AGENTS.md` and, if present, `AGENTS.local.md` before editing. Follow the
STAR-suite branch policy: raw work on `feature/*`, then merge through
`dev-release` before `master`.

## Goal

Add an optional STAR-suite guide-calling mode that complements the existing
CR-compatible GMM caller:

```text
CR-compatible mode: GMM thresholds, Cell Ranger-style crispr_analysis outputs.
Ambient-FDR mode: transparent q-values from empty-droplet guide background,
                  tunable by users and usable directly by QC/MuData.
```

This belongs in STAR-suite, not only Multiomics Suite, because it is useful for
ordinary Perturb-seq and because STAR-suite already owns guide count generation
and `outs/crispr_analysis/`.

## Scope And Non-Goals

In scope:

- add an ambient-FDR guide caller for Perturb-seq guide counts;
- keep current Cell Ranger-compatible GMM caller unchanged and available;
- let integrated STAR multi-feature paths run `gmm`, `ambient-fdr`, `both`, or
  `none`;
- let `star_feature_call --call-only` run ambient-FDR from existing MEX outputs;
- emit machine-readable QC for FDR tuning.

Out of scope for the first implementation:

- changing guide assignment/counting logic;
- changing EmptyDrops or filtered-cell selection;
- replacing Cell Ranger-compatible `crispr_analysis` files;
- materializing dense cell-by-guide q-value matrices;
- computing q-values for rejected/non-cell barcodes.

The raw MEX participates in ambient estimation only. The call universe is the
final filtered cell barcode set.


## Why This Matters

The current automatic STAR perturb-seq path generates Cell Ranger-compatible
`outs/crispr_analysis/` through GMM calling after filtered MEX generation. That
is the right default for compatibility, but it has two limitations:

- the decision threshold is model-derived and not directly interpretable as an
  error rate;
- changing the cutoff generally means re-calling or choosing a different
  thresholding method.

Ambient-FDR guide calling should instead store the statistic needed for any
cutoff:

```text
call at FDR alpha = guide_qvalue <= alpha
```

The QC report and downstream MuData can then expose an FDR slider without
re-running guide assignment.

## Current STAR-Suite Surfaces

Relevant code:

- `core/features/process_features/src/call_features.c`
- `core/features/process_features/include/call_features.h`
- `core/features/process_features/src/call_features_main.c`
- `core/legacy/source/PfMultiProcess.cpp`
- `core/legacy/source/CrMultiProcess.cpp`
- `core/legacy/source/star_feature_call.cpp`
- `core/legacy/source/Parameters.cpp`
- `core/legacy/source/Parameters.h`

Current automatic path:

1. `PfMultiProcess.cpp` writes `outs/raw_feature_bc_matrix/`.
2. It writes `outs/filtered_feature_bc_matrix/` on the final GEX cell set.
3. If any `CRISPR Guide Capture` feature library is present, it filters the
   **filtered** combined MEX to CRISPR features.
4. It writes a temporary CRISPR-only MEX.
5. It calls `cf_process_mex_dir_gmm()`.
6. Outputs land under `outs/crispr_analysis/`.

That is sufficient for GMM compatibility. Ambient-FDR additionally needs raw
guide counts from non-cell barcodes to estimate the ambient null.

The existing GMM caller is already post-cell-calling: it reads the filtered MEX,
filters that matrix to `CRISPR Guide Capture`, and emits calls only for filtered
cell barcodes. Ambient-FDR should preserve that boundary.

## Required Inputs

Ambient-FDR must consume:

- raw CRISPR guide count MEX with all observed barcodes;
- final called-cell barcode set;
- guide feature names/ids;
- optional default FDR threshold, default `0.01`;
- optional minimum guide UMI count floor, default `1`.

In the integrated STAR path, use:

```text
ambient source:   outs/raw_feature_bc_matrix filtered to CRISPR Guide Capture
call/qvalue set:  outs/filtered_feature_bc_matrix/barcodes.tsv(.gz)
output location:  outs/crispr_analysis/ambient_fdr/
```

Raw non-cell barcodes are for background estimation only. Do not compute or
store guide q-values for rejected barcodes.

In `star_feature_call --call-only`, require either:

```text
--guide-caller ambient-fdr --raw-mex-dir RAW --filtered-barcodes CELLS
```

or a documented equivalent.

## Statistical Contract

For guide `g` and final called cell `c`:

```text
ambient_rate[g] = raw guide UMIs in non-cell barcodes for g /
                  raw guide UMIs in non-cell barcodes over all guides

cell_depth[c] = total guide UMIs in c

null[c,g] = Poisson(ambient_rate[g] * cell_depth[c])

pvalue[c,g] = P(X >= observed_count[c,g] | null[c,g])
```

Apply Benjamini-Hochberg across the full called-cell by guide universe:

```text
total_tests = n_called_cells * n_guides
```

Cells/guides with zero observed counts have `pvalue = 1` and should not be
materialized. This matters for scaling.

## Scaling Contract

Do **not** write a dense cell x guide q-value matrix. Genome-scale perturb
screens can have enough cells and guides that dense q-values are not acceptable,
and nearly all cell-guide pairs have zero observed counts.

Default output must be sparse:

```text
guide_qvalues.mtx(.gz)       sparse cell x guide q-values for observed entries
guide_counts.mtx(.gz)        optional sparse cell x guide counts, if not already available
guide_ambient_rates.tsv
guide_fdr_calls_per_cell.csv
guide_fdr_summary.json
```

Only entries with observed count > 0 need q-values stored. Missing entries imply
`qvalue = 1`.

## Output Contract

Keep existing compatibility outputs unchanged:

```text
outs/crispr_analysis/protospacer_calls_per_cell.csv
outs/crispr_analysis/protospacer_calls_summary.csv
outs/crispr_analysis/protospacer_umi_thresholds.csv
outs/crispr_analysis/protospacer_umi_thresholds.json
```

Add ambient-FDR outputs under:

```text
outs/crispr_analysis/ambient_fdr/
```

Required files:

```text
guide_fdr_calls_per_cell.csv
guide_fdr_summary.json
guide_ambient_rates.tsv
guide_qvalues.mtx
guide_qvalues_barcodes.tsv
guide_qvalues_features.tsv
```

Recommended `guide_fdr_calls_per_cell.csv` columns:

```text
cell_barcode
num_features
feature_call
num_umis
min_qvalue
num_features_at_default_fdr
call_status          # none|singlet|multiplet
default_fdr
caller               # ambient-fdr
```

For compatibility with existing downstream code, `feature_call` should use
`None` for no call and `|`-joined guide names for multiplets.

## CLI/Parameter Proposal

For integrated STAR:

```text
--crGuideCaller auto|gmm|ambient-fdr|both|none
--crGuideFdr 0.01
--crGuideFdrMinUmi 1
--crGuideFdrEmitQvalues sparse|none
```

Recommended defaults:

- current CR-compat profiles: `--crGuideCaller auto`
- new perturb-QC profiles: `--crGuideCaller auto` or explicit `both`
- Multiomics/CAT-ATAC profiles: `--crGuideCaller auto` where raw+filtered guide MEX are available
- default FDR: `0.01`
- default q-value storage: `sparse`

For `star_feature_call`:

```text
--guide-caller auto|gmm|ambient-fdr|both|dominant
--guide-fdr 0.01
--guide-fdr-min-umi 1
--raw-mex-dir DIR
--filtered-barcodes FILE
```

`--compat-perturb` in `auto` mode writes GMM compatibility outputs and adds
ambient-FDR QC when raw guide MEX are available. Explicit `gmm` skips the
ambient-FDR sidecar.

## Implementation Order

1. Add a shared ambient-FDR implementation in `process_features`.
   - Read raw and filtered MEX.
   - Filter both matrices to `CRISPR Guide Capture`.
   - Use raw barcodes not present in the filtered cell set as the ambient pool.
   - Compute q-values only for observed entries in filtered cells.
   - Write sparse q-value output and per-cell calls.
2. Expose the caller through `call_features` / `star_feature_call`.
   - Add CLI options for `--guide-caller`, `--guide-fdr`,
     `--guide-fdr-min-umi`, `--raw-mex-dir`, and `--filtered-barcodes` or the
     equivalent existing MEX controls.
   - Keep `--compat-perturb` behavior GMM-only unless explicitly overridden.
3. Integrate with STAR `PfMultiProcess.cpp` and `CrMultiProcess.cpp`.
   - Pass both `raw_feature_bc_matrix` and `filtered_feature_bc_matrix`.
   - Run the selected caller(s) after filtered MEX is written.
   - Preserve current GMM outputs and locations.
4. Add parameters in `Parameters.*` and `parametersDefault`.
   - Default current CR-compatible profiles to GMM.
   - Use `both` only for new perturb-QC or explicitly requested profiles.
5. Add focused tests.
   - Tiny synthetic ambient-FDR unit/regression test.
   - GMM no-regression test when `gmm` is selected.
   - CLI parsing tests for invalid/missing raw MEX in call-only mode.
   - Optional CAT-ATAC smoke after guide-arm integration is available.

## Suggested API Shape

Keep the exact names consistent with local style, but aim for an API like:

```c
typedef struct cf_ambient_fdr_config {
    double fdr_threshold;
    int min_umi;
    int emit_sparse_qvalues;
} cf_ambient_fdr_config;

int cf_process_mex_dir_ambient_fdr(const char *raw_mex_dir,
                                   const char *filtered_mex_dir,
                                   const char *output_dir,
                                   const cf_ambient_fdr_config *config);
```

If `star_feature_call` already has a filtered-barcode extraction path that is
cleaner than passing a full filtered MEX, use it, but keep the internal
implementation centered on:

```text
ambient counts from raw non-cell barcodes
observed guide counts from filtered cells
```

Do not create a second divergent statistical implementation in STAR core.
Process-features should own the shared caller.

## QC Integration

STAR-suite should emit a machine-readable QC JSON even before a full HTML UI is
added:

```text
guide_fdr_qc.json
```

Include:

- assignment rate at default FDR;
- singlet/multiplet/none counts;
- assignment rate grid over common FDR thresholds;
- top guide counts and ambient rates;
- warnings for weak ambient support;
- counts of cells absent from guide MEX, if any.

Multiomics Suite can render the interactive slider from these files and from the
sparse q-value matrix.

## Fallbacks And Warnings

If raw non-cell guide counts are unavailable or empty:

- do not silently fabricate an ambient null;
- log a clear warning;
- keep GMM outputs if requested;
- skip ambient-FDR outputs or write `guide_fdr_summary.json` with
  `status = "skipped_no_ambient"`.

If all ambient guide counts are zero:

- use the same skipped behavior unless an explicit pseudocount option is added.

If filtered barcodes are unavailable:

- in integrated STAR this should not happen after GEX cell calling;
- in `star_feature_call`, require `--filtered-barcodes` for ambient-FDR.

## Validation Plan

### Unit / Synthetic

Create a tiny raw guide MEX:

- 3 guides;
- 5 cell barcodes;
- 20 empty barcodes;
- one clear singlet;
- one multiplet;
- one ambient-like false positive;
- one no-call.

Validate:

- ambient rates sum to 1;
- q-values are monotonic under BH;
- zero entries are not materialized and imply q=1;
- calls change predictably as `--crGuideFdr` changes.

### Existing Perturb-Seq

Run on A375/CR-compatible guide data:

- confirm existing GMM files are byte-stable when `--crGuideCaller gmm`;
- run `auto` or `both` and compare ambient-FDR calls to GMM as secondary output;
- report assignment-rate/FDR curve.

### CAT-ATAC

Run after guide-arm support lands:

- compare count matrix before calls;
- compare ambient-FDR ranked guide evidence to `DMSO1_calls.rds`;
- emit q-value layer for Multiomics Suite MuData/QC.

## Implementation Notes

- Prefer implementing the statistical core in `process_features` so
  `call_features`, `star_feature_call`, `PfMultiProcess`, and `CrMultiProcess`
  can all share it.
- Reuse htslib `kfunc` or existing numeric utilities for Poisson tails if that
  is cleaner than hand-rolling a fragile CDF.
- BH should be sparse-aware: sort only observed-entry p-values, but use
  `n_called_cells * n_guides` as the number of tests.
- Keep memory O(nnz + n_cells + n_guides), not O(n_cells * n_guides).
- Preserve existing `crispr_analysis/` files and names for compatibility.
- Use sparse Matrix Market orientation consistently with existing MEX writers.
- Log the number of raw barcodes, filtered cells, ambient barcodes, guides,
  observed filtered guide entries, and calls at the default FDR.

## Done Criteria

- New mode is documented.
- Existing GMM compatibility remains default and regression-stable.
- Ambient-FDR outputs are produced automatically when requested.
- STAR logs clearly state which guide caller(s) ran.
- Sparse q-value output is accepted by downstream Multiomics Suite.
- Tests cover tiny synthetic, A375 no-regression, and CAT-ATAC once available.
