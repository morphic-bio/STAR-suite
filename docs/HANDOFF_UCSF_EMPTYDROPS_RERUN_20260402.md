# Handoff: Apply Fixed EmptyDrops Routine to UCSF Production Set (2026-04-02)

## Status

In progress on `feature-emptyDrops`. The guarded ambient fix is now in
`core/features/libscrna/src/OrdMagStage.cpp` but not yet merged to `master`.

## Objective

Re-run the EmptyDrops cell filtering step on the 16 UCSF production samples
using the fixed ambient cell selection logic from the `feature-emptyDrops`
branch. The current production run already completed with `--soloCellFilter
EmptyDrops_CR`, but it used the old `libscrna` ambient pool which is known to
under-rescue tail cells compared to the legacy rank-window approach.

## Background: The Ambient Pool Bug

The `libscrna` EmptyDrops path constructs its ambient profile from the
`SimpleED` low-UMI cell set — typically a very large pool of barcodes below the
knee. The legacy STAR backend instead uses a fixed rank window (e.g.
`45,000–89,999`). On the MSK 30polyKO dataset, this difference alone accounts
for nearly all of the parity gap vs CellRanger 9:

- current `libscrna` ambient: `207,413` cells → `30,411` final cells (Jaccard vs CR9: `0.943`)
- legacy rank-window ambient: `45,000` cells → `32,303` final cells (Jaccard vs CR9: `0.998`)

Full analysis in:

- `/mnt/pikachu/STAR-suite/docs/MSK_AMBIENT_RESCUE_ISOLATION_20260402.md`

## What Was Fixed

The `feature-emptyDrops` branch contains two categories of changes:

### 1. Committed fix (already on master)

- **Commit `b39a2e5d`**: Fix `soloCellFiltering` EmptyDrops flag parsing in
  `ParametersSolo.cpp`. Ensures `--soloEmptyDropsLegacy`, `--soloEmptyDropsLegacyKnee`,
  and `--soloEmptyDropsMode` are parsed at the correct point in initialization.

### 2. Pending changes on `feature-emptyDrops`

These need to be committed and merged before the rerun.

**`core/features/libscrna/src/OrdMagStage.cpp`**
- Guarded minimum for `minAmbientCells` in both `runCRSimpleFilterBootstrap()`
  and `runCRSimpleFilter()`:
  - Old: `(nCB >= 1000) ? (nCB / 10) : min(100, nCB)`
  - New: for `nCB >= 1000`, use `min(nCB/10, max(5000, nCB/50))` — at least 2%
    or 5000 cells, capped at 10%. This prevents the fallback from swallowing
    the entire rank window on large datasets.

**`core/features/libscrna/tools/scrna_simpleed.cpp`**
- New CLI flags: `--lower-testing-bound`, `--ambient-umi-max`,
  `--apply-bh-correction`, `--use-bootstrap`, `--direct-ed-surface`,
  `--use-legacy-rank-ambient`, `--use-guarded-rank-ambient`,
  `--ambient-fallback-min-abs`, `--ambient-fallback-min-frac`
- Two new execution paths:
  - `run_direct_ed_surface()` — direct EmptyDrops with custom ambient
  - `run_simpleed_custom_ambient()` — SimpleED with configurable ambient source
- These were used for the MSK ambient swap experiments and will be the
  entry point for the UCSF rerun.

**`core/legacy/source/SoloFeature_emptyDrops_CR.cpp`**
- Debug audit file output at end of `emptyDrops_CR()`:
  - `EmptyDrops/backend_debug_summary.json`
  - `EmptyDrops/backend_debug_candidates.tsv`
  - `EmptyDrops/backend_debug_ambient.tsv`

**`core/legacy/source/SoloFeature_emptyDrops_libscrna.cpp`**
- Matching debug audit file output for the `libscrna` backend.

## UCSF Production Run Root

Authoritative STAR outputs to reprocess:

- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009`

Per-sample raw matrix location (example for EBs1_1):

- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/samples/EBs1_1/run/Solo.out/GeneFull/raw/`

Per-sample filtered barcode location:

- `.../Solo.out/GeneFull/filtered/barcodes.tsv`

## Sample Set

All 16 samples:

- `EBs1_1`, `EBs1_2`, `EBs1_3`, `EBs1_4`, `EBs1_5`
- `EBs2_1`, `EBs2_2`, `EBs2_3`, `EBs2_4`, `EBs2_5`
- `iPSC1_1`, `iPSC1_2`, `iPSC1_3`
- `iPSC2_1`, `iPSC2_2`, `iPSC2_3`

## Production STAR Parameters (for reference)

From the production `RUN_COMMAND.sh`:

```
--soloCellFilter EmptyDrops_CR
--soloMultiMappers Unique
--soloCrMultimapRescue yes
--soloStrand Forward
--soloFeatures GeneFull Velocyto
--soloCBwhitelist .../3M-february-2018_TRU.txt
```

## Recommended Procedure

### Step 1: Commit and merge the emptyDrops changes

```bash
git checkout feature-emptyDrops
git add core/features/libscrna/src/OrdMagStage.cpp \
        core/features/libscrna/tools/scrna_simpleed.cpp \
        core/legacy/source/SoloFeature_emptyDrops_CR.cpp \
        core/legacy/source/SoloFeature_emptyDrops_libscrna.cpp
git commit -m "EmptyDrops: guarded ambient minimum, debug audit, and custom ambient CLI"
git checkout master
git merge feature-emptyDrops
```

### Step 2: Rebuild STAR

```bash
make -C core/legacy/source -j$(nproc) STAR
```

Verify the binary picks up the OrdMag guard change:

```bash
# Quick check: grep for the new constant in the binary
strings core/legacy/source/STAR | grep -i "guarded" || echo "check OrdMagStage.cpp compiled in"
```

### Step 3: Run EmptyDrops refiltering on each sample

For each sample, the existing raw matrix does NOT need to be regenerated. The
STAR alignment, UMI deduplication, and raw count matrix are unchanged. Only the
cell-calling step needs to rerun.

The recommended approach is to use the `scrna_simpleed` tool with the
`--use-guarded-rank-ambient` flag against each sample's existing raw MEX:

```bash
PROD_ROOT="/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009"

for sample in EBs1_{1..5} EBs2_{1..5} iPSC1_{1..3} iPSC2_{1..3}; do
    RAW="${PROD_ROOT}/samples/${sample}/run/Solo.out/GeneFull/raw"
    OUTDIR="${PROD_ROOT}/samples/${sample}/emptydrops_guarded_rerun"
    mkdir -p "${OUTDIR}"

    core/features/libscrna/tools/scrna_simpleed \
        --matrix "${RAW}/matrix.mtx" \
        --barcodes "${RAW}/barcodes.tsv" \
        --out "${OUTDIR}/barcodes.tsv" \
        --out-dir "${OUTDIR}" \
        --mode full \
        --use-guarded-rank-ambient \
        --ambient-fallback-min-abs 5000 \
        --ambient-fallback-min-frac 0.02 \
        --expected-cells 3000 \
        --sim-n 10000 \
        --fdr 0.01 \
        --use-fdr-gate \
        2>&1 | tee "${OUTDIR}/rerun.log"
done
```

**Important**: Run samples sequentially, not in parallel, to avoid resource
contention that could affect reproducibility.

### Step 4: Compare old vs new filtered barcodes

For each sample, compute:

1. Old filtered cell count (from existing `filtered/barcodes.tsv`)
2. New filtered cell count (from `emptydrops_guarded_rerun/barcodes.tsv`)
3. Jaccard overlap between old and new
4. If CR9 reference is available, Jaccard vs CR9

Collect into a summary table.

### Step 5: Update downstream if cell sets changed materially

If the guarded ambient rerun produces materially different cell sets, the
downstream pipeline will need to be re-run for affected samples:

- `scripts/run_scrna_downstream_gene_full_velocyto.sh`
- CellBender denoising (remote GPU)
- Adaptive QC refiltering

See `docs/HANDOFF_UCSF_ADAPTIVE_DOWNSTREAM_REPAIR_20260331.md` for the
downstream repair procedure.

### Step 6: Record results

Write a summary document at:

- `docs/UCSF_EMPTYDROPS_GUARDED_RERUN_RESULTS_20260402.md`

Include: per-sample cell count comparison table, Jaccard metrics, any samples
where the change was negligible vs significant, and whether downstream re-run
is needed.

## Expected Outcome

Based on the MSK experiments, the guarded ambient policy should:

- Preserve the legacy rank-window ambient selection when the window is large
  enough (which it will be for these 14K+ cell UCSF samples)
- Produce cell sets closer to CR9 parity
- Not change OrdMag-only results on smaller datasets

The UCSF samples have 14K–20K cells each, well above the `indMin=45,000`
threshold where the ambient window is meaningful. The fix should have a real
effect here.

## Key Files

| File | Purpose |
| --- | --- |
| `core/features/libscrna/src/OrdMagStage.cpp` | Guarded `minAmbientCells` fix |
| `core/features/libscrna/tools/scrna_simpleed.cpp` | CLI tool with `--use-guarded-rank-ambient` |
| `core/legacy/source/SoloFeature_emptyDrops_CR.cpp` | Legacy backend debug audit |
| `core/legacy/source/SoloFeature_emptyDrops_libscrna.cpp` | libscrna backend debug audit |
| `core/legacy/source/ParametersSolo.cpp` | Flag parsing fix (already merged) |
| `docs/MSK_AMBIENT_RESCUE_ISOLATION_20260402.md` | MSK experiments proving the ambient pool is the root cause |
| `docs/HANDOFF_UCSF_ADAPTIVE_DOWNSTREAM_REPAIR_20260331.md` | Downstream repair procedure if cell sets change |

## Caveats

1. The `scrna_simpleed` tool parameters above (`--expected-cells 3000`,
   `--sim-n 10000`, `--fdr 0.01`) are starting-point defaults. Cross-check
   against the production `RUN_COMMAND.sh` for each sample to ensure consistency.

2. The production run used `--soloCellFilter EmptyDrops_CR` which triggers the
   integrated STAR EmptyDrops path. The tool-based rerun isolates just the
   cell-calling step. Verify that the simple-cell threshold from the tool
   matches what STAR originally computed by checking
   `Solo.out/GeneFull/UMIperCellSorted.txt`.

3. Some UCSF samples still have old fixed-threshold downstream QC (not adaptive).
   See the sample split in `docs/HANDOFF_UCSF_DATA_ORGANIZATION_AND_RELEASE_NOTES_20260331.md`.
   The EmptyDrops rerun should happen BEFORE downstream repair so that the
   adaptive QC step uses the corrected cell set.

4. The debug audit files (`backend_debug_*.tsv`) are generated by the in-STAR
   path only. The standalone tool does not emit them. If you need per-candidate
   p-value audit for the UCSF rerun, capture them from the tool's stdout/log.
