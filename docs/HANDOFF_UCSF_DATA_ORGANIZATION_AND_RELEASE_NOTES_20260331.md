# Handoff: UCSF Data Organization And Release Notes (2026-03-31)

## Status

The corrected UCSF production rerun is complete at:

- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009`

Completion markers:

- STAR/batch log: `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/STAR_ONLY_RELAUNCH_FIX_20260331_072949.log`
  - ends with `All requested samples complete`
- remote CellBender watcher log: `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/REMOTE_SCAN_RELAUNCH_FIX_20260331_072949.log`
  - ends with `All expected samples completed`
  - ends with `Watcher exiting`

All 16 samples now have:

- `RUN_COMPLETE.ok`
- `downstream_genefull_velocyto_cellbender/final_counts.h5ad`
- `downstream_genefull_velocyto_cellbender/cellbender/cellbender_counts.h5`

## What This Handoff Is For

The next agent should do two things:

1. Organize the UCSF production outputs into a clean, publication/release-ready
   structure with clear provenance.
2. Write UCSF release notes using the same overall format/style we used for the
   FLEX benchmark summary.

Reference FLEX format:

- `docs/HANDOFF_FLEX_PIPELINE_BENCHMARK_SUMMARY_20260323.md`

Relevant UCSF run docs already present:

- `docs/PAPER_UCSF_CORRECTED_PRODUCTION_WORKFLOW_20260329.md`
- `docs/HANDOFF_UCSF_PRODUCTION_RUN_20260329.md`
- `docs/HANDOFF_UCSF_ADAPTIVE_DOWNSTREAM_REPAIR_20260331.md`

## Authoritative Run Root

Use this run root as the authoritative source of truth:

- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009`

Earlier UCSF production roots exist but are superseded and should not be used
for release-note numbers:

- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_172809`
- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_172913`
- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_173242`
- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305`
- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_20260329_175052`

## Key Provenance / Surface

The final corrected run uses:

- `24` STAR threads
- `STAR.release`
- `--soloFeatures GeneFull Velocyto`
- `--outSAMtype BAM Unsorted`
- Y-removal enabled:
  - `--emitNoYBAM yes`
  - `--emitYNoYFastq yes`
- GEX/Solo whitelist fixed to `TRU`
- CR whitelist fixed to translated `NXT`
- remote `rsync` CellBender fanout on `10.159.4.53`

The whitelist fix is critical. The broken production attempt used the wrong
Solo whitelist surface and should not be summarized as the production result.

## Sample Set

Sample count: `16`

Samples:

- `EBs1_1`
- `EBs1_2`
- `EBs1_3`
- `EBs1_4`
- `EBs1_5`
- `EBs2_1`
- `EBs2_2`
- `EBs2_3`
- `EBs2_4`
- `EBs2_5`
- `iPSC1_1`
- `iPSC1_2`
- `iPSC1_3`
- `iPSC2_1`
- `iPSC2_2`
- `iPSC2_3`

## Important Remaining Data-Quality Caveat

The run finished, but the downstream QC surface is mixed because the adaptive
QC default was fixed partway through the run.

Already documented in:

- `docs/HANDOFF_UCSF_ADAPTIVE_DOWNSTREAM_REPAIR_20260331.md`

Current adaptive-QC artifact split:

Adaptive artifacts present:

- `EBs1_1` (repair pilot)
- `EBs2_5`
- `iPSC1_1`
- `iPSC1_2`
- `iPSC1_3`
- `iPSC2_1`
- `iPSC2_2`
- `iPSC2_3`

Likely still on the old fixed-threshold downstream surface and should be
treated as repair targets before final release-note cell-count summaries are
frozen:

- `EBs1_2`
- `EBs1_3`
- `EBs1_4`
- `EBs1_5`
- `EBs2_1`
- `EBs2_2`
- `EBs2_3`
- `EBs2_4`

This matters because the old downstream path used:

- `min_genes=200`
- `max_genes=2500`
- `mt_pct_cutoff=5`

and can massively under-call filtered cells for these UCSF data.

## Recommended Immediate Next Steps

### 1. Freeze provenance before reorganizing

Do not move or rename the authoritative run root until the final organization
plan is written down.

Capture or preserve:

- `LAUNCH_COMMAND.txt`
- `RUN_WRAPPER_COMMAND.sh`
- `STAR_ONLY_RELAUNCH_FIX_20260331_072949.log`
- `REMOTE_SCAN_RELAUNCH_FIX_20260331_072949.log`
- per-sample `RUN_COMMAND.sh`
- per-sample `RUN_MANIFEST.txt`

### 2. Repair/adjudicate downstream QC before publication numbers

Before writing final release-note summary counts, verify whether the non-adaptive
subset above has been repaired. If not, repair those samples first using the
flow in:

- `docs/HANDOFF_UCSF_ADAPTIVE_DOWNSTREAM_REPAIR_20260331.md`

Do not mix pre-repair and post-repair filtered-cell numbers in a single summary
table without labeling them clearly.

### 3. Organize the final UCSF deliverables

At minimum, the next agent should produce a clean inventory of:

- per-sample STAR root
- per-sample downstream root
- per-sample `final_counts.h5ad`
- per-sample CellBender outputs
- per-sample CRISPR outputs
- per-sample raw/filtered Velocyto packaged outputs
- Globus destination root used for BAM/YFASTQ transfer

### 4. Write UCSF release notes in FLEX-summary style

Use the same high-level format as:

- `docs/HANDOFF_FLEX_PIPELINE_BENCHMARK_SUMMARY_20260323.md`

Recommended UCSF release-note structure:

1. Title / date / branch / dataset / machine
2. Executive summary
3. Run surface and key corrections
4. Final sample inventory
5. Output inventory
6. Notable fixes that materially changed the result
7. Remaining caveats / follow-up items
8. Links to authoritative artifacts and handoffs

## UCSF Release Notes: Content Guidance

### Executive summary should mention

- corrected production rerun completed successfully
- GEX/Solo whitelist regression was fixed (`TRU` for GEX, `NXT` for CR)
- Velocyto regression was fixed earlier in code
- Y-removal outputs were enabled
- remote GPU CellBender fanout completed successfully

### Key corrections to call out explicitly

These should appear as a dedicated section, not buried in prose:

- Velocyto/BAM lifetime fix (`packedReadInfo` early clear regression)
- UCSF production wrapper whitelist fix
- downstream adaptive-QC default fix
- remote CellBender offload / watcher path

### Things to avoid in release notes

- Do not summarize superseded UCSF run roots as if they were final.
- Do not present old fixed-threshold filtered-cell counts as final unless they
  have been repaired or clearly labeled.
- Do not imply that local `Aligned.out_Y.bam` / `Aligned.out_noY.bam` absence
  means Y-removal was skipped; many large artifacts were transferred to Globus
  and cleaned locally after transfer.

## Concrete Files The Next Agent Should Inspect

Run-wide:

- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/LAUNCH.log`
- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/LAUNCH_COMMAND.txt`
- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/STAR_ONLY_RELAUNCH_FIX_20260331_072949.log`
- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/REMOTE_SCAN_RELAUNCH_FIX_20260331_072949.log`

Per-sample examples:

- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/samples/EBs1_1/RUN_COMMAND.sh`
- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/samples/EBs1_1/RUN_MANIFEST.txt`
- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/samples/EBs1_1/run/outs/crispr_analysis/protospacer_calls_summary.csv`
- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/samples/EBs1_1/downstream_genefull_velocyto_cellbender/final_counts.h5ad`

Adaptive-repair pilot:

- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/samples/EBs1_1/downstream_genefull_velocyto_cellbender/adaptive_qc_threshold.json`
- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009/samples/EBs1_1/downstream_genefull_velocyto_cellbender/adaptive_repair_20260331_085330.log`

## Suggested Deliverables For The Next Agent

1. A UCSF organization/inventory doc that states exactly where the final
   deliverables live, including any Globus-only large artifacts.
2. A UCSF release-notes doc in FLEX-summary format.
3. If needed, an addendum describing which samples were repaired from the old
   fixed-threshold QC surface and which were produced directly under adaptive QC.

## Final Warning

The production run itself is done. The remaining risk is not “run completion”
but “summary correctness.” The next agent should assume the run root is final,
but the filtered-cell release numbers may still require adaptive-QC cleanup for
the mid-run sample subset listed above.
