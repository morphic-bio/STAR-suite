# Handoff: UCSF Production Run

Date: 2026-03-29  
Repo: `/mnt/pikachu/STAR-suite`  
Branch: `production/ucsf-run-20260329`

## Authoritative Live Run

Use this run root, not the earlier aborted launch attempts:

- run root: `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305`

Earlier launcher roots exist but are not authoritative:

- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_172809`
- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_172913`
- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_173242`

## Launch Files

- wrapper log: [`LAUNCH.log`](/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/LAUNCH.log)
- launch command: [`LAUNCH_COMMAND.txt`](/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/LAUNCH_COMMAND.txt)
- launcher script: [`launch_setsid.sh`](/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/launch_setsid.sh)
- wrapper command emitted by the paper wrapper: [`RUN_WRAPPER_COMMAND.sh`](/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/RUN_WRAPPER_COMMAND.sh)
- detached session pid: [`SESSION_PID.txt`](/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/SESSION_PID.txt)

## Requested Runtime Surface

The production run was launched with:

- `24` threads
- `STAR.release`
- `--soloFeatures GeneFull Velocyto`
- no `Gene`
- `--outSAMtype BAM Unsorted`
- Y-removal outputs enabled
- downstream enabled
- CellBender enabled

The first sample command confirms that exact surface:

- [`samples/EBs1_1/RUN_COMMAND.sh`](/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/samples/EBs1_1/RUN_COMMAND.sh)

## Important Ordering Change

The batch runner was changed before this launch so that downstream does **not** run sample-by-sample during the STAR pass.

Current behavior:

1. Run STAR for each requested sample.
2. Package Velocyto outputs and finalize per-sample STAR artifacts.
3. After the STAR loop completes, start downstream across the completed sample set.

That behavior comes from [`run_ucsf_perturb_yremove_batch.sh`](/mnt/pikachu/STAR-suite/scripts/run_ucsf_perturb_yremove_batch.sh).

## Initial State Observed

At launch:

- [`LAUNCH.log`](/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/LAUNCH.log) shows wrapper startup and `Starting sample EBs1_1`
- [`samples/EBs1_1/RUN_MANIFEST.txt`](/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/samples/EBs1_1/RUN_MANIFEST.txt) exists
- [`samples/EBs1_1/pf_multi_config.csv`](/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/samples/EBs1_1/pf_multi_config.csv) exists
- STAR temp/log files exist:
  - [`samples/EBs1_1/run/Log.out`](/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/samples/EBs1_1/run/Log.out)
  - [`samples/EBs1_1/run/Log.progress.out`](/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/samples/EBs1_1/run/Log.progress.out)

## Monitoring

Use these checks:

1. Wrapper log:
   `tail -n 50 /mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/LAUNCH.log`
2. Current sample files:
   `find /mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/samples -maxdepth 2 -type f | sort | tail`
3. STAR progress for the active sample:
   `tail -n 50 /mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/samples/<sample>/run/Log.progress.out`
4. Completed samples:
   `find /mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/samples -name RUN_COMPLETE.ok | wc -l`
5. Downstream should only begin after the STAR loop finishes:
   `find /mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_20260329_174305/samples -maxdepth 2 -type d -name 'downstream_*'`

## What To Verify When Following Up

1. STAR outputs are present for each sample:
   - `Aligned.out_Y.bam`
   - `Aligned.out_noY.bam`
   - `y_separated/*.fastq.gz`
2. Velocyto packaged outputs exist:
   - `outs/raw_velocyto_feature_bc_matrix`
   - `outs/filtered_velocyto_feature_bc_matrix`
   - `outs/velocyto_feature_bc_matrix_manifest.json`
3. Downstream outputs appear only after the STAR phase completes.
4. If CellBender falls back on small/degenerate matrices, confirm the wrapper records the fallback rather than failing the whole batch.

## If Intervention Is Needed

- Do not edit STAR code on this branch while the run is active.
- If the run must be restarted, preserve this run root for forensics and launch a new timestamped root.
- If a sample fails, capture:
  - `RUN_COMMAND.sh`
  - `RUN_MANIFEST.txt`
  - `run/Log.out`
  - `run/Log.progress.out`
  - any downstream summary or failure marker
