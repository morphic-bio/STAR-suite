# Paper UCSF Corrected Production Workflow

This records the UCSF production workflow used for the corrected symlinked
dataset surface and the paper-facing scripts that drive it.

## Binary Pin

- Git tag: `v0.90.1`
- Workflow binary: `core/legacy/source/STAR.release`
- `STAR.release` was refreshed from the local `v0.90.1` source tree so the
  saved release path exposes the current CR-compat flags used by this workflow
- Repo state to reproduce: this checkout is at the `v0.90.1` tag commit
  (`7bb13eb39a27fd43026ef6bc1717e7ce07638466`)

## Canonical Input Surface

- Corrected UCSF root: `/mnt/pikachu/ucsf-perturb-seq-corrected`
- Per-sample layout:
  - `GEX/`
  - `guides/`
  - `pf_multi_config.csv`

The production wrapper defaults to the corrected root and follows symlinks when
discovering FASTQs.

## Canonical Production Scripts

- `scripts/paper/run_ucsf_corrected_production_workflow.sh`
  - paper-owned UCSF entrypoint
  - defaults to:
    - `STAR.release`
    - corrected UCSF root
    - `24` threads
    - downstream h5ad generation
    - CellBender enabled
- `scripts/run_ucsf_perturb_yremove_batch.sh`
  - per-sample STAR runner
  - emits Y/noY BAMs and Y/noY FASTQs
  - packages raw/filtered velocyto MEX
  - writes per-sample transfer helpers and manifests
- `scripts/run_scrna_downstream_gene_full_velocyto.sh`
  - builds `counts.h5ad`, filtered h5ads, feature-library h5ads
  - runs CellBender on raw-backed `unfiltered_counts.h5ad`
  - writes `final_counts.h5ad`
- `tests/run_ucsf_corrected_production_100k_smoke.sh`
  - 100K per-library smoke on corrected `EBs2_2`
  - validates the full path through Y-removal, packaged velocyto, downstream
    h5ads, and CellBender
  - accepts either a real `cellbender_counts.h5` or a
    `CELLBENDER_FAILED.txt` marker when the sparse smoke matrix is too small for
    CellBender prior estimation
- `tests/run_ucsf_velocyto_cbub_smoke.sh`
  - host-local 100K corrected `EBs2_2` smoke for the shared `packedReadInfo` lifetime
  - validates non-zero Velocyto plus CB/UB tag injection in unsorted noY BAM on
    the same production-style surface
  - uses protected UCSF inputs/staging and its outputs must not be redistributed

## STAR Surface

The UCSF production runner now uses the current benchmark parameter surface for
mapping and CR-assignment, with one production-specific adjustment:
`--soloInlineHashMode` is not used here because this workflow must emit BAMs and
run Velocyto in the same pass.

- `--clipAdapterType CellRanger4`
- `--clip3pPolyG yes`
- `--alignEndsType Local`
- `--chimSegmentMin 1000000`
- `--soloType CB_UMI_Simple`
- `--soloBarcodeReadLength 0`
- `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts`
- `--soloUMIfiltering MultiGeneUMI_CR`
- `--soloUMIdedup 1MM_CR`
- `--soloMultiMappers Unique`
- `--soloCellFilter EmptyDrops_CR`
- `--soloCbUbRequireTogether no`
- `--soloStrand Forward`
- `--soloFeatures GeneFull Velocyto`
- no `Gene`
- `--soloCrGexFeature genefull`
- `--soloCrMultimapRescue yes`
- `--crChemistry auto`
- `--crOutputChemistry TRU`
- `--crAssignFeatureOffset -1`
- `--crAssignConsumerThreads -1`
- `--crAssignSearchThreads 1`
- `--crAssignSkipQcOutputs 1`
- `--dynamicThreadInterface 1`
- `--dynamicThreadConstMapPermits <threads>`
- `--dynamicThreadTelemetry 1`

Y-removal is always enabled in this production path:

- `--outSAMtype BAM Unsorted`
- `--emitNoYBAM yes`
- `--emitYNoYFastq yes`
- `--emitYNoYFastqCompression gz`

The production runner also explicitly unsets `STAR_SOLO_NONFLEX_HASH_BRIDGE`
before invoking STAR so the saved workflow does not inherit benchmark-only
environment settings.

## Commands

100K smoke on corrected `EBs2_2`:

```bash
tests/run_ucsf_corrected_production_100k_smoke.sh
tests/run_ucsf_velocyto_cbub_smoke.sh
```

Full corrected UCSF production run:

```bash
scripts/paper/run_ucsf_corrected_production_workflow.sh \
  --all-samples \
  --out-root /mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_<stamp>
```

With Globus upload:

```bash
scripts/paper/run_ucsf_corrected_production_workflow.sh \
  --all-samples \
  --globus-src-endpoint <pikachu_collection_id> \
  --globus-dst-endpoint <destination_collection_id> \
  --globus-dst-root /UCSF-perturb/<run_name>
```

## Recorded Outputs

Per sample the workflow should produce:

- `run/Aligned.out_Y.bam`
- `run/Aligned.out_noY.bam`
- `run/y_separated/*.fastq.gz`
- `run/outs/raw_velocyto_feature_bc_matrix/`
- `run/outs/filtered_velocyto_feature_bc_matrix/`
- `downstream_genefull_velocyto_cellbender/counts.h5ad`
- `downstream_genefull_velocyto_cellbender/unfiltered_counts.h5ad`
- `downstream_genefull_velocyto_cellbender/filtered_counts.h5ad`
- `downstream_genefull_velocyto_cellbender/default_singlet_filtered_counts.h5ad`
- `downstream_genefull_velocyto_cellbender/final_counts.h5ad`
- `downstream_genefull_velocyto_cellbender/cellbender/cellbender_counts.h5`
  or `downstream_genefull_velocyto_cellbender/cellbender/CELLBENDER_FAILED.txt`
- `RUN_COMMAND.sh`
- `RUN_MANIFEST.txt`
- `RUN_GLOBUS_TRANSFER.sh`
- `globus_batch.tsv`

The top-level wrapper also writes `RUN_WRAPPER_COMMAND.sh` for paper audit and
exact reruns.
