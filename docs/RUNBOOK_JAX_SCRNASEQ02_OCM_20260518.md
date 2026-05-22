# JAX scRNAseq02 OCM Runbook

Date: 2026-05-18
Status: production-ready after clean STAR rebuild and the focused OCM tests.
The STAR core OCM multi-compatibility writer is implemented with native
run-level and per-sample Velocyto MEX materialization.

## Dataset Inventory

Primary dataset:

```text
/mnt/pikachu/JAX_scRNAseq02
```

Metadata mirror:

```text
/mnt/pikachu/JAX_scRNAseq02_metadata
```

The dataset is 10x GEM-X Single Cell 3' v4 OCM from KOLF2.2J-derived lines.
The raw download is complete for six library sets: 24 FASTQs total, with two
lanes per library (`L007`, `L008`) and paired `R1`/`R2` files. Total raw size is
about 1.2 TB.

Use `/mnt/pikachu/JAX_scRNAseq02/readme.md` and the workbook's `Library
preparation` sheet for OCM sample design. The workbook `Sequence file` sheet
appears to have shifted `INPUT_LIBRARY_PREPARATION_ID` values in the rows
previewed, so do not derive library-to-FASTQ mapping from that column until it
is corrected. Derive each library set from the FASTQ stem prefix.

| Library set | FASTQ stem | Pool | Raw size | OCM mapping |
| --- | --- | --- | ---: | --- |
| `25E32-L3` | `25E32-L3_GT25-03394_ACCTCGAGCT-ATCGAACACA_S44` | Day 4 pool 1 | 201.4 GB | `OB1=GCM1-Day-4`, `OB2=GRHL1-Day-4`, `OB3=OVOL1-Day-4`, `OB4=WT-PrS-20pct-Day-4` |
| `25E32-L4` | `25E32-L4_GT25-03395_CGAAGTATAC-CTCCAAGTTC_S45` | Day 4 pool 2 | 219.0 GB | `OB1|OB2=EPAS1-Day-4`, `OB3|OB4=WT-PrS-3pct-Day-4` |
| `25E34-L3` | `25E34-L3_GT25-03396_GCACTGAGAA-TTCACGCATA_S40` | Day 5 pool 3 | 217.1 GB | `OB1=GCM1-Day-5`, `OB2=GRHL1-Day-5`, `OB3=OVOL1-Day-5`, `OB4=WT-PrS-20pct-Day-5` |
| `25E34-L4` | `25E34-L4_GT25-03397_GCTACAAAGC-AGGGCACGTG_S41` | Day 5 pool 4 | 189.4 GB | `OB1=ISL1-Day-5`, `OB2=EPAS1-Day-5`, `OB3=WT-PrS-3pct-Day-5`, `OB4=WT-ExM-Day-5` |
| `25E35-L3` | `25E35-L3_GT25-03398_CGCTGAAATC-GCAGACACCT_S42` | Day 6 pool 5 | 197.9 GB | `OB1=GCM1-Day-6`, `OB2=GRHL1-Day-6`, `OB3=OVOL1-Day-6`, `OB4=WT-PrS-20pct-Day-6` |
| `25E35-L4` | `25E35-L4_GT25-03399_GAGCAAGGGC-CCAAGTCAAT_S43` | Day 6 pool 6 | 207.6 GB | `OB1=ISL1-Day-6`, `OB2=EPAS1-Day-6`, `OB3=WT-PrS-3pct-Day-6`, `OB4=WT-ExM-Day-6` |

The Cell Ranger log bundle in `cellranger-logs/` is for `25E32-L3` only. It is
useful as a validation oracle: Cell Ranger 9.0.1 reported `Single Cell 3' v4
(polyA) OCM`, 1,317,343,894 read pairs, 15,088 filtered cells, and
`cells_per_tag.json` counts of `OB1=4114`, `OB2=3830`, `OB3=4010`, `OB4=3134`.
The bundled MRI archive also contains Cell Ranger multi output layout entries
for `outs/multi/` and `outs/per_sample_outs/<sample>/count/`, so we can test
STAR's compatibility writer against the expected directory and filename
structure, not just against counts.

Important caveat: the Cell Ranger oracle `config.csv` has `include-introns=false`,
so its count matrix is the exonic `Gene` surface. The STAR-suite production
surface remains `GeneFull Velocyto` to match the MSK/UCSF downstream path. Treat
the oracle as a multi-layout and OCM cell-assignment smoke target, not as a
strict GeneFull MEX count-parity target.

## Reference And Barcode Policy

Use the same GEX reference surface as MSK/UCSF:

```text
/storage/autoindex_110_44/bulk_index
```

Use the staged GEM-X 3' v4 OCM TRU whitelist:

```text
/storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_TRU.txt
```

Because this is KOLF2.2J-derived material, run the UCSF-style Y-removal path:

```text
--outSAMtype BAM Unsorted
--emitNoYBAM yes
--emitYNoYFastq yes
--emitYNoYFastqCompression gz
```

Also keep the UCSF/MSK expression surface:

```text
--soloFeatures GeneFull Velocyto
--soloCellFilter None
--soloCrGexFeature genefull
--soloCrMultimapRescue yes
```

The OCM production route is split-before-ED. STAR counts on the effective
`CB16+OCM_TAG8` barcode axis and the native OCM materializer then runs the
CR-compatible EmptyDrops implementation separately per OCM biological sample.

## Preferred Boundary

Run STAR and OCM materialization locally on `pikachu`. Move only post-MEX
downstream work to the GPU host.

Local work:

1. FASTQ preflight and manifest generation.
2. Per-library STAR GeneFull/Velocyto run with Y-removal.
3. Velocyto MEX materialization.
4. STAR core OCM multi-compatibility writer materializes raw/filtered GeneFull
   and raw/filtered Velocyto matrices into per-sample run directories,
   including per-sample CR-compatible EmptyDrops when no pool filtered MEX is
   present.
5. Globus transfer of generated large files, followed by local cleanup after
   Globus success.

Remote GPU work:

1. Stage each per-sample post-MEX run directory with `rsync`.
2. Run `scripts/run_scrna_downstream_gene_full_velocyto.sh` with adaptive QC and
   CellBender GPU enabled.
3. Copy the downstream h5ad/QC directory back to the local sample directory.
4. Remove the remote staging directory unless debugging requires `--keep-remote`.

This boundary avoids requiring STAR, genome indices, or OCM FASTQs on the GPU
instance.

## STAR Command Shape

Create one local run directory per library set:

```text
/mnt/pikachu/JAX_scRNAseq02_processed/<run_name>/samples/<library_set>/run/
```

For each library, pass comma-separated `R2` cDNA files first and comma-separated
`R1` barcode/UMI files second:

```bash
STAR_BIN=/mnt/pikachu/STAR-suite/core/legacy/source/STAR
GENOME_DIR=/storage/autoindex_110_44/bulk_index
WL=/storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_TRU.txt
CONFIG_CSV=/path/to/cellranger_multi_config.csv

"${STAR_BIN}" \
  --runThreadN 16 \
  --dynamicThreadInterface 1 \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${R2_FILES}" "${R1_FILES}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${RUN_DIR}/" \
  --outTmpDir "${SAMPLE_DIR}/tmp" \
  --outSAMtype BAM Unsorted \
  --emitNoYBAM yes \
  --emitYNoYFastq yes \
  --emitYNoYFastqCompression gz \
  --clipAdapterType CellRanger4 \
  --clip3pPolyG yes \
  --alignEndsType Local \
  --chimSegmentMin 1000000 \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 \
  --soloCBlen 16 \
  --soloUMIstart 17 \
  --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloCBwhitelist "${WL}" \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloInlineCBCorrection yes \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloMultiMappers Unique \
  --soloCellFilter None \
  --soloCbUbRequireTogether no \
  --soloStrand Forward \
  --soloFeatures GeneFull Velocyto \
  --soloCrGexFeature genefull \
  --soloCrMultimapRescue yes \
  --soloInlineHashMode no \
  --ocmMultiEnable auto \
  --ocmMultiConfig "${CONFIG_CSV}" \
  --ocmMultiBarcodeMode flex \
  --ocmMultiOutputCompat cellranger
```

Use `--soloInlineHashMode no` for this Y-removal + Velocyto MEX path. The
non-Flex direct hash bridge (`STAR_SOLO_NONFLEX_HASH_BRIDGE=1` with
`--soloInlineHashMode yes`) is only valid for the no-BAM benchmark surface and
does not write the legacy `Solo.out/GeneFull` / `Solo.out/Velocyto` MEX trees
that the downstream wrapper consumes.

Pin the binary before production. If source changes or branches are switched,
perform a clean rebuild before debugging any failure, but do not rebuild the
binary while another production run is using it.

Use the Cell Ranger multi config that matches the library being processed. For
the oracle smoke this is `/mnt/pikachu/JAX_scRNAseq02/cellranger-logs/config.csv`;
for the full production set use the per-library config generated from the JAX
metadata workbook so `sample_id` and `ocm_barcode_ids` match the pool design.

## OCM Multi-Compatibility Writer

STAR now materializes Velocyto MEX internally as part of the post-Solo output
path, producing `outs/raw_velocyto_feature_bc_matrix` and
`outs/filtered_velocyto_feature_bc_matrix` without the legacy
`prepare_velocyto_mex.py` helper.

The STAR core OCM multi-compatibility writer then splits the pool-level
matrices by OCM overhang and sample design. This lives with the STAR-suite
CR-compatible output surface, not as a Python production path. The generic OCM policy is in
`docs/RUNBOOK_SCRNA_OCM_CR_COMPAT.md`.

OCM barcode classification:

| Overhang | OCM ID |
| --- | --- |
| `GT` | `OB1` |
| `CA` | `OB2` |
| `TC` | `OB3` |
| `AG` | `OB4` |

The writer should preserve STAR-native outputs and also write a Cell
Ranger-multi-compatible structure:

```text
outs/multi/count/raw_feature_bc_matrix/
outs/multi/multiplexing_analysis/cells_per_tag.json
outs/per_sample_outs/<sample_id>/count/sample_raw_feature_bc_matrix/
outs/per_sample_outs/<sample_id>/count/sample_filtered_feature_bc_matrix/
outs/per_sample_outs/<sample_id>/count/sample_filtered_barcodes.csv
```

For STAR downstream compatibility, also expose or mirror the per-sample
matrices under the directory shape expected by
`scripts/run_scrna_downstream_gene_full_velocyto.sh`:

```text
samples/<ocm_sample>/
  run/
    outs/raw_feature_bc_matrix/
    outs/filtered_feature_bc_matrix/
    outs/raw_velocyto_feature_bc_matrix/
    outs/filtered_velocyto_feature_bc_matrix/
    outs/multiplexing_analysis/cells_per_tag.json
```

For `OB1|OB2` and `OB3|OB4` samples, union the matching OCM barcode sets before
subsetting matrices.

Native implementation details for the per-OCM-sample MEX materializer are
tracked in
`docs/RUNBOOK_SCRNA_OCM_MULTI_MEX_MATERIALIZER_IMPLEMENTATION_20260519.md`.

## Oracle Smoke Harness

Do not rebuild STAR while another production run is using the binary. The smoke
harness uses an existing binary and is safe to prepare while production mapping
or downstream handoff continues:

```bash
scripts/run_jax_scrnaseq02_ocm_oracle_smoke.sh \
  --out-root /mnt/pikachu/JAX_scRNAseq02_processed/ocm_oracle_smoke_<stamp>
```

By default this only stages a 2,000,000 read-pair downsample of the `25E32-L3`
oracle library and writes `RUN_STAR.sh`. It does not launch STAR.

When the production STAR binary is safe to update, run the smoke with a cleanly
rebuilt STAR binary that includes the native OCM materializer:

```bash
scripts/run_jax_scrnaseq02_ocm_oracle_smoke.sh \
  --out-root /mnt/pikachu/JAX_scRNAseq02_processed/ocm_oracle_smoke_<stamp> \
  --run-star \
  --validate
```

The smoke harness now passes these STAR flags and Velocyto low-memory
environment defaults:

```text
--ocmMultiEnable auto
--ocmMultiConfig /mnt/pikachu/JAX_scRNAseq02/cellranger-logs/config.csv
--ocmMultiBarcodeMode flex
STAR_VELOCYTO_LOW_MEM=1
STAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS=8192
STAR_VELOCYTO_UMI_RESERVE_CAP=32
STAR_SOLO_BINARY_SPOOL=1
MALLOC_ARENA_MAX=2
MALLOC_TRIM_THRESHOLD_=131072
```

Override the environment before invoking the harness if a different spill budget
is needed. Override or extend STAR flags after `--`.

The validator can also be run independently against an existing smoke run:

```bash
scripts/validate_jax_scrnaseq02_ocm_oracle.py \
  --star-run-dir /mnt/pikachu/JAX_scRNAseq02_processed/ocm_oracle_smoke_<stamp>/samples/25E32-L3/run \
  --oracle-dir /mnt/pikachu/JAX_scRNAseq02/cellranger-logs
```

The smoke validator checks:

- `outs/multi/count/raw_feature_bc_matrix/`;
- `outs/multi/multiplexing_analysis/cells_per_tag.json`;
- `outs/per_sample_outs/<sample>/count/sample_raw_feature_bc_matrix/`;
- `outs/per_sample_outs/<sample>/count/sample_filtered_feature_bc_matrix/`;
- `outs/per_sample_outs/<sample>/count/sample_filtered_barcodes.csv`;
- per-sample filtered barcode count consistency;
- OCM tag proportions against the full `25E32-L3` Cell Ranger oracle as a
  warning by default;
- Cell Ranger oracle barcode recall in the matching STAR OCM tags;
- optional oracle-precision overlap against the full oracle cells.

The validator intentionally does not compare matrix values against the Cell
Ranger MEX because the oracle is `Gene` while the planned STAR output is
`GeneFull`. For the same reason, full-depth STAR `GeneFull` outputs with native
per-sample OCM EmptyDrops may call additional cells relative to the exonic Cell
Ranger oracle. The default pass/fail check is therefore oracle recall, not
equality of tag proportions. Use `--strict-tag-proportion-delta` only for a
same-feature-surface comparison where extra STAR-called cells should be treated
as a failure. If exact Cell Ranger count parity is needed later, run a separate
STAR comparison surface with `--soloFeatures Gene` and keep that out of the
production downstream path.

If the 2M smoke passes and we still want a deeper compatibility check, reuse
the same harness for full-depth `25E32-L3` structure/cell-assignment validation
by staging symlinks to the complete FASTQs rather than copying another large
fixture:

```bash
scripts/run_jax_scrnaseq02_ocm_oracle_smoke.sh \
  --out-root /mnt/pikachu/JAX_scRNAseq02_processed/ocm_oracle_full_<stamp> \
  --full-fastqs \
  --run-star \
  --validate
```

## Remote GPU Downstream

Use the existing remote downstream runner once each per-sample post-MEX run
directory exists:

```bash
scripts/run_remote_scrna_downstream_rsync.sh \
  --sample-dir "${SAMPLE_DIR}" \
  --remote-host 10.159.4.53 \
  --remote-root /home/lhhung/jax_scrnaseq02_remote_downstream \
  --output-name downstream_genefull_velocyto_cellbender_remote \
  --run-cellbender \
  --cellbender-gpu \
  --adaptive-filter
```

The script stages only MEX-level inputs and local helper scripts to the GPU
host, runs CellBender/downstream h5ad construction there, rsyncs the output
directory back, and removes the remote staging directory by default.

## Globus Transfer And Cleanup

Use the same destination endpoint currently used for JAX large-file handoffs:

```text
61fb8b9a-9b52-456e-928c-30c0fb0140bf
```

Proposed destination root:

```text
/JAX_scRNAseq02_processed/large_files
```

Transfer generated large artifacts:

```text
run/Aligned.out_Y.bam
run/Aligned.out_noY.bam
run/y_separated/*.fastq.gz
optional STAR trim/QC files if enabled
```

If this JAX handoff should be self-contained, also include the source FASTQs
under a `raw/` subdirectory for each library set. Mark those entries as raw
inputs in the transfer inventory and never delete them automatically after
transfer. Automatic cleanup is only for generated BAM/YFASTQ artifacts.

After the Globus task succeeds, delete generated local BAM/YFASTQ artifacts and
the local STAR temp directory. Preserve:

```text
source FASTQs under /mnt/pikachu/JAX_scRNAseq02/raw
STAR logs
MEX outputs
OCM split matrices
downstream h5ad/QC outputs
Globus batch manifests, task IDs, and cleanup markers
```

## Validation Gates

Before a new production launch:

1. Generate a manifest from FASTQ stems and the OCM design; verify exactly six
   library sets and 22 OCM biological samples.
2. Preflight one `R1` per library against the May-2023 GEM-X TRU whitelist.
3. Dry-run STAR commands and confirm the reference, whitelist, `GeneFull
   Velocyto`, `--soloCellFilter None`, native OCM per-sample EmptyDrops, and
   Y-removal flags.
4. Prepare the 2M `25E32-L3` oracle smoke harness without launching STAR while
   another production mapping run is active from the same checkout.
5. Once the STAR binary is free to update, run the 2M `25E32-L3` smoke through
   STAR, native Velocyto packaging, the OCM multi writer, and the oracle
   validator.
6. Compare the `25E32-L3` smoke against the bundled Cell Ranger log oracle:
   filtered total near 15,088 and `cells_per_tag.json` counts close to
   `4114/3830/4010/3134` for `OB1/OB2/OB3/OB4` on the full run. For the 2M
   smoke, require the same structure and reasonable OCM tag proportions rather
   than exact full-depth cell counts.
7. Compare the Cell Ranger multi compatibility layout against
   `25E32-L3_Day4-pool-1.mri.tgz`, especially `outs/multi/count/`,
   `outs/multi/multiplexing_analysis/`, and
   `outs/per_sample_outs/<sample>/count/sample_*_feature_bc_matrix/`.
8. If the 2M smoke passes and more confidence is needed, run the full
   `25E32-L3` oracle set for structure, OCM tag proportions, and cell-overlap
   validation. Do not require GeneFull-vs-Cell Ranger MEX count parity because
   the oracle is exonic `Gene`.
9. Confirm generated large files are deleted only after successful Globus task
   completion.

Production launch should use a background wrapper with per-library completion
markers, per-sample remote downstream logs, and a Globus state table. The
wrapper should not remove source FASTQs.
