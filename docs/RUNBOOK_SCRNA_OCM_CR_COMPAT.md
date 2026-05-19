# scRNA-seq OCM CR-Compatibility Runbook

Date: 2026-05-15

## Policy

OCM support is a small extension to the existing STAR-suite scRNA-seq
CR-compatible GEX path. It is not a new STARsolo-only workflow.

Default OCM processing must use the STAR-suite CR-compatible routines:

- STAR-suite `core/legacy/source/STAR`, after a clean rebuild when source
  changed or branches were switched.
- `--soloCellFilter EmptyDrops_CR`, which uses the libscrna backend by default.
- `GeneFull` GEX counting for Cell Ranger compatibility.
- the same multimapper, UMI, barcode, adapter, poly-G, and dynamic-thread flags
  used by the MSK/UCSF CR-compatible scRNA-seq surfaces.

Do not use vanilla STARsolo cell calling for OCM unless the user explicitly asks
for a comparison run. In particular, do not use these as the production OCM
path:

- `--soloCellFilter CellRanger2.2 ...`
- standalone `--runMode soloCellFiltering`
- ad-hoc STARsolo-only filtered MEX generation followed by OCM splitting

The OCM-specific work is only:

1. parse Cell Ranger multi-style OCM config/sample metadata;
2. use the correct GEM-X whitelist family;
3. demultiplex raw and filtered GEX matrices by OCM overhang;
4. write Cell Ranger-style per-sample outputs and downstream per-sample
   GeneFull/Velocyto mirrors.

The code-level implementation plan for item 4 lives in
`docs/RUNBOOK_SCRNA_OCM_MULTI_MEX_MATERIALIZER_IMPLEMENTATION_20260519.md`.

## Reference And Whitelist

Use the same reference surface as the MSK production set:

```text
/storage/autoindex_110_44/bulk_index
```

This is the MSK benchmark/production STAR index for GRCh38 2024-A-compatible
GEX processing. Do not use the ad-hoc STAR-Spatial index generated during the
earlier OCM investigation.

For GEM-X 3' v4 OCM GEX, use the staged May-2023 TRU whitelist:

```text
/storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_TRU.txt
```

Do not point STAR at the gzipped Cell Ranger whitelist directly. Use the staged
plain-text whitelist above.

## Public 10x OCM Fixture

The downloaded public 10x human OCM fixture is external to the repo:

```text
/mnt/pikachu/star-spatial/10x/ocm/human_20k_gemx_ocm/
```

Relevant inputs:

```text
config.csv
cells_per_tag.json
count_raw_feature_bc_matrix.h5
fastqs/20k_Human_Donor1-4_PBMC_3p_gem-x_multiplex_GEX_fastqs/
Donor1_official_filtered_barcodes.csv
Donor2_official_filtered_barcodes.csv
Donor3_official_filtered_barcodes.csv
Donor4_official_filtered_barcodes.csv
```

The official H5 and per-sample barcode CSVs are validation oracles. FASTQs and
generated outputs must stay external and untracked.

## Command Shape

The OCM CR-compatible run should follow the current STAR-suite GEX-only modern
surface, with only the whitelist and FASTQ selection changed for GEM-X OCM.

Example skeleton:

```bash
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR

export STAR_SOLO_NONFLEX_HASH_BRIDGE=1

core/legacy/source/STAR \
  --runThreadN 32 \
  --dynamicThreadInterface 1 \
  --genomeDir /storage/autoindex_110_44/bulk_index \
  --readFilesIn "${R2_FILES}" "${R1_FILES}" \
  --outFileNamePrefix "${OUTDIR}/" \
  --outTmpDir "${OUTDIR}/tmp" \
  --outSAMtype BAM Unsorted \
  --emitNoYBAM yes \
  --emitYNoYFastq yes \
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
  --soloCBwhitelist /storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_TRU.txt \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloMultiMappers Unique \
  --soloCellFilter EmptyDrops_CR \
  --soloCbUbRequireTogether no \
  --soloStrand Forward \
  --soloFeatures GeneFull Velocyto \
  --soloCrGexFeature genefull \
  --soloCrMultimapRescue yes \
  --soloInlineHashMode no \
  --ocmMultiEnable yes \
  --ocmMultiConfig /mnt/pikachu/JAX_scRNAseq02/cellranger-logs/config.csv \
  --ocmMultiOutputCompat cellranger
```

Notes:

- If using a source checkout, clean rebuild before crash or parity debugging.
- Use `GeneFull Velocyto` for production OCM runs so downstream h5ad generation
  receives both expression counts and raw/filtered Velocyto layers.
- Use `--soloInlineHashMode no` on the production BAM/Y-removal surface. The
  inline hash no-BAM surface is useful for benchmarks, but it is not the OCM
  production path.
- `--ocmMultiEnable yes` materializes OCM outputs natively after Solo completes.
  `--ocmMultiConfig` may be omitted only when `--pfMultiConfig` points at the
  same Cell Ranger multi config with `[samples]`.
- `--ocmMultiOutputCompat cellranger` is the only implemented layout mode.
- Do not add `--soloCellFilter CellRanger2.2`; that is the vanilla STARsolo
  path and is not the STAR-suite standard.
- Do not add feature-library flags for public OCM GEX-only fixture runs. OCM is
  encoded in the cell barcode sequence, not as a Guide Capture feature library.

## OCM Demultiplexing

After STAR-suite writes raw and filtered GEX matrices:

```text
Solo.out/GeneFull/raw/
Solo.out/GeneFull/filtered/
```

`--ocmMultiEnable yes` splits both matrices into per-sample outputs using the
corrected cell barcode. Strip any `-1` suffix only for classification; preserve
the output barcode namespace used by STAR-suite.

OCM assignment is based on bases 8-9 of the 16 bp cell barcode:

| Overhang | OCM ID |
| --- | --- |
| `GT` | `OB1` |
| `CA` | `OB2` |
| `TC` | `OB3` |
| `AG` | `OB4` |

Map `OB1`-`OB4` to sample IDs from the Cell Ranger multi config `[samples]`
section. Do not hardcode donor names.

Expected output layout:

```text
outs/
  multi/count/raw_feature_bc_matrix/
  multi/multiplexing_analysis/cells_per_tag.json
  per_sample_outs/<sample_id>/count/
    sample_raw_feature_bc_matrix/
    sample_filtered_feature_bc_matrix/
    sample_filtered_barcodes.csv
samples/<sample_id>/run/outs/
  raw_feature_bc_matrix/
  filtered_feature_bc_matrix/
  raw_velocyto_feature_bc_matrix/
  filtered_velocyto_feature_bc_matrix/
  multiplexing_analysis/cells_per_tag.json
```

Velocyto per-sample mirrors are also native. STAR writes run-level raw/filtered
Velocyto MEX through `VelocytoMexWriter`, then the OCM materializer maps
GeneFull sample columns to Velocyto columns by barcode key before subsetting
`spliced`, `unspliced`, `ambiguous`, and total matrices. Do not use
`prepare_velocyto_mex.py` for new OCM production runs.

## Validation Plan

1. Preflight FASTQs:
   - resolve all configured FASTQ IDs and lanes;
   - verify R1 first 16 bp matches the May-2023 GEM-X TRU whitelist family;
   - verify no February-2018 whitelist is accidentally selected.

2. Dry-run STAR command:
   - confirm `/storage/autoindex_110_44/bulk_index`;
   - confirm `/storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_TRU.txt`;
   - confirm `EmptyDrops_CR`, `GeneFull`, `soloCrGexFeature genefull`,
     `soloCrMultimapRescue yes`, and no `CellRanger2.2`.

3. Full FASTQ-to-MEX run:
   - write to a fresh external output directory;
   - treat wrapper-written completion summaries or successful STAR exit plus
     final logs as completion signals.

4. OCM split:
   - split raw and filtered matrices by OCM overhang through
     `--ocmMultiEnable yes`;
   - write `cells_per_tag.json` from the filtered split.
   - verify `raw_velocyto_feature_bc_matrix` and
     `filtered_velocyto_feature_bc_matrix` exist for each sample.

5. Compare against official 10x oracles:
   - total raw and filtered UMI counts;
   - filtered barcode counts per donor;
   - barcode overlap after normalizing optional `-1` suffixes;
   - confirm STAR-only barcodes are not enriched in the wrong OCM tag.

The earlier vanilla STARsolo comparison produced Donor4 as a strict subset of
the official calls. That result should be treated as a negative control showing
why the STAR-suite CR-compatible cell caller is required.

## Implementation Checklist

- Add an OCM config parser for Cell Ranger multi CSV sections.
- Add a FASTQ resolver that supports Illumina sample number between FASTQ ID and
  lane, for example `<fastq_id>_S5_L001_R1_001.fastq.gz`.
- Add a GEM-X OCM preflight using the staged May-2023 TRU whitelist.
- Add a command builder that reuses the STAR-suite CR-compatible GEX-only
  surface above.
- Add a native post-Solo OCM materializer for raw and filtered MEX.
- Add native per-sample Velocyto MEX mirrors.
- Add tests for config parsing, OCM overhang classification, FASTQ resolution,
  MEX splitting, path-safe sample IDs, and Velocyto barcode-order mapping.
- Add an external 10x fixture validation entry to `tests/ARTIFACTS.md` after
  the corrected CR-compatible run is complete.
