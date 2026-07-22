# scRNA-seq OCM Multi-MEX Materializer Implementation Runbook

Date: 2026-05-19
Status: implemented in STAR-suite (`OcmMultiMaterialize`, `OcmMultiConfig`,
`VelocytoMexWriter`, `--ocmMultiEnable`, `--ocmMultiConfig`,
`--ocmMultiBarcodeMode`, `--ocmMultiOutputCompat`). Unit tests:
`tests/test_ocm_{config_parser,barcode_classifier,mex_materializer_tiny,sample_id_validation,velocyto_barcode_shuffle}.sh`.
Run the 2M oracle smoke with `scripts/run_jax_scrnaseq02_ocm_oracle_smoke.sh
--run-star --validate` (rebuild STAR first) to confirm validator `PASS`.
Production OCM runs should use `--ocmMultiBarcodeMode flex`; `posthoc` is kept
for historical materialization comparisons only.

## Goal

Add a native STAR-suite OCM materializer that takes a completed OCM GEX run and
writes Cell Ranger multi-compatible outputs with one raw and filtered MEX per
OCM biological sample.

This is not a new mapper. Mapping, CR-compatible UMI collapse, Y-removal, and
Velocyto counting remain on the existing STAR-suite path. In the current OCM
production mode, CR-compatible EmptyDrops is applied by the native OCM
materializer per biological sample after the raw OCM split:

```text
pool-level GeneFull raw/filtered MEX
  -> OCM tag classification from corrected CB16 or effective CB16+TAG8 barcode
  -> Cell Ranger multi-compatible outs/multi and outs/per_sample_outs
  -> STAR downstream per-sample MEX mirrors, including Velocyto
```

## Verification

After rebuilding STAR (`make -C core/legacy/source STAR`), run:

```bash
scripts/run_jax_scrnaseq02_ocm_oracle_smoke.sh \
  --out-root /mnt/pikachu/JAX_scRNAseq02_processed/ocm_oracle_smoke_<stamp> \
  --run-star \
  --validate
```

To materialize OCM outputs on an existing completed pool run without remapping:

```bash
export OCM_TEST_RUN_DIR=/path/to/samples/<library>/run
export OCM_TEST_CONFIG=/path/to/cellranger-logs/config.csv
core/legacy/source/ocm_multi_unit_tests materialize
```

## Non-Goals

- Do not use `--soloCellFilter CellRanger2.2` or standalone
  `--runMode soloCellFiltering`.
- Do not implement OCM as a Guide Capture or feature-barcode assignment path.
  OCM is encoded in the cell barcode sequence.
- Do not require Gene-level count parity against the Cell Ranger oracle for
  the production surface. The oracle uses `include-introns=false`; production
  STAR-suite uses `GeneFull`.
- Do not use the non-Flex direct hash bridge for this BAM/Y-removal path.
  `STAR_SOLO_NONFLEX_HASH_BRIDGE=1` with `--soloInlineHashMode yes` is the
  no-BAM benchmark surface and does not write the MEX trees consumed here.

## Implemented Native Flags

OCM materialization is controlled by explicit flags so normal `pfMultiConfig`
feature workflows are not changed by surprise:

```text
--ocmMultiEnable yes|no|auto
--ocmMultiConfig <Cell Ranger multi config.csv>
--ocmMultiBarcodeMode posthoc|flex
--ocmMultiOutputCompat cellranger
```

- `--ocmMultiEnable no` by default.
- `--ocmMultiBarcodeMode posthoc` is the default historical materializer
  behavior: correct/count on CB16 first, then split materialized matrices by
  OCM.
- `--ocmMultiBarcodeMode flex` derives an internal `CB16+OCM_TAG8` effective
  barcode before barcode correction/counting and strips TAG8 from
  Cell Ranger-compatible output labels.
- `--ocmMultiEnable yes` requires `--ocmMultiConfig` or a reusable
  `--pfMultiConfig`.
- `--ocmMultiEnable auto` is accepted, but still requires an explicit config
  path today; it does not auto-discover OCM metadata without a config.
- If `--ocmMultiConfig` is unset and `--pfMultiConfig` is set, allow
  `--ocmMultiEnable yes` to reuse that config path, but keep this fallback
  explicit in the log.
- `--ocmMultiOutputCompat cellranger` is the only implemented compatibility
  mode and writes `outs/multi`, `outs/per_sample_outs`, and downstream
  `samples/<sample_id>/run/outs` mirrors.

The writer should use the active STAR-suite GEX feature selection:

```text
--soloCrGexFeature genefull
```

For JAX scRNAseq02, production/smoke commands should keep:

```text
--soloFeatures GeneFull Velocyto
--soloCellFilter None
--soloInlineCBCorrection yes
--soloInlineHashMode no
--outSAMtype BAM Unsorted
--emitNoYBAM yes
--emitYNoYFastq yes
```

## Config Parsing

Extend the existing Cell Ranger multi config parser rather than adding an
ad-hoc CSV reader.

Add a sample entry model with at least:

```text
sample_id
ocm_barcode_ids
description
```

The JAX oracle config has:

```text
[samples]
sample_id,ocm_barcode_ids,description
GCM1-Day-4,OB1,iPSCs
GRHL1-Day-4,OB2,iPSCs
OVOL1-Day-4,OB3,iPSCs
WT-PrS-20pct-Day-4,OB4,iPSCs
```

The full JAX dataset also needs pipe-union samples, for example:

```text
EPAS1-Day-4,OB1|OB2,...
WT-PrS-3pct-Day-4,OB3|OB4,...
```

Validation rules:

- Every sample must have a non-empty `sample_id`.
- `sample_id` must be path-safe: alphanumeric plus `-`, `_`, and `%`, no `/`,
  no `\`, no `..`, and no leading `.` or `-`.
- Every `ocm_barcode_ids` token must be one of `OB1`, `OB2`, `OB3`, `OB4`.
- Duplicate `sample_id` is fatal.
- Duplicate tag assignment across samples is allowed only for intentional
  union samples if the config expresses that design; log it clearly.
- Preserve config sample order for deterministic output.

## OCM Classification

In `posthoc` mode, classify corrected 16 bp cell barcodes by bases 8-9 using
1-based indexing. In C++ this is `barcode.substr(7, 2)` after stripping an
optional `-1` suffix for classification only.

In `flex` mode, STAR appends an OCM TAG8 suffix to the raw CB16 before
correction/counting. The materializer classifies TAG8 first, then falls back to
the CB16 overhang for historical fixtures.

| Bases 8-9 | OCM ID | Internal TAG8 |
| --- | --- | --- |
| `GT` | `OB1` | `GTGTGTGT` |
| `CA` | `OB2` | `CACACACA` |
| `TC` | `OB3` | `TCTCTCTC` |
| `AG` | `OB4` | `AGAGAGAG` |

Important compatibility details:

- Preserve the original output barcode string, including `-1` if present; in
  `flex` mode strip the internal TAG8 from CR-compatible output barcodes.
- Use the stripped barcode only for overhang classification and joins.
- Unknown overhangs are excluded from per-tag outputs and counted in a summary.
- The Cell Ranger oracle `cells_per_tag.json` preserves `-1` suffixes.
- The Cell Ranger oracle `sample_filtered_barcodes.csv` has no header and rows
  are `GRCh38,<barcode>`.

## Input Matrices

Minimum GEX inputs:

```text
Solo.out/GeneFull/raw/barcodes.tsv
Solo.out/GeneFull/raw/features.tsv
Solo.out/GeneFull/raw/matrix.mtx
Solo.out/GeneFull/filtered/barcodes.tsv
Solo.out/GeneFull/filtered/features.tsv
Solo.out/GeneFull/filtered/matrix.mtx
```

Velocyto inputs when `Velocyto` is requested:

```text
Solo.out/Velocyto/raw/barcodes.tsv
Solo.out/Velocyto/raw/features.tsv
Solo.out/Velocyto/raw/spliced.mtx
Solo.out/Velocyto/raw/unspliced.mtx
Solo.out/Velocyto/raw/ambiguous.mtx
Solo.out/Velocyto/filtered/barcodes.tsv
Solo.out/Velocyto/filtered/features.tsv
Solo.out/Velocyto/filtered/spliced.mtx
Solo.out/Velocyto/filtered/unspliced.mtx
Solo.out/Velocyto/filtered/ambiguous.mtx
```

The first implementation can read these on-disk MEX files after Solo output.
That is acceptable because this materialization boundary is after mapping and
cell calling. A later optimization may use in-memory matrices, but it must
produce byte-stable output structure first.

Run-level Velocyto `outs/` materialization is now a STAR core post-Solo step
implemented by `VelocytoMexWriter`, not the legacy `prepare_velocyto_mex.py`
helper. The OCM materializer reuses that internal writer for per-sample Velocyto
mirrors. Per-sample Velocyto subsetting is barcode-key based, not positional:
GeneFull column indices are mapped to Velocyto column indices with normalized
barcode keys before subsetting `spliced`, `unspliced`, `ambiguous`, and total
matrices. This allows the Velocyto raw/filtered barcode order to differ from
the GeneFull barcode order without silently swapping counts between OCM tags.

## Output Layout

Write Cell Ranger multi-compatible GEX outputs:

```text
outs/
  multi/
    count/
      raw_feature_bc_matrix/
        matrix.mtx.gz
        barcodes.tsv.gz
        features.tsv.gz
    multiplexing_analysis/
      cells_per_tag.json
  per_sample_outs/
    <sample_id>/
      count/
        sample_raw_feature_bc_matrix/
          matrix.mtx.gz
          barcodes.tsv.gz
          features.tsv.gz
        sample_filtered_feature_bc_matrix/
          matrix.mtx.gz
          barcodes.tsv.gz
          features.tsv.gz
        sample_filtered_barcodes.csv
```

Also write STAR downstream mirrors so the existing downstream runner can treat
each OCM biological sample like a normal sample:

```text
samples/<sample_id>/
  run/
    outs/
      raw_feature_bc_matrix/
      filtered_feature_bc_matrix/
      raw_velocyto_feature_bc_matrix/
      filtered_velocyto_feature_bc_matrix/
      multiplexing_analysis/cells_per_tag.json
```

For union samples, subset by the union of all declared OCM IDs.

## Materialization Algorithm

In the historical/posthoc path, the pool run already has raw and filtered MEX.
In the current split-before-ED production path, STAR writes raw MEX only and the
OCM materializer runs per-sample CR-compatible EmptyDrops after streaming raw
per-sample matrices.

1. Read raw GeneFull MEX, and read filtered GeneFull MEX only if the pool
   filtered tree exists.
2. Validate raw and filtered feature axes match exactly.
3. Build `tag -> raw column indices` from raw barcodes.
4. Build `tag -> filtered column indices` from filtered barcodes.
5. Write `outs/multi/count/raw_feature_bc_matrix` as the full raw GeneFull
   pool matrix, not just OCM-classified columns.
6. Write `cells_per_tag.json` from filtered barcode groups:

   ```json
   {
     "OB1": ["AAACCCTGTAAGCGCG-1"],
     "OB2": ["AAACCCGCAACTAGAC-1"],
     "OB3": ["AAACCATTCACCTGGG-1"],
     "OB4": ["AAACCAAAGCATTGAT-1"]
   }
   ```

7. For each configured sample:
   - resolve its tag set, including pipe unions;
   - subset raw GeneFull columns by raw tag indices;
   - subset filtered GeneFull columns by filtered tag indices, or run
     per-sample EmptyDrops from the raw per-sample MEX in split-before-ED mode;
   - write per-sample raw and filtered MEX;
   - write `sample_filtered_barcodes.csv` as `GRCh38,<barcode>` rows.
8. If Velocyto is present:
   - validate Velocyto barcodes can be joined to the GeneFull barcode namespace;
   - map each GeneFull sample column to the corresponding Velocyto column by
     normalized barcode key;
   - subset `spliced`, `unspliced`, `ambiguous`, and total matrices by the
     mapped per-sample Velocyto columns;
   - write per-sample downstream Velocyto MEX mirrors.
9. Write a machine-readable summary:

   ```text
   outs/multi/multiplexing_analysis/ocm_materialization_summary.json
   ```

   Include per-tag raw/filtered counts, unknown overhang counts, per-sample
   raw/filtered counts, config path, GEX feature surface, and STAR version.

## C++ Integration Points

Preferred new files:

```text
core/legacy/source/OcmMultiConfig.h
core/legacy/source/OcmMultiConfig.cpp
core/legacy/source/OcmMultiMaterialize.h
core/legacy/source/OcmMultiMaterialize.cpp
```

Implementation should reuse existing primitives where possible:

- config parsing patterns from `PfMultiConfig`;
- MEX loading/subsetting logic from `CrMultiMerge`;
- MEX writing via `MexWriter` or the existing MEX utilities;
- STAR output path and logging conventions from `PfMultiProcess`.

Hook point:

- Run after Solo writes `Solo.out/GeneFull` and `Solo.out/Velocyto`.
- Run before remote downstream handoff.
- Treat materialization failure as fatal when `--ocmMultiEnable yes`.
- Treat materialization failure as a warning only if `--ocmMultiEnable auto`
  and no OCM config/sample section is present.

Avoid reusing the old feature-assignment path for OCM. OCM does not need
`assignBarcodes`, feature references, or guide calling.

## Tests

Add unit-level tests:

```text
tests/test_ocm_config_parser.sh
tests/test_ocm_barcode_classifier.sh
tests/test_ocm_mex_materializer_tiny.sh
tests/test_ocm_sample_id_validation.sh
tests/test_ocm_velocyto_barcode_shuffle.sh
```

Tiny fixture requirements:

- four barcodes, one per OCM tag;
- one union sample;
- raw matrix with all four barcodes;
- filtered matrix with a subset;
- expected `cells_per_tag.json`;
- expected `sample_filtered_barcodes.csv` with `GRCh38,<barcode>` rows.

The tiny test must compare:

- MEX dimensions;
- barcode order;
- feature axis identity;
- matrix column subsetting;
- gzip file presence for Cell Ranger-compatible outputs;
- native Velocyto run-level and per-sample layer outputs, including
  `spliced.mtx.gz`, `unspliced.mtx.gz`, and `ambiguous.mtx.gz`.
- order-independent Velocyto subsetting, using a shuffled Velocyto barcode
  fixture where old positional indexing would assign the wrong OCM count.
- path-safe sample id validation for production configs.

## Smoke Test

Use the existing 2M oracle harness:

```bash
scripts/run_jax_scrnaseq02_ocm_oracle_smoke.sh \
  --out-root /mnt/pikachu/JAX_scRNAseq02_processed/ocm_oracle_smoke_<stamp> \
  --run-star \
  --validate \
  -- \
  --ocmMultiEnable yes \
  --ocmMultiConfig /mnt/pikachu/JAX_scRNAseq02/cellranger-logs/config.csv \
  --ocmMultiBarcodeMode flex
```

Expected pre-implementation result:

- STAR succeeds.
- GeneFull and Velocyto MEX exist.
- Validator fails only on missing `outs/multi` and `outs/per_sample_outs`.

Expected post-implementation result:

- Validator status is `PASS`.
- `outs/multi/count/raw_feature_bc_matrix` exists.
- `outs/multi/multiplexing_analysis/cells_per_tag.json` exists.
- Every config sample has raw and filtered per-sample MEX.
- Every config sample has `sample_filtered_barcodes.csv`.

For the 2M downsample, do not require full-depth Cell Ranger cell counts.
Require structure and reasonable OCM tag proportions. Full-depth `25E32-L3`
can be run later with `--full-fastqs` if deeper confidence is needed.

## Production Gate

Full JAX scRNAseq02 production is allowed only after:

1. Tiny materializer tests pass.
2. The Velocyto barcode-order shuffle and sample-id validation tests pass.
3. The 2M `25E32-L3` oracle smoke passes the validator.
4. Y-removal side outputs are still produced.
5. Velocyto per-sample mirrors exist for downstream CellBender/h5ad work.
6. The implementation is documented in `tests/ARTIFACTS.md` with the smoke
   output location.
7. The runbook command in `docs/RUNBOOK_JAX_SCRNASEQ02_OCM_20260518.md`
   matches the implemented flags.

For production, keep the native writer on the critical path. Do not route new
OCM runs through Python-only materialization unless debugging requires a
side-by-side comparison surface. The native path should write per-sample
GeneFull and Velocyto mirrors directly under `samples/<sample_id>/run/outs/`
so the existing downstream h5ad/CellBender wrapper can treat every OCM
biological sample like a normal single-cell sample.
