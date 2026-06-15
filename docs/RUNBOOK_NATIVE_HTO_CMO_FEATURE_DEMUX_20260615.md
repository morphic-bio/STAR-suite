# Native HTO/CMO Feature Demultiplexing

Date: 2026-06-15

## Goal

Implement native STAR-suite demultiplexing for hashtag / cell-multiplexing
feature libraries:

```text
raw feature FASTQs -> feature-count MEX -> hash assignments -> STAR-suite outputs
```

This belongs in STAR-suite because HTO/CMO assignment is part of the raw
feature-barcode processing contract. Multiomics Suite may package the resulting
`hash` modality, but it should not own the caller for production runs.

## Background

DOGMA-HIV currently has a Multiomics Suite post-MEX compatibility bridge:
existing `Antibody Capture` rows named `hashtag1` through `hashtag6` are split
out of a merged feature MEX and classified after alignment. That is useful when
the only available artifact is a processed matrix, but it is not the production
architecture.

Production HTO/CMO should be implemented like the native feature-library arms:

- feature reads are assigned by `assignBarcodes`;
- deduplicated counts are emitted as 10x-style MEX;
- demux calls are generated from the same counts and cell barcode namespace;
- pf-multi records provenance and permit telemetry;
- Multiomics Suite consumes STAR-suite outputs without reclassifying MEX rows.

This is different from OCM effective-barcode GEX demultiplexing. OCM changes
the GEX barcode axis before correction and EmptyDrops. HTO/CMO here is a
feature-library sample-tag assay: count hash features per cell, then classify
cells as singlet, doublet/multiplet, or negative.

## Scope

In scope:

- `assignBarcodes` support for `Multiplexing Capture`, `HTO`, and `CMO`
  feature libraries.
- pf-multi routing and provenance for hash feature arms.
- First-class STAR-suite outputs:
  - hash count MEX;
  - optional protein-without-hash MEX when ADT and hash rows share a feature
    reference;
  - `hash_demux_assignments.tsv`;
  - `singlet_barcodes.tsv`, `doublet_barcodes.tsv`, `negative_barcodes.tsv`;
  - `hash_demux_summary.json`;
  - demux parameters and feature provenance.
- Focused synthetic tests plus a DOGMA-HIV smoke using existing raw ADT FASTQs.

Out of scope:

- Exact Seurat `MULTIseqDemux` reproduction as the first implementation.
- CLR/DSB normalization or downstream MuData assembly.
- Re-running the HIV DNA/RNA assignment arms.
- Changing OCM effective-barcode GEX behavior.

## Input Contract

### Feature reference

Accept hash libraries with standard feature reference columns:

| column | requirement |
|---|---|
| `id` | MEX feature ID and assignment label |
| `name` | display name; defaults to `id` if absent in legacy paths |
| `sequence` | HTO/CMO tag sequence |
| `feature_type` | `Multiplexing Capture`, `HTO`, `CMO`, or compatible configured alias |
| `pattern` | optional offset/layout override, as in existing feature matching |

For DOGMA-HIV, the current feature reference contains ADT rows plus hashtag
rows in one ADT/HTO library. STAR-suite should support both layouts:

1. hash-only library, e.g. all rows `Multiplexing Capture`;
2. mixed ADT + hash feature reference, with hash rows selected by feature type
   or explicit demux selector.

### pfMultiConfig

Hash-only library row:

```csv
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_library_id,star_feature_ref,star_whitelist
/path/to/hto_fastqs,S1,Multiplexing Capture,Multiplexing Capture,TRU,hto_s1,/path/to/hto_feature_ref.csv,/path/to/whitelist.txt
```

Mixed ADT + HTO library row:

```csv
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_library_id,star_feature_ref,star_whitelist,star_hash_feature_selector
/path/to/adt_fastqs,S1,Antibody Capture,Protein,TRU,adt_s1,/path/to/adt_hto_feature_ref.csv,/path/to/whitelist.txt,feature_type:HTO
```

Recommended selector forms:

| selector | meaning |
|---|---|
| empty | if `feature_types` is hash-like, all rows are hash rows |
| `feature_type:Multiplexing Capture` | select rows by feature type |
| `feature_type:HTO` | select rows by feature type |
| `id_prefix:hashtag` | select DOGMA-style hashtag rows |
| `name_regex:<regex>` | explicit fallback for legacy references |

Selector support can be staged. Minimal DOGMA support only needs hash-like
feature types plus `id_prefix:hashtag`.

## STAR CLI Surface

Standalone `assignBarcodes` should gain hash demux flags while preserving the
existing ADT MEX mode:

```bash
core/features/process_features/assignBarcodes \
  --whitelist /path/to/cells_or_whitelist.txt \
  --featurelist /path/to/adt_hto_feature_ref.csv \
  --directory /path/to/out \
  --output-mode adt_mex \
  --hash-demux yes \
  --hash-feature-selector id_prefix:hashtag \
  --hash-demux-method ratio \
  --hash-min-total 3 \
  --hash-min-top 3 \
  --hash-min-ratio 2.0 \
  /path/to/feature_fastqs \
  -b 16 -u 12
```

Suggested new options:

| option | default | meaning |
|---|---:|---|
| `--hash-demux yes/no/auto` | `auto` | run demux when hash features are present |
| `--hash-feature-selector SPEC` | empty | select hash rows in a mixed feature reference |
| `--hash-demux-method ratio` | `ratio` | first implementation method |
| `--hash-min-total N` | `3` | minimum total hash UMIs for non-negative calls |
| `--hash-min-top N` | `3` | minimum top hash UMI count for non-negative calls |
| `--hash-min-ratio X` | `2.0` | singlet if top / second is at least X |

pf-multi should expose the same controls through optional `star_` columns:

| pf-multi column | maps to |
|---|---|
| `star_hash_demux` | `--hash-demux` |
| `star_hash_feature_selector` | `--hash-feature-selector` |
| `star_hash_demux_method` | `--hash-demux-method` |
| `star_hash_min_total` | `--hash-min-total` |
| `star_hash_min_top` | `--hash-min-top` |
| `star_hash_min_ratio` | `--hash-min-ratio` |

## Output Contract

For a hash-capable feature library under:

```text
cr_assign/<feature_types>/<star_library_id>/<sample>/
```

emit:

| path | description |
|---|---|
| `hash/barcodes.tsv.gz` | barcodes on the assignment namespace |
| `hash/features.tsv.gz` | selected hash features with feature type `Multiplexing Capture` |
| `hash/matrix.mtx.gz` | hash feature x barcode count matrix |
| `protein/barcodes.tsv.gz` | optional protein-only MEX for mixed ADT+HTO libraries |
| `protein/features.tsv.gz` | ADT rows with hash rows removed |
| `protein/matrix.mtx.gz` | protein feature x barcode count matrix |
| `hash_demux_assignments.tsv` | per-barcode demux table |
| `singlet_barcodes.tsv` | singlet barcodes in output namespace |
| `doublet_barcodes.tsv` | doublet/multiplet barcodes |
| `negative_barcodes.tsv` | low/no-hash barcodes |
| `hash_demux_summary.json` | counts, thresholds, method, feature IDs, provenance |
| `hash_demux_command.txt` | exact effective demux arguments |

`hash_demux_assignments.tsv` columns:

```text
barcode
hash_assignment
hash_classification
hash_total_umis
hash_top_feature
hash_top_count
hash_second_feature
hash_second_count
hash_top_ratio
```

Classification labels:

| label | rule |
|---|---|
| `negative` | total hash UMIs < `hash_min_total` or top count < `hash_min_top` |
| `singlet` | top / second >= `hash_min_ratio`; second == 0 and top > 0 is infinite |
| `doublet` | otherwise |

Use `doublet` for the file/class label. Downstream QC may display this as
`multiplet` if matching other feature-call vocabulary.

## Merge Semantics

For CR-compatible merged outputs:

- hash MEX rows should merge as `Multiplexing Capture`;
- protein rows should exclude hash rows when the library was mixed;
- raw and filtered merged feature matrices should keep the same barcode
  namespace rules as existing pf-multi arms;
- `pf_library_provenance.tsv` should include hash demux fields:
  - `hash_demux_enabled`;
  - `hash_demux_method`;
  - `hash_feature_selector`;
  - `hash_min_total`;
  - `hash_min_top`;
  - `hash_min_ratio`;
  - `n_hash_features`;
  - `n_hash_singlet`;
  - `n_hash_doublet`;
  - `n_hash_negative`.

The filtered merged MEX should remain governed by the main GEX/cell-calling
surface. Singlet-only filtering is a downstream option, not the default merged
STAR output.

## Implementation Plan

1. Refactor the ADT MEX writer so row selection can write multiple MEX outputs
   from one deduplicated feature x barcode count object.
2. Add a hash feature selector parser in `process_features`:
   - match by feature type;
   - match by ID prefix;
   - optionally match by regex in a later pass.
3. Add a small hash demux classifier over the selected hash count matrix.
   Keep it independent from Seurat/MULTIseq internals for the first version.
4. Extend `assignBarcodes` options and `pf_adt_mex_config` / related structs
   with hash demux settings.
5. Extend pf-multi config parsing and `PfMultiAssign::runAssignBarcodes` to
   pass hash settings per library.
6. Write output provenance and telemetry.
7. Add tests before broad refactoring:
   - synthetic mixed ADT+HTO fixture;
   - hash-only fixture;
   - pf-multi config parse test for `star_hash_*` columns;
   - pf-multi smoke asserting per-library hash outputs and merged
     `Multiplexing Capture` rows.
8. Run a DOGMA-HIV smoke from existing ADT FASTQs and compare:
   - 165 protein rows;
   - 6 hash rows;
   - singlet/doublet/negative counts close to the Multiomics post-MEX bridge
     when using the same ratio thresholds.

## Validation Gates

Focused tests:

```bash
make -C core/features/process_features assignBarcodes
bash core/features/process_features/tests/test_adt_mex.sh
bash core/features/process_features/tests/test_hash_demux_mex.sh
bash tests/multi_feature/test_hto_cmo_hash_demux_arm.sh
bash tests/multi_feature/test_multi_feature_config.sh
```

DOGMA-HIV smoke:

```bash
bash tests/multi_feature/test_hiv_dogma_hto_demux_smoke.sh
```

Expected DOGMA-HIV compatibility target with ratio method (`min_total=3`,
`min_top=3`, `min_ratio=2.0`) on the current YW8 STAR-filtered barcode surface:

```text
hash_features=6
protein_features=165
singlet=8386
doublet=633
negative=15
```

If native STAR-suite assignment uses a different barcode universe than the
post-MEX bridge, compare on the shared finalized barcode set and document the
difference.

## Documentation Updates Required

Update these docs when implementation lands:

- `docs/feature_barcodes.md`
- `docs/RUNBOOK_PROCESS_FEATURES_ADT_MEX.md`
- `docs/PF_MULTI_TABLE_BACKED_FEATURE_ARM.md` only to clarify that HTO/CMO
  should be FASTQ-backed when raw reads are available, not table-backed.
- Multiomics Suite DOGMA-HIV inventory with the native STAR-suite output path.

## Acceptance Criteria

- Hash demux works from raw feature FASTQs without post-MEX Python
  reconstruction.
- Mixed ADT+HTO libraries emit both protein-only and hash-only MEX outputs.
- Hash assignment tables are attached to the same barcode namespace as the
  feature MEX.
- pf-multi provenance records method, thresholds, selector, and classification
  counts.
- Multiomics Suite can build a DOGMA-HIV MuData by consuming STAR-suite native
  `hash` outputs and optional singlet barcode list.
- Existing ADT/protein, CRISPR, LARRY/custom, and table-backed pf-multi tests
  remain green.
