# pf-multi table-backed feature arms

STAR pf-multi supports feature libraries whose assignments are already resolved
to per-barcode counts. Set `star_input_format=table` in the `[libraries]`
section of `--pfMultiConfig` and put the table path in the existing `fastqs`
column.

This mode is intended for cases where the feature signal is not encoded as a
FASTQ barcode read that `assignBarcodes` should search directly. Examples
include viral-state calls, lineage tables, or other per-cell feature assignments
produced by an upstream caller. The table-backed arm still runs inside pf-multi:
it uses the same per-library feature reference, barcode namespace validation,
dynamic FEATURE permit telemetry, per-library MEX output, provenance manifests,
and merge into `outs/raw_feature_bc_matrix` and
`outs/filtered_feature_bc_matrix`.

**HTO/CMO hash demultiplexing** is not a table-backed arm. When raw hashtag or
cell-multiplexing feature FASTQs are available, use the normal FASTQ-backed
`assignBarcodes --output-mode adt_mex` path with optional `star_hash_*` columns
(see `docs/RUNBOOK_NATIVE_HTO_CMO_FEATURE_DEMUX_20260615.md`). Do not reconstruct
hash assignments from a merged feature MEX when raw reads exist.

## Config

Minimal non-GEX table-backed library row:

```csv
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_library_id,star_feature_ref,star_whitelist,star_input_format
/path/to/counts.tsv,S1,Custom,Custom,TRU,viral_state,/path/to/viral_features.csv,/path/to/3M-february-2018_TRU.txt,table
```

Rules:

- `star_input_format` may be empty, `fastq`, or `table`; empty means `fastq`.
- `star_feature_ref` is required for table-backed non-GEX libraries.
- split-read feature-layout columns are incompatible with `star_input_format=table`.
- relative table paths are resolved against the pfMultiConfig directory.
- table-backed rows bypass FASTQ discovery and chemistry autodetect; use an
  explicit `star_chemistry`/whitelist that matches the table barcode namespace.

## Table Schema

Input tables are headered TSV or CSV files with these required columns:

```text
barcode feature_id count
```

The importer accepts `feature_id` or `featureid` as the feature column name.
Counts must be non-negative integers. Zero-count rows are skipped. Duplicate
`barcode,feature_id` pairs are collapsed by summing counts.

Barcodes are matched against the normalized assignment whitelist. GEM-well
suffixes such as `-1` are accepted and stripped when the whitelist is unsuffixed.
The table import telemetry reports:

- `barcode_namespace_input`: namespace shape observed in table rows
  (`SUFFIXED`, `UNSUFFIXED`, `MIXED`, or `UNKNOWN`)
- `barcode_namespace_output`: namespace shape emitted in the per-library MEX
- `table_rows_suffix_normalized`: rows whose suffix was normalized during match

## Permits And Provenance

When `--dynamicThreadInterface 1` is enabled, table import waits for the dynamic
permit interface before processing rows. It acquires and releases FEATURE-domain
permit chunks during import so table-backed work is visible to the same permit
system as FASTQ-backed feature assignment.

Per-library outputs include:

- `table_feature_import.api_run.txt`
- `pf_library_provenance.tsv`
- per-library `matrix.mtx`, `features.tsv`, and `barcodes.tsv`

Useful telemetry keys:

- `table_rows_read`
- `table_rows_retained`
- `table_rows_rejected_barcode`
- `table_rows_rejected_feature`
- `table_rows_rejected_count`
- `table_duplicate_pairs_collapsed`
- `table_feature_permit_acquires`
- `dynamicPermitDelta.feature.acquires`

## Validation

Focused tests:

```bash
bash tests/multi_feature/test_table_feature_import.sh
bash tests/multi_feature/test_multi_feature_config.sh
bash tests/multi_feature/test_table_gex_pf_multi_star_smoke.sh
```

The fast STAR smoke runs GEX plus a table-backed Custom arm and asserts nonzero
table counts in the per-library MEX, merged raw MEX, and merged filtered MEX.

The HIV DOGMA four-arm smoke exercises the proof-of-concept composition:

```text
GEX + ATAC + protein/ADT FASTQ + HIV table-backed Custom arm
```

Run:

```bash
bash tests/multi_feature/test_hiv_dogma_four_arm_table_smoke.sh
```

That smoke materializes `HIV_DNA`/`HIV_RNA` table rows from the local DOGMA
assignment files and asserts nonzero HIV counts in the per-library, raw merged,
and filtered merged MEX outputs.
