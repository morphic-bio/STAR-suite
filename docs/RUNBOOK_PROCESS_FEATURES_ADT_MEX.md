# process_features ADT / Protein MEX Output

STAR-suite owns raw ADT (antibody-derived tag) quantification for Multiomics Suite.
This path reuses the existing `assignBarcodes` / `process_features` feature-matching
engine and emits a standard 10x-style protein MEX directory.

## Scope

- Input: feature-barcode FASTQs + antibody feature reference CSV
- Output: deduplicated protein MEX (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`)
- Not in scope: CLR normalization, denoising, or Multiomics orchestration

Multiomics Suite consumes the output as `protein.mex_dir` in a four-factor manifest.
Normalization (`protein.normalization = clr`) is applied downstream in Multiomics Suite.

## Read layout (10x feature barcode)

| Read | Content |
|------|---------|
| R1 (`_R1_`) | cell barcode + UMI |
| R2 (`_R2_`) | antibody tag sequence |

Defaults match standard 10x Feature Barcode libraries: 16 bp CB + 12 bp UMI on R1,
feature tag on R2 at offset 0 unless a `pattern` column overrides offset detection.

Matching uses the same Hamming-tolerant `process_features` search as CRISPR/ADT
libraries in CR-multi mode.

## Feature reference CSV

Required columns:

| Column | Description |
|--------|-------------|
| `id` | Feature ID (MEX row 1) |
| `name` | Display name (MEX row 2) |
| `sequence` | Antibody tag sequence |
| `feature_type` | Should be `Antibody Capture` for protein/ADT |

Optional columns (`target_gene`, `clone`, `vendor`, `isotype_control`, …) are
preserved in the copied snapshot `feature_reference.csv` written beside the MEX.

## Command

```bash
core/features/process_features/assignBarcodes \
  --whitelist /path/to/gex_whitelist.txt \
  --featurelist /path/to/protein_feature_ref.csv \
  --directory /path/to/protein_out \
  --output-mode adt_mex \
  --skip_empty_drops \
  /path/to/adt_fastq_dir \
  -b 16 -u 12
```

`--adt-mex` is an alias for `--output-mode adt_mex`.

ADT mode defaults:

- Skips EmptyDrops (barcodes stay in the GEX whitelist namespace)
- Skips histogram/heatmap QC artifacts
- Writes gzipped 10x MEX plus provenance files

Use `--filtered_barcodes` with `--source_namespace` / `--target_namespace` when
restricting to a GEX-called cell set while preserving barcode namespace parity.

## Output contract

Directory (one sample subdirectory per input FASTQ folder) contains:

| File | Description |
|------|-------------|
| `barcodes.tsv.gz` | Cell barcodes (GEX namespace) |
| `features.tsv.gz` | `id`, `name`, `Antibody Capture` |
| `matrix.mtx.gz` | Feature × barcode UMI counts (deduplicated) |
| `feature_reference.csv` | Snapshot of input feature ref |
| `protein_quant_summary.json` | Mode, ref fingerprint, layout, counts |
| `protein_quant_command.txt` | Command provenance |

Legacy assignBarcodes artifacts (`barcodes.txt`, `matrix.mtx`, `stats.txt`, …)
are still written for debugging.

## Multiomics Suite manifest fields

| Manifest field | Source |
|----------------|--------|
| `protein.mex_dir` | Output sample directory with gz MEX |
| `protein.feature_ref` | `feature_reference.csv` in that directory |
| `protein.normalization` | Set to `clr` in Multiomics Suite (not computed here) |

## Test

```bash
core/features/process_features/tests/test_adt_mex.sh
```
