# Feature Barcode Tools (process_features vendor)

This module vendors the `process_features` toolchain into STAR-suite for early perturb-seq testing.
It builds standalone binaries and does not modify STAR alignment output.

## Build

```bash
make feature-barcodes-tools
```

Binaries are produced in `core/features/feature_barcodes/`:
- `assignBarcodes`
- `demux_fastq`
- `demux_bam`

## Quick Usage (assignBarcodes)

```bash
core/features/feature_barcodes/assignBarcodes \
  --whitelist /path/to/whitelist.txt \
  --featurelist /path/to/feature_ref.csv \
  --directory /path/to/output \
  /path/to/fastq_dir
```

## MEX Stub Conversion

`assignBarcodes` outputs `features.txt` and `barcodes.txt`. To create 10x-style TSVs:

```bash
tools/feature_barcodes/assignbarcodes_mex_stub.py \
  --assign-out /path/to/output \
  --feature-csv /path/to/feature_ref.csv
```

This writes `features.tsv` and `barcodes.tsv` alongside the existing `matrix.mtx`.

## Notes
- This is a tool-level vendor intended for early parity testing.
- STAR output integration will follow once the MEX path is validated.

---

## CRISPR Feature Calling (CR-Compat Mode)

When running STAR with `--crMultiConfig` and CRISPR Guide Capture features, STAR automatically runs GMM-based feature calling after EmptyDrops filtering.

### Parameter: `--crMinUmi N`

**Default:** 10 (CR-compatible for CRISPR guides)

Controls the minimum UMI threshold for feature calling. Adjust based on assay type:

| Assay Type | Recommended Value | Rationale |
|------------|-------------------|-----------|
| **CRISPR Guide Capture** | 10 (default) | Guides have variable capture efficiency and noise |
| **Lineage Barcodes** | 2-3 | Stable features with minimal noise - lower threshold captures more signal |
| **FLEX Probes** | 10 | Similar characteristics to CRISPR guides |

### Example

```bash
# CRISPR Guide Capture (default)
STAR ... --crMultiConfig config.csv

# Lineage Barcodes (lower threshold)
STAR ... --crMultiConfig config.csv --crMinUmi 3
```

### Output

CR-compat mode produces `outs/crispr_analysis/`:
- `protospacer_calls_per_cell.csv` - Per-cell feature assignments
- `protospacer_calls_summary.csv` - Calling statistics
- `protospacer_umi_thresholds.csv` - GMM-derived UMI thresholds

See `tests/crispr_feature_calling_comparison_report.md` for validation details.
