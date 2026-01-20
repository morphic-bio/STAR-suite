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
