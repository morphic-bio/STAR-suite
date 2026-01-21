# Cell Ranger Multi Input (Smoke Pipeline)

This tool parses Cell Ranger `multi` configs and runs the feature-barcode pipeline to produce
10x-style MEX outputs suitable for early STAR-suite parity testing.

## Usage

```bash
python3 tools/cr_multi/cr_multi.py \
  --multi-config /path/to/multi_config.csv \
  --feature-ref /path/to/feature_reference.csv \
  --whitelist /path/to/whitelist.txt \
  --assign-barcodes core/features/feature_barcodes/assignBarcodes \
  --mex-stub tools/feature_barcodes/assignbarcodes_mex_stub.py \
  --gex-mex /path/to/gex_mex_dir \
  --gex-filter-type "Gene Expression" \
  --outdir /path/to/output
```

The script expects the `gex-mex` directory to contain `matrix.mtx(.gz)`,
`barcodes.tsv(.gz)`, and `features.tsv(.gz)`.

## Notes

- This is a smoke-path implementation for feature MEX integration.
- GEX processing via STAR is still separate; provide GEX MEX from STAR outputs when ready.
