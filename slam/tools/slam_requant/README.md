# slam_requant

Re-quantify SLAM outputs from a STAR dump without re-alignment.

## Build
```
make
```

## Usage (minimal)
```
./slam_requant \
  --dump dump.bin \
  --out out_prefix \
  --slamSnpMaskIn mask.bed.gz
```

## Key options
- `--dump <path>`: STAR dump created with `--slamDumpBinary`.
- `--out <prefix>`: output prefix for `SlamQuant.out` and diagnostics.
- `--slamSnpMaskIn <bed.gz>`: apply SNP mask during replay (optional).
- `--trim5p/--trim3p`: manual trims (optional).
- `--autoTrim variance`: compute trims from dump using variance method.
- `--trimScope first|per-file`: auto-trim scope (default: first).
- `--strandness none|sense|antisense`: drop opposite strand reads (default: none).
- `--slamWeightMode dump|alignments|uniform`: use stored weights or recompute (default: dump).
- `--slamWeightFile <path>`: keyed weight sidecar file (overrides `--slamWeightMode`).
- `--slamWeightMatch auto|order|key`: matching strategy for weight file (default: auto).
- `--slamQcReport <prefix>`: write QC JSON + HTML (optional).

## Outputs
- `<out>.SlamQuant.out`
- `<out>.SlamQuant.out.diagnostics`
- `<out>.SlamQuant.out.transitions.tsv`
- `<out>.SlamQuant.out.mismatches.tsv`
- `<out>.SlamQuant.out.mismatchdetails.tsv`
- QC JSON/HTML if `--slamQcReport` is provided
