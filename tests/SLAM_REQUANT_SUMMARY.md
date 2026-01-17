## Summary: External SLAM Re-Quant (Dump + slam_requant)

This summarizes the implementation of the external re-quant workflow described
in `plans/external_SLAM_plan.md`.

### STAR: Dump Writer
New parameters:
- `--slamDumpBinary <path>`: write binary dump for external re-quant
- `--slamDumpMaxReads <N>`: cap number of dumped reads (default 1,000,000)
- `--slamDumpWeights <path>`: write keyed weight sidecar for external re-quant
- `--slamDumpWeightsMode dump|vbGene`: choose weight source for sidecar

The dump stores:
- gene IDs + names
- chromosome names + starts (for mask mapping)
- errorRate + convRate (defaults for requant)
- per-read payload (positions, bases, qual, gene IDs, weight, intronic flag)

### Tool: `tools/slam_requant`
CLI:
```
./tools/slam_requant/slam_requant \
  --dump <dump.bin> \
  --out <prefix> \
  --slamSnpMaskIn <mask.bed.gz> \
  --trim5p N --trim3p N \
  --autoTrim variance \
  --trimScope first|per-file \
  --strandness none|sense|antisense \
  --slamWeightMode dump|alignments|uniform \
  --slamWeightFile <weights.bin> \
  --slamWeightMatch auto|order|key \
  --slamQcReport <prefix>
```

Outputs:
- `<prefix>SlamQuant.out`
- `<prefix>SlamQuant.out.diagnostics`
- `<prefix>SlamQuant.out.transitions.tsv`
- `<prefix>SlamQuant.out.mismatches.tsv`
- `<prefix>SlamQuant.out.mismatchdetails.tsv`
- QC JSON/HTML if `--slamQcReport` is set

### Parity Test
`tests/run_slam_requant_parity.sh`:
- runs STAR on the 100k fixture with `--slamDumpBinary`
- runs `slam_requant` on the dump
- compares STAR vs requant with `tests/slam/compare_star_outputs.py`

### Weight Sidecar Test
`tests/test_slam_requant_weights.sh`:
- runs STAR with `--slamDumpBinary` + `--slamDumpWeights`
- reverses the weight sidecar order to force key matching
- runs `slam_requant --slamWeightFile --slamWeightMatch key`
- compares STAR vs requant output
