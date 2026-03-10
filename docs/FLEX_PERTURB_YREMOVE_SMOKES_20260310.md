# Flex And Perturb Y-Removal Smokes

Date: March 10, 2026

This note records the new end-to-end Y-removal smokes for single-cell modes:

- [`tests/run_flex_yremove_smoke.sh`](/mnt/pikachu/STAR-suite/tests/run_flex_yremove_smoke.sh)
- [`tests/run_perturb_yremove_smoke.sh`](/mnt/pikachu/STAR-suite/tests/run_perturb_yremove_smoke.sh)

## Purpose

These smokes compare integrated STAR Y-removal against the standalone
`remove_y_reads` splitter on small real fixtures:

- Flex: downsampled `SC2300771`
- Perturb: UCSF `100k` `pfMultiConfig` fixture (`iPSC2_1_AALG2`)

## Key Semantics Locked

Single-cell Y/noY FASTQ emission now behaves as:

- the aligned cDNA read stream is emitted in `y_separated/`
- the corresponding barcode read is emitted to the same `Y` or `noY` side
- emitted FASTQ headers strip the Illumina comment suffix after the read name

Because of that, the smoke comparisons are defined as:

- compare emitted cDNA `R2` and barcode `R1` FASTQs
- ignore header comments
- require exact record order
- require identical normalized read names, sequence, and quality

The standalone reference is produced by:

- running `remove_y_reads` on the same cDNA `R2` FASTQs
- aggregating per-input outputs in the original input order

## Additional Guardrails

Flex smoke also forces:

- `--soloFlexOutputPrefix <outdir>/per_sample`

This avoids the historical default:

- `soloFlexOutputPrefix -`

Without an explicit override, Flex writes per-sample outputs under a top-level
`-` directory, and this repo already tracks legacy files there.

## Validation Results

Passing artifacts:

- Flex:
  [`tests/flex_yremove_smoke_output_20260310_082732/SUMMARY.txt`](/mnt/pikachu/STAR-suite/tests/flex_yremove_smoke_output_20260310_082732/SUMMARY.txt)
- Perturb:
  [`tests/perturb_yremove_smoke_output_20260310_082723/SUMMARY.txt`](/mnt/pikachu/STAR-suite/tests/perturb_yremove_smoke_output_20260310_082723/SUMMARY.txt)

Observed counts:

- Flex:
  - Y BAM reads: `499`
  - emitted cDNA/barcode FASTQ parity vs `remove_y_reads`: pass
- Perturb:
  - Y BAM reads: `57`
  - emitted GEX cDNA/barcode FASTQ parity vs `remove_y_reads`: pass
  - `outs/crispr_analysis/` present: pass
