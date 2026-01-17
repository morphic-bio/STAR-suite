# Gold Standard Test Data

This directory contains baseline MEX output files for validating the STAR-Flex pipeline.

## Contents

- `raw/` - Raw MEX output (barcodes.tsv, features.tsv, matrix.mtx)
- `per_sample/` - Per-sample filtered MEX output for 4 sample tags (BC004, BC006, BC007, BC008)
  - Each contains `Gene/filtered/` with barcodes.tsv, features.tsv, matrix.mtx
  - `flexfilter_summary.tsv` - Summary statistics

## Source

Generated from 100K downsampled SC2300771 dataset using the reference STAR-Flex implementation.

## Usage

The `run_flex_multisample_test.sh` script compares test output against these files to verify correctness.

