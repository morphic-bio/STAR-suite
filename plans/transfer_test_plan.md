# Transfer Test Plan (STAR-suite)

## Goals
- Validate that Bulk/core, Flex, and SLAM test scripts run from STAR-suite paths.
- Confirm parity metrics and key outputs match prior STAR-Flex baselines.
- Record any deviations and identify whether they are data-path or code issues.

## Preflight
- Build binaries/tools:
  - `make core` (core/legacy STAR binary)
  - `make flex` (Flex tools)
  - `make slam` (SLAM tools)
- Ensure data dependencies exist (examples):
  - `/storage/flex_filtered_reference/star_index`
  - `/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt`
  - `/storage/downsampled_100K/...`
  - `/storage/SLAM-Seq-prod-compare-20260109/input/...`
  - `/storage/autoindex_110_44/bulk_index`
- Optional tools/paths:
  - `GEDI_BIN` if GEDI comparisons are required.
- Output hygiene:
  - Use per-run work dirs under `tests/` or `/storage/` as defined in each script.

## Core/Bulk Tests
1. Smoke/basic
   - `tests/run_baseline_test.sh`
   - `tests/run_solo_smoke.sh`
2. Y-removal / BAM split
   - `tests/run_ychrom_fastq_se_test.sh`
   - `tests/run_ychrom_fastq_r3_ignore_test.sh`
   - `tests/run_ychrom_fastq_edge_cases.sh`
   - `tests/run_ychrom_bam_split_test.sh`
   - `tests/run_ychrom_bulk_pe_test.sh`
   - `tests/run_ychrom_bulk_pe_samtools_sort_test.sh`
   - `tests/run_ychrom_regression_test.sh`
3. Extended (long)
   - `tests/run_pe_preprod_validation.sh`

## Flex Tests
1. Smoke/basic
   - `tests/run_flex_inline_test.sh`
   - `tests/run_flex_with_tags_test.sh`
2. Validation/coverage
   - `tests/run_flex_cbub_validation_test.sh`
   - `tests/run_flex_multisample_test.sh`

## SLAM Tests
1. Unit/smoke
   - `tests/run_slam_unit_tests.sh`
   - `tests/run_slam_solver_test.sh`
   - `tests/test_slam_requant_weights.sh`
   - `tests/run_slam_requant_parity.sh`
2. Fixture parity
   - `tests/run_slam_fixture_parity.sh`
   - `tests/run_slam_end_to_end_fixture.sh`
3. Optional full E2E (GEDI)
   - `tests/run_slam_end_to_end.sh`
   - `tests/run_slam_end_to_end_vb_commonmask.sh`
   - `tests/run_slam_end_to_end_gedi_parity.sh`
4. SNP mask checks
   - `tests/test_binomial_snp_mask.sh`
   - `tests/run_star_snp_mask_sweep.sh`

## Validation Criteria
- Core/Flex: expected output files exist and match prior run summaries.
- SLAM: correlation thresholds and parity checks meet existing docs.
- Record any failures with log excerpts and input paths used.

## Notes
- Historical documentation files in `tests/` still reference STAR-Flex paths; this does not affect script execution.
- If any scripts reference external data not present, mark the test as blocked and document required inputs.
