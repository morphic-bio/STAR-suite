# STAR Suite v1.9.5 Release Notes

Corrected publication, 2026-09-20: this release replaces the incomplete
2026-09-15 v1.9.5 publication at commit `c95c57d`. Its tag and binary assets
have been rebuilt to include integrated LARRY feature calling. If you
downloaded v1.9.5 before this correction, download the new assets and verify
them against the updated `SHA256SUMS.download` file; the old package filenames
are reused, but their contents and checksums differ.

STAR Suite 1.9.5 extends the 1.9.4 half-probe Flex work to Visium HD spatial
data. The fused spatial route uses the shared half-probe classifier without a
genome index or alignment fallback, while retaining candidate-specific
molecule resolution and the strict, soft-expected, hard, and gated-hard
matrices at 2, 8, and 16 micrometers. The older spatial alignment route is
retired; Chromium Flex legacy compatibility remains available.

An opt-in `star_feature_caller=dominant` setting in `--pfMultiConfig` now calls
an individual non-GEX feature library during the same STAR run. For MSK LARRY,
the default assignment rule is the production top-count rule:
`top_count > second_count`. An exact tie remains unassigned. Optional
`star_feature_call_min_umi` and `star_feature_call_min_ratio` columns support
stricter policies, including 2 UMIs and 2:1. Calls are written to
`outs/feature_analysis/<library_id>/feature_calls.csv`. CRISPR GMM calling is
unchanged, and the new caller is disabled unless requested.

The integrated branch also includes FASTQ mate-token pairing and portable
process-features build fixes. `STAR --version` reports `1.9.5`; upstream STAR
remains `2.7.11b`, and genome-index compatibility remains `2.7.4a` (legacy
`2.7.1a`). Existing indexes do not need rebuilding.

## Validation

- Spatial development fixtures passed ordinary gzip, native BGZF, and forced
  spill parity, malformed-pair rejection, probe ambiguity, and matrix checks.
  See `docs/RELEASE_1_9_5.md` for fixture details and known scope.
- In the MSK 30KO ES three-library SSD run, the integrated LARRY caller
  produced 32,898 GeneFull filtered cells in 1,801 seconds with 32 threads.
  The GeneFull MEX was byte-identical to the unmodified 1.9.4 control. After
  barcode-order alignment, the LARRY MEX counts were identical too.
- On 26,521 cells shared with the released ES H5AD, LARRY assignments matched
  for every cell: 26,017 identical calls and 504 jointly unassigned.
- The merged candidate passed a clean build, 33 configuration checks, the
  dominant-caller unit test, and a three-library integration fixture whose
  call CSV matched the standalone caller byte-for-byte.

These are development validation results. Existing paper benchmark numbers
retain their measured 1.9.4 identity; no Cell Ranger comparison was rerun for
this feature.
