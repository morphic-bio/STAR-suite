# STAR Suite 1.9.5

Visium HD Flex now uses the shared half-probe classifier through fused FASTQ
input. The spatial path retains raw R1 candidate coordinates and UMIs, resolves
molecules per candidate, and emits strict, soft-expected, hard and gated-hard
matrices at 2, 8 and 16 micrometers. It loads no genome and performs no alignment.
The previous spatial alignment route is rejected, including `--flexLegacy yes`.
Chromium Flex legacy compatibility remains available.

The implementation uses the existing spatial decoder, resolver and bounded spill
engine. It adds deterministic lane/record identities and checks complete paired
transactions. Plain FASTQ, ordinary gzip and native BGZF share the same feature
policy; malformed or truncated pairs fail rather than producing a shorter slide.
See [the spatial processing contract](SPATIAL_FLEX_HALF_PROBE.md).

## Development validation

The candidate built from `0c8d01a2ad4f7fd3500bb70ff3885e793be32cae` passed:

- The ordinary Flex route suite: 27 checks, including matrix equality with the
  public 1.9.4 artifact and cell-caller scheduling checks.
- Spatial input-format and spill parity, transaction invariants, feature-axis
  validation, synthetic probe ambiguity and malformed-pair rejection tests.
- CRC 100K ordinary gzip, native BGZF and forced spill: all 36 matrix components
  identical between formats and memory modes.
- Spatial GEX 100K: all 36 components byte-identical to public 1.9.4.
- Both SPATCH 100K probe reference arms: spatial barcode-decoding accounting
  unchanged from the archived route.

Feature counts can change because half-probe assignment replaces alignment
fallback. In the development fixtures, hard molecules changed as follows:

| Fixture | Archived 1.9.4 spatial route | New half-probe route |
| --- | ---: | ---: |
| CRC, 2020-A | 79,074 | 77,595 |
| SPATCH, 2020-A | 90,045 | 88,969 |
| SPATCH, 2024-A | 88,950 | 89,522 |

These are correctness fixtures, not full-slide performance or biological
benchmarks. Full-slide paper results require the accepted public artifact and
separate sealed runs. Existing README benchmark numbers and bundled evidence
retain their original 1.9.4 identities.

## Integrated feature calling

The corrected 1.9.5 branch adds opt-in dominant calls for an individual
non-GEX library in a three-library `--pfMultiConfig` run. Setting
`star_feature_caller=dominant` on the LARRY row writes
`outs/feature_analysis/<library_id>/feature_calls.csv` during STAR finalization.
The default rule matches the existing MSK production integration:
`top_count > second_count`; tied top counts are not assigned. Optional
`star_feature_call_min_umi` and `star_feature_call_min_ratio` columns allow a
stricter policy, such as 2 UMIs and a 2:1 ratio. The CRISPR GMM caller and
runs without this setting are unchanged.

On the MSK 30KO ES three-library SSD run, the v1.9.4-based development build
produced 32,898 GeneFull filtered cells in 1,801 seconds with 32 threads,
without BAM or Velocyto. The filtered GeneFull MEX was byte-identical to the
unmodified v1.9.4 P02 control. The LARRY filtered MEX had identical counts
after aligning barcode order. Among 26,521 cells shared with the released ES
H5AD, all LARRY calls agreed: 26,017 identical assignments and 504 jointly
unassigned. These are development validation results, not a new paper benchmark.

The merged branch passed a clean build, 33 multi-feature configuration checks,
the dominant-caller unit test, and a three-library integration fixture whose
call CSV matched the standalone caller byte-for-byte.
