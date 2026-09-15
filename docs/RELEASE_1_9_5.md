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
