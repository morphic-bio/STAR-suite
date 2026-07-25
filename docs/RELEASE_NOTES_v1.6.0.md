# STAR Suite v1.6.0 Release Notes

Date: 2026-07-25

`v1.6.0` adds an opt-in integrated Visium HD 3-prime GEX path and fixes four
Flex alignment edge cases found during full-slide Visium validation. Release
artifacts use `v1.6.0`, Debian packages use `1.6.0-1`, and `STAR --version`
plus the molecule-first companion tools report `1.6.0`. The upstream STAR base
remains `2.7.11b`; genome-index compatibility remains `2.7.4a`.

## Integrated Visium HD GEX

- `--soloSpatialGexIntegrated yes` joins direct R1 spatial candidates to gene
  evidence after modern GeneFull and CR-compatible multimap resolution, before
  barcode correction and UMI collapse.
- Spatial barcode windows containing N use the bounded edit-distance fallback
  and normal uniqueness/ambiguity rules. UMIs containing N remain invalid and
  are retained in accounting without creating molecules.
- Strict, soft expected-count, hard, and gated-hard products are materialized
  atomically at 2, 8, and 16 micrometer resolutions.
- Bounded accumulation and downstream contribution/matrix spools preserve the
  accepted MEX bytes while avoiding the prior full-slide in-memory boundary.
- Annotated rescue remains the modern default. Explicit compatibility evidence
  is available for controlled comparisons; legacy GX/UR-derived gene policy is
  not used as a substitute for the integrated modern path.
- Diagnostic feature sidecars are optional and are intended as an overflow or
  debugging aid rather than the production data path.

## Flex alignment fixes

- Skip expected-cell validation only when Flex processing is explicitly
  disabled.
- Avoid legacy Solo replay in inline skip-processing mode.
- Route no-sample hash-cache misses to alignment instead of dropping them.
- Account for hash-kept reads that do not receive a barcode assignment.

These behaviors remain gated behind their existing Flex/spatial options. The
normal scRNA-seq and STARsolo paths are unchanged by default.

## Validation

- Clean Chromap-enabled STAR build and embedded-parameter regeneration.
- Focused spatial GEX, spill, sidecar, R1 decoder, multigene UMI, and CR rescue
  units; molecule-first native, frozen-reference, bounded-materializer, and
  production-adapter tests.
- Flex skip-processing parameter/inline tests, hash-to-alignment fallback,
  barcode-free keep bookkeeping, and checksum-pinned 800K decision replay.
- Standard Solo, CR-compatible GEX/CRISPR, Flex, Y/no-Y, SLAM, TranscriptVB,
  and normal-scRNA sidecar-off regression coverage.
- CRC Flex 100K deterministic gate: accepted open candidate surface, all eight
  resolver artifacts, and policy matrices reproduced exactly.
- Ovarian GEX 100K annotated and compatibility spill gates reproduced their
  accepted 36-component MEX and policy-summary hashes exactly.

Implementation and acceptance details are recorded in
`docs/RUNBOOK_VISIUM_HD_GEX_DOWNSTREAM_SPOOL_20260725.md` and
`docs/RUNBOOK_VISIUM_HD_GEX_IN_MEMORY_1MM_CR_20260724.md`.
