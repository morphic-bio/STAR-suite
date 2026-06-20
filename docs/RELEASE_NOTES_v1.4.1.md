# STAR Suite v1.4.1 Release Notes

Date: 2026-06-20

`v1.4.1` is a hotfix release for Flex/Solo build hygiene and the public Flex
tiny smoke wrapper. The release artifact version is `v1.4.1`, and
`STAR --version` reports `1.4.1`. The upstream STAR base remains available
through `STAR --upstream-version` (`2.7.11b`), and the genome index
compatibility string remains available through `STAR --genome-compat-version`
(`2.7.4a`).

## Fixes

- Include VPATH-built top-level Flex sources in the core dependency scan so
  shared header/layout changes rebuild the corresponding Flex/Solo objects.
- Preserve the wrapper-level `--out-samtype none` contract through
  `scripts/run_flex_cr_config.sh`, allowing no-BAM Flex smoke runs to validate
  that no alignment BAM is emitted.
- Add `bam-unsorted` validation for the tiny public Flex binary smoke wrapper.
- Make release validation pass the expected binary version derived from the tag,
  avoiding stale hard-coded `STAR --version` checks in tag-triggered releases.

## Validation

- Clean rebuilt STAR from `core/legacy/source` before reproducing the Flex case.
- Passed direct Flex no-BAM run using the cached tiny fixture.
- Passed `tests/run_flex_tiny_public_binary_smoke.sh --out-samtype none`.
- Passed `tests/run_flex_tiny_public_binary_smoke.sh --out-samtype bam-unsorted`.
