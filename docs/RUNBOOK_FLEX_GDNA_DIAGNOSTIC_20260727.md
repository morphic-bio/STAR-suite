# Flex gDNA Diagnostic: Withdrawn

The optional diagnostic introduced on 2026-07-27 was removed from the
corrected v1.9.5.a release on 2026-09-27 at the repository owner's request.

The previous runbook explicitly recorded Cell Ranger source-derived
implementation and numerical test fixtures. That provenance conflicts with
the project's documentation-only clean-room boundary. The estimator, its
metadata loader, per-cell/gene diagnostic collection, reports and derived
tests have been deleted. No replacement has been implemented.

- Flex expression counting, UMI correction and cell calling are unchanged by
  this removal.
- Existing STAR cache formats and packed count decoding remain compatible.
  Retained region bits are only compatibility data, not an estimator.
- No new `gdna_metrics.json`, `flex_gdna_library.json` or
  `flex_gdna_summary.tsv` is produced.
- `--soloFlexGdna auto/no` are deprecated no-ops. Required `yes` and explicit
  `--soloFlexGdnaProbeSet` paths are rejected.
- Historical output directories and Git history are not silently rewritten.
  Old gDNA reports in reused directories are stale; use fresh output directories.
- Historical release notes and benchmark logs describe older binaries, not
  an available feature. Do not use them as an implementation recipe.

Regression coverage checks removal, retired-option handling, packed count
compatibility and existing Flex count/cache behavior. See
`tests/test_flex_gdna_removed.py` and
`core/legacy/test/test_flex_probe_region.cpp`.
