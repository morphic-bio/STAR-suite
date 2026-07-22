# STAR Suite v1.0.1 Release Notes

Date: 2026-05-28

`v1.0.1` is a hotfix release for CR-compatible perturb-seq feature assignment
determinism. The release artifact version is `1.0.1`, and `STAR --version`
reports `1.0.1`. The upstream STAR base remains available through
`STAR --upstream-version` (`2.7.11b`), and the genome index compatibility
string remains available through `STAR --genome-compat-version` (`2.7.4a`).

## Hotfix Scope

- Serialize the first `process_features` feature-bootstrap window to one
  consumer thread before enabling multi-consumer feature assignment. This keeps
  the bootstrap-learned feature-offset/search parameters deterministic when
  `--dynamicThreadInterface 1` and auto-sized CR assignment consumers are used.
- Fix the single-thread feature matcher used by `--crAssignSearchThreads 1` so
  non-exact matches report the selected match position. This restores barcode
  correction parity with parallel feature search.

## Validation

- Clean rebuilt `process_features` and STAR from `master`.
- Re-ran the 16-thread CR-compatible CRISPR smoke on the SSD fixture with:
  `--dynamicThreadInterface 1 --crAssignConsumerThreads -1
  --crAssignSearchThreads 1`.
- Confirmed byte-identical SHA matches versus the single-consumer baseline for
  raw and filtered GeneFull MEX matrices and CRISPR call outputs.
