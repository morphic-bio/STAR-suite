# STAR Suite v1.0.2 Release Notes

Date: 2026-05-29

`v1.0.2` is a hotfix release for the non-Flex Solo direct inline-hash bridge
collapse path. The release artifact version is `1.0.2`, and `STAR --version`
reports `1.0.2`. The upstream STAR base remains available through
`STAR --upstream-version` (`2.7.11b`), and the genome index compatibility
string remains available through `STAR --genome-compat-version` (`2.7.4a`).

## Hotfix Scope

- Combine duplicate exact UMI rows for each `(cell barcode, gene, UMI)` before
  `1MM_CR` correction in the non-Flex direct bridge collapse path.
- Preserve overflow checking when duplicate counts are merged.
- Book a production-scale repeat-run regression for exact raw GeneFull MEX
  identity under the same FASTQ input, same binary, and same direct bridge
  production flags. This is needed because tiny smoke fixtures do not reliably
  expose worker-local multithreading reduction bugs.

## Known Follow-Up

This release fixes one real bridge-collapse contract violation. It does not
claim that every residual schedule-sensitive difference in large production
direct-bridge runs has been eliminated. The next correctness item is a
repeat-run determinism regression that compares raw and filtered MEX exactly
on a large enough perturb fixture.

## Validation

- Clean rebuilt STAR from the hotfix branch.
- Confirmed `STAR --version` reports `1.0.2`.
- Confirmed `star_feature_call -v` reports `1.0.2`.
- Prior A375 diagnostic runs showed the exact-UMI collapse fix reduced the
  raw GeneFull matrix delta on the direct bridge path, while exposing the need
  for the broader repeat-run regression above.
