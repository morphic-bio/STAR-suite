# STAR Suite v1.4.4 Release Notes

Date: 2026-07-17

`v1.4.4` is a patch release for feature-barcode namespace remapping in
integrated AnnData/MuData exports. The release artifact version is `v1.4.4`,
Debian packages use `1.4.4-1`, and `STAR --version` reports `1.4.4`. The
upstream STAR base remains `2.7.11b`; genome-index compatibility remains
`2.7.4a`.

## Fixes

- Translate process_features top-feature indices back to the source feature
  namespace before annotating integrated feature-library outputs.
- Preserve blank calls for genuinely ambiguous multi-feature rows instead of
  assigning the first remapped feature.
- Add a regression fixture covering one-based source indices, two-feature
  ambiguous rows, and the resulting AnnData annotations.

## Workflow hardening

- Add a typed, deterministic UCSF per-sample Slurm workflow contract.
- Make the UCSF filtering combination stage deterministic and cover it with
  workflow-schema and execution-contract tests.

## Validation

- Clean STAR core build and version check.
- Feature namespace remapping regression.
- UCSF Slurm workflow contract tests.
- Release-script version-default audit.
