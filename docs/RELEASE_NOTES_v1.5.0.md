# STAR Suite v1.5.0 Release Notes

Date: 2026-07-17

`v1.5.0` adds an opt-in, native molecule-first barcode assignment reference
while retaining the v1.4.4 decoder and ambiguity resolver as the default
compatibility baseline. The release artifact version is `v1.5.0`, Debian
packages use `1.5.0-1`, and both `STAR --version` and
`molecule_first_resolver --version` report `1.5.0`. The upstream STAR base
remains `2.7.11b`; genome-index compatibility remains `2.7.4a`.

## Molecule-first assignment

- Preserve every finite sequence candidate through deterministic
  read-to-clique-to-UMI accounting.
- Combine per-read sequence/Phred log likelihoods with PCR-inclusive exact-read
  frequency as a prior applied once per read clique, after which posterior
  mass participates in candidate-specific UMI correction and occupancy.
- Emit strict, soft expected-count, ungated hard, and gated-hard products in a
  single invocation without allowing any policy to overwrite another.
- Enforce the invariant that ambiguous molecules may refine conservative broad
  support but cannot introduce a candidate absent from sequence evidence.
- Keep spatial, expression, cell-type, neighborhood, graph, and GPU priors out
  of the v1.5.0 reference model.

The release surface is the companion `molecule_first_resolver` executable. It
accepts the documented normalized candidate ledger and is invoked before cell
or spatial-bin calling. Omitting this explicit stage preserves existing STAR,
STARsolo, and Flex behavior.

## Distribution

- Build `molecule_first_resolver` through `make flex-tools` and the default
  suite tool target.
- Install it in Docker and Debian packages.
- Include it beside STAR in static binary tarballs and each compatibility
  bundle variant.
- Require matching executable versions during tarball, installer-bundle, and
  Debian install/uninstall validation.

## Validation

- Native unit, malformed-input, sanitizer, and installed-binary tests.
- Exact conformance on the tracked small ledger against the frozen Python
  reference.
- Visium HD 3-prime, Visium HD Flex, Chromium 3-prime, and Chromium Flex 100K
  conformance: all integer product residuals were zero; soft-mass residuals
  were zero for the sealed spatial summaries and below `6.1e-8` for Chromium.
- Clean STAR build plus STARsolo, Flex on/off, Flex configuration, and 100K
  bridge determinism regressions.
- Byte-identical default-off STARsolo and Flex products versus v1.4.4.
- Ubuntu 22.04 release-build and runtime-install validation.
- Pull request Tier A Docker CI.

The full commands, fixture hashes, artifact paths, exclusions, and numeric
results are recorded in
`docs/RUNBOOK_MOLECULE_FIRST_BARCODE_ASSIGNMENT_20260717.md` and
`docs/MOLECULE_FIRST_BARCODE_ASSIGNMENT_IMPLEMENTATION_20260717.md`.
