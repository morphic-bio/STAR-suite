# Molecule-first barcode assignment: implementation and validation record

Date: 2026-07-17

Branch: `agent/molecule-first-barcode-assignment`

Base: `v1.4.4` (`20efad80bd71ec9004b73c2f003fa8f3decdf20e`)

This record accompanies
`RUNBOOK_MOLECULE_FIRST_BARCODE_ASSIGNMENT_20260717.md`. It documents the
implementation, the scripts used to validate it, and the immutable inputs and
outputs needed to reproduce the paper-facing conformance claims.

## Implementation surface

- `flex/source/solo/MoleculeFirstResolver.{h,cpp}`: deterministic C++ reference
  library for read cliques, Bayesian candidate posterior, candidate-specific
  UMI correction, expected occupancy, and integer policy products.
- `flex/tools/molecule_first_resolver/`: command-line release surface, native
  unit tests, and fail-closed CLI tests.
- `tests/molecule_first/compare_reference.py`: small-ledger comparison against
  the frozen Python model in `visium-hd-processing`.
- `tests/molecule_first/run_100k_conformance.py`: normalizes the four frozen
  100K ledgers, runs the native resolver, compares summaries, and records input
  hashes and numeric residuals.
- `tests/run_molecule_first_native_smoke.sh`: self-contained installed-binary
  regression included in Tier A container CI.
- Top-level, Docker, Debian, tarball, compatibility-bundle, and release
  install/runtime checks now build and require `molecule_first_resolver`.

The existing STAR, STARsolo, and Flex ambiguity behavior remains the default.
The new resolver is invoked explicitly on the normalized candidate ledger and
runs before a downstream cell or spatial-bin call.

## Native and reference tests

Commands:

```bash
make -C flex/tools/molecule_first_resolver clean test
bash tests/run_molecule_first_reference_conformance.sh
MOLECULE_FIRST_RESOLVER="$PWD/flex/tools/molecule_first_resolver/molecule_first_resolver" \
  tests/run_molecule_first_native_smoke.sh
```

All passed. The native unit binary was also compiled and run with
AddressSanitizer and UndefinedBehaviorSanitizer using `-Wall -Wextra -Werror`.
The tracked small ledger matched the frozen Python reference for clique
membership, candidate sets, posterior, UMI roots, and policy products.

## Frozen 100K conformance

Output root:

```text
/mnt/pikachu/star-spatial/runs/20260717_star_suite_molecule_first_integration_v1
```

Each successful directory contains the normalized ledger, native product
tables, and `conformance.json`. The JSON records source paths and SHA-256
digests, expected and native summaries, residuals, and pass status.

| Assay | Candidate reads | Read cliques | Integer residuals | Soft-mass residual |
| --- | ---: | ---: | --- | ---: |
| Visium HD 3-prime | 54,542 | 41,808 | all zero | 0 |
| Visium HD Flex | 89,872 | 78,494 | all zero | 0 |
| Chromium 3-prime | 96,535 | 87,840 | all zero | at most `6.04e-8` |
| Chromium Flex | 90,869 | 88,727 | all zero | at most `1.98e-8` |

The spatial runs compare against sealed summaries. The Chromium runs compare
against the live frozen Python reference at processing commit
`ae598fbf14051972c2f1ac7bdbdf705f288913a0`. Floating residuals are
summation-order effects and are below the declared tolerance.

Successful reports:

```text
visium_hd_3prime_100k/conformance.json
visium_hd_flex_100k_v2/conformance.json
chromium_3prime_100k/conformance.json
chromium_flex_100k/conformance.json
```

The earlier `visium_hd_flex_100k/` attempt used a training-split prior rather
than the frozen full human exact-read prior and correctly failed its gated
count comparison. It is not a successful artifact and is excluded from all
claims. `visium_hd_flex_100k_v2/` is the valid run.

## Default-off STAR regressions

The following were run after a clean legacy STAR rebuild:

```bash
bash tests/run_solo_smoke.sh
FLEX_FIXTURE_ROOT=/mnt/pikachu/STAR-suite/reference/tests/100K \
  bash tests/flex_smoke/run_flex_smoke.sh
OUT_BASE=/tmp/star_molecule_first_flex_config_smoke \
  bash tests/run_flex_cr_config_smoke.sh
READ_MAP_NUMBER=100000 THREADS_A=4 THREADS_B=4 \
  OUT_ROOT=/tmp/star_molecule_first_bridge_determinism \
  bash tests/run_solo_bridge_determinism_harness.sh --force
```

All applicable tests passed. Bridge source stages, post-collapse stages, and
raw/filtered MEX products were identical between the two repeated runs;
`first_different_source_stage=none` and `first_different_stage=none`.

The feature branch and v1.4.4 baseline raw STARsolo MEX files were
byte-identical:

| File | SHA-256 |
| --- | --- |
| `matrix.mtx` | `f51084929dd000e648a50b46cc8447b2d682d7f7a54d31568ae6beaf78387a24` |
| `barcodes.tsv` | `3357896364378946c3ab3f4e015413c517e02a22ac4a53a0c53f99999786427a` |
| `features.tsv` | `87a968c8cdca623b1d1575a87e6f5cad6414e95a577c59a049361618926d3cb2` |

The Flex-off raw MEX files were also byte-identical:

| File | SHA-256 |
| --- | --- |
| `matrix.mtx` | `865dad59e899298546dd42a5fd71d9d93674ed8baad93d9043ae77c7c0a3ca33` |
| `barcodes.tsv` | `74c2c9d186f4ffb3ea86fcb35d0f7a0de91bc0aee3b88e35d84bad89bbf7f25f` |
| `features.tsv` | `4eee35c6d94a78aec8629e8e1cf179dc6ed3e5d76795304f29dfb82fc375df1c` |

Recursive comparisons of both the Flex-off and Flex-on `Solo.out` trees found
no difference between v1.4.4 and the feature branch.

`tests/run_cbub_regression_test.sh` could not run because
`/storage/SC2300771_filtered_2M/sample_whitelist.tsv` is absent on this host.
This is an explicit fixture skip, not a pass. The 100K bridge determinism and
Flex regressions cover the available CB/UB paths.

## Release validation

- v1.4.4 PR Tier A: GitHub Actions run `29560701174`, passed.
- v1.4.4 artifact publication: GitHub Actions run `29561096150`, passed.
- v1.4.4 release: `https://github.com/morphic-bio/STAR-suite/releases/tag/v1.4.4`.
- Integration Tier A: GitHub Actions run `29563199754`, passed.
- Integration PR #5 merged as
  `b03930adff707f97c009d4d6116a677ff38168dc`.
- Release PR #6 Tier A: GitHub Actions run `29563878968`, passed.
- Release PR #6 merged as
  `fac38ca6485b3b3308792e114f88b6c943e7541b`; annotated tag `v1.5.0` peels to
  that exact merge.
- v1.5.0 artifact publication: GitHub Actions run `29564366223`, passed.
- Published release:
  `https://github.com/morphic-bio/STAR-suite/releases/tag/v1.5.0`.

Published primary artifact SHA-256 digests:

| Artifact | SHA-256 |
| --- | --- |
| glibc 2.34 tarball | `f04c1f7609a6f395a9b304390e2432ee88978b8401ceed9e451fab96ce855c59` |
| glibc 2.39 tarball | `fb7744fdd46070dda97be62f1ae4dfee5e6c6a55c466088c97a57867ea00737f` |
| compatibility installer | `14128256ed39fefd9c3076aebc0b216a48a6039b720b9cb6b9c3a2ccc3e98bac` |
| Debian binary package | `f0fbbaa425e1fe6b872f7b03f4fe67b5b5338aff3ebf2bcb1e2b017bab199630` |

The committed feature snapshot was also built with
`star-suite-release-build:ubuntu22` and validated with
`star-suite-runtime-check:ubuntu22`. The installed STAR and resolver both
reported `1.4.4`; the installed resolver passed its native smoke. The local
validation tarball SHA-256 was
`bcf4ca39a91dcfcf9f3ed4fbfe0e249d6532217208617252305a6fc88cbae87e`.
Installer conflict preflight and idempotent reinstall behavior passed.

The native-host `STARstatic` link is not used as a release verdict because
this host's system HTSlib lacks `sam_hdr_destroy`; release artifacts are built
and runtime-tested in the repository's pinned Ubuntu release containers.
