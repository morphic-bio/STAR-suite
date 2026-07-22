# Binary Test Runbook (Ubuntu 22.04 / 24.04)

## Goal

Define a reproducible binary-validation matrix for installed STAR-suite release
artifacts, using:

- clean Docker containers
- Ubuntu `22.04` and `24.04`
- public downloads where possible
- existing local artifacts where public-download coverage does not exist yet

The immediate target is installed-binary testing from release artifacts
(tarball or installer bundle), not developer-tree testing.

## Scope

Required workflow families:

- PE bulk
- scRNA GEX
- Perturb-seq
- FLEX
- SLAM

Required output/mode coverage:

- default no-BAM
- BAM `SortedByCoordinate`
- BAM `Unsorted`
- Y-removal where applicable
- adapter trimming when applicable
- Salmon / TranscriptVB where applicable

## Baseline Assumptions

- Container baselines are Ubuntu `22.04` and `24.04`.
  Note: Ubuntu `24.02` is not a real release; use `24.04`.
- Binary validation should execute against an installed release artifact inside
  the container, not the repo build tree.
- Public/download-based tests are the preferred release-gate surface.
- Artifact-backed tests are still valuable, but should remain a secondary lane
  until they are portable enough for container execution.

## Current Release Coverage

Current release smoke validation is narrow.

- Runtime/install validation already exists for tarballs and installer bundles
  on Ubuntu `22.04` and `24.04` via:
  - [release.yml](/mnt/pikachu/STAR-suite/.github/workflows/release.yml)
  - [run_tarball_runtime_container_check.sh](/mnt/pikachu/STAR-suite/scripts/release/run_tarball_runtime_container_check.sh)
  - [run_installer_bundle_container_check.sh](/mnt/pikachu/STAR-suite/scripts/release/run_installer_bundle_container_check.sh)
- Installed-binary smoke execution currently runs only on Ubuntu `24.04` via
  [container_check_release_smokes.sh](/mnt/pikachu/STAR-suite/scripts/release/docker/container_check_release_smokes.sh):
  - [run_solo_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_solo_smoke.sh)
  - [test_snp_mask_build_smoke.sh](/mnt/pikachu/STAR-suite/tests/slam/test_snp_mask_build_smoke.sh)
  - [run_flex_tiny_public_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_flex_tiny_public_smoke.sh)

That is not yet a real workflow-parity matrix.

## Design Principles

1. Test the installed binary, not the in-tree binary.
2. Keep the default lane no-BAM unless BAM is the feature under test.
3. Prefer public-download fixtures for release gating.
4. Reuse existing scripts first; add thin wrappers before inventing new
   end-to-end harnesses.
5. Split the matrix into:
   - Tier A: public/download-based, container-safe, release-gate candidates
   - Tier B: artifact-backed, host-dependent, extended parity coverage
6. Do not require a single mega-test to cover every feature. Use one primary
   no-BAM parity/smoke per workflow, then targeted extension tests for BAM,
   Y-removal, trimming, or VB.

## Proposed Matrix

| Workflow | Primary binary test | Data source | Default no-BAM | Sorted BAM | Unsorted BAM | Y-removal | Trim/adapters | Salmon/VB | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PE bulk | [run_public_bulk_pe_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_public_bulk_pe_smoke.sh) | public download | no | yes | no | yes | yes | yes | reusable now |
| scRNA GEX | new installed-binary wrapper around [run_default_bundle_scrna_fixture.sh](/mnt/pikachu/STAR-suite/tests/run_default_bundle_scrna_fixture.sh) or [run_solo_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_solo_smoke.sh) | repo fixture | needs new no-BAM variant | possible | possible | no | no | no | gap |
| Perturb | [run_a375_public_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_a375_public_smoke.sh) | public download | yes for 5' | yes | no | no | chemistry-specific | no | reusable now |
| FLEX | new no-BAM variant of [run_flex_tiny_public_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_flex_tiny_public_smoke.sh) | public download | needs change | possible | possible | maybe | no | no | gap |
| SLAM | new installed-binary wrapper around [run_slam_end_to_end_fixture.sh](/mnt/pikachu/STAR-suite/tests/run_slam_end_to_end_fixture.sh) plus current [test_snp_mask_build_smoke.sh](/mnt/pikachu/STAR-suite/tests/slam/test_snp_mask_build_smoke.sh) | existing artifact today | needs new no-BAM path | current script uses sorted BAM | no | no | fixed trims currently | no | gap |

## Workflow Notes

### 1. PE bulk

Best immediate binary parity candidate:
[run_public_bulk_pe_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_public_bulk_pe_smoke.sh)

Why:

- public downloadable fixture
- exercises adapter trimming
- exercises Y-removal split outputs
- exercises TranscriptVB
- exercises Salmon comparison
- already validates PE-oriented outputs

Limitation:

- currently writes sorted BAM by design
- not a default no-BAM lane

Plan:

- Keep this as the PE extension test.
- Add a separate no-BAM PE binary smoke later if we want a default no-BAM PE
  lane; this is lower priority because the PE-specific value here is the
  trim/Y/VB path.

### 2. scRNA GEX

Current scripts:

- [run_default_bundle_scrna_fixture.sh](/mnt/pikachu/STAR-suite/tests/run_default_bundle_scrna_fixture.sh)
- [run_solo_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_solo_smoke.sh)

Current problem:

- neither is a parity test
- [run_solo_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_solo_smoke.sh) always
  writes sorted BAM
- [run_default_bundle_scrna_fixture.sh](/mnt/pikachu/STAR-suite/tests/run_default_bundle_scrna_fixture.sh)
  is good for default-bundle validation but is not a real installed-binary
  workflow parity surface

Plan:

- Create a binary-oriented GEX smoke that:
  - runs on a tiny synthetic or public small fixture
  - defaults to `--outSAMtype None`
  - verifies `Solo.out/GeneFull/raw` and filtered outputs
  - can optionally toggle sorted and unsorted BAM for CB/UB validation

Recommended target script:

- `tests/run_gex_binary_smoke.sh`

### 3. Perturb-seq

Best immediate binary parity candidate:
[run_a375_public_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_a375_public_smoke.sh)

Why:

- public downloadable fixture
- already supports auto-download
- already defaults to no-BAM on the `5'` path
- exercises CR-compat multi-library behavior
- validates combined GEX + CRISPR outputs

Limitation:

- current script is a smoke, not a deep parity benchmark
- unsorted BAM is not covered

Plan:

- Use this as the Tier A perturb binary test now.
- Add an artifact-backed Tier B parity wrapper later if exact CR9 metrics need
  to be checked from installed binaries.

### 4. FLEX

Current public candidate:
[run_flex_tiny_public_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_flex_tiny_public_smoke.sh)

Why it is close:

- public/download-based tiny dataset
- builds the tiny filtered reference and STAR index
- runs the real Flex wrapper:
  [run_flex_cr_config.sh](/mnt/pikachu/STAR-suite/scripts/run_flex_cr_config.sh)

Current problem:

- it currently requires `Aligned.out.bam`
- it is not aligned with the desired default no-BAM release surface

Plan:

- Refactor or add a sibling script:
  - `tests/run_flex_tiny_public_binary_smoke.sh`
- Default behavior:
  - `--outSAMtype None`
  - validate per-sample filtered matrices and manifest/log outputs
- Optional extension flags:
  - `--write-bam sorted`
  - `--write-bam unsorted`

Artifact-backed extension:

- [run_flex_cr_config_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_flex_cr_config_smoke.sh)
  remains useful as Tier B coverage on existing local Flex artifacts.

### 5. SLAM

Current state:

- [test_snp_mask_build_smoke.sh](/mnt/pikachu/STAR-suite/tests/slam/test_snp_mask_build_smoke.sh)
  is already in the release smoke lane, but it is only a build/smoke test.
- [run_slam_end_to_end_fixture.sh](/mnt/pikachu/STAR-suite/tests/run_slam_end_to_end_fixture.sh)
  is the real workflow test, but it is artifact-backed and currently BAM-based.

Current problem:

- no public-download, container-stable SLAM E2E fixture exists in the current
  repo workflow
- the existing E2E script expects sorted BAM outputs and GEDI comparison assets

Plan:

- Keep `test_snp_mask_build_smoke.sh` in Tier A immediately.
- Treat [run_slam_end_to_end_fixture.sh](/mnt/pikachu/STAR-suite/tests/run_slam_end_to_end_fixture.sh)
  as Tier B until a portable downloaded fixture exists.
- Later, add a smaller binary-friendly SLAM smoke that validates:
  - `SlamQuant.out`
  - non-empty mutation profile
  - optional sorted BAM mode only

## Extension Coverage

These should not be forced into the primary no-BAM workflow scripts. They
should be separate focused checks.

### BAM sorted / unsorted

Recommended existing scripts:

- sorted non-Flex CB/UB:
  [run_sorted_bam_cbub_nonflex_test.sh](/mnt/pikachu/STAR-suite/tests/run_sorted_bam_cbub_nonflex_test.sh)
- unsorted non-Flex/Flex-style CB/UB coverage:
  [run_unsorted_bam_cbub_test.sh](/mnt/pikachu/STAR-suite/tests/run_unsorted_bam_cbub_test.sh)

Current limitation:

- both are artifact-backed
- they are not container-portable release-gate tests today

Plan:

- Keep them as Tier B initially.
- Longer term, build tiny downloadable or synthetic CB/UB fixtures so the same
  checks can run against installed binaries in containers.

### Y-removal

Recommended existing scripts:

- bulk PE Y split:
  [run_ychrom_bulk_pe_test.sh](/mnt/pikachu/STAR-suite/tests/run_ychrom_bulk_pe_test.sh)
- Solo Y split:
  [run_solo_yremove_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_solo_yremove_smoke.sh)

Plan:

- PE Y-removal is already covered best by
  [run_public_bulk_pe_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_public_bulk_pe_smoke.sh).
- Keep the stronger artifact-backed Y-removal scripts as Tier B until a
  portable GEX/Solo Y-removal fixture exists.

### Adapter trimming / trim QC

Recommended scripts:

- [run_public_bulk_pe_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_public_bulk_pe_smoke.sh)
- [run_trim_qc_replay_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_trim_qc_replay_smoke.sh)

Plan:

- PE bulk remains the main end-to-end adapter/trim parity lane.
- `trim_qc_replay` stays as a focused binary/tool smoke.

### Salmon / TranscriptVB

Primary coverage:

- [run_public_bulk_pe_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_public_bulk_pe_smoke.sh)

This is sufficient for the first binary matrix because it already validates the
main PE + TranscriptVB + Salmon surface.

## Tier Definition

### Tier A: release-gate candidates

These should be runnable inside a clean Ubuntu `22.04` or `24.04` container
from an installed release artifact.

Initial Tier A target set:

- scRNA GEX binary smoke, no-BAM
- perturb binary smoke:
  [run_a375_public_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_a375_public_smoke.sh)
- Flex tiny public binary smoke, no-BAM
- PE public bulk smoke:
  [run_public_bulk_pe_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_public_bulk_pe_smoke.sh)
- SLAM binary smoke:
  [test_snp_mask_build_smoke.sh](/mnt/pikachu/STAR-suite/tests/slam/test_snp_mask_build_smoke.sh)

### Tier B: extended binary parity / artifact-backed

- [run_flex_cr_config_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_flex_cr_config_smoke.sh)
- [run_slam_end_to_end_fixture.sh](/mnt/pikachu/STAR-suite/tests/run_slam_end_to_end_fixture.sh)
- [run_sorted_bam_cbub_nonflex_test.sh](/mnt/pikachu/STAR-suite/tests/run_sorted_bam_cbub_nonflex_test.sh)
- [run_unsorted_bam_cbub_test.sh](/mnt/pikachu/STAR-suite/tests/run_unsorted_bam_cbub_test.sh)
- [run_solo_yremove_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_solo_yremove_smoke.sh)
- [run_ychrom_bulk_pe_test.sh](/mnt/pikachu/STAR-suite/tests/run_ychrom_bulk_pe_test.sh)

## Container Dependency Gaps

The current runtime-check image in
[Dockerfile.runtime-check](/mnt/pikachu/STAR-suite/scripts/release/docker/Dockerfile.runtime-check)
is sufficient for the existing narrow smoke lane, but not yet for the proposed
workflow matrix.

Known gaps:

- [run_public_bulk_pe_smoke.sh](/mnt/pikachu/STAR-suite/tests/run_public_bulk_pe_smoke.sh)
  requires `salmon`
- [download_a375_5prime_fixture.sh](/mnt/pikachu/STAR-suite/scripts/download_a375_5prime_fixture.sh)
  requires `curl`
- [download_a375_3prime_fixture.sh](/mnt/pikachu/STAR-suite/scripts/download_a375_3prime_fixture.sh)
  requires `curl`

Recommended first-pass runtime image additions:

- `curl`
- `salmon`

`tar`, `python3`, `git`, `git-lfs`, and `samtools` are already accounted for or
available from the current base image plus repo tooling.

## Container Execution Model

The existing release container wrappers are the right foundation:

- [run_release_smokes_container_check.sh](/mnt/pikachu/STAR-suite/scripts/release/run_release_smokes_container_check.sh)
- [container_check_release_smokes.sh](/mnt/pikachu/STAR-suite/scripts/release/docker/container_check_release_smokes.sh)

Recommended expansion:

1. Keep install/runtime validation exactly as it is.
2. Replace the hard-coded three-test smoke block with a profile-driven runner.
3. Add explicit smoke profiles:
   - `core`
   - `binary-workflows-tier-a`
   - `binary-workflows-tier-b`
4. Run Tier A on both Ubuntu `22.04` and `24.04`.
5. Keep Tier B manual, scheduled, or self-hosted until fixtures are portable.

## Implementation Order

### Phase 1: establish binary workflow lane

1. Add this runbook.
2. Create no-BAM binary wrappers where missing:
   - `tests/run_gex_binary_smoke.sh`
   - `tests/run_flex_tiny_public_binary_smoke.sh`
3. Update
   [container_check_release_smokes.sh](/mnt/pikachu/STAR-suite/scripts/release/docker/container_check_release_smokes.sh)
   to support a workflow profile.
4. Run Tier A in Ubuntu `22.04` and `24.04`.

### Phase 2: add BAM-mode extensions

1. Add sorted/unsorted toggles to the new GEX and Flex binary wrappers.
2. Add installed-binary BAM extension checks as a separate profile.
3. Keep artifact-backed CB/UB tests as the reference oracle until portable
   fixtures exist.

### Phase 3: add artifact-backed parity lane

1. Wrap existing artifact-backed tests so they can consume `STAR_BIN` from an
   installed tarball or bundle.
2. Document required host datasets in
   [tests/ARTIFACTS.md](/mnt/pikachu/STAR-suite/tests/ARTIFACTS.md).
3. Run these on self-hosted or manual jobs only.

## Immediate Next Tasks

1. Build the Tier A no-BAM wrappers for GEX and FLEX.
2. Convert release smoke execution from a hard-coded list to a named profile.
3. Add Ubuntu `22.04` smoke execution, not just runtime validation.
4. Keep PE bulk and perturb public smokes as-is for the first pass.
5. Defer full SLAM parity and portable BAM/Y-removal parity until after the
   first containerized workflow lane is stable.

## Recommended Acceptance Criteria

For the first binary-test milestone, require all of the following from an
installed release artifact in Ubuntu `22.04` and `24.04`:

- `STAR --version` passes
- scRNA GEX no-BAM binary smoke passes
- perturb public no-BAM smoke passes
- Flex tiny public no-BAM smoke passes
- PE bulk trim/Y/VB/Salmon smoke passes
- SLAM SNP-mask smoke passes

The second milestone adds:

- sorted BAM extension checks
- unsorted BAM extension checks
- artifact-backed Y-removal and CB/UB parity checks
- artifact-backed SLAM E2E parity

## Summary

The repo already has most of the ingredients, but not yet the right
composition for release-quality binary workflow testing. The short path is:

- reuse public PE bulk and A375 perturb scripts now
- refactor Flex into a no-BAM public binary smoke
- add a real no-BAM GEX binary smoke
- keep SLAM E2E and BAM/Y-removal checks as Tier B until portable fixtures
  exist
- run Tier A against installed artifacts on Ubuntu `22.04` and `24.04`
