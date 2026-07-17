# STAR Suite Binary Distribution Plan

## Goal

Ship discoverable, installable Ubuntu binaries for STAR Suite with low-friction install and clear provenance.

## Implementation Status (2026-02-14)

- Phase 1 (Packaging Foundation): **implemented (initial)**
  - Completed:
    - Debian packaging skeleton added under `debian/`
      (`control`, `rules`, `changelog`, `copyright`,
      `source/format`, `docs`).
    - Release helper scripts added:
      - `scripts/release/build_static_tarball.sh`
      - `scripts/release/build_source_package.sh`
      - `scripts/release/build_deb_binary_package.sh`
      - `scripts/release/create_checksums.sh`
    - Clean-source Debian binary package build validated on Ubuntu toolchain
      (`dpkg-buildpackage -b` via `scripts/release/build_deb_binary_package.sh`).
    - Post-build package payload validation passed from extracted `.deb`:
      - `STAR --version`
      - `tests/slam/test_snp_mask_build_smoke.sh`
      - `tests/run_solo_smoke.sh`
    - Automated containerized install/uninstall validator added:
      - `scripts/release/validate_deb_install_uninstall.sh`
      - validates `apt install /path/to/star-suite*.deb`, binary presence,
        `STAR --version`, and `apt purge` cleanup on Ubuntu 24.04.
    - Packaging/build fixes required for clean export builds:
      - `Depend.list` generation hardened for `parametersDefault.xxd`.
      - explicit `libscrna` build rule added to core Makefile.
      - `remove_y_reads` Makefile now preserves distro flags while linking
        against vendored `htslib` and `zlib`.
  - Pending:
    - optional validation on additional host OS versions (for example Ubuntu 22.04)
    - optional package split decisions
- Phase 2 (Signing and Provenance): **implemented (initial)**
  - Completed:
    - source-package signing script:
      - `scripts/release/sign_source_package.sh`
    - release workflow support for optional source signing from CI secrets:
      - `DEB_SIGNING_KEY` (armored private key)
      - `DEB_SIGNING_KEY_ID` (key fingerprint or key ID)
      - `DEB_SIGNING_KEY_PASSPHRASE` (optional passphrase)
  - Pending:
    - production key provisioning/rotation and launchpad credential wiring
- Phase 3 (PPA Publishing): **pending**
- Phase 4 (GitHub Releases): **in progress**
  - Completed:
    - release workflow scaffold for artifact publication
    - release workflow builds compatibility tarballs + installer bundle + Debian binary package + Debian source package
  - Pending:
    - production signing and release key setup
    - stable release note templates/badges
- Phase 5 (CI/CD Automation): **implemented (initial)**
  - Completed:
    - `.github/workflows/ci-pr.yml`
    - `.github/workflows/ci-dev-release.yml`
    - `.github/workflows/ci-master.yml`
    - `.github/workflows/release.yml`
  - Pending:
    - Tier B fixture-backed runs in self-hosted/scheduled path
    - optional automated PPA upload gate

## Checkpoint Validation (2026-02-14)

- Compatibility tarball build:
  - `scripts/release/build_static_tarball.sh --version v0.0.0-test --compat-label glibc234 --glibc-baseline 2.34`
  - output verified: `STAR-suite-v0.0.0-test-linux-amd64-glibc234.tar.gz`
  - release workflow now builds two amd64 tarballs:
    - `glibc234` with `scripts/release/docker/Dockerfile.build-ubuntu22` for broader Linux compatibility
    - `glibc239` with `scripts/release/docker/Dockerfile.build-ubuntu24` for current systems
- Installer bundle build:
  - `scripts/release/build_installer_bundle.sh --version v0.0.0-test --tarball <glibc234> --tarball <glibc239>`
  - output verified: `STAR-suite-v0.0.0-test-linux-amd64-installer.tar.gz`
- Container runtime validation:
  - runtime checks are executed with `scripts/release/docker/Dockerfile.runtime-check`
  - host wrappers:
    - `scripts/release/run_tarball_runtime_container_check.sh`
    - `scripts/release/run_installer_bundle_container_check.sh`
  - release workflow writes per-container manifests under `dist/release/runtime-manifests/`
  - validated on March 10, 2026 in clean Ubuntu 22.04 and 24.04 containers:
    - `glibc234` tarball:
      - max referenced glibc symbol: `GLIBC_2.34`
      - Ubuntu 22.04 runtime packages: `libc6`, `zlib1g`, `libssl3`
    - `glibc239` tarball:
      - max referenced glibc symbol: `GLIBC_2.38`
      - Ubuntu 24.04 runtime packages: `libc6`, `zlib1g`, `libssl3t64`
    - installer bundle selected `glibc234` on Ubuntu 22.04 and `glibc239` on Ubuntu 24.04
- Dynamic Debian binary packaging:
  - script implemented and validated:
    - build: `scripts/release/build_deb_binary_package.sh --out-dir /tmp/star-release-test/deb`
    - extract: `dpkg-deb -x /tmp/star-release-test/deb/star-suite_*_amd64.deb /tmp/star-deb-root`
    - runtime checks:
      - `/tmp/star-deb-root/usr/bin/STAR --version`
      - `STAR_BIN=/tmp/star-deb-root/usr/bin/STAR tests/slam/test_snp_mask_build_smoke.sh`
      - `STAR_BIN=/tmp/star-deb-root/usr/bin/STAR tests/run_solo_smoke.sh`
  - install/uninstall validation script:
    - `scripts/release/validate_deb_install_uninstall.sh --deb /tmp/star-release-test/deb/star-suite_*_amd64.deb`
    - validates package install and purge lifecycle in a clean `ubuntu:24.04` container
- Debian source packaging:
  - script implemented with clean `git archive` snapshot:
    `scripts/release/build_source_package.sh`
  - optional signing script:
    - `scripts/release/sign_source_package.sh --in-dir dist/release/source --key-id <KEY_ID>`
- Compatibility tarball runtime check:
  - extracted binary reports expected version: `2.7.11b`
- Tier A smoke tests against extracted compatibility tarball binary (`STAR_BIN=<extracted>/bin/STAR`):
  - `tests/slam/test_snp_mask_build_smoke.sh`: passed
  - `tests/run_solo_smoke.sh`: passed

## Published Release Smoke (2026-03-11)

- Release tested:
  - `v0.50.3`
  - GitHub prerelease assets downloaded from the published release page, not from local build outputs
- Download verification:
  - `SHA256SUMS` checked against the published installer bundle and both published compatibility tarballs
- External-user installer checks:
  - host install on a `glibc 2.35` machine selected `glibc234` and `STAR --version` returned `2.7.11b`
  - manual direct install from the published `glibc234` tarball also returned `2.7.11b`
  - clean `ubuntu:24.04` container install from the published installer bundle selected `glibc239` and `STAR --version` returned `2.7.11b`
- This confirms the published release assets, checksums, and installer selection behavior match the intended compatibility model outside the repo build environment

## Production Release Prep (2026-05-22)

- Candidate tag: `v1.0.0`
- Release notes: `docs/RELEASE_NOTES_v1.0.0.md`
- Packaging metadata: `debian/changelog` entry `1.0.0-1`
- Scope added since `v0.94.0`: native OCM composite-barcode production mode,
  native per-sample OCM MEX/Velocyto materialization, updated JAX OCM and MSK
  40KO production recipes, adaptive MT downstream packaging, and the current
  multiome post-MEX handoff docs.
- Version note: release artifacts and Debian packages use STAR Suite `1.0.0`;
  `STAR --version` reports `1.0.0`. Upstream provenance is reported by
  `STAR --upstream-version` (`2.7.11b`), and genome index compatibility is
  reported by `STAR --genome-compat-version` (`2.7.4a`).
- Local caveat: rebuild release binaries only from an idle production checkout
  or from CI so no running wrapper observes a partially rebuilt
  `core/legacy/source/STAR`.

## Hotfix Release Prep (2026-05-28)

- Candidate tag: `v1.0.1`
- Release notes: `docs/RELEASE_NOTES_v1.0.1.md`
- Packaging metadata: `debian/changelog` entry `1.0.1-1`
- Scope added since `v1.0.0`: deterministic `process_features`
  feature-bootstrap handling for CR-compatible dynamic multi-consumer runs and
  corrected single-thread feature-match position reporting for
  `--crAssignSearchThreads 1`.
- Version note: release artifacts and Debian packages use STAR Suite `1.0.1`;
  `STAR --version` reports `1.0.1`. Upstream provenance and genome index
  compatibility remain `2.7.11b` and `2.7.4a`.

## Feature Release Prep (2026-05-29)

- Candidate tag: `v1.1.0`
- Release notes: `docs/RELEASE_NOTES_v1.1.0.md`
- Packaging metadata: `debian/changelog` entry `1.1.0-1`
- Scope added since `v1.0.3`: native CBQ/BINSEQ input support through the
  shared C++ reader, STAR read-buffer adapter, ordered FASTQ/FASTQ.gz-to-CBQ
  encoder, STAR mapper/STARsolo/OCM/Flex/SLAM/process_features smoke coverage,
  and PE-SLAM CBQ divergence testing.
- Version note: release artifacts and Debian packages use STAR Suite `1.1.0`;
  `STAR --version` reports `1.1.0`. Upstream provenance and genome index
  compatibility remain `2.7.11b` and `2.7.4a`.

## Minor Feature-Label Refresh (2026-06-01)

- Candidate tag: `v1.3.0b`
- Release notes: `docs/RELEASE_NOTES_v1.3.0b.md`
- Packaging metadata: `debian/changelog` entry `1.3.0b-1`
- Scope added since `v1.3.0`: explicit STAR-Flex Fixed RNA CBQ MCP recipe,
  CBQ-native multiome recipe inputs for GEX plus libchromap ATAC, OCM CBQ
  Y/noY recipe metadata, and recipe defaults that avoid FLEX `zcat` unless
  explicitly requested.
- Version note: release artifacts and Debian packages use STAR Suite
  `1.3.0b`; `STAR --version` reports `1.3.0b`. Upstream provenance and genome
  index compatibility remain `2.7.11b` and `2.7.4a`.
- Hosted static tarballs and Debian packages are built through the portable
  no-Chromap target so GitHub release artifacts do not depend on a sibling
  Chromap-suite checkout; local production multiome builds still use the
  Chromap-enabled `make core` default.

## Feature Release Prep (2026-06-01)

- Candidate tag: `v1.3.0`
- Release notes: `docs/RELEASE_NOTES_v1.3.0.md`
- Packaging metadata: `debian/changelog` entry `1.3.0-1`
- Scope added since `v1.1.0`: parallel CBQ integration across STAR core and
  Flex, native process_features CBQ range-reader coverage, native libchromap
  multiome ATAC CBQ input, and PBMC multiome GEX+ATAC CBQ benchmark
  documentation with ATAC Y/noY BAM validation.
- Version note: release artifacts and Debian packages use STAR Suite `1.3.0`;
  `STAR --version` reports `1.3.0`. Upstream provenance and genome index
  compatibility remain `2.7.11b` and `2.7.4a`.

## Patch Release v1.4.4 (2026-07-17)

- Release tag: `v1.4.4`
- Release notes: `docs/RELEASE_NOTES_v1.4.4.md`
- Packaging metadata: `debian/changelog` entry `1.4.4-1`
- Scope added since `v1.4.3`: corrected feature-barcode namespace remapping in
  integrated AnnData/MuData exports and deterministic UCSF per-sample Slurm
  workflow contracts.
- Version note: release artifacts use tag `v1.4.4`, Debian packages use STAR
  Suite `1.4.4`, and `STAR --version` reports `1.4.4`. Upstream provenance and
  genome index compatibility remain `2.7.11b` and `2.7.4a`.

## Patch Release v1.4.2 (2026-06-29)

- Release tag: `v1.4.2`
- Release notes: `docs/RELEASE_NOTES_v1.4.2.md`
- Packaging metadata: `debian/changelog` entry `1.4.2-1`
- Scope added since `v1.4.1`: TranscriptVB Salmon parity improvements,
  STAR-generated transcriptome GC/effective-length sidecar integration, OpenMP
  linkage for the GC-bias helper, PE bulk benchmark production-mode
  documentation, and SE/PE decompressed byte-parity regression coverage for the
  samtools-style BAM sorter.
- Version note: release artifacts use tag `v1.4.2`, Debian packages use STAR
  Suite `1.4.2`, and `STAR --version` reports `1.4.2`. Upstream provenance and
  genome index compatibility remain `2.7.11b` and `2.7.4a`.

## Hotfix Release Prep (2026-06-20)

- Candidate tag: `v1.4.1`
- Release notes: `docs/RELEASE_NOTES_v1.4.1.md`
- Packaging metadata: `debian/changelog` entry `1.4.1-1`
- Scope added since `v1.3.1`: Flex/Solo dependency tracking now covers
  VPATH-built top-level Flex sources so incremental builds do not reuse stale
  objects after shared header/layout changes. The tiny public Flex smoke wrapper
  now validates both no-BAM and unsorted-BAM output modes.
- Version note: release artifacts use tag `v1.4.1`, Debian packages use STAR
  Suite `1.4.1`, and `STAR --version` reports `1.4.1`. Upstream provenance and
  genome index compatibility remain `2.7.11b` and `2.7.4a`.

## Hotfix Release Prep (2026-05-29)

- Candidate tag: `v1.0.3`
- Release notes: `docs/RELEASE_NOTES_v1.0.3.md`
- Packaging metadata: `debian/changelog` entry `1.0.3-1`
- Scope added since `v1.0.2`: non-Flex Solo bridge determinism is hardened by
  canonicalizing worker-local bridge hashes before collapse, aggregating
  ambiguous cell-barcode evidence deterministically, and using the resolved
  bridge result for read-accounting statistics.
- Version note: release artifacts and Debian packages use STAR Suite `1.0.3`;
  `STAR --version` reports `1.0.3`. Upstream provenance and genome index
  compatibility remain `2.7.11b` and `2.7.4a`.

- Candidate tag: `v1.0.2`
- Release notes: `docs/RELEASE_NOTES_v1.0.2.md`
- Packaging metadata: `debian/changelog` entry `1.0.2-1`
- Scope added since `v1.0.1`: exact duplicate UMI rows are combined before
  1MM_CR correction on the non-Flex Solo direct inline-hash bridge collapse
  path. This prevents duplicate exact UMIs split across worker-local bridge
  hashes from entering the CR correction routine as independent observations.
- Version note: release artifacts and Debian packages use STAR Suite `1.0.2`;
  `STAR --version` reports `1.0.2`. Upstream provenance and genome index
  compatibility remain `2.7.11b` and `2.7.4a`.

## Distribution Strategy

Use two channels in parallel:

1. **Launchpad PPA** for Ubuntu-native apt discovery and updates.
2. **GitHub Releases** for project-native downloads, release notes, and provenance artifacts.
3. **Container registry (multi-arch)** for reproducible runtime images.

## Scope

- Primary package: `star-suite`
- Initial target OS:
  - Ubuntu 24.04
  - Ubuntu 22.04 compatibility target via release tarball validation
- Build type:
  - Dynamic binaries (default runtime path)
  - Compatibility tarballs plus installer bundle (portable fallback path)

## Artifact Matrix

Per release tag (`v*`), publish both binary styles:

- Dynamic:
  - `star-suite_<version>_amd64.deb`
  - `star-suite_<version>_arm64.deb`
- Compatibility tarballs:
  - `STAR-suite-<version>-linux-amd64-glibc234.tar.gz`
  - `STAR-suite-<version>-linux-amd64-glibc239.tar.gz`
- Installer bundle:
  - `STAR-suite-<version>-linux-amd64-installer.tar.gz`
- Containers:
  - `biodepot/star-suite:<version>` (multi-arch manifest)
  - `biodepot/star-suite:latest` (stable pointer only)

Notes:

- `.deb` remains the preferred Ubuntu installation path.
- The installer bundle is the preferred non-`.deb` installation path.
- Direct compatibility tarballs are for manual installs and HPC/container edge cases.
- `v0.*` tags are prereleases for testing. They publish release artifacts and a versioned container tag, but do not move `latest`.
- `v1.0.0` is the first production release and may move the `latest` container
  tag when release CI succeeds.
- Native `build_static_tarball.sh` runs still inherit the host toolchain/glibc baseline.
- The CI/release path publishes multiple tarballs because Linux may reject a binary built for a newer runtime environment before STAR Suite starts.
- The installer bundle checks the host environment and selects the highest compatible bundled binary automatically.
- Release artifacts now include runtime-manifest text files captured from clean Ubuntu 22.04 and 24.04 containers.
- Those manifests record the resolved runtime libraries and Ubuntu/Debian package names seen by the validator.

## OS Compatibility Note

For non-technical users:

- If a downloaded Linux binary does not run, that is usually the operating system rejecting a binary built for a newer system environment.
- STAR Suite itself is not crashing and is not asking the user to install extra libraries manually.
- The release fix is to ship multiple Linux binaries built for different compatibility levels, plus an installer bundle that picks the right one.

For technical users:

- The tarballs are only partially static.
- Current validated runtime SONAMEs are:
  - `libc.so.6`
  - `libm.so.6`
  - `libz.so.1`
  - `libcrypto.so.3`
- On Ubuntu 22.04 the validated package set is:
  - `libc6`
  - `zlib1g`
  - `libssl3`
- On Ubuntu 24.04 the validated package set is:
  - `libc6`
  - `zlib1g`
  - `libssl3t64`
- The `glibc239` asset label tracks the Ubuntu 24.04 build environment.
- The March 10, 2026 runtime manifest showed that this build currently references `GLIBC_2.38`, so the label is conservative.
- The exact resolved package names are recorded in the release `runtime-manifests/` artifacts for each container baseline.

## Phase 1: Packaging Foundation

1. Add Debian packaging metadata under `debian/`:
   - `control`, `rules`, `changelog`, `copyright`
   - `install`, `docs`, and post-install metadata as needed
2. Define package contents:
   - Core binaries (`STAR`, `star_feature_call`, `slam_requant`, `pileup_snp`, `flexfilter`, etc.)
   - Optional helper scripts (if stable and versioned)
3. Ensure install paths:
   - Binaries under `/usr/bin` or `/usr/lib/star-suite/bin` with symlinks
   - Docs under `/usr/share/doc/star-suite`
4. Validate package builds locally:
   - `dpkg-buildpackage` in clean Ubuntu 24.04 environment
   - Confirm `dpkg -i` works on fresh system

## Phase 2: Signing and Provenance

1. Create release GPG key (team-owned, documented rotation owner).
2. Sign source and binary artifacts.
3. Publish checksums (`sha256sum`) per release.
4. Include commit hash and build metadata in release notes.

## Phase 3: PPA Publishing (Discoverability)

1. Create Launchpad project/PPA:
   - Suggested name: `ppa:biodepot/star-suite`
2. Upload source package for each release tag.
3. Validate install flow:
   - `add-apt-repository`, `apt update`, `apt install star-suite`
4. Add PPA install instructions to top-level `README.md`.

## Phase 4: GitHub Releases (Discoverability + Auditability)

1. Publish each release tag with:
   - `.deb` artifacts
   - checksums
   - short install instructions
   - compatibility notes (tested Ubuntu versions)
2. Add badges/links in top-level `README.md`:
   - Latest release
   - PPA install path

## Phase 5: CI/CD Automation (GitHub Actions)

1. Add `ci-pr.yml`:
   - Trigger: `pull_request`
   - Run fast checks only (`amd64` build + Tier A smoke), no publish.
2. Add `ci-dev-release.yml`:
   - Trigger: push to `dev-release` and `dev-release-*`
   - Run integration checks and optional publish to `dev-release-<sha>` /
     `dev-release-latest`.
3. Add `ci-master.yml`:
   - Trigger: push to `master`
   - Build/push multi-arch runtime images (`amd64`, `arm64`).
4. Add `release.yml`:
   - Trigger: tags `v*`
   - Build and publish:
     - dynamic artifacts (`.deb` path, preferably via Launchpad source upload)
     - compatibility tarball + installer bundle artifacts
     - multi-arch image tags
   - Publish checksums and release notes.
5. Add PPA upload step in release workflow:
   - upload signed source package
   - Launchpad builds native architecture `.deb` packages.

## Quality Gates

- Package installs on clean Ubuntu 24.04 with no manual dependency fixing.
- `STAR --version` and core smoke tests pass from installed binaries.
- Static binary tarball runs on clean target host (documented dependency expectations).
- `amd64` and `arm64` release artifacts are both generated and verifiable.
- Release notes contain:
  - STAR Suite version
  - commit hash
  - dependency list
  - known limitations
- PPA and GitHub release both updated for every stable tag.

## Risks and Mitigations

- **Dependency drift across Ubuntu versions**
  - Mitigation: pin build matrix to explicit Ubuntu versions and run install tests.
- **Packaging stale binaries instead of fresh build artifacts**
  - Mitigation: enforce clean build in CI and verify artifact checksums.
- **Signing key management risk**
  - Mitigation: documented ownership, rotation procedure, and backup policy.

## Deliverables

- `debian/` packaging files in repo
- Build/install docs in top-level `README.md`
- Launchpad PPA configured and documented
- GitHub Releases with signed dynamic + static artifacts and checksums
- Multi-arch container images published per release
- CI workflows for repeatable build/test/publish paths

## Immediate Next Steps

1. Add release signing and checksum verification policy (GPG key ownership + rotation).
2. Configure Docker Hub and release secrets in GitHub repository settings.
3. Add README install section (`apt` via PPA + manual `.deb` + installer bundle).
4. Add optional PPA upload step (guarded by release environment approval).
5. Validate install/uninstall behavior on Ubuntu 22.04 as secondary target.
