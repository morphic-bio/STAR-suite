# STAR-suite Binary Distribution Plan

## Goal

Ship discoverable, installable Ubuntu binaries for STAR-suite with low-friction install and clear provenance.

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
    - release workflow builds static tarball + Debian binary package + Debian source package
  - Pending:
    - production signing and release key setup
    - stable release note templates/badges
- Phase 5 (CI/CD Automation): **implemented (initial)**
  - Completed:
    - `.github/workflows/ci-pr.yml`
    - `.github/workflows/ci-dev.yml`
    - `.github/workflows/ci-master.yml`
    - `.github/workflows/release.yml`
  - Pending:
    - Tier B fixture-backed runs in self-hosted/scheduled path
    - optional automated PPA upload gate

## Checkpoint Validation (2026-02-14)

- Static release artifact build:
  - `scripts/release/build_static_tarball.sh --version v0.0.0-test`
  - output verified: `STAR-static-v0.0.0-test-linux-amd64.tar.gz`
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
- Static binary runtime check:
  - extracted binary reports expected version: `2.7.11b`
- Tier A smoke tests against extracted static binary (`STAR_BIN=<extracted>/bin/STAR`):
  - `tests/slam/test_snp_mask_build_smoke.sh`: passed
  - `tests/run_solo_smoke.sh`: passed

## Distribution Strategy

Use two channels in parallel:

1. **Launchpad PPA** for Ubuntu-native apt discovery and updates.
2. **GitHub Releases** for project-native downloads, release notes, and provenance artifacts.
3. **Container registry (multi-arch)** for reproducible runtime images.

## Scope

- Primary package: `star-suite`
- Initial target OS:
  - Ubuntu 24.04 (required)
  - Ubuntu 22.04 (optional follow-up)
- Build type:
  - Dynamic binaries (default runtime path)
  - Static binaries (portable fallback path, like upstream STAR distribution model)

## Artifact Matrix

Per release tag (`v*`), publish both binary styles:

- Dynamic:
  - `star-suite_<version>_amd64.deb`
  - `star-suite_<version>_arm64.deb`
- Static:
  - `STAR-static-<version>-linux-amd64.tar.gz`
  - `STAR-static-<version>-linux-arm64.tar.gz`
- Containers:
  - `biodepot/star-suite:<version>` (multi-arch manifest)
  - `biodepot/star-suite:latest` (stable pointer only)

Notes:

- `.deb` remains the preferred Ubuntu installation path.
- Static tarballs are for minimal-host-dependency deployments and HPC/container edge cases.
- `v0.*` tags are prereleases for testing. They publish release artifacts and a versioned container tag, but do not move `latest`.

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
2. Add `ci-dev.yml`:
   - Trigger: push to `dev`
   - Run integration checks and optional publish to `dev-<sha>`.
3. Add `ci-master.yml`:
   - Trigger: push to `master`
   - Build/push multi-arch runtime images (`amd64`, `arm64`).
4. Add `release.yml`:
   - Trigger: tags `v*`
   - Build and publish:
     - dynamic artifacts (`.deb` path, preferably via Launchpad source upload)
     - static tarball artifacts
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
  - STAR-suite version
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
3. Add README install section (`apt` via PPA + manual `.deb` + static tarball).
4. Add optional PPA upload step (guarded by release environment approval).
5. Validate install/uninstall behavior on Ubuntu 22.04 as secondary target.
