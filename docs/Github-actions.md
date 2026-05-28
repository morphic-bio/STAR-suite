# GitHub Actions CI/CD Policy

This document defines when STAR Suite CI should run, what it should execute,
and when artifacts are published.

## Goals

- Keep PR feedback fast.
- Avoid expensive publishes on every push.
- Publish reproducible multi-arch artifacts for stable branches/tags.
- Reuse the Docker multi-stage build already in this repo.

## Trigger Matrix

| Event | Branch/Tag | Build Scope | Tests | Publish |
|---|---|---|---|---|
| Pull Request | any | Docker build sanity (`amd64`) | Tier A | No |
| Push | `dev-release`, `dev-release-*` | integration build (`amd64`) | Tier A + selected fixture-free checks | Optional `dev-release-<sha>` image |
| Push | `master` | full multi-arch build (`amd64`,`arm64`) | required checks | Yes (`latest`, `master-<sha>`) |
| Tag Push | `v*` | release multi-arch build | release gates | Yes (version tags + release assets + PPA source upload) |
| Scheduled (optional) | nightly | drift detection build | extended checks | Optional `nightly` tag |

Branch and tag policy:

- `master` is production/stable. Stable release tags (`vX.Y.Z`) are cut from
  `master` after release-candidate validation.
- `dev-release` is the moving advanced-user preview branch for the next release.
  It should contain only work that is intended for the next production release,
  not raw experimental development.
- `dev-release-vX.Y.Z` is the optional version-scoped release-candidate branch.
  Use it when a candidate needs multiple fixes without moving the general
  `dev-release` pointer unpredictably.
- `vX.Y.Z-rcN` tags are immutable prerelease tags cut from `dev-release` or
  `dev-release-vX.Y.Z`. Use these for reproducible advanced-user testing and
  bug reports.
- After an RC is accepted, merge the version branch to `master` with
  `git merge --no-ff`, then tag `vX.Y.Z`.

Tag policy:

- Tags matching `v0.*` are treated as prereleases.
- Tags containing `-alpha`, `-beta`, `-rc`, or `-dev` are treated as prereleases.
- Prerelease tags publish versioned container images and GitHub prereleases, but do not move Docker `latest`.
- Stable tags (`v1.0`, `v1.1.0`, etc.) publish versioned images and update Docker `latest`.

## Path Filters (Implemented)

`ci-pr.yml`, `ci-dev-release.yml`, and `ci-master.yml` are path-scoped and trigger only when at least one of these changes:

- `.github/workflows/ci-pr.yml`
- `.github/workflows/ci-dev-release.yml`
- `.github/workflows/ci-master.yml`
- `.github/workflows/release.yml`
- `docker/**`
- `scripts/docker/**`
- `scripts/release/**`
- `core/**`
- `flex/**`
- `slam/**`
- `tests/**`
- `debian/**`
- `Makefile`
- `README.md`
- `docs/Github-actions.md`
- `docs/Star-binary-distribution.md`

`release.yml` is tag-driven (`v*`) and intentionally does not use path filters.

## Test Tiers

- Tier A:
  - self-contained smoke tests
  - public-data-only or synthetic surfaces
  - current Docker set: `tests/run_solo_smoke.sh`,
    `tests/slam/test_snp_mask_build_smoke.sh`,
    `tests/run_flex_tiny_public_smoke.sh`
  - should run on PRs and all protected branch pipelines
- Tier B:
  - fixture-backed tests (for example `/storage` data)
  - run on self-hosted/manual/scheduled workflows only unless fixtures are
    explicitly available in CI

## Container Publishing Policy

- Image name: `biodepot/star-suite`
- `dev-release` pushes:
  - optional `biodepot/star-suite:dev-release-<sha>`
  - optional `biodepot/star-suite:dev-release-latest`
- `master` pushes:
  - `biodepot/star-suite:latest`
  - `biodepot/star-suite:master-<sha>`
- `v*` tags:
  - `biodepot/star-suite:vX.Y.Z`
  - update `latest` only for stable tags

## Binary Distribution Policy

For Ubuntu discoverability and package updates:

- Build/sign source package from CI release workflow.
- Upload source package to Launchpad PPA.
- Let Launchpad produce architecture-native binaries (`amd64`, `arm64`).
- Also publish release artifacts and checksums in GitHub Releases.

## Workflow Structure (Recommended)

- `ci-pr.yml`
  - trigger: pull_request
  - checks: lint + `amd64` build + Tier A
- `ci-dev-release.yml`
  - trigger: push to `dev-release` and `dev-release-*`
  - checks: integration build + Tier A (+ optional lightweight extras)
  - publish: optional `dev-release-<sha>` and `dev-release-latest`
- `ci-master.yml`
  - trigger: push to `master`
  - checks: required checks
  - publish: multi-arch stable images
- `release.yml`
  - trigger: tags `v*`
  - checks: release gates
  - publish: versioned images, GitHub release artifacts/checksums, PPA source
    upload
  - prerelease logic:
    - `v0.*` and `-alpha|-beta|-rc|-dev` tags are published as prereleases
    - prereleases do not retag Docker `latest`
  - Linux tarball builds:
    - build with Dockerfiles under `scripts/release/docker/`
    - publish an `ubuntu:22.04`-based `glibc234` tarball for broader compatibility
    - publish an `ubuntu:24.04`-based `glibc239` tarball for current systems
    - publish an installer bundle that auto-selects the correct tarball for the host
  - includes:
    - Docker container runtime/version validation for release tarballs on Ubuntu 22.04 and 24.04
    - Docker container installer-bundle selection/version validation on Ubuntu 22.04 and 24.04
    - Docker container smoke execution against installed release artifacts on Ubuntu 24.04
    - runtime dependency manifests captured from the validation containers and uploaded with release artifacts
    - Debian `.deb` install/uninstall validation in clean `ubuntu:24.04` container
    - optional Debian source package signing when signing secrets are configured

## Implemented Workflow Files

- `.github/workflows/ci-pr.yml`
- `.github/workflows/ci-dev-release.yml`
- `.github/workflows/ci-master.yml`
- `.github/workflows/release.yml`

## Local Make Helpers

The top-level Makefile includes release-branch helpers. They create local
branches/tags only; pushing remains a human-controlled step.

```bash
# Show the dev-release workflow commands.
make dev-release-help

# Create local dev-release and version-scoped candidate branches.
make dev-release-branch RELEASE_VERSION=1.1.0
# Optionally set DEV_RELEASE_BASE=<branch-or-sha> to seed the shared branch.

# Run the local build/preflight checks expected before pushing a candidate.
make dev-release-check

# After committing on dev-release-v1.1.0, create an immutable RC tag.
make dev-release-tag RELEASE_VERSION=1.1.0 RC=1
```

## Branch Protection and Release Gates

- Require PR checks before merge to protected branches.
- Require `CI Dev Release` before promoting from `dev-release-vX.Y.Z` to
  `master`.
- Restrict stable release tags to commits reachable from `master`.
- Restrict prerelease tags (`vX.Y.Z-rcN`) to commits reachable from
  `dev-release` or `dev-release-vX.Y.Z`.
- Use immutable tags for release images (`vX.Y.Z`) and keep checksums for all
  published assets.

## Required Secrets (Release Workflows)

- Container registry credentials (for push).
- Debian source-signing key material:
  - `DEB_SIGNING_KEY` (ASCII-armored private key)
  - `DEB_SIGNING_KEY_ID` (fingerprint or key ID)
  - `DEB_SIGNING_KEY_PASSPHRASE` (optional, for encrypted key)
- Launchpad credentials or upload token (PPA source upload).

## Release Validation and Signing Commands

- Local install/uninstall validation of generated `.deb`:
  - `scripts/release/validate_deb_install_uninstall.sh --deb <path-to-deb>`
- Local signing of source package artifacts:
  - `scripts/release/sign_source_package.sh --in-dir dist/release/source --key-id <KEY_ID>`

## Notes

- This policy is intentionally conservative on publish events to control CI
  cost and registry churn.
- It can be tightened further (for example, only publish `latest` on tags) if
  release discipline becomes stricter.
