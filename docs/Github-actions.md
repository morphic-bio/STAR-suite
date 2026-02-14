# GitHub Actions CI/CD Policy

This document defines when STAR-suite CI should run, what it should execute,
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
| Push | `dev` | integration build (`amd64`) | Tier A + selected fixture-free checks | Optional `dev-<sha>` image |
| Push | `master` | full multi-arch build (`amd64`,`arm64`) | required checks | Yes (`latest`, `master-<sha>`) |
| Tag Push | `v*` | release multi-arch build | release gates | Yes (version tags + release assets + PPA source upload) |
| Scheduled (optional) | nightly | drift detection build | extended checks | Optional `nightly` tag |

## Test Tiers

- Tier A:
  - self-contained smoke tests
  - should run on PRs and all protected branch pipelines
- Tier B:
  - fixture-backed tests (for example `/storage` data)
  - run on self-hosted/manual/scheduled workflows only unless fixtures are
    explicitly available in CI

## Container Publishing Policy

- Image name: `biodepot/star-suite`
- `dev` pushes:
  - optional `biodepot/star-suite:dev-<sha>`
- `master` pushes:
  - `biodepot/star-suite:latest`
  - `biodepot/star-suite:master-<sha>`
- `v*` tags:
  - `biodepot/star-suite:vX.Y.Z`
  - optionally update `latest` from tagged release only

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
- `ci-dev.yml`
  - trigger: push to `dev`
  - checks: integration build + Tier A (+ optional lightweight extras)
  - publish: optional `dev-<sha>`
- `ci-master.yml`
  - trigger: push to `master`
  - checks: required checks
  - publish: multi-arch stable images
- `release.yml`
  - trigger: tags `v*`
  - checks: release gates
  - publish: versioned images, GitHub release artifacts/checksums, PPA source
    upload
  - includes:
    - Debian `.deb` install/uninstall validation in clean `ubuntu:24.04` container
    - optional Debian source package signing when signing secrets are configured

## Implemented Workflow Files

- `.github/workflows/ci-pr.yml`
- `.github/workflows/ci-dev.yml`
- `.github/workflows/ci-master.yml`
- `.github/workflows/release.yml`

## Branch Protection and Release Gates

- Require PR checks before merge to protected branches.
- Restrict release workflow to tags created from `master`.
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
