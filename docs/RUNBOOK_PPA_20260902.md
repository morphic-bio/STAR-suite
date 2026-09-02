# Runbook: Launchpad PPA for STAR Suite

Date: 2026-09-02. Goal: `sudo add-apt-repository ppa:<owner>/star-suite &&
sudo apt install star-suite` works on stock Ubuntu — and becomes the opening
line of the cloud 320k reproducibility demonstration. Launchpad builds the
package from source on Canonical's builders, which is itself a third-party
attestation that the suite builds clean on stock Ubuntu.

## One-time setup (operator-side; ~30 min)

1. Launchpad account (an organizational owner such as `biodepot` keeps the PPA
   address institution-named).
2. Generate or choose a GPG key; push it to `keyserver.ubuntu.com`; add the
   fingerprint to the Launchpad account and complete its email verification.
3. Create the PPA on the account page: name `star-suite`, default settings
   (amd64 builds are on by default).

## Per-release upload (v1.8.1 example; run AFTER the release tag exists)

From a clean checkout of the release tag on an Ubuntu machine with
`devscripts`, `debhelper`, and `dput` installed:

```bash
git checkout v1.8.1
# Re-version for the PPA and target the series: noble (24.04).
# Version sorts BELOW a future official 1.8.1-1 upload, which is intended.
dch --newversion 1.8.1-0ppa1~noble1 --distribution noble \
    "PPA build of the v1.8.1 release."
# Build the SOURCE package with the repo's own script (runs dpkg-buildpackage -S -sa):
scripts/release/build_source_package.sh --out-dir ../ppa-out
# Sign with the Launchpad-registered key, then upload:
debsign -k <FINGERPRINT> ../ppa-out/star-suite_1.8.1-0ppa1~noble1_source.changes
dput ppa:<owner>/star-suite ../ppa-out/star-suite_1.8.1-0ppa1~noble1_source.changes
```

Launchpad emails an acceptance, builds within the hour, and publishes ~15 min
after a green build. For a 22.04 build as well, repeat with
`--newversion 1.8.1-0ppa1~jammy1 --distribution jammy` (the glibc234 release
variant demonstrates the codebase already builds there), or use Launchpad's
"Copy packages" with rebuild.

## Verify (clean container)

```bash
docker run --rm -it ubuntu:24.04 bash -c '
  apt-get update -qq && apt-get install -y -qq software-properties-common
  add-apt-repository -y ppa:<owner>/star-suite && apt-get install -y star-suite
  STAR --version'
```

## Cloud demo wiring

`docs/benchmarks/cloud_320k/run_cloud_320k_matrix.sh` accepts
`INSTALL_MODE=ppa PPA=ppa:<owner>/star-suite` and then installs the suite via
`add-apt-repository` + `apt install star-suite` instead of downloading a
release asset. The manuscript's worked example quotes whichever mode actually
ran; the PPA mode is the stronger sentence.

## Notes

- Launchpad requires source-only uploads; the CI-built `.deb` cannot be
  uploaded directly. The repo's `build_source_package.sh` produces exactly the
  right artifacts.
- If the Launchpad build fails where CI's container checks passed, the usual
  cause is a missing entry in `Build-Depends` (Launchpad chroots are more
  minimal than the check containers). Fix in `debian/control`, bump `~noble2`.
- `debian/copyright` is currently thin; fine for a PPA, but expand per-component
  (vendored htslib, opal, PCG, khash, kseq) before any Debian Med submission.
