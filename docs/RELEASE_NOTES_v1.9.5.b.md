# STAR Suite v1.9.5.b

Release candidate: portable Launchpad and native Multiome execution.

- Launchpad uses checkout/config-relative locations and explicit user data and
  reference paths. Workstation dataset defaults are removed from shipped forms.
- The bundled 10x Multiome recipe has explicit RNA R1/R2 and ATAC R1/R2/R3 roles,
  lane validation, managed execution, logs, cancellation, and run records with
  executable hashes. Namespaced official recipe routes work in the UI.
- `STAR --build-features` reports whether the executable includes Chromap.
  Multiome checks this before starting and rejects incompatible binaries.
- Tarballs and Debian packages include the UI/server and `star-suite-launchpad`.
  `--setup` installs Python dependencies into a user-owned environment. Browser
  assets are bundled for use without a CDN.
- Release artifact builds depend on relocation, input-validation, process-control
  and Chromium browser tests, alongside the existing build gates.

The portable STAR binaries still omit Chromap. Multiome requires a
Chromap-enabled STAR source build plus `star_multiome_atac_peak_mex`; downloading
the standalone Chromap executable does not enable ATAC in an existing STAR.
See [the Multiome Launchpad guide](LAUNCHPAD_MULTIOME.md) for build and launch
commands. This maintenance release does not change Multiomics Suite ownership
or claim new biological parity results.
