# STAR Suite v1.0.0 Release Notes

Date: 2026-05-22

`v1.0.0` is the first production STAR Suite release. The release artifact
version is `1.0.0`, and `STAR --version` now reports the suite version
`1.0.0`. The upstream STAR base is available through `STAR --upstream-version`
(`2.7.11b`), and the genome index compatibility string is available through
`STAR --genome-compat-version` (`2.7.4a`).

## Production Surfaces

- CR-compatible perturb-seq production path for UCSF/MSK-style GEX + feature
  libraries, including native Velocyto `outs/` MEX packaging and adaptive
  downstream MT/n_genes QC.
- Native OCM scRNA-seq production path using
  `--ocmMultiBarcodeMode flex`, where STAR counts on the effective
  `CB16+OCM_TAG8` barcode before correction, UMI collapse, EmptyDrops, and
  Velocyto.
- Cell Ranger multi-compatible OCM materialization under `outs/multi` and
  `outs/per_sample_outs`, plus downstream per-sample mirrors under
  `samples/<sample>/run/outs`.
- JAX scRNAseq02 OCM recipes for smoke, full production, Y/noY side outputs,
  remote post-MEX downstream processing, and Globus large-file handoff.
- MSK 40KO production wrapper and manifest for the February-2018 TRU/NXT
  GEX/PolyIII/LARRY layout.
- Multiome production docs for the preferred local STAR/Chromap MEX boundary
  followed by remote post-MEX h5ad/h5mu construction.

## Velocyto

- `VelocytoMexWriter` is the current production packaging path. New production
  runs should not invoke `scripts/prepare_velocyto_mex.py`; that script remains
  a legacy repair/backfill helper for old STAR outputs.
- OCM per-sample Velocyto materialization maps GeneFull columns to Velocyto
  columns by normalized barcode key, avoiding positional assumptions when
  barcode order differs between matrices.
- The delivered velocity surface is still the legacy STAR Velocyto policy.
  The count-preserving `CRCompat` policy is documented as a future
  implementation target and is not enabled in this release.

## Documentation And Recipes

- Added and updated production runbooks for OCM, JAX scRNAseq02, MSK 40KO, and
  STAR Suite binary distribution.
- Updated `scripts/README.md` and top-level README entries for OCM and native
  Velocyto production behavior.
- Added MCP/Launchpad recipe schemas for JAX OCM smoke and MSK 40KO production
  wrapper usage.

## Provenance

- Record the final tagged commit hash in the GitHub release body and checksum
  manifest when `v1.0.0` is cut.
- Release CI should build compatibility tarballs, installer bundle, Debian
  packages, source packages, checksums, runtime manifests, and container images
  from the final tagged `v1.0.0` commit.
