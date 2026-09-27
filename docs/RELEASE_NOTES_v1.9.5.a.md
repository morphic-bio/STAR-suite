# STAR Suite v1.9.5.a Release Notes

This maintenance release stays in the 1.9.5 series. It does not replace or move
the original `v1.9.5` tag, and does not relabel or rerun the paper benchmarks.
`STAR --version` reports `1.9.5.a`; Debian source packaging uses `1.9.5.a-1`.
The upstream STAR version remains `2.7.11b` and genome-index compatibility
remains `2.7.4a` (legacy `2.7.1a`). Existing indexes do not need rebuilding.

## Source Compilation

- Discover matching external HTSlib headers and libraries with `pkg-config`,
  including installations outside system include paths. Dependency generation
  now honors `CPPFLAGS` and `CXXFLAGSextra` and replaces `Depend.list` atomically.
- Diagnose missing external HTSlib or Chromap prerequisites before scanning
  the full source tree, instead of failing with an unexplained missing
  `htslib/khash.h`. Never mix bundled old HTSlib headers with an external ABI.
- `make core` remains Chromap-enabled. Set `CHROMAP_SUITE_DIR` when Chromap
  Suite is not a sibling checkout; install its prerequisites as documented in
  [compilation instructions](https://github.com/morphic-bio/STAR-suite/blob/v1.9.5.a/docs/compile_instructions.md).
  Explicit `make core-portable` uses bundled HTSlib for RNA/Flex/SLAM, including
  STARsolo poly-G trimming, without Chromap.
- Correct installer instructions: the installer archive contains precompiled
  binaries and `install.sh`; it is not a source package and needs no `make`.

## Conventional Solo Counts

Restore returned molecule counts to the output matrix for `1MM_All`,
`1MM_CR`, `1MM_Directional`, and `1MM_Directional_UMItools` on the conventional
counting path. Previously, alignment could succeed while these counts were
discarded. The default `1MM_All` configuration was affected.

The production `1MM_CR` + `MultiGeneUMI_CR` counting path is separate and is
unchanged. The PBMC paper benchmark used that CR-compatible configuration,
not the affected defaults. `Exact` and `NoDedup` already retained their counts.
If an affected conventional configuration produced empty matrices, regenerate
its counts with the fixed release; this fix cannot recover discarded counts
from those matrices alone.

## Withdrawn Flex Diagnostic

The optional Flex gDNA estimator, its metadata collection, QC reports and
source-derived numerical tests have been removed. Its recorded implementation
provenance did not meet this project's documentation-only clean-room boundary.
There is no replacement estimator in this release.

Expression counting and cell-calling logic are unchanged by this removal.
Existing cache formats remain readable; the retained packed-region utilities
only preserve STAR cache/count compatibility. No gDNA JSON/TSV is emitted.
Legacy `--soloFlexGdna auto/no` values are compatibility no-ops; `yes` and
explicit `--soloFlexGdnaProbeSet` paths now fail with a removal message.

The original, unfinished `v1.9.5.a` tag was repointed with the repository
owner's explicit approval after its publication was canceled. Users who
fetched the earlier tag must refresh it; verify `STAR --source-revision`
against the final release commit. The original `v1.9.5` tag is unchanged.

## Regression Coverage

- Removal guards run in Tier A and release-artifact smokes. Cache-format and
  packed-count tests preserve legacy decoding. On the eight-lane JAX 100K
  fixture, all 121 non-diagnostic outputs are byte-identical to the pre-removal
  candidate; new runs emit no gDNA reports.

- The existing Solo smoke requires three primary alignments and exactly two
  molecules at a known gene/barcode coordinate. It rejects the released
  v1.9.5 default-Solo output, which had an empty matrix, and passes the fix.
- An analytical GEX-only test checks all supported UMI methods, raw/filtered
  Gene and GeneFull counts, known cell identities, intronic reads and poly-G
  trimming. It runs in Tier A and release-artifact smoke checks.
- Nine HTSlib build tests exercise dependency discovery, an indirect khash
  include compiled to an object, missing headers/libraries, custom include
  flags, portable bundled headers and atomic dependency-file generation.
- Clean, independent partial-Make targets gate PR, development, master and
  release workflows. The full Chromap build is also available as an explicit
  local acceptance test with an external checkout.
- The fixture-backed PBMC 100K gate tests default Solo and CR-compatible
  GEX-only profiles with/without BAM. Existing perturb coverage is retained.

Before the version bump, the fixed 100K PBMC default profile produced 47,358
Gene UMIs and 70,527 GeneFull UMIs, versus zero before the fix. CR-compatible
counts remained exactly unchanged at 47,203 and 70,270, respectively, with
identical cell sets and with/without-BAM parity. These are shallow regression
fixtures, not revised biological cell estimates or paper timing benchmarks.
See the [regression runbook](https://github.com/morphic-bio/STAR-suite/blob/v1.9.5.a/docs/RUNBOOK_SCRNA_GEX_100K_REGRESSION.md).
