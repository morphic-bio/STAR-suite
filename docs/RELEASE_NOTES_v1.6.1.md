# STAR Suite v1.6.1 Release Notes

Date: 2026-07-28

`v1.6.1` is the immutable replacement release for the `v1.6.0` packaging
defect. The public `v1.6.0` source tree was byte-for-byte identical to the
final private integration tree and already contained every scientific change
listed below. Its Git-free release builds, however, did not embed the source
revision consumed by spatial spill finalization. Consequently, the official
`v1.6.0` artifacts could pass all-memory 100K checks but could not finalize a
full-slide spill run.

The initial patch changed provenance, packaging, and release validation. The
final candidate also restores the native spatial Flex producer/consumer path
as the production default and reduces hard-only materialization bookkeeping.
Those execution changes retain the accepted hard-count and hard-MEX policy;
FASTQ and indexed-CBQ parity gates are required before release. The old
`v1.6.0` assets remain immutable and should not be used for integrated spatial
runs that can spill.

The release artifact version is `v1.6.1`, Debian packages use `1.6.1-1`, and
both `STAR --version` and `molecule_first_resolver --version` report `1.6.1`.
The upstream STAR base remains `2.7.11b`; genome-index compatibility remains
`2.7.4a`.

## Packaging and immutable source identity

- Embed the exact full 40-character release commit in STAR builds made from
  Git-free static-tarball, installer, Debian/source-package, and Docker
  contexts.
- Add `STAR --source-revision` and record the same value in STAR logs, spatial
  sidecars, spill records, run summaries, and release metadata.
- Reject integrated spatial and sidecar runs before opening FASTQs when the
  binary lacks an exact immutable source revision. Release construction also
  rejects missing, short, sentinel, or malformed revisions.
- Assert equality among the expected tag commit, packaged metadata, installed
  binary, and runtime manifest for both glibc compatibility tarballs, the
  installer bundle, Debian packages, and the Docker build.
- Publish this curated document as the GitHub release body. The automatic
  `v1.6.0` body listed only two pull requests and omitted most direct and
  private-folded work.

## Integrated Visium HD 3-prime GEX

- Enable the path explicitly with `--soloSpatialGexIntegrated yes` and the
  spatial barcode/coordinate contract; it remains off by default.
- Decode R1 spatial candidates while STAR maps R2, then join them to
  post-alignment GeneFull evidence before barcode correction and UMI collapse.
- Retain finite coordinate candidates through read-clique formation and
  candidate-specific 1MM_CR correction, then materialize strict, soft
  expected-count, hard, and gated-hard matrices at 2, 8, and 16 micrometers.
- Route N-containing spatial barcodes through bounded alignment/DP fallback.
  Invalid UMIs remain excluded from molecule creation but present in audit
  accounting.
- Keep the feature sidecar optional as a diagnostic surface rather than a
  production interchange requirement.

## Full-slide bounded execution

- Spill candidate accumulation when configured and use the bounded downstream
  contribution spool for clique formation, coordinate-local correction,
  multigene reconciliation, and matrix materialization.
- Preserve deterministic record formats, ordering, integrity checks,
  admission guards, and all-memory/spill parity.
- Stabilize coarse soft-matrix reduction order so thread, spill, and
  materialization choices do not change released products.

## Gene-evidence policy

- Keep `compatibility` as the shipped/default CR-compatible exon-first rescue
  policy when CR multimap rescue is enabled.
- Retain `annotated` as an explicit deterministic, score-first retained-GTF
  alternative. Under that policy, unannotated equal-best alignments remain
  decoy evidence instead of silently promoting an annotated alternative.
- Leave ordinary STARsolo behavior unchanged unless CR multimap rescue is
  explicitly enabled.

## Native Visium HD Flex

- Enable the native path with `--soloSpatialFlexIntegrated yes`; products are
  written under `SpatialFlex.out`.
- Use `--flexPipeline auto` in production. FASTQ inputs select fused per-lane
  reader/routers with separate spatial cache-hit consumers and alignment-miss
  workers; indexed CBQ inputs select fully fused range readers.
- Retain `--flexPipeline no` only as an explicitly authorized diagnostic
  compatibility mode. The earlier native-spatial requirement for this setting
  was a temporary debugging restriction and was not a production policy.
- Resolve direct H0/H1 probe hits and alignment fallback into retained spatial
  candidate families, then delay candidate-specific UMI correction until
  molecule resolution.
- Emit the same strict, soft expected-count, hard, and gated-hard policies at
  2, 8, and 16 micrometers through the same bounded spill engine used by GEX.
- Require no BAM, GX/UR text bridge, external decoder/resolver, or external
  materializer.
- Use liberal `--soloMapqMode off` behavior by default for ordinary and spatial
  Flex; explicit MAPQ filtering remains available.
- Honor explicit `--soloRunFlexFilter no`, resolve R2 features before CB/sample
  eligibility, and report exact H0/H1/alignment terminal-route accounting.

## Flex skip, fallback, and gDNA fixes

- Do not require expected-cell settings under `--soloSkipProcessing yes`, and
  avoid legacy Solo replay in inline skip-processing mode.
- Route true no-sample cache misses to alignment while retaining and accounting
  valid hash hits that lack a barcode assignment.
- Compute the 10x-style Flex gDNA diagnostic from final filtered
  barcode/gene/UMI families after barcode and UMI correction, without requiring
  BAM output.
- Carry spliced/unspliced probe-region evidence through H0/H1 cache keeps,
  alignment fallback, ambiguous-barcode resolution, and UMI re-keying without
  changing molecule keys or allocating per-read sidecars.
- Add FH01SEQ1 cache v3 without changing the 24-byte record size. STAR and the
  replay utilities retain v1/v2 read compatibility; older caches map normally
  but report the diagnostic unavailable.
- Emit per-sample Cell Ranger-compatible JSON metrics and a library/per-sample
  audit TSV. `--soloFlexGdna auto|yes|no` is matrix-inert and inert outside
  Flex.

## Molecule-first and production-adapter improvements since v1.5

- Stabilize tolerance-level ties, stream feature-sorted inputs, and bound
  Visium HD materialization memory.
- Fix H1 matching for complex barcode whitelists and support `oneExact` gating
  in bridge collapse.
- Package `molecule_first_resolver`, `molecule_first_bam_ledger`, and
  `molecule_first_materialize` together with the post-collapse production
  adapters.

## Validation contract

- The public `v1.6.0` tag tree and final private integration tree have the same
  Git tree identity; the v1.6.1 scientific inclusion audit therefore starts
  from the complete folded source, not a partial re-integration.
- Focused spatial, downstream-spool, resolver, Flex hash/fallback, Flex gDNA,
  CR rescue, complex-whitelist, default-off isolation, and release rails must
  pass from a clean candidate checkout.
- Installed packaged binaries must force and complete both integrated GEX and
  native Flex spill/materialization paths, with their summaries carrying the
  exact release revision.
- An independent clean checkout rebuild verifies version, source revision,
  retained feature inventory, runtime behavior, and artifact manifests before
  the stable tag is published.

Implementation and safety contracts remain documented in
`docs/RUNBOOK_VISIUM_HD_GEX_FEATURE_SIDECAR_20260722.md`,
`docs/RUNBOOK_VISIUM_HD_GEX_IN_MEMORY_1MM_CR_20260724.md`,
`docs/RUNBOOK_VISIUM_HD_GEX_DOWNSTREAM_SPOOL_20260725.md`,
`docs/RUNBOOK_VISIUM_HD_NATIVE_FLEX_SPATIAL_FAMILIES_20260726.md`,
`docs/RUNBOOK_FLEX_GDNA_DIAGNOSTIC_20260727.md`, and
`docs/RUNBOOK_POST_COLLAPSE_PRODUCTION_INTEGRATION_20260717.md`.
