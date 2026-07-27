# STAR Suite v1.6.0 Release Notes

Date: 2026-07-25

`v1.6.0` promotes the opt-in integrated Visium HD 3-prime GEX pipeline and
its bounded molecule-first materializer. It also fixes Flex-as-aligner paths
that skip Solo processing or fall back from the feature hash to alignment.
The release artifact version is `v1.6.0`, Debian packages use `1.6.0-1`, and
both `STAR --version` and `molecule_first_resolver --version` report `1.6.0`.
The upstream STAR base remains `2.7.11b`; genome-index compatibility remains
`2.7.4a`.

## Integrated Visium HD GEX

- Enable the new path explicitly with `--soloSpatialGexIntegrated yes` and a
  spatial barcode/coordinate contract. It is off by default.
- Decode R1 spatial coordinate candidates while STAR maps R2, then join those
  candidates to post-alignment GeneFull evidence before barcode correction and
  UMI collapse.
- Resolve read cliques and candidate-specific UMIs before materializing
  strict, soft expected-count, hard, and gated-hard spatial products at 2, 8,
  and 16 micrometer scales.
- Route spatial barcodes containing `N` through the bounded alignment/DP
  fallback under the ordinary uniqueness rules. UMIs containing non-ACGT
  bases remain invalid for molecule creation and are retained in accounting.
- Keep the feature sidecar optional for diagnostics or overflow handling; it
  is not the production interchange requirement.

## Memory-bounded full-slide processing

- Spill candidate accumulation when configured, then use a bounded downstream
  spool for clique formation, coordinate-local correction, multigene
  reconciliation, and matrix materialization.
- Preserve exact in-memory/spooled parity with versioned binary records,
  deterministic ordering, integrity checks, and fail-closed admission gates.
- Stabilize floating-point reduction order for coarse soft matrices so thread,
  spill, and materialization choices do not change released products.

## Gene evidence and compatibility policy

- Make annotated CR-compatible rescue deterministic and score-first while
  retaining unannotated alignments as decoy evidence rather than silently
  promoting them.
- Keep `annotated` as the normal evidence policy. The `compatibility` evidence
  mode is explicit and experimental; it is not implied by GX/UR tags and does
  not change ordinary STARsolo defaults.
- Preserve the accepted deterministic annotation/MECOM rescue behavior already
  present on the stable line.

## Flex behavior and fixes

- Use the liberal `--soloMapqMode off` policy by default in ordinary and
  spatial Flex. MAPQ filtering remains available as an explicit request.
- Resolve hash-miss R2 feature evidence before CB/sample eligibility, then feed
  the chosen gene and raw UMI into the existing deferred family correction.
  This keeps feature assignment independent of whether R1 later forms an
  eligible single-cell or spatial family.
- Do not require expected-cell settings when `--soloSkipProcessing yes`.
- Avoid legacy Solo replay in inline skip-processing mode.
- Route true no-sample hash-cache misses to alignment while keeping valid H0
  hash hits, including reads without a barcode assignment, on the hash path.
- Account hash-kept reads even when no barcode assignment exists.

## Isolation and validation

- Integrated spatial GEX is gated by its explicit enable flag and contract;
  sidecar and spill settings are inert on ordinary scRNA runs.
- Normal scRNA sidecar-off outputs reproduced their tracked byte-level golden
  files. Solo, GEX, CR multi, Flex, CB/UB, Y/no-Y, SLAM, TranscriptVB, and
  tximport regression rails passed or cleanly reported unavailable optional
  dependencies.
- Core spatial, spill, downstream-spool, CR multimap, 20,000-case UMI parity,
  molecule-first resolver/materializer, N-barcode, and concurrency tests
  passed.
- The fresh human CRC 100K gate reproduced 109,852 candidate rows, 79,644
  cliques, 66,845 strict molecules, 79,144 hard molecules, and 71,262
  gated-hard molecules, with exact accepted ledgers, resolver artifacts, and
  policy matrices.
- The final MAPQ-off CRC 100K gate assigned 3,875 of 6,359 hash misses. Its
  ordinary and native-spatial arms covered all 6,359 reads with byte-equivalent
  terminal resolver fields and zero per-read decision mismatches.
- The fresh ovarian 100K compatibility/spool gate reproduced 58,395 strict,
  71,830 hard, and 63,105 gated-hard molecules; all 36 MEX components were
  byte-identical to the accepted compatibility oracle.
- The retained full 212.6-million-read CRC validation reproduced the open
  candidate surface and all eight resolver artifacts byte-for-byte while
  reducing materialization to 10 minutes 49 seconds.
- The CR-compatible CRISPR fixture passed with 1,043 single-feature cells and
  35 multi-feature cells; pinned SLAM parity passed after manifest/index
  verification.

Implementation, safety contracts, commands, and artifact provenance are in
`docs/RUNBOOK_VISIUM_HD_GEX_IN_MEMORY_1MM_CR_20260724.md`,
`docs/RUNBOOK_VISIUM_HD_GEX_DOWNSTREAM_SPOOL_20260725.md`, and the associated
validation records under `docs/`.
