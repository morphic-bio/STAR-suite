# STAR Suite v1.12.2 Release Notes

Date: 2026-10-27

This patch makes spatial R1 barcode decoding faster without changing its
results. It adds commit `2c7159e` ("Reuse invariant work in spatial barcode
decoding") to v1.12.1.

The decoder previously repeated two calculations that do not change within a
loop:

- The N-base target scan prepared an edit-distance bit pattern for every
  candidate target. It now prepares the observed span once per scan and
  compares it with each target. Levenshtein distance is symmetric, and the
  64-base truncation, length filter and literal treatment of `N` are unchanged.
- The raw-candidate alignment likelihood computed the base-quality error terms
  in every dynamic-program cell. They depend only on the observed base and its
  quality, so they are now computed once per row. Each cell keeps the original
  operations and order, so likelihoods, including ties, are bit-identical.

Candidate enumeration, thresholds, ordering, policies and output manifests are
unchanged. No other module changes.

On full Visium HD slides with matched builds on one host, decoding a
colorectal slide took 38.7% less wall time and an ovarian slide 27.7% less
(one matched pair each). Gains elsewhere depend on how much of a run is spent
decoding barcodes.

Validation:

- A new regression, `make -C core/legacy/source test_SpatialR1DecoderHotpaths`,
  is part of `make test`. It checks 151,072 likelihoods bit-for-bit against the
  previous recurrence, 12,675 edit distances against an independent scalar
  oracle, and 800 target-scan and cache cases. It also passes under UBSan.
- All core test targets, the synthetic Solo smoke, all eight scRNA
  deduplication modes and the spatial half-probe suite pass. The half-probe
  suite's five reader, thread and spill variants agree in all 36 MEX files.
- On the full colorectal and ovarian slides, all 36 MEX components of each
  run match the frozen v1.9.5 reference checksums byte for byte, as the
  matched v1.12.1 control does. All retained summary fields agree with the
  control, apart from source revision, memory budget and timers.

`STAR --version` reports `1.12.2`. Debian source packaging uses `1.12.2-1`;
Ubuntu packages use `1.12.2-1~ubuntu22.04.1` and `1.12.2-1~ubuntu24.04.1`.
Upstream STAR remains `2.7.11b`, genome-index compatibility remains `2.7.4a`
and legacy compatibility remains `2.7.1a`. Existing indexes need no rebuild.
