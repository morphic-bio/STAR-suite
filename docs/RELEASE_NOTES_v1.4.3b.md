# STAR Suite v1.4.3b Release Notes

Date: 2026-08-09

`v1.4.3b` is a paper-era maintenance release. It preserves the `v1.4.3`
source and benchmark line while backporting the scoped fix for ordinary Fastx
STARsolo thread scaling and startup memory.

The paper benchmarks did not use the affected plain-STARsolo path. Their
reported results therefore remain attributed to `v1.4.3` and are not replaced
or restated by this release. This maintenance tag exists so reviewers testing
STAR Suite as a general drop-in STARsolo executable receive the corrected
behavior and can see the deviation from the paper tag explicitly.

The release artifact version is `v1.4.3b`, Debian packages use `1.4.3b-1`, and
`STAR --version` reports `1.4.3b`. The upstream STAR base remains `2.7.11b`;
genome-index compatibility remains `2.7.4a`.

## Fixes

- Restore STAR's direct chunk reader for ordinary `--readFilesType Fastx`
  input. The generic record adapter had parsed and reconstructed each record
  while STAR held the global input mutex, flattening mapping throughput as
  `--runThreadN` increased.
- Construct the expanded one-mismatch `CbCorrector` only for inline modes that
  require it. Plain STARsolo uses its established whitelist-matching path and
  did not consume that multi-gigabyte precompute.
- Preserve Fastx contract validation and plain/gzip, explicit decompression,
  multi-lane, manifest, single-end, and Y/no-Y behavior.
- This release does not change the separate Visium/spatial barcode decoders or
  introduce persistent correction-table caches.

## Validation

- Clean-build Fastx contract, standalone parser, STARsolo, explicit inline
  correction, and single-end Y/no-Y tests passed.
- Plain STARsolo raw and filtered GeneFull matrices were byte-identical before
  and after the fix and across tested thread counts.
- In upstream-compatible mode, fixed outputs were byte-identical to official
  upstream STAR 2.7.11b.
- The full diagnosis, benchmark design, timing table, and artifacts are in
  `docs/SOLO_FASTX_SCALING_FIX_20260809.md`.

The same fix is included in the current production line beginning with STAR
Suite `v1.7.1`.
