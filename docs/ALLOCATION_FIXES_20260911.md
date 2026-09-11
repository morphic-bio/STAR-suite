# Allocation fixes — 2026-09-11

Source audit: `/mnt/pikachu/STAR-suite/docs/AUDIT_ALLOCATION_PATTERNS_20260911.md`.
Development branch: `fix-larry-feature-regression-20260911`.
LARRY search corrections were committed first as `7ecd7c5`; the earlier scheduler
and PF correctness work is `e2bb279`. Full MSK paper P02 remains pending.

## Whitelist batch (audit items 2–4)

- CbCorrector retains packed lookup tables and whitelist size. The two Bayesian
  consumers use Solo's existing canonical strings; the corrector no longer copies
  millions of strings. Its tables remain independent of the constructor argument's lifetime.
- Default GEX normalization borrows the matching Solo packed hash. It verifies
  barcode length and the indexed canonical sequence. Other whitelist inputs keep
  their existing loader. A375 logs 3,686,400 reused keys and zero copied keys.
- The existing full whitelist scan records whether an output mapping can exist.
  MEX conversion receives no map path when the normalized whitelist is known to
  have none. Explicit maps and two-column translation mappings are retained.
  A mapping after the first row remains detectable; this is not a first-line guess.

Clean build and legacy CbCorrector differential ASan/UBSan tests pass, including
input destruction after construction. Namespace tests cover borrowed packed lookup,
wrong-length/N rejection, one/two/mixed-column metadata and overlapping namespaces.
The old metadata test linked unused CBQ dependencies; its harness now builds the
normalizer with section garbage collection. Failure logs are preserved.

Full A375, fixed four decoder workers, auto input detection, same inputs and policy:

| Binary | Whole wall | Peak RSS (KiB) |
|---|---:|---:|
| LARRY-corrected control | 146.22 s | 41,816,724 |
| Whitelist batch | 142.72 s | 41,559,480 |

All six matrices and three guide-call tables match exactly. The observed reduction
is 3.50 s (2.4%) and 251 MiB RSS; these are single executions, not replicated estimates.
Frozen STAR SHA256: `f28c9037fb02730b9492079baa8e0131135bba7dbbc4ebc23dfd24400573527d`.

Artifacts: `/home/lhhung/pf_larry_regression_20260911/allocation_audit/`:
`whitelist_build2`, `whitelist_cb_unit`, `whitelist_namespace_test2` (helper pass),
`whitelist_namespace_remaining` (metadata/overlap pass), `a375_whitelist`, and
`whitelist_comparison.json`. Earlier failed build/test harness attempts are separately
marked and never used as timing controls.

## UMI counter batch (audit item 1)

Single-feature barcode–UMI counters now occupy the pooled record, with a separate
visited flag. A general hash is allocated only upon observing a second feature.
Gather copies values into destination-owned storage, so source pools can be freed
independently. This removes the common per-UMI table while retaining the existing
hash representation for rare multi-feature cases; a shared rare-counter arena is
not implemented.

The differential fixture exercises cross-thread duplicates, UMI-neighbor connected
components, 64-way competing features, ties and count/stringency thresholds. Old
and new ledgers match exactly; AddressSanitizer also passes. The actual 200,000-pair
LARRY diagnostic retains every matrix entry (159,027 UMIs), with auto gzip detection.
Runtime is 15.97 s versus 16.57 s; 159,058 of 159,704 counter records use inline storage.

Full A375 preserves all six matrices and three guide tables:

| Metric | Whitelist control | UMI counter change |
|---|---:|---:|
| Whole wall | 142.72 s | 142.78 s |
| Peak RSS (KiB) | 41,559,480 | 40,010,100 |
| PF read/assign/join | 13.124 s | 9.984 s |
| PF thread hash merge/cleanup | 9.717 s | 1.775 s |
| PF sample cleanup | 2.164 s | 0.182 s |

Of 3,333,784 gathered barcode–UMIs, 3,332,717 are inline and only 1,067 need the
general table (99.968% inline). Peak RSS falls 1.48 GiB. The shortened feature work
overlaps GEX, so this run demonstrates no whole-job speedup. Do not add overlapping
phase savings to predict total runtime.

Artifacts: `umi_baseline`, `umi_probe_before`, `umi_probe_after`, `umi_probe_asan`,
`umi_build`, `umi_asan_build`, `umi_1000`, `umi_larry_comparison.json`, `a375_umi`,
`umi_a375_comparison.json`. Frozen STAR SHA256:
`75baa872a4b4f67d6b933ccc62195c5c4f5d72b32e06e47bb5605299bad21167`.
Reproduce isolated comparisons with `tests/run_umi_storage_probe.py`, supplying
preserved headers/libraries and `--expected-ledger`; controls need not be rerun.

## Bridge storage batch (audit items 7–8)

Immediate barcode read counters now use an owning khash with the existing packed
64-bit count values. Copy, move, reserve, clear, gather and snapshot semantics are
covered by sanitizer tests against a standard-map oracle. Final bridge rows use
one append buffer per worker and a stored slice per barcode, then copy into the
same final CSR order. No cell or gene ordering policy changes.

Full A375 reports 350,678 integer read-counter keys and 294,633 barcode slices in
32 worker buffers (15,134,378 output slots). All six matrices and three guide
CSV tables exactly match the UMI-counter control. Whole wall is 140.73 s versus
142.78 s; peak RSS is 39,941,536 KiB versus 40,010,100 KiB (67 MiB lower).
These are single-run observations. No snapshot format changed; its integer rows
may have a different hash iteration order, with the same key/value contents.

Artifacts: `bridge_build`, `bridge_unit`, `a375_bridge`, `bridge_a375_comparison.json`.
Frozen STAR SHA256: `9652d78f0b8a1ac88f926250fd38224c72f5ed2265437ccd4d6b2a5c3f300f44`.

## Optional inline correction (audit item 5)

InlineCBCorrection reuses its exact khash for N-path membership and stores
ambiguous variant ranges in khash with one shared parent-index array. A counting
pass and replay preserve the old per-variant candidate order, including duplicate
whitelist handling. The public ambiguity/evidence payloads are unchanged.

Before/after and ASan+UBSan ledgers match byte for byte across exact, H1, N,
ambiguous and random queries, a 48-parent collision, quality-based resolution and
shard merge. Ledger SHA256:
`4ef1b8e64659637fb4f22d1e752eef2679973aa290a952dde0fd6e05efaf5e32`.
A clean STAR build passes (SHA256
`64890391d9fc26841a18c7c8e16a5a2c0a79fcfce64f23fded64b38176cda834`).
This optional path is not enabled by the A375 control; no full-dataset speed or
memory claim is made for it. Existing collision-size/fanout diagnostics are retained.

Artifacts: `inline_baseline`, `inline_probe_before`, `inline_probe_after`,
`inline_probe_asan`, `inline_build`. Reproduction driver:
`tests/run_inline_cb_storage_probe.py` with a preserved `--expected-ledger`.

## Occupancy grouping (audit item 10)

Observed-tag occupancy packs each CB16+TAG8 into 48 bits, sorts/deduplicates the
keys and counts adjacent GEM groups. Only rejected GEMs become output strings.
The general Monte Carlo path sorts indices into its existing barcode strings,
replacing two nested maps with one array. It still handles non-DNA strings and
variable tag lengths. Removed indices now have deterministic ascending order;
callers consume them as a set. No sixteen-tag limit was introduced.

The Poisson calculation, fitted mean, lambda, percentile cutoff, simulations and
filtering rules are unchanged. Ninety observed-fit comparisons (including up to
80 tags per GEM), sixteen Monte Carlo/fallback comparisons, duplicates, invalid
barcodes/percentiles and variable-length tags match the preserved implementation.
ASan/UBSan and clean libscrna/core builds pass. Ledger SHA256:
`4f84865f447e38747b75747ba2fe464821b0305bf469ffc7b2dbab2cc99972ec`.
Existing observed-GEM diagnostics remain available. No full Flex performance
measurement was made for this storage-only change.

Artifacts: `occupancy_baseline`, `occupancy_before`, `occupancy_after`,
`occupancy_asan`, `occupancy_build`. Frozen STAR SHA256:
`be1800e803de55d20448a7d827113c93991e0c71e548b95c9bb63616466ed2bb`.
Reproduction: `tests/run_occupancy_storage_probe.py`.

## File-input MEX storage (audit item 11)

The file loader counts entries per cell, allocates one CSR data buffer, and
rewinds/fills it using cell cursors. It preserves each cell's original entry order,
duplicate coordinates, explicit zeros and existing header-count tolerance. It
checks 32-bit sparse-offset capacity and validates bounds/counts on replay.
The loader logs cells, actual stored entries and the single sparse buffer.

Before/after and ASan/UBSan fixtures cover unsorted/duplicate entries, empty input,
header/axis mismatches, zeros and out-of-bounds indices. The actual LARRY MEX also
matches exactly: 35,626 columns, 105,709 entries. Its complete loaded-data ledger
SHA256 is `ea49dc131765a5d4e8ee9a9881d4a0fb483cd21c1de3f1b003039a5c814d06b0`.
A clean core/libflex build passes; frozen STAR SHA256:
`8dc6529af4a64ef14430fc0d2ce5baea0aea7977e1f65fc0f4378036f232f55d`.
This file-input path is distinct from normal in-process mapping; no whole-job
speedup is inferred from it.

Artifacts: `mex_baseline`, `mex_inputs`, `mex_before`, `mex_after`, `mex_asan`,
`mex_real_inputs`, `mex_real_manifest.json`, `mex_real_compare`, `mex_build`.
Reproduction: `tests/run_flex_mex_storage_probe.py` with the preserved fixtures.

## Merge/import barcode indices (audit item 12)

The merge and table-import indices use khash keys that borrow pointer/length
views of existing strings. Strings stay stable for the entire lookup lifetime;
no packed-DNA assumption excludes suffixes, N bases or other valid identifiers.
Table import finishes vector growth, then compacts first occurrences in place
before storing views of the retained slots. This also handles short strings whose
characters live inside the string object. Output duplicate detection reuses the
existing barcode sort, with deterministic diagnostic ordering.

Sanitizer tests cover long/short/empty/NUL-containing keys, prefix queries,
duplicates, compaction, clear and capacity overflow. Table-import tests pass for
CSV/TSV, suffix normalization, rejections, duplicate pair collapse and stable
first-occurrence order. The test harness can retain all fixtures and results.

A clean build and full A375 pass. All six matrices and three guide CSVs match the
bridge control. Whole wall: 139.75 s versus 140.73 s; RSS: 39,911,908 KiB versus
39,941,536 KiB. These small differences are single-run observations. The merge
logs 3,686,400 keys and zero copied strings. The intervening optional-inline,
occupancy and file-MEX changes are also present in this binary; those conditional
paths are not exercised by this A375 configuration.

Artifacts: `merge_baseline`, `barcode_view_unit`, `merge_table_test`, `merge_build`,
`a375_merge`, `merge_a375_comparison.json`. Frozen STAR SHA256:
`74e4716f1dc2673ff21b59a70bcc0988d597f11503e02e2cf0a743d72fcbcab3`.

## Remaining work, in order

Review conditional items 6 and 9 separately with their relevant fixtures;
   preserve iteration-dependent outputs and do not infer Flex-scale savings from A375.

This file records completed changes separately from pending audit directions.
