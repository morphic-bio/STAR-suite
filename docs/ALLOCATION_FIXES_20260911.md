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

## Remaining work, in order

Review conditional items 6 and 9–12 separately with their relevant fixtures;
   preserve iteration-dependent outputs and do not infer Flex-scale savings from A375.

This file records completed changes separately from pending audit directions.
