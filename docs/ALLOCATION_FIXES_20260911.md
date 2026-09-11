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

## Remaining work, in order

1. Audit item 1: remove the per-UMI table from the single-feature case; retain general
   counters for actual multi-feature UMIs. Validate gather, connected components,
   ties and source cleanup against the preserved library before benchmarking.
2. Items 7–8: bounded integer read counters and flat bridge output rows.
3. Review conditional items 5–6 and 9–12 separately with their relevant fixtures;
   preserve iteration-dependent outputs and do not infer Flex-scale savings from A375.

This file records completed changes separately from pending audit directions.
