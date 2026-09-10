# STAR Suite v1.9.1 Release Notes

Date: 2026-09-10

STAR Suite 1.9.1 brings the recent execution and memory improvements to
scRNA-seq, perturb-seq, bulk RNA-seq, and SLAM-seq. It also removes repeated
OCM matrix parsing by routing counts once and calling samples directly from
compact count records. It includes the merged cross-module and OCM changes
through `a04b8a58b45ec522107bf1428b4433c80bbc5cf7`.

`STAR --version` reports `1.9.1`. Debian source packaging uses `1.9.1-1`;
Ubuntu packages use `1.9.1-1~ubuntu22.04.1` and
`1.9.1-1~ubuntu24.04.1`. Upstream STAR remains `2.7.11b`, genome-index
compatibility remains `2.7.4a`, and legacy compatibility remains `2.7.1a`.
Existing indexes do not need rebuilding.

## Input and feature assignment

- Native BGZF feeds ordered windows through the established STAR parser,
  extending support beyond fused FLEX to ordinary alignment, scRNA-seq and
  SLAM input paths. Pairing, qualities, lane transitions, read limits and
  supported output modes retain their existing semantics.
- `process_features` supports native BGZF with independent mate block
  boundaries and leased record batches. Eligible direct workers avoid
  repeated copying; chemistry detection, three-stream input and
  order-sensitive rescue retain the queued path.
- Decoders, mappers and feature workers share bounded CPU permits. Corrupt
  input, mismatched mates, truncated records and worker failures are checked.
- Core input uses `--readFilesBgzfMode`; feature input uses
  `--crAssignBgzfMode`. An explicit `--readFilesCommand` retains precedence.

## Shared cell calling and OCM

- Shared scRNA/perturb EmptyDrops reuses the actual OrdMag result for
  diagnostics, borrows split or strided sparse arrays, and uses the available
  post-mapping thread budget for Monte Carlo simulation. Logical random
  streams, model settings and call decisions are preserved in accepted tests.
- OCM admits independent sample preparations within shared CPU and estimated
  matrix-memory budgets. Oversized samples run alone.
- Native OCM calling scans the pooled raw MEX once, routes compact count
  records to samples and tag unions, prepares caller arrays directly, and
  writes required sample matrices after calling. Generated sample MEX files
  are no longer read back into the caller.
- Count retention uses at most one quarter of `--ocmCellCallMaxMemory`
  (default 1 GiB); excess records spill through buffered binary storage.
  The setting is an admission estimate, not a hard process RSS limit.
- When filtered calls already exist, OCM scans each raw/filtered source once
  to route all sample outputs. Required directory, matrix and downstream
  contracts remain intact. This path still uses temporary text bodies.
- Axis errors, malformed/truncated counts, overflow and gzip failures are
  rejected. Coordinator-owned directory creation avoids concurrent logging
  through STAR's shared stream.

## SLAM and bulk quantification

- SLAM caches invariant probability tables, reuses fitted gene results
  across output writers, and fits independent genes within the thread
  budget. Histogram, rate or model changes invalidate cached results.
  Standalone `slam_requant` also supports parallel fitting.
- `--quantVBComponentParallel 1` enables component execution without the
  large accumulation buffer. It remains opt-in and preserves global
  convergence and GC-update scheduling.
- **Native BGZF is disabled for TranscriptVB online model learning.**
  Automatic input selection falls back to the established reader; forced
  range mode rejects the configuration. Parallel reader scheduling changed
  online fragment-length/EC evidence during validation. General bulk
  alignment/GeneCounts BGZF support remains available.

## Validation and measured effects

These are single, serialized **pre-version-bump implementation measurements**,
with workload-specific controls. They are not new full-workflow paper
benchmarks or measurements of the final packaged 1.9.1 binaries.

| Workload | Control → changed implementation | Interpretation |
| --- | --- | --- |
| Full A375 feature library, 9,748,584 pairs | 28.2722 → 18.0532 s | 36.1% less wall time; exact counts. Feature library only, not GEX + guide end-to-end. |
| Saved A375 GEX caller, 100K simulations, MC1 → MC8 | 24.0404 → 3.68962 s; RSS 380,400 → 296,876 KiB | Exact 1,188 calls and diagnostic decisions; caller-only. |
| SLAM fitting, 1,500 histograms | 0.377186 → 0.093858 s | Exact fits; approximately 4× kernel speedup. |
| Bulk saved-EC engine | 69.61 → 68.52 s | 1.6% engine gain; removes a 14,464,320-byte buffer; exact quantification and iterations. |
| OCM native routing, matched 100 MB budget | 28.7536 → 8.5893 s | 70.1% less wall time; 95 files exact; RSS increases about 15.2 MiB. |
| OCM existing-call materialization, same budget | 27.5632 → 9.4468 s | 65.7% less wall time; 75 files exact; RSS and temporary writes increase. |

The OCM measurements use real A375 counts with **synthetic sample tags**;
native calling uses 1,000 simulations. They do not establish a biological
OCM or production 100K-simulation speedup. Larger admission budgets change
the time/memory tradeoff and are not matched-budget comparisons.

Accepted regression evidence covers integrated perturb counting, guarded
bulk quantification, SLAM blank/treatment and requantification, scRNA input
layouts, FLEX reader ownership, OCM retained/spilled routing and tag unions,
malformed input, and sanitizer checks for the count store. See the
[cross-module report](https://github.com/morphic-bio/STAR-suite/blob/v1.9.1/docs/benchmarks/CROSS_MODULE_PERFORMANCE_20260910.md)
and [OCM report](https://github.com/morphic-bio/STAR-suite/blob/v1.9.1/docs/benchmarks/OCM_DIRECT_COUNTS_20260910.md)
for exact scope, source identities and limitations.

The completed full 320K FLEX timings remain explicitly labeled **1.9.0**.
New paper benchmark executions should pin `v1.9.1`, record its resolved
source commit and binary hash, and measure whole-workflow time and memory.
This release does not change the validated cell-calling estimator, feature
policy, rescue thresholds or statistical targets. No Cell Ranger source
was inspected for these implementations.

## Distribution

Release packaging provides amd64 glibc234/glibc239 tarballs, a compatibility
installer bundle, Ubuntu 22.04/24.04 Debian binaries, Debian source packaging,
runtime manifests and SHA-256 checksums. Hosted tarball/Debian builds retain
the established portable no-Chromap configuration; local production source
builds retain the Chromap-enabled default. `STAR --source-revision` reports
the exact tagged source embedded in each official binary.
