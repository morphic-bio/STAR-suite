# STAR Suite v1.8.2 Release Notes

Date: 2026-09-04

`v1.8.2` is a correctness, robustness, and performance release for the fully
fused STAR-Flex path. It fixes progress when residual alignment is enabled,
restores STARsolo-compatible ambiguous-cell-barcode posterior resolution, and
pipelines barcode-bucket ordering through producer-local sorted runs and
asynchronous hierarchical consolidation.

The release artifact version is `v1.8.2`, the Debian source package uses
`1.8.2-1`, and `STAR --version` reports `1.8.2`. GitHub Releases provides
binary packages built and validated on Ubuntu 22.04 and Ubuntu 24.04 as
`1.8.2-1~ubuntu22.04.1` and `1.8.2-1~ubuntu24.04.1`. The upstream STAR base
remains `2.7.11b`; genome-index compatibility remains `2.7.4a`. Existing
indexes do not need to be rebuilt.

## Fully fused Flex correctness and progress

- Keep residual-alignment queues making progress when producers fill the
  bounded queue. Producers help drain alignment work until they can enqueue,
  avoiding the startup deadlock that appeared when hash misses exceeded queue
  capacity.
- Validate paired BGZF read-name stems on every record rather than assuming
  that equal ordinals imply a valid pair. Truncated, filtered, or independently
  reordered mates now fail at their first mismatch.
- Resolve ambiguous Flex cell barcodes with the established STARsolo
  posterior: candidate base quality, exact-barcode read-count priors, and the
  configured posterior threshold. The fused bucket path retains every
  observation assigned by that decision and rejects malformed evidence.

## Pipelined sorted bucket runs

- Sort each immutable producer tranche locally before publication, allowing
  tranche ordering to overlap continued mapping by other producers.
- Consolidate same-level sorted runs asynchronously in bounded 64-way batches
  while input is still arriving. Once the input is closed, a deterministic
  tournament-tree merge emits the globally ordered bucket stream; final output
  cannot safely be emitted earlier because a later tranche may contain the
  lowest remaining key.
- Preserve sorted-run boundaries for RAM, spill, and automatic RAM-to-spill
  storage. The spill payload and validation format remain checksummed and the
  final matrices remain deterministic.
- On the full JAX CBQ hash-only control (2.011B read pairs, 32 threads, local
  NVMe), this release reduced wall time from 3:44.57 for the corrected baseline
  to 3:38.23. The bucket loop fell from 14.52 to 8.48 seconds, with
  byte-identical per-sample outputs. The small wall benefit is expected because
  this dataset has relatively few surviving ambiguities; the run-based design
  is principally a more robust scaling path for larger tails.

## Other fused-tail improvements

- Avoid the per-read sample-table mutex after each sample token's one-time
  registration.
- Accumulate ambiguous observations in one striped store and parallelize the
  producer fan-in and teardown work.
- Inline packed-record accessors used on the bucket-ordering hot path and keep
  detailed stage timing for reproducible performance diagnosis.
- Retain the full-budget FlexFilter Monte Carlo allocation used by the
  corrected JAX benchmark. Output parity remains the release gate.

## Ubuntu packages and portable binaries

- Publish two dynamic `.deb` packages from the same tagged source: an Ubuntu
  22.04 build and an Ubuntu 24.04 build. The suffix records the build and test
  baseline; it does not select a different STAR algorithm or output contract.
- Keep C11 atomic declarations confined to C translation units so the packaged
  `star_feature_call` companion builds under Ubuntu 22.04's GCC 11 as well as
  Ubuntu 24.04's newer compiler.
- Validate the Ubuntu 22.04 package on both Ubuntu 22.04 and 24.04, and validate
  the Ubuntu 24.04 package on Ubuntu 24.04. A binary built against a newer
  system runtime is not promised to run on an older distribution.
- The release continues to provide `glibc234` and `glibc239` compatibility
  tarballs plus an installer bundle that selects the compatible binary. Use
  the installer bundle when the host distribution is not known in advance.

## Cloud benchmark tooling

- Harden the public 320K cloud matrix around blank instance-store detection,
  archive extraction, fail-closed completion markers, durable progress copies,
  staged reference assets, and portable CBQ encoder construction.
- Record scratch high-water marks and preserve the commands and manifests used
  for the matched FASTQ/BGZF/CBQ controls.
