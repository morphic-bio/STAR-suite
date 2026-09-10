# PF native BGZF — phase 1

Date: 2026-09-10. Source baseline: `9e146eef20080108cc075a977938757f90221dbc`.
Status: phase 1 implemented and validated in the isolated performance branch; full PF comparison complete.

## Implementation

- Reuse `BgzfRangeReader` / `BgzfBlockReader` through `pf_bgzf_input`, a C ABI
  linked into standalone PF without STAR or a genome dependency.
- Detect regular BGZF inputs by their headers. Select native or legacy input
  per lane; require every stream in a native lane to be BGZF. Forced range
  rejects unsupported lanes before processing. Split-read synthesis currently
  logs a gzip fallback; forced native range rejects it.
- Validate mate ordinals and name tokens (terminal `/1`, `/2`, `/3` ignored),
  preserving independent compressed block boundaries. Return borrowed spans;
  copy them into the existing ring before releasing the decoder leases.
- Fail on input/CRC/name/count/capacity errors and suppress final count export.
  Canonical decoder behavior allowing a missing terminal EOF marker is retained.
- Distribute the requested inflater count across lanes and streams. Reuse PF's
  existing FEATURE permit hooks for inflate work. Release partial assignment
  permit batches when input queues are empty, avoiding producer starvation.
- Preserve standalone failures from earlier child batches even when a later sample succeeds.
- Catch worker exceptions in the shared reader and return input errors instead
  of terminating the process. Return held permits if inflation throws.

## Controls and accounting

Standalone: `--readFilesBgzfMode auto|off|range`, `--bgzfReaderThreads N`,
`--bgzfCrcCheck 0|1`. Zero workers means synchronous inflation on producer
threads, preserving the default thread footprint. Explicit workers are
additional to existing producers/consumers; the log prints all counts.

API: `pf_config_set_bgzf_input(config, mode, threads, crc)`.

Embedded: `--crAssignBgzfMode auto|off|range` selects PF input independently of
GEX. Existing `bgzfReaderThreads` and `bgzfCrcCheck` supply decoder settings.
With dynamic scheduling, zero requests derive the inflater count from
`runThreadN`; inflate and assignment use the shared FEATURE pool. The chemistry
probe runs before that pool is active and uses synchronous inflation when
threads were automatic. Existing auto-sized consumer/search policy is retained.
Producer/parse coordination remains outside the compute permits in this first
phase; phase 2 must account for it while removing the ring handoff.

## Validation and evidence

Evidence root:
`/mnt/pikachu/star_suite_paper/analysis/cross_module_performance_20260910/`.

- `00_inventory`: frozen source identity, clean build logs, Chromap-enabled STAR.
- `01_pf_bgzf_tests`: 21 distinct arms, exact ordered name/sequence/quality
  digests; two/three streams, unequal blocks, spanning records, CRLF, missing
  EOF, malformed records, oversized fields, CRC, mate mismatch, early close,
  single-permit backpressure, mixed lanes and forced-mode failure.
- `03_pf_cli_gates`: six successful directory CLI arms comparing baseline and
  native BGZF for triplet reads, ADT export and chemistry detection.
- `04_final_validation`: final clean build, shared-reader injected exception,
  single-permit paired/triplet decode and assignment parity.
- `02_a375_pf_100k`: exactly 200,000 paired feature reads (100K per lane), source
  and decoded/BGZF hashes, frozen baseline and native API outputs.
- `05_a375_pf_full`: full input conversion/manifest and completed benchmark evidence.
- `06_cli_multisample_failure`: early failed sample followed by a successful sample; parent correctly returns failure and only the successful sample exports counts.

Integer count comparisons resolve matrix indices against barcode/feature names.
Existing PF hash iteration can change barcode column order; bytewise matrix
comparison alone is invalid. Assignment statistics must also match. The raw
ADT matrix header uses `real` and its compressed export uses `integer` in the
baseline too; numeric content and axes are the comparison surfaces.

Small A375 result: 8,423 barcodes, 11 features, 8,581 nonzero entries,
146,979 UMIs, exact counts/axes/assignment diagnostics. Whole process 2.0695 s
baseline vs 2.0689 s native, RSS 210,688 KiB both. No speed or memory reduction
is established by this small test. Conversion is outside timing. No BAM,
sidecar, alignment reference load, GEX processing or cell calling is involved.

A frozen standalone explicit-lane-list invocation crashed before read handling;
its failed log is retained and not counted as a benchmark. The API lane-list
path and standalone directory path both passed. Investigate the older CLI
failure separately if that interface is needed.

## Full A375 feature-library result

All **9,748,584 paired feature reads**, across two lanes, were processed in
serialized runs using the same BGZF files, eight-CPU affinity, four assignment
consumers, search threads = 1, queue = 8,192 lines and six shared work permits.
The native reader used four inflater workers distributed across both lanes.
Inputs, reference hashes, exact argv, environment and completion records are
in `05_a375_pf_full`.

| Metric | Frozen PF baseline | Native BGZF |
| --- | ---: | ---: |
| Whole-process wall time | 28.2722 s | 25.5696 s |
| Reported read-processing phase | 22.25 s | 21.26 s |
| Maximum RSS | 1,554,892 KiB | 1,593,580 KiB |
| UMIs | 3,218,159 | 3,218,159 |
| Barcodes | 110,772 | 110,772 |
| Features | 11 | 11 |
| Occupied matrix entries | 124,871 | 124,871 |

Counts, named axes and assignment diagnostics match exactly. This single
comparison shows **9.6% lower whole-process wall time**, with **37.8 MiB more
peak RSS**. The reported ingestion/assignment phase improved 4.4%; do not
attribute the whole-process difference solely to decompression. This is the
full feature library, not full perturb GEX processing or a cell-calling run.
The six work permits include inflate and assignment; the two producer/parse
threads are separately accounted for within the eight-CPU affinity.

## Next steps and limits

Preserve the completed full PF baseline for later phases.
Next connect leased batches to the existing worker-owned PF kernel and remove
per-record ring locking/copies. Preserve gates for legacy rescue and online
feature-offset bootstrap, which the existing direct API does not support.
Do not simply bypass those gates. Whole perturb GEX+PF overlap has not been
benchmarked; core STAR has been clean-built with Chromap enabled. General GEX,
bulk, scRNA/OCM and SLAM work remains in later phases of the parent runbook.
No code was inspected from Cell Ranger. No release, merge, push or cloud job
is part of this phase's validation.
