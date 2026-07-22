# Runbook: process_features Native CBQ Input

Date: 2026-05-28

Status: Phase 1 implemented and smoke-tested. The current CBQ path for
`process_features` is still a harness/API surface, not a first-class production
CLI input mode.

## Goal

Run feature-barcode assignment from CBQ without materializing FASTQ:

```text
shared CBQ parser -> input-format gate -> PF decoded line queue -> PF consumer
```

The production surface should feed barcode and feature reads directly into
`process_features`. Temporary FASTQ is only for oracle tests and debugging.

## Current State

The useful pieces already exist:

- `core/legacy/source/input/CbqInputModule.{h,cpp}` reads CBQ natively and emits
  `CbqReadBatchView`.
- `core/features/process_features/include/pf_api.h` exposes streaming
  `pf_process_records_begin()` / `pf_process_record_views()` /
  `pf_process_records_end()` APIs for non-FASTQ providers.
- `core/legacy/source/input/cbq_pf_adapter_harness.cpp` streams
  `CbqReadBatchView` records into the PF decoded input gate.
- `tests/run_cbq_pf_adapter_smoke.sh` compares gzipped FASTQ input with native
  CBQ input and checks MEX/count parity.

The gaps:

- CBQ is not a first-class `process_features` CLI/workflow input.
- The production CLI still lacks an explicit `--input_format cbq` gate.
- Multi-sample/multi-lane CBQ manifest conventions for PF are not documented.
- The current CBQ gate decodes to the same ASCII line shape consumed by the
  existing PF queue. This is an in-memory adapter, not disk FASTQ
  materialization. A direct coded-sequence fast path remains optional future
  work.

## Reader Policy

`process_features` does not use the STAR aligner, so it should not be forced
through STAR mapper internals. A PF-specific CBQ input module is reasonable.

However, it should not become a second CBQ parser. Keep one source of truth for
CBQ file/block/index parsing, zstd decompression, two-bit decode, `N`
restoration, quality/header handling, and lane iteration.

Recommended policy:

- Share the low-level CBQ parser.
- Add a PF-local consumer that maps decoded CBQ batches to PF barcode/feature
  input.
- If needed, move the CBQ parser into a small linkable common component so
  `process_features` can depend on it without depending on STAR alignment
  code.
- Do not duplicate CBQ format parsing in `process_features` unless the common
  parser becomes impossible to link cleanly.

This gives PF its own input surface while avoiding divergent CBQ semantics.

## Target Architecture

Add a PF-native input-format gate:

```text
if input_format == FASTQ:
  existing FASTQ/gzip producers -> decoded line queue
else if input_format == CBQ:
  shared CBQ parser -> CBQ producer/adapter -> decoded line queue

decoded line queue -> existing PF consumer threads -> assignment/counting/finalize
```

Minimum batch record fields:

- barcode sequence span: CB + UMI read, usually R1;
- barcode quality span, with explicit missing-quality handling;
- feature sequence span, usually R2;
- feature quality span, optional;
- optional second feature sequence/quality span for dual-orientation libraries;
- lane/sample ordinal for logging and future multisample accounting.

The shared CBQ reader keeps STAR-facing sequence data in packed two-bit storage.
The PF compatibility adapter may materialize bounded ASCII sequence spans,
because current `process_features` consumers still expect FASTQ-style
sequence/quality line buffers. A later optimization can add a 2-bit PF path if
profiling shows that sequence decoding, not I/O or feature search, is limiting.

TODO: reduce the decoded-record interface surface after the production PF CBQ
gate is stable. CBQ is the immediate driver, but this is a general PF input
architecture cleanup: callers should not need to understand FASTQ-line newline
conventions, default-quality synthesis, or multiple record view shapes. Hiding
those details behind a smaller decoded-record adapter will make all non-FASTQ
producers more robust and will let PF change its internal sequence format
without widening the compatibility burden.

## Phase 1: Make The Existing Adapter Streaming

Goal: remove all-record accumulation without changing output.

Status: complete. The CBQ harness streams one `CbqReadBatchView` at a time into
the PF streaming API, and that API feeds the existing decoded line queue and
consumer thread path.

Tasks:

- Add a batch API beside `pf_process_records()`, for example
  `pf_process_record_batch()` plus explicit `pf_begin_sample()` /
  `pf_finish_sample()`, or one higher-level streaming API.
- Reuse the existing downstream PF consumer path so FASTQ, array records, and
  CBQ batches share assignment/counting behavior after input decoding.
- Preserve current locking around the process-global PF runtime state.
- Keep `pf_process_records()` as a compatibility wrapper that feeds one or more
  batches into the new streaming implementation.
- Update `cbq_pf_adapter_harness.cpp` to process one `CbqReadBatchView` at a
  time instead of building a full `std::vector<OwnedPfRecord>`.

Exit criteria:

- `tests/run_cbq_pf_adapter_smoke.sh` still passes.
- Peak memory scales with batch size, not total reads.
- FASTQ API output is unchanged.
- CBQ records are not materialized as disk FASTQ; they are copied into the
  existing bounded in-memory decoded line queue.

## Phase 2: Add A PF CBQ CLI/Workflow Input

Goal: make CBQ a production `process_features` input mode.

Tasks:

- Add explicit CBQ input flags to the PF-facing tool surface. Suggested narrow
  first surface:

```text
--input_format cbq
--readFilesIn feature_reads.cbq
```

- Validate that paired CBQ records provide barcode and feature segments in the
  configured mate roles.
- Keep the existing FASTQ flags and discovery behavior unchanged.
- For multi-sample input, either require a manifest or mirror the existing
  sample grouping explicitly; do not infer sample names from CBQ paths until
  that convention is documented.
- Wire the production command to the streaming batch API, not to FASTQ
  materialization.

Exit criteria:

- A production PF command can process CBQ directly.
- The command reports clear errors for wrong mate count, missing qualities when
  qualities are required, empty sequences, and mismatched sample metadata.
- Existing FASTQ commands are unchanged.

## Phase 3: Regression Tests

Add tests before broad use:

- Tiny synthetic parity: current `tests/run_cbq_pf_adapter_smoke.sh`, updated
  to exercise the production CBQ CLI path.
- Compression coverage: default-compressed and level-0 CBQ.
- Fixture-scale parity: one downsampled production feature-barcode sample,
  comparing MEX/count outputs against gzipped FASTQ.
- Negative tests: single-end CBQ passed to paired PF mode, missing mate,
  malformed quality length, and feature/barcode role mismatch.
- Register the production case in
  `tests/run_cbq_e2e_module_regression.sh` and the production regression
  manifest once stable.

## Phase 4: Optional 2-Bit Fast Path

Only add this after correctness and production wiring are done.

The CBQ sequence columns and PF feature matching both have compact sequence
representations, so there may be a real speedup from avoiding ASCII expansion
for feature reads. Treat this as an optimization, not the initial production
contract.

Current code check:

- Exact/common-length feature matching already uses integer sequence hashes
  (`seq_encode_64_fixed`, `seq_encode_128_fixed`, `feature_hamming_le1/le2`).
- `feature_lookup_code()` is currently disabled, so there is no clean coded
  feature-lookup API yet.
- The consumer still needs sequence spans for offset scans, anchor scans,
  `N` expansion, barcode/UMI correction, and recording the matched feature
  sequence. A direct coded path should be added as a targeted optimization
  after the production CBQ CLI gate is stable.
- Before adding a coded fast path, narrow the ASCII adapter surface so the
  default route remains easy to reason about and new input producers do not
  duplicate newline/quality/defaulting details.

Guardrails:

- Preserve exact handling of `N` positions and PF `max_feature_n` /
  `max_barcode_n` behavior.
- Keep barcode correction and UMI extraction identical to the ASCII path.
- Require parity tests before and after enabling the fast path.
- Make the fast path optional until it has production-scale timing evidence.

## Recommendation

Use a PF-local CBQ input module, but build it on the shared CBQ parser. That is
the right compromise: PF gets a direct reader/consumer appropriate for a
standalone feature caller, while STAR-suite keeps one authoritative
implementation of the CBQ format.
