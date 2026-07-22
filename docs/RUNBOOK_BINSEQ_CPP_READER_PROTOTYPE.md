# BINSEQ C++ Reader And STAR Adapter Prototype Runbook

## Goal

Build the first native C++ reader for ARC BINSEQ CBQ input and a direct adapter
into STAR read buffers as a correctness prototype. The reader emits a cheap
CBQ-native interchange view and must reconstruct the same input-contract state
produced from the source FASTQ files used to create the CBQ:

```text
lane_index	read_ordinal	read_filter	read_name	mate	mate_len	seq_sha256	qual_sha256
```

This phase proved reader correctness and byte-identical STAR internal input
state before wiring CBQ into production mapper routing. Production STAR routing
is now handled in Phase 6 of `docs/RUNBOOK_BINSEQ_INPUT_CONTRACT.md`.

## Scope

- Target CBQ files (`CBQFILE` version 1).
- Support single-end and paired-end records stored in one CBQ file per lane.
- Decode zstd-compressed CBQ columns natively from C++.
- Decode ARC two-bit sequence words and restore `N` positions from the
  serialized Elias-Fano side column.
- Preserve STAR-suite input-contract behavior for read names, Illumina filter
  parsing, mate synchronization, lane indexing, qualities, and optional
  source-order semantics.
- Fill STAR's internal per-mate read buffers directly from `CbqReadBatchView`
  spans without constructing synthetic FASTQ text.
- Keep the existing bqtools-backed `binseq_probe_harness` for comparison and
  upstream fixture probing during the transition.

## Non-Goals

- Do not depend on Rust headers, Rust FFI, or bqtools in the native reader.
- Do not implement random access through the CBQ index in this phase.

## Implementation Shape

1. Add a standalone `CbqInputModule` implementing `InputModule`.
2. Stream one CBQ lane at a time:
   - read the 64-byte file header;
   - read each 96-byte block header;
   - stop at the `CBQINDEX` section;
   - decode the column payloads in block order.
3. Decompress columns with libzstd loaded at runtime through `dlopen`, so the
   prototype does not require zstd development headers.
4. Decode columns into a borrowed `CbqReadBatchView`:
   - sequence lengths and header lengths are little-endian `u64` arrays;
   - sequence words hold 32 bases per `u64`, low bits first,
     `A/C/G/T = 0/1/2/3`;
   - `N` positions remain as a compact per-segment side view over the decoded
     Elias-Fano positions;
   - sequence views point at the packed two-bit CBQ block storage, while quality
     spans point into the decoded quality column;
   - the first mate's header drives `read_name`, `read_name_extra`, and
     `read_filter`, matching `FastxInputModule`.
5. Keep `next_record()` as a compatibility adapter from `CbqReadBatchView` to
   the older owned-string `InputRecord` contract.
6. Add `cbq_reader_harness` to consume `CbqReadBatchView` directly and emit the
   same TSV or debug FASTQ dump as the existing FASTX/BINSEQ probe harnesses.
7. Add `CbqStarAdapter` and `ReadAlign::loadCbqReadView()` to fill STAR's
   existing `readNameMates`, `Read0`, `Qual0`, numeric `Read1`, length,
   filter, lane, and file-type fields directly from the CBQ view.
8. Add `cbq_star_adapter_harness` to compare the direct adapter state against a
   reference path that emulates STAR's current FASTQ/readLoad buffer state.
9. Add app-specific adapters for non-STAR mapper surfaces:
   - `pf_process_records()` accepts in-memory barcode/feature records for
     process_features so CBQ feature libraries do not need temporary FASTQ.
   - `CbqChromapAdapter` materializes synchronized R1/R2/barcode FASTQs from
     paired-read and barcode CBQ sources because Chromap's current contract is
     still path-based.

The CBQ interchange view is block-backed and borrowed: no per-read string
ownership, no full-block ASCII sequence expansion, no synthetic FASTQ text, and
no FASTQ reparse. The batch exposes shared backing storage; consumers that queue
a batch beyond the next reader call must keep that backing alive. ASCII sequence
materialization is a compatibility helper for probes and non-STAR adapters, not
the production STAR mapper path.

## Acceptance Tests

Run:

```bash
make -C core/legacy/source fastx-input-harness cbq-reader-harness cbq-star-adapter-harness
BQTOOLS=/path/to/bqtools tests/run_cbq_cpp_reader_smoke.sh
```

The smoke must pass for:

- paired FASTQ -> CBQ at default compression;
- paired FASTQ -> CBQ at zstd level 0;
- single-end FASTQ -> CBQ at default compression;
- single-end FASTQ -> CBQ at zstd level 0;
- direct `--readFilesIn`;
- manifest input (`CBQ<TAB>-<TAB>ReadGroup`).

The comparison is order-independent because the input contract treats
`read_ordinal` as module-emitted order unless a source plan explicitly
guarantees source order. Record content, mate lengths, read names, filter
status, sequence hashes, and quality hashes must match exactly.

The STAR-adapter comparison is byte-identical and order-preserving for the
paired and single-end synthetic fixtures. It dumps STAR-style internal fields
from both the direct CBQ adapter and the reference FASTQ/readLoad emulation and
requires `cmp` equality.

Additional app-adapter smokes:

```bash
make -C core/legacy/source cbq-pf-adapter-harness cbq-chromap-adapter-harness
BQTOOLS=/path/to/bqtools tests/run_cbq_pf_adapter_smoke.sh
BQTOOLS=/path/to/bqtools tests/run_cbq_chromap_adapter_smoke.sh
```

The process_features smoke compares gzipped FASTQ reference output against the
CBQ in-memory record API for stable count/MEX files. The Chromap smoke verifies
FASTQ payload parity for the materialized Chromap contract files and runs a
tiny Chromap mapping check when `CHROMAP_BIN` is available.

## Transition Policy

- Keep the standalone harnesses as the regression oracle for native reader and
  adapter changes.
- Production `STAR --readFilesType Binseq PE|SE` now consumes the same
  `CbqReadBatchView` through the STAR adapter instead of materializing
  synthetic FASTQ.
- The next integration step is app-specific validation for Solo, Flex, SLAM,
  Chromap, feature-calling, and any other FASTQ-consuming surfaces that need
  CBQ support.
