# BINSEQ C++ Reader Prototype Runbook

## Goal

Build the first native C++ reader for ARC BINSEQ CBQ input as a correctness
prototype. The reader must reconstruct the same input-contract state produced
from the source FASTQ files used to create the CBQ:

```text
lane_index	read_ordinal	read_filter	read_name	mate	mate_len	seq_sha256	qual_sha256
```

This phase proves reader correctness before optimizing for direct STAR mapper
handoff.

## Scope

- Target CBQ files (`CBQFILE` version 1).
- Support single-end and paired-end records stored in one CBQ file per lane.
- Decode zstd-compressed CBQ columns natively from C++.
- Decode ARC two-bit sequence words and restore `N` positions from the
  serialized Elias-Fano side column.
- Preserve STAR-suite input-contract behavior for read names, Illumina filter
  parsing, mate synchronization, lane indexing, qualities, and optional
  source-order semantics.
- Keep the existing bqtools-backed `binseq_probe_harness` for comparison and
  upstream fixture probing during the transition.

## Non-Goals

- Do not integrate CBQ into the production `STAR` binary yet.
- Do not map CBQ fields directly into STAR's internal read buffers yet.
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
4. Decode columns into a plain `InputRecord`:
   - sequence lengths and header lengths are little-endian `u64` arrays;
   - sequence words hold 32 bases per `u64`, low bits first,
     `A/C/G/T = 0/1/2/3`;
   - `N` positions are restored after two-bit decode;
   - qualities are sliced with sequence offsets;
   - the first mate's header drives `read_name`, `read_name_extra`, and
     `read_filter`, matching `FastxInputModule`.
5. Add `cbq_reader_harness` to emit the same TSV or debug FASTQ dump as the
   existing FASTX/BINSEQ probe harnesses.

## Acceptance Tests

Run:

```bash
make -C core/legacy/source fastx-input-harness cbq-reader-harness
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

## Transition Policy

- Keep new code behind the standalone harness until the parity test is stable.
- Use this harness as the regression oracle before adding a production
  `STAR --readFilesType Binseq` surface.
- The next implementation step is an optimized adapter that fills STAR read
  buffers directly from CBQ block slices without constructing intermediate
  FASTQ text.
