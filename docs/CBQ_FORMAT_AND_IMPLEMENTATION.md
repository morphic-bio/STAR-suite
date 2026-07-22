# CBQ Format And Implementation

This document describes the CBQ/BINSEQ format subset supported by STAR Suite
and how the native C++ reader is wired into STAR Suite modules. It is an
implementation reference for STAR Suite `v1.1.0`, not a complete external ARC
BINSEQ specification.

User-facing command examples live in `docs/EXPERIMENTAL_BINSEQ_INPUT.md`.
Ordered FASTQ-to-CBQ encoding policy lives in
`docs/RUNBOOK_CBQ_ORDERED_ENCODER.md`.

## Supported Scope

STAR Suite supports CBQ file version 1 through:

- `STAR --readFilesType Binseq PE --readFilesIn sample.cbq`
- `STAR --readFilesType Binseq SE --readFilesIn sample.cbq`
- manifest input with one CBQ path in column 1 and `-` in column 2 for paired
  CBQ records
- multi-lane and batch-mode STAR inputs, where each lane is one CBQ file
- STAR mapper, STARsolo, OCM, Flex, SLAM, and process_features adapter
  surfaces

FASTQ and FASTQ.gz paths are unchanged. Production CBQ support means the
consumer receives an in-memory decoded view or a direct adapter into its native
read buffers. It does not mean writing temporary FASTQ files, except for the
current Chromap compatibility adapter.

The generic input contract does not require source-order preservation. The
native `cbq_ordered_encoder` does preserve FASTQ source order and should be
used for parity tests and production workflows where deterministic downstream
boundaries depend on read order.

## File Layout

All integer fields below are little-endian. CBQ columns are zstd-compressed
payloads unless their stored length is zero.

### File Header

The file begins with a 64-byte header:

| Offset | Size | Field |
|---:|---:|---|
| 0 | 7 | magic bytes `CBQFILE` |
| 7 | 1 | version, currently `1` |
| 8 | 8 | presence flags |
| 16 | 8 | compression level recorded by the encoder |
| 24 | 8 | virtual block size requested by the encoder |
| 32 | 32 | reserved, currently zero-filled by the STAR Suite encoder |

Presence flags used by the C++ implementation:

| Bit | Name | Meaning |
|---:|---|---|
| 0 | `PRESENCE_PAIRED` | one record contains two sequences |
| 1 | `PRESENCE_QUALITIES` | quality column is present |
| 2 | `PRESENCE_HEADERS` | header column is present |
| 3 | `PRESENCE_FLAGS` | optional flags column is present |

The STAR Suite ordered encoder always writes headers and qualities. It sets
`PRESENCE_PAIRED` for paired input. The reader can parse a flags column, but
STAR Suite `v1.1.0` integrations do not consume record flags.

### Blocks

After the file header, the reader streams block headers and block payloads
until it sees the `CBQINDEX` marker.

Each block header is 96 bytes:

| Offset | Size | Field |
|---:|---:|---|
| 0 | 3 | magic bytes `BLK` |
| 3 | 1 | block version, currently `1` |
| 4 | 4 | reserved/filler, ignored by the reader |
| 8 | 8 | compressed byte length of `seq_len` |
| 16 | 8 | compressed byte length of `header_len` |
| 24 | 8 | compressed byte length of `npos` |
| 32 | 8 | compressed byte length of `seq` |
| 40 | 8 | compressed byte length of `flags` |
| 48 | 8 | compressed byte length of `headers` |
| 56 | 8 | compressed byte length of `qual` |
| 64 | 8 | total bases in the block, `nuclen` |
| 72 | 8 | decompressed byte length of the `npos` Elias-Fano payload |
| 80 | 8 | number of records in the block |
| 88 | 8 | number of sequences in the block |

The payload columns immediately follow the header in this exact order:

```text
z_seq_len
z_header_len
z_npos
z_seq
z_flags
z_headers
z_qual
```

There are no per-column length prefixes in the payload. The compressed lengths
in the block header define the column boundaries.

### Columns

`seq_len`
: zstd-compressed `uint64_t` length per sequence. The expected decompressed
  size is `num_sequences * 8`.

`header_len`
: zstd-compressed `uint64_t` header length per sequence when headers are
  present. If headers are absent, the reader synthesizes empty header lengths.

`npos`
: zstd-compressed sucds-compatible Elias-Fano payload for global ambiguous
  base positions in the concatenated block sequence. Ambiguous non-ACGT input
  bases are stored as two-bit code `0` and restored as `N` through this
  position list.

`seq`
: zstd-compressed two-bit packed sequence words for all sequences in the block
  concatenated together. Each little-endian `uint64_t` stores 32 bases with
  `A=0`, `C=1`, `G=2`, `T=3`. The expected decompressed size is
  `ceil(nuclen / 32) * 8`.

`flags`
: optional zstd-compressed record flags. The STAR Suite reader validates the
  column length when the presence flag is set, but current STAR Suite
  integrations do not use the values.

`headers`
: zstd-compressed concatenated FASTQ header payloads. The ordered encoder
  removes the leading `@` before storage. The reader accepts payloads with or
  without a leading `@` or `>` when parsing names.

`qual`
: zstd-compressed concatenated quality strings, one quality string per
  sequence, in the same sequence order as `seq_len`. If qualities are absent,
  the reader currently synthesizes `A` qualities as a compatibility fallback.

### Records And Mates

For single-end CBQ, one record has one sequence. For paired CBQ, one record has
two sequences in record order:

```text
record 1: segment 0, segment 1
record 2: segment 0, segment 1
...
```

`num_sequences` must equal `num_records * mate_count`; the reader rejects a
block if this does not match the requested `SE` or `PE` mode.

Segment roles are assigned by position: segment 0 is `Read1`, segment 1 is
`Read2`. Some STAR Suite workflows intentionally encode mates in command order
instead of biological R1/R2 order. For example, Flex and OCM recipes can encode
cDNA R2 before barcode R1 because that is the order expected by the STAR
command surface.

### Index

The ordered encoder writes a trailing `CBQINDEX` section:

```text
CBQINDEX
uint64_t uncompressed_index_bytes
uint64_t compressed_index_bytes
zstd-compressed index payload
uint64_t compressed_index_bytes
CBQINDEX
```

The uncompressed index payload stores `(block_offset, cumulative_records)` as
little-endian `uint64_t` pairs for each block. The STAR Suite reader currently
streams blocks and stops when it sees `CBQINDEX`; it does not use the index for
random access.

## Ordered Encoder

The reference encoder is implemented in
`core/legacy/source/input/cbq_ordered_encoder.cpp`. It uses zlib for gzipped
FASTQ input and loads libzstd dynamically for column compression.

Build and run it with:

```bash
make -C core/legacy/source cbq-ordered-encoder

core/legacy/source/cbq_ordered_encoder \
  --readFilesIn R1.fastq.gz R2.fastq.gz \
  --outFile lane.cbq
```

Single-end input uses one FASTQ path:

```bash
core/legacy/source/cbq_ordered_encoder \
  --readFilesIn reads.fastq.gz \
  --outFile reads.cbq
```

The encoder supports plain FASTQ and `.gz` FASTQ. Paired read names are
validated by default after stripping the leading `@`, whitespace suffixes, and
trailing `/1` or `/2`. Use `--no-validate-read-names` only for a known
non-standard fixture.

The encoder writes records sequentially in input order. It does not sort,
parallelize output publication, or allow later blocks to be written before
earlier blocks. The default virtual block size is 1 MiB. The compression level
is passed to zstd and recorded in the file header.

## Native Reader

The shared reader is implemented in:

- `core/legacy/source/input/CbqInputModule.h`
- `core/legacy/source/input/CbqInputModule.cpp`

The reader also loads libzstd dynamically. A runtime host must provide
`libzstd.so.1` or `libzstd.so`.

`make_cbq_input_source_plan()` creates an `InputSourcePlan` with:

- `format = SourceFormat::Binseq`
- `module_name = "Cbq"`
- one CBQ path per lane
- `mate_count` of 1 or 2
- `preserves_source_order = true` for the STAR Suite native reader
- no helper stream and no internal gzip path

`CbqInputModule::next_batch()` emits a borrowed `CbqReadBatchView`. Pointers in
the view remain valid only while the returned `backing` shared pointer is
alive. Consumers that queue decoded batches must keep `batch.backing`.

Important view structs:

- `CbqReadBatchView`: one decoded block of records plus backing storage
- `CbqReadView`: read name, read-name suffix, lane, ordinals, filter, segments
- `CbqSegmentView`: segment role, sequence view, quality span, packed sequence
  view, original length
- `CbqPackedSequenceView`: pointer to packed sequence bytes, base offset,
  length, and `N` positions for that segment

The batch view keeps sequence in packed two-bit block storage. ASCII sequence
is materialized only through explicit helpers or compatibility adapters.

`CbqInputModule::next_record()` is a compatibility path that materializes
`InputRecord` strings. It is useful for harnesses and older adapter surfaces,
but it is not the preferred high-throughput STAR integration path.

## Adapter Surfaces

### STAR Mapper, STARsolo, OCM, Flex, And SLAM

STAR production routing uses:

- `core/legacy/source/input/CbqStarAdapter.h`
- `core/legacy/source/input/CbqStarAdapter.cpp`
- `core/legacy/source/ReadAlign_loadCbqReadView.cpp`
- `core/legacy/source/ReadAlignChunk_processChunks.cpp`
- `core/legacy/source/ReadAlignChunk_mapChunk.cpp`

The adapter fills the same STAR read buffers used after FASTQ parsing. Packed
CBQ bases are translated through byte lookup tables into STAR numeric bases and
ASCII buffers where legacy STAR code still expects ASCII. `N` positions are
applied after two-bit decode so STAR sees base code `4` and ASCII `N` at the
same positions as the original FASTQ.

This keeps the existing optimized FASTQ path intact while adding a CBQ gate
inside the STAR chunk-processing path.

### process_features

The current process_features CBQ support uses the shared reader and fills the
existing decoded feature-record path. It does not write temporary FASTQ files.
The compatibility surface may materialize ASCII sequence because existing
feature matching code still consumes sequence strings in places. The format
parser remains shared; process_features should not grow a second CBQ parser.

Relevant files:

- `core/legacy/source/PfMultiAssign.cpp`
- `core/legacy/source/input/cbq_pf_adapter_harness.cpp`
- `docs/RUNBOOK_PROCESS_FEATURES_CBQ_NATIVE.md`

### Chromap

STAR still keeps the historical Chromap compatibility adapter as a test oracle,
but the production libchromap contract can now pass CBQ paths directly to
Chromap-suite's native CBQ reader:

```text
--chromapAtacInputFormat cbq
--chromapAtacReadPairCbq lane1.reads.cbq,lane2.reads.cbq
--chromapAtacBarcodeCbq lane1.barcodes.cbq,lane2.barcodes.cbq
```

This path avoids temporary FASTQ materialization and shares the existing
Chromap ATAC mapping, BAM/fragments, sidecar, and peak-MEX surfaces.

Relevant files:

- `core/legacy/source/input/CbqChromapAdapter.h`
- `core/legacy/source/input/CbqChromapAdapter.cpp`
- `core/legacy/source/input/cbq_chromap_adapter_harness.cpp`
- `core/features/libchromap_contract/include/star_chromap_contract.h`
- `core/features/libchromap_contract/src/star_chromap_contract.cpp`
- `core/legacy/source/star_chromap_orchestration.cpp`

## Command Surface

Direct paired CBQ:

```bash
STAR \
  --readFilesType Binseq PE \
  --readFilesIn sample.cbq \
  ...
```

Direct single-end CBQ:

```bash
STAR \
  --readFilesType Binseq SE \
  --readFilesIn sample.cbq \
  ...
```

Manifest paired CBQ:

```bash
printf '%s\t-\tID:sample1\n' sample.cbq > binseq_manifest.tsv

STAR \
  --readFilesType Binseq PE \
  --readFilesManifest binseq_manifest.tsv \
  ...
```

Multiple lanes are represented as one CBQ path per lane in the existing STAR
comma-separated input conventions or manifest rows, depending on the workflow.

## Validation Matrix

Core smoke and regression coverage:

- `tests/run_cbq_ordered_encoder_smoke.sh`
- `tests/run_cbq_cpp_reader_smoke.sh`
- `tests/run_cbq_star_input_smoke.sh`
- `tests/run_cbq_solo_e2e_smoke.sh`
- `tests/run_cbq_pf_adapter_smoke.sh`
- `tests/run_cbq_chromap_adapter_smoke.sh`
- `tests/run_cbq_ocm_composite_smoke.sh`
- `tests/run_cbq_flex_tiny_public_smoke.sh`
- `tests/run_slam_cbq_divergence_harness.sh`
- `tests/run_cbq_e2e_module_regression.sh`
- `tests/production_module_regression_manifest.tsv`

The release gate for STAR Suite `v1.1.0` included the downsampled CBQ E2E
module suite, STAR mapper and STARsolo parity smokes, process_features and
Chromap adapter smokes, OCM/Flex production-shaped smokes, and SLAM divergence
harnesses with pre-NTR parity.

## Limitations And Follow-Ups

- Only CBQ file version 1 is supported.
- `--readFilesCommand` is not supported with `--readFilesType Binseq`.
- Y/noY FASTQ emission remains FASTQ-only until non-FASTQ emission has a
  separate contract.
- SLAM per-file skipping is rejected for BINSEQ input.
- Chromap ATAC CBQ requires Chromap-suite with native CBQ support and uses
  split sources: one paired-read CBQ plus one barcode CBQ per lane.
- The STAR Suite reader streams blocks and does not use the CBQ index for
  random access.
- External CBQ encoders may not preserve source order. Use
  `cbq_ordered_encoder` when source-order parity matters.
- The reader accepts optional flags at the column-validation level, but STAR
  Suite consumers do not currently interpret record flags.
