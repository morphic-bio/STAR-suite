# STAR Suite v1.8.3 Release Notes

Date: 2026-09-05

`v1.8.3` is a performance and robustness release for fully fused STAR-Flex
input and Matrix Market output. It reduces per-record copying and producer
stalls for BGZF FASTQ, streamlines the existing plain-gzip path, reads the
needed CBQ bases directly from packed records, and writes raw and per-sample
MEX matrices in parallel without changing their bytes.

The release artifact version is `v1.8.3`, the Debian source package uses
`1.8.3-1`, and `STAR --version` reports `1.8.3`. GitHub Releases provides
binary packages built and validated on Ubuntu 22.04 and Ubuntu 24.04 as
`1.8.3-1~ubuntu22.04.1` and `1.8.3-1~ubuntu24.04.1`. The upstream STAR base
remains `2.7.11b`; genome-index compatibility remains `2.7.4a`. Existing
indexes do not need to be rebuilt.

## Fused FASTQ input

- Pass fixed-capacity FASTQ records through the fused pipeline without
  constructing transient strings for fields that are not consumed.
- Lease BGZF sequence buffers across fused batches and publish completions
  through a bounded ring, allowing inflater workers to continue while
  completed batches are being processed.
- Balance BGZF inflater work by compressed stream size and retain fine-grained
  work claims. The reader continues to discover each member from its inline
  `BC/BSIZE` field; it does not require a prebuilt index or sidecar.
- Preserve exact record-ordinal pairing and mate-name validation. Cross-member
  FASTQ lines remain bounded and malformed overlong records fail explicitly.
- Apply the compatible fixed-buffer parsing improvements to ordinary gzip
  FASTQ while retaining its established serial zlib stream. Plain gzip is not
  treated as seekable BGZF and is not parallelized.

On the full JAX BGZF FASTQ dataset (2.011B pairs, 32 threads, cold local NVMe),
the retained BGZF reader changes reduced wall time from 269.43 to 238.69
seconds, an 11.4% reduction, with identical pipeline counters, all 16 sample
outputs, and the same output-manifest digest. On the matching ordinary-gzip
input, the optimized median was 487.63 seconds versus 497.82 seconds for the
control, a 2.0% reduction with identical outputs.

## Packed CBQ processing

- Extract the fused no-alignment probe and cell-barcode bases directly from
  CBQ packed sequence storage instead of materializing ASCII sequences and
  converting them back into encoded lookup keys.
- Load hash tables in the packed input's native base order. FASTQ lookup keeps
  its existing order, avoiding a per-read bit-order conversion in either
  input mode.
- Retain the simplified bounded CBQ decoder. Full JAX measurements found its
  wall-time effect neutral, but the implementation makes prefix, full-byte,
  and trailing-base bounds explicit and preserves exact output parity.

The full JAX CBQ control remained effectively unchanged at 217.66 seconds for
the decoder candidate versus the established 218.23-second result. The CBQ
changes are retained for simpler bounds handling and to avoid redundant ASCII
materialization in the fused no-alignment path, not as a claimed wall-time
speedup.

## Parallel Matrix Market output

- Write matrix entries through bounded 262,144-entry chunks. A parallel sizing
  pass computes exact byte offsets, then workers format disjoint mapped ranges;
  no whole-matrix temporary string or unbounded mapping is required.
- Generate per-sample matrices concurrently while partitioning the configured
  STAR thread budget between sample and matrix workers. Raw fused MEX output
  may use the full budget.
- Preserve row ordering, column ordering, headers, dimensions, and every output
  byte across thread counts. Existing single-thread `MexWriter` entry points
  remain available to callers.

On the full JAX CBQ run, the 99,266,019-entry raw matrix was written in 0.332
seconds and all 16 per-sample matrices in 0.393 seconds using 32 threads.
`processRecords` fell from 47.8708 to 38.7911 seconds (19.0%), and total wall
time fell from 3:29.40 to 3:20.73. Peak RSS remained essentially unchanged
(41,794,120 versus 41,761,232 kB), and all 70 compared raw and per-sample files
were SHA-256 byte-identical. Because MEX output follows fused ingest, this
improvement applies equally to supported gzip, BGZF, and CBQ Flex runs.

## Compatibility

- The optimizations are confined to the fully fused STAR-Flex input and output
  paths. Classic STARsolo, bulk RNA-seq, scRNA-seq ingest, and other non-Flex
  behavior are unchanged.
- No htslib code or format is changed, and no external dependency is added.
- BGZF CRC checking and all existing reader-mode controls keep their prior
  defaults and behavior.
