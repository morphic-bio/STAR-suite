# TODO: Optimize CBQ conversion and add FASTQ-side caching

Status: deferred until after the preprint benchmarks. The preprint production
path remains direct FASTQ input. Do not replace a direct-FASTQ benchmark with a
CBQ run unless the operator explicitly requests it.

## Motivation

CBQ remains a useful reusable input format, but the current standalone
`cbq_ordered_encoder` is single-threaded. On the full SPATCh Flex slide it
encoded 290,000,000 of 676,879,511 read pairs in 22:25 before the operator
cancelled the experiment. A separate conversion pass therefore costs enough
to obscure the first-run advantage even when subsequent CBQ processing is
fast.

Most required machinery is already present:

- `core/legacy/source/input/CbqWriter.{h,cpp}` implements the CBQ format writer.
- `core/legacy/source/input/CbqYNoYWriter.{h,cpp}` provides an ordered sidecar
  writer precedent.
- The fused Flex FASTQ producer/consumer path already owns decompressed and
  parsed paired records before hash triage and alignment-fallback routing.

## Work order

1. Optimize the standalone converter before promoting CBQ in benchmarks.
   Separate block encoding from ordered output, encode independent blocks in
   parallel, preserve deterministic input order, and avoid per-record string
   allocation/copying. Benchmark conversion wall time, CPU utilization,
   compressed size, and byte/record parity on 100K before a full slide.
2. Add an opt-in FASTQ `process + cache` mode. Tee each already-decompressed,
   parsed read pair into an asynchronous CBQ block encoder while the same
   record continues through Flex processing. Do not decompress or parse FASTQ
   twice.
3. Use a bounded queue for backpressure, parallel block encoders, and a single
   ordered file/index writer. Pass batch/block ownership from the FASTQ reader
   instead of copying each read through synchronous `CbqWriter::add_record`.
4. Emit one CBQ per original input lane. This preserves lane/source order and
   allows multi-lane conversion and later reads to remain parallel.
5. Write to temporary files and publish atomically only after all blocks and
   the final index are complete. Record record counts, input identities,
   format version, checksum, compression settings, stage timings, and failure
   status. A partial file must never be accepted as a valid cache.
6. Keep caching disabled by default. If explicitly requested and cache writing
   fails, fail the requested cache operation rather than silently reporting a
   complete reusable artifact.

## Acceptance gates

- Standalone and fused-cache modes reproduce the same ordered CBQ records and
  the same STAR hard MEX as direct FASTQ on a deterministic 100K fixture.
- Run the authorized full-slide benchmark only after the 100K gates pass.
- Report direct FASTQ time, standalone encode time plus CBQ-processing time,
  and fused FASTQ-processing-plus-cache time separately.
- Demonstrate bounded memory and deterministic output across thread counts.
- Do not launch an identical successful dataset run without new explicit user
  authorization; failed or cancelled attempts remain in provenance.

## Cancelled baseline

The cancelled SPATCh standalone-conversion attempt is retained as incomplete
provenance at:

`/storage/star-spatial-ssd-benchmark-20260729/optimization/spatch_hard_packed_cbq_fallthrough_v1`

It stopped during encoding with `exit_status=130`; STAR never ran. The partial
CBQ is not a valid input artifact.
