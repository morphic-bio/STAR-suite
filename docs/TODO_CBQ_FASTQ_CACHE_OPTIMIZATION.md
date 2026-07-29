# TODO: Parallel FASTQ admission and durable CBQ caching

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

## General input-library direction

Do not make this optimization specific to Flex or to STAR alignment. Build a
FASTQ-aware input library beneath those consumers. Its primary contract is a
set of independent, synchronized record producers rather than one ordered
decompressed byte stream:

```text
plan = fqgzip_plan(R1, R2, shard_count)
handle[i] = fqgzip_open_shard(plan, i)
while (fqgzip_next_pair(handle[i], order_key, ...)) consume(...)
```

Each handle owns a non-overlapping half-open record interval and emits the
R1/R2 record with a deterministic source-order key. Do not require a dense
global record ordinal before the shards have been consumed: use the
lexicographic key `(lane, anchor_ordinal, local_record)` during processing and,
if a consumer requires dense ordinals, normalize it after prefix-summing the
final per-shard record counts. Anchor targets and order must not depend on
worker scheduling. The planning backend is selected from file contents rather
than filename suffixes:

- plain FASTQ: probe approximate byte offsets and advance to validated record
  boundaries;
- ordinary gzip FASTQ without an index: use `rapidgzip` speculative DEFLATE
  recovery, then validate FASTQ framing and synchronize mates by canonical
  read identifier;
- BGZF or an existing compatible auxiliary index: use exact checkpoints as a
  fast path;
- indexed CBQ: bypass FASTQ parsing and open existing logical record ranges.

Expose three surfaces over this engine:

1. an ordered stdout mode for compatibility with programs that currently call
   `gzip`, `pigz`, or `zcat`;
2. a multi-stream process mode that publishes a shard manifest for callers
   capable of parallel ingestion;
3. a zero-copy C++ record API and a small stable C ABI for in-process callers.

The ordered pipe can accelerate decompression, but necessarily serializes
record delivery again. The multi-stream and callable forms preserve the full
parallel parsing and delivery benefit.

## Index-free paired-anchor discovery

For ordinary paired `fastq.gz`, choose approximate fractional compressed
offsets independently in R1 and R2; equal compressed offsets are not expected
because mate compression ratios differ. At each probe:

1. Ask `rapidgzip` to recover a valid DEFLATE restart state and decode a small
   window from each mate.
2. Scan for complete FASTQ records. Validate the four-line structure, header
   marker, third-line `+`, sequence/quality length equality, sequence alphabet,
   and quality range.
3. Canonicalize read identifiers and find an R1/R2 intersection. Do not assume
   that compressed position alone identifies the same record.
4. Require a configurable run of consecutive matching pair identifiers (start
   with 256-1024) before accepting the boundary. Abort rather than silently
   repair a missing, duplicated, reordered, truncated, or corrupt mate.
5. Retain the recoverable compressed state plus the decoded byte/record skip
   to the confirmed pair. An anchor is not merely an uncompressed ordinal.
6. Use adjacent confirmed anchors to define exact half-open shard intervals.
   Consume any already-decoded records from the probe buffer, then continue
   from the validated restart state. Assert mate identity for every emitted
   pair and attach the deterministic source-order key.

Discover anchors independently and in parallel. Once adjacent anchors are
known, processing a shard must not rescan either file from its beginning.
Persisting the discovered anchors is optional; the durable output of the first
pass is CBQ.

This design builds on, and must cite rather than claim, generic speculative
gzip recovery in `pugz` and `rapidgzip`. The closest semantic comparison is
`mim`, which creates a small FASTQ-aware auxiliary index and supports
paired-file synchronization. The proposed distinction to evaluate is
index-free, on-demand paired synchronization coupled to immediate processing
and durable CBQ publication. Do not make a novelty claim until a complete
literature audit is recorded.

Primary comparison starting points:

- `rapidgzip`: https://github.com/mxmlnkn/rapidgzip
- `pugz`: https://arxiv.org/abs/1905.07224
- `mim`: https://pmc.ncbi.nlm.nih.gov/articles/PMC12697355/

## First-pass process-and-retain architecture

The user-visible pipeline is:

```text
ordinary paired fastq.gz
    -> synchronized parallel record producers
    -> application sink + ordered asynchronous CBQ sink
```

On the first pass, downstream work begins immediately while the same decoded
records are encoded into independent CBQ blocks. Completed blocks may be
produced out of order, but publication and the trailing index must preserve
exact source order. On successful finalization, atomically publish the CBQ as
the durable, indexed representation for future analysis. A conversion-only
command is the same engine with the application sink disabled.

Bind each CBQ to the source FASTQ identities and record in its manifest the
source hashes, exact mate record counts, canonical paired-name validation,
lane order, format and encoder versions, compression settings, block index,
and timing. Never publish a partial CBQ as reusable. The existing CBQ contract
already preserves source order, paired association, headers, qualities, and
`A/C/G/T/N` sequence content, so this work should extend the admission and
writer architecture rather than invent another durable format.

## Work order

1. Prototype the generalized shard planner/reader on plain paired FASTQ, then
   on ordinary gzip through the `rapidgzip` C++ library. Keep it independent
   of STAR mapping and CBQ encoding.
2. Optimize the standalone converter before promoting CBQ in benchmarks.
   Feed its workers from the generalized shard reader, separate block encoding
   from ordered output, encode independent blocks in parallel, preserve
   deterministic input order, and avoid per-record string allocation/copying.
   Benchmark conversion wall time, CPU utilization, compressed size, and
   byte/record parity on 100K before a full slide.
3. Add an opt-in FASTQ `process + cache` mode. Tee each already-decompressed,
   parsed read pair into an asynchronous CBQ block encoder while the same
   record continues through Flex processing. Do not decompress or parse FASTQ
   twice.
4. Use a bounded queue for backpressure, parallel block encoders, and a single
   ordered file/index writer. Pass batch/block ownership from the FASTQ reader
   instead of copying each read through synchronous `CbqWriter::add_record`.
5. Emit one CBQ per original input lane. This preserves lane/source order and
   allows multi-lane conversion and later reads to remain parallel.
6. Write to temporary files and publish atomically only after all blocks and
   the final index are complete. Record record counts, input identities,
   format version, checksum, compression settings, stage timings, and failure
   status. A partial file must never be accepted as a valid cache.
7. Keep caching disabled by default. If explicitly requested and cache writing
   fails, fail the requested cache operation rather than silently reporting a
   complete reusable artifact.

## Acceptance gates

- Standalone and fused-cache modes reproduce the same ordered CBQ records and
  the same STAR hard MEX as direct FASTQ on a deterministic 100K fixture.
- Shard concatenation reproduces the exact source record order with no gaps or
  overlaps, and every paired record identifier matches. Injected truncation,
  corruption, dropped mates, duplicate identifiers, and reordering fail
  closed.
- Plain FASTQ, ordinary gzip at several compressor/level combinations, BGZF,
  indexed and index-free rapidgzip paths, and single-/paired-end inputs are
  covered. Benchmark against `zcat`, `pigz`, rapidgzip alone, `pugz`, and
  `mim` where licensing and supported inputs permit.
- Run the authorized full-slide benchmark only after the 100K gates pass.
- Report direct FASTQ time, standalone encode time plus CBQ-processing time,
  and fused FASTQ-processing-plus-cache time separately.
- Demonstrate bounded memory and deterministic output across thread counts.
- Do not launch an identical successful dataset run without new explicit user
  authorization; failed or cancelled attempts remain in provenance.

## Possible independent software paper

If the generalized library is exact and its speedup survives diverse FASTQs,
evaluate it as a separate bioinformatics software/application paper. STAR
Suite should be one demanding real-world consumer and parallel CBQ conversion
the clean format-conversion consumer; neither should be required to use the
library. Report anchor-discovery overhead, false-anchor rejection, shard-count
scaling, first-pass time, retained-CBQ reuse time, CPU, memory, I/O, compressed
size, and failure behavior. Frame any contribution around FASTQ-aware,
index-free paired synchronization plus simultaneous processing and durable
publication—not around speculative DEFLATE recovery itself.

## Cancelled baseline

The cancelled SPATCh standalone-conversion attempt is retained as incomplete
provenance at:

`/storage/star-spatial-ssd-benchmark-20260729/optimization/spatch_hard_packed_cbq_fallthrough_v1`

It stopped during encoding with `exit_status=130`; STAR never ran. The partial
CBQ is not a valid input artifact.
