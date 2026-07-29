# TODO: In-place sharding of gzipped FASTQ files for increased performance

Status: deferred until after the preprint benchmarks. The preprint production
path remains the current direct FASTQ input. Do not substitute a new input
implementation in a benchmark unless the operator explicitly requests it.

## Headline and contribution

The headlining method is in-place sharding of gzipped FASTQ files for increased
performance.

An ordinary unindexed `fastq.gz` is a serial input bottleneck for otherwise
parallel bioinformatics programs.

The goal is to expose one existing FASTQ as multiple independent, exact
FASTQ-record producers during its first use, without requiring a preparatory
indexing, conversion, or recompression pass. Paired inputs must be synchronized
by record identity, not by assuming that equal compressed offsets represent
equal record positions.

This FASTQ sharding method is the contribution. It must be usable by arbitrary
callers and must not depend on STAR, Flex, spatial assays, or CBQ. Optional
outputs such as CBQ are downstream conveniences, not part of the novelty
claim.

Relevant existing machinery includes:

- `rapidgzip` provides general speculative recovery and parallel decompression
  for ordinary gzip streams;
- The fused Flex FASTQ producer/consumer path already owns decompressed and
  parsed paired records before application-specific routing;
- existing callback-based FASTQ readers provide practical adapter surfaces for
  other bioinformatics programs.

## General sharding-library direction

Do not make this optimization specific to Flex or to STAR alignment. Build a
FASTQ-aware input library beneath those consumers. Its primary contract is a
set of independent, synchronized record producers rather than one ordered
decompressed byte stream:

```text
plan = fastq_shard_plan(R1, R2, shard_count)
handle[i] = fastq_shard_open(plan, i)
while (fastq_shard_next_pair(handle[i], order_key, ...)) consume(...)
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
2. a multi-stream mode that publishes synchronized shard streams and their
   manifest for callers capable of parallel ingestion;
3. zero-copy C++ and small stable C APIs that open independent shard handles
   over the source file in place.

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

Discover anchors independently; discovery may be serial or parallel. Once
adjacent anchors are known, processing a shard must not rescan either file from
its beginning.
Persisting the discovered anchors is optional. Correct and efficient first-pass
operation must not require any durable output.

This design builds on, and must cite rather than claim, generic speculative
gzip recovery in `pugz` and `rapidgzip`. The closest semantic comparison is
`mim`, which creates a small FASTQ-aware auxiliary index and supports
paired-file synchronization. The proposed distinction to evaluate is
index-free, on-demand paired synchronization that directly supplies independent
record producers to an application. Do not make a novelty claim until a
complete literature audit is recorded.

Primary comparison starting points:

- `rapidgzip`: https://github.com/mxmlnkn/rapidgzip
- `pugz`: https://arxiv.org/abs/1905.07224
- `mim`: https://pmc.ncbi.nlm.nih.gov/articles/PMC12697355/

## Streamed and in-place consumption

The two primary user-visible forms end at the application:

```text
ordinary paired fastq.gz
    -> synchronized parallel record producers
    -> synchronized shard streams -> arbitrary process workers
    or
    -> independent in-place shard handles -> arbitrary in-process workers
```

Both forms use the same validated shard plan. Stream consumers receive
ordinary FASTQ records over multiple ordered shard streams. In-place consumers
open the source file through the library and iterate only their assigned record
interval without an intermediate FASTQ materialization.

## Optional downstream sinks

The library should allow optional sinks to observe the same decoded records.
CBQ is a convenient existing durable sink: it avoids inventing another format
and lets later runs bypass gzip recovery and FASTQ parsing. It is an afterthought
to the sharding method and must not be required or presented as its purpose.

If requested, encode independent CBQ blocks asynchronously while the application
processes the first pass. Completed blocks may be produced out of order, but
publication and the trailing index must preserve exact source order. Atomically
publish the CBQ only after successful finalization. A conversion-only command
is the same engine with the application sink disabled.

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
2. Expose independent record handles, an ordered compatibility stream, a
   multi-stream process surface, a C++ API, and a small C ABI. Validate adapter
   use with at least one consumer outside STAR Suite.
3. Benchmark the input engine alone before integrating it into an aligner.
   Measure anchor discovery, decompression, FASTQ parsing, mate validation,
   delivery throughput, thread scaling, CPU, memory, and I/O.
4. Integrate the same library with STAR Suite without changing mapping or
   counting policy. Preserve deterministic source order independently of worker
   scheduling and shard count.
5. As optional follow-up work, optimize the standalone CBQ converter.
   Feed its workers from the generalized shard reader, separate block encoding
   from ordered output, encode independent blocks in parallel, preserve
   deterministic input order, and avoid per-record string allocation/copying.
   Benchmark conversion wall time, CPU utilization, compressed size, and
   byte/record parity on 100K before a full slide.
6. Add an opt-in FASTQ `process + cache` mode. Tee each already-decompressed,
   parsed read pair into an asynchronous CBQ block encoder while the same
   record continues through application processing. Do not decompress or parse FASTQ
   twice.
7. Use a bounded queue for backpressure, parallel block encoders, and a single
   ordered file/index writer. Pass batch/block ownership from the FASTQ reader
   instead of copying each read through synchronous `CbqWriter::add_record`.
8. Emit one CBQ per original input lane. This preserves lane/source order and
   allows multi-lane conversion and later reads to remain parallel.
9. Write to temporary files and publish atomically only after all blocks and
   the final index are complete. Record record counts, input identities,
   format version, checksum, compression settings, stage timings, and failure
   status. A partial file must never be accepted as a valid cache.
10. Keep CBQ output disabled by default. If explicitly requested and cache writing
    fails, fail the requested cache operation rather than silently reporting a
    complete reusable artifact.

## Acceptance gates

- Shard concatenation reproduces the exact source record order with no gaps or
  overlaps, and every paired record identifier matches. Injected truncation,
  corruption, dropped mates, duplicate identifiers, and reordering fail
  closed.
- Plain FASTQ, ordinary gzip at several compressor/level combinations, BGZF,
  indexed and index-free rapidgzip paths, and single-/paired-end inputs are
  covered. Benchmark against `zcat`, `pigz`, rapidgzip alone, `pugz`, and
  `mim` where licensing and supported inputs permit.
- The same library demonstrates correct integration with multiple independent
  bioinformatics consumers; STAR Suite is the demanding end-to-end example,
  not the only usable host.
- Optional standalone and fused CBQ modes reproduce the same ordered CBQ
  records and the same STAR results as direct FASTQ.
- Run the authorized full-slide benchmark only after the 100K gates pass.
- Report the input-only result before any optional CBQ timings. Report direct
  FASTQ time, standalone encode time plus CBQ-processing time, and fused
  FASTQ-processing-plus-cache time separately.
- Demonstrate bounded memory and deterministic output across thread counts.
- Do not launch an identical successful dataset run without new explicit user
  authorization; failed or cancelled attempts remain in provenance.

## Possible independent software paper

If the generalized library is exact and its speedup survives diverse FASTQs,
evaluate it as a separate bioinformatics software/application paper. STAR
Suite should be one demanding real-world consumer, not a required host. Report
anchor-discovery overhead, false-anchor rejection, shard-count scaling,
first-pass time, CPU, memory, I/O, and failure behavior. Frame the contribution
around FASTQ-aware, index-free paired synchronization and independent producer
streams/handles—not around speculative DEFLATE recovery itself. The title and
abstract should lead with in-place sharding of gzipped FASTQ files for increased
performance. "In-place" distinguishes the method from approaches
that first create new files or convert the input into a serial stream. Present
optional CBQ output only as a convenient demonstration that the one unavoidable
gzip pass can also leave an existing durable indexed representation.

## CBQ afterthought: cancelled conversion baseline

The cancelled SPATCh standalone-conversion attempt is retained as incomplete
provenance at:

`/storage/star-spatial-ssd-benchmark-20260729/optimization/spatch_hard_packed_cbq_fallthrough_v1`

It stopped during encoding with `exit_status=130`; STAR never ran. The partial
CBQ is not a valid input artifact.
