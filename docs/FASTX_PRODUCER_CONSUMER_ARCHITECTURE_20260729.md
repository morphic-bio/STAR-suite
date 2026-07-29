# Ordered FASTX producer/consumer architecture

Status: implemented and synthetic-parity tested on 2026-07-29. No biological
dataset or full-slide benchmark was launched as part of this implementation.

## Purpose

Standard STAR mapping workers historically acquired one global input mutex,
decompressed and parsed their next FASTQ chunk, released the mutex, and then
mapped it. Mapping overlapped input, but gzip decode and FASTQ parsing remained
serial and could leave mapping threads idle.

`input/FastxProducerPool.{h,cpp}` moves source reading ahead of that mutex. A
bounded set of producers independently opens lanes, decompresses and parses
records, and publishes owned batches. STAR mapping workers remain the
consumers. The existing chunk builder still assigns global read ordinals and
constructs STAR's internal FASTQ-equivalent buffers, so mapping and counting
policy are unchanged.

This is a general `InputModule`; it contains no GEX, Flex, spatial, alignment,
or counting logic. Integrated spatial GEX is the first automatic production
consumer. Other FASTX callers can opt in while their regression surfaces are
audited.

## Ordering and correctness contract

- Producers may decompress lanes concurrently.
- Delivery is deterministic and lane-major: all records from lane 0, then lane
  1, and so on, in source order within each lane.
- Bounded per-lane queues provide backpressure; later lanes cannot consume
  unbounded memory while an earlier lane is mapped.
- A producer parse/open failure is propagated to the consumer and wakes all
  blocked producers. Partial input is never reported as a successful end.
- Closing a partially consumed input wakes and joins producers, including
  producers blocked on full queues.
- The pool delegates source-plan validation and record parsing to
  `FastxInputModule`, so pooled and sequential paths share FASTA/FASTQ,
  internal-gzip, helper-command, manifest, quality, and lane semantics.

The gzip line reader now uses a 1 MiB zlib input buffer and 64 KiB line blocks
instead of calling `gzgetc` once per byte. Lines larger than the block are
assembled without truncation.

## User controls

`--readFilesFastxProducerConsumer auto|off|on`

- `auto` (default): enable for `--soloSpatialGexIntegrated yes` and retain the
  sequential reader elsewhere.
- `off`: force the sequential compatibility/diagnostic reader.
- `on`: enable the ordered producer pool for any FASTX mapping run.

`--readFilesFastxProducerThreads N`

- `0` (default): `min(4, lane_count)`.
- positive values: requested producer count, bounded by the lane count.

Each producer currently reads all mates for one claimed lane. This preserves
paired synchronization. It does not yet split one gzip member into independent
in-place shards; a single-lane input therefore has one decompression producer.
That later capability plugs in below the same `InputModule` boundary.

## Verification

The following tests are deterministic and use only synthetic data:

- `tests/run_fastx_input_harness_smoke.sh`
  - sequential versus pooled plain and gzip records;
  - single-, paired-, and three-mate inputs;
  - multiple lanes and manifests;
  - exact lane, ordinal, sequence, and quality digest parity;
  - early close while bounded queues are full;
  - producer-side truncated FASTQ failure without deadlock.
- `tests/run_fastx_contract_star_smoke.sh`
  - full STAR input-to-alignment execution with two mapping workers;
  - sequential versus pooled input-count and mapping-category parity;
  - exact header-stripped SAM alignment parity;
  - production module/thread-count logging.

Clean build gate:

```bash
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR fastx-input-harness
bash tests/run_fastx_input_harness_smoke.sh
bash tests/run_fastx_contract_star_smoke.sh
```

Before enabling `auto` for another STAR Suite pipeline, run its existing
synthetic and authorized dataset parity gates with `on` versus `off`. An
identical successful dataset run still requires explicit user authorization.

## Precursor to in-place gzip FASTQ indexing

The future in-place reader should replace a lane producer with multiple exact,
synchronized shard producers while retaining the same bounded ownership,
failure propagation, and consumer boundary. The ordered STAR adapter can merge
shards by deterministic source key. Consumers that do not require serial
delivery can open the independent shard handles directly.

The primary external demonstration should be `bwa-mem2`:

1. present synchronized shard streams or callable shard handles as ordinary
   paired FASTQ records;
2. give each `bwa-mem2` worker a disjoint half-open record interval;
3. merge SAM/BAM by deterministic source key where ordered output is required;
4. prove read-name, flag, reference position, CIGAR, MAPQ, and optional-tag
   identity against ordinary serial FASTQ input, allowing only documented
   differences caused by `bwa-mem2` itself;
5. report input-only and end-to-end scaling separately.

That demo is valuable precisely because it is outside STAR Suite: it tests the
claim that in-place gzipped-FASTQ indexing and processing is a reusable input
facility rather than a STAR-specific optimization. CBQ remains only an
optional durable sink after the first decode.
