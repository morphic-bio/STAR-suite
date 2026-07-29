# Correction: native spatial Flex input pipeline

Date: 2026-07-29

## Finding

The accepted native spatial Flex commands used `--flexPipeline no`. Source
inspection confirmed that this bypassed the Flex producer/consumer system and
returned execution to STAR's ordinary chunk reader. The option originated as
a temporary implementation/debugging guard. It was not requested by the user,
was not justified as a production policy, and should have been called out in
the run provenance.

The resulting matrices remain valid compatibility artifacts because their
scientific inputs and downstream hard-count policy were accepted. Their
timings must be labeled **diagnostic compatibility-path timings**, not
production Flex-pipeline timings.

## Correction

Native spatial Flex now permits the producer/consumer path and treats
`--flexPipeline auto` as production:

- paired FASTQ: one reader/router per lane, separate raw-R1 spatial cache-hit
  consumers, and the remaining workers for alignment misses;
- indexed paired CBQ: fully fused indexed-range workers that claim independent
  ranges before draining alignment misses.

Every route retains the same zero-based input ordinal. H0/H1 keeps and hash
denies decode and complete raw R1 in the spatial accumulator on cache-consumer
threads; true cache misses decode raw R1 on the alignment worker and complete
it with the existing Flex alignment resolver.

For multi-lane FASTQ, the ordinal is packed as lane-major `(lane, local read)`
rather than allocated by a cross-lane race. Numeric sorting therefore matches
STAR's original lane-concatenated input order without a preliminary full FASTQ
scan. This is required because spatial clique partitioning intentionally uses
the original input order as its final deterministic tie break.

Fully fused input with alignment enabled reserves half of the worker threads
for immediate alignment-queue draining. The other half claim CBQ ranges (or
whole lanes) until input is exhausted and then switch to alignment. This
producer budget is mandatory: if every worker claims an input range, all
workers can block on a full alignment-miss queue before any reaches the
role-switch phase. Per-lane counters are atomic because indexed ranges from one
lane are processed concurrently.

Indexed CBQ planning deliberately creates up to 64 work-stealing ranges per
producer, with a minimum range size of 8,192 records. A first full-slide CBQ
attempt used only one balanced range per producer. Near the end of CRC lane
L008, one range entered a spatial-barcode fallback-heavy imaging tile while all
other producers had exhausted their ranges; throughput fell to about 139 reads
per second and the alignment queue was empty. The incomplete attempt is
preserved as rejected at
`/storage/star-spatial-ssd-benchmark-20260729/optimization/crc_hard_packed_cbq_v1`.
Oversubscribed ranges retain the same deterministic global source ordinals and
let idle producers claim work beyond a locally expensive tile.

`--flexPipeline no` is retained only to diagnose compatibility with historical
artifacts. It must not be added to a production recipe, benchmark, or relaunch
without explicit user authorization for that exact run.

## Required validation

For each input format, run one authorized 100K correctness gate before the
authorized full slide:

1. require exact hard MEX parity with the accepted compatibility artifact;
2. require exact paired-read and terminal-route conservation;
3. record binary SHA-256, immutable source revision, complete argv, inputs,
   input-format preparation, topology log, timing, and output hashes;
4. stop before the full slide on any unexplained difference.

Do not perform an independent repeat of any successful definition unless the
user explicitly authorizes that duplicate. Failed or cancelled attempts may be
relaunched and must remain in provenance.
