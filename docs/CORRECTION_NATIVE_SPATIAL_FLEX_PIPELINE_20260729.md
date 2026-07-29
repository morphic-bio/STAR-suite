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

Indexed CBQ planning caps work-stealing ranges at 8,192 records. A first
full-slide CBQ attempt used only one balanced range per producer. Near the end
of CRC lane L008, one range entered a spatial-barcode fallback-heavy imaging
tile while all other producers had exhausted their ranges; throughput fell to
about 139 reads per second and the alignment queue was empty. A second attempt
used 1,026 ranges, but its final 207,446-record task isolated the same expensive
tile and reproduced the serialized tail. The incomplete attempts are preserved
as rejected at
`/storage/star-spatial-ssd-benchmark-20260729/optimization/crc_hard_packed_cbq_v1`
and
`/storage/star-spatial-ssd-benchmark-20260729/optimization/crc_hard_packed_cbq_v2`.
The 8,192-record ceiling distributes that tile across the producer pool while
retaining deterministic global source ordinals.

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

## Accepted validation results

The final implementation is source revision
`5e8f3619eb50d5fb8f8ce3a2497ad0b714770524`; its STAR binary SHA-256 is
`ab8e8c6a5080da3674d687d09cb175c6063555d5f72508a2b94ef0739343f387`.

- The SPATCH ovarian Flex 100K CBQ correctness gate completed at
  `/storage/star-spatial-ssd-benchmark-20260729/optimization/spatch_100k_hard_packed_cbq_v4`.
  It decoded 100,000 reads, produced 89,675 read cliques and 89,623 hard
  molecules, and its complete hard MEX is byte-identical to the accepted
  comparator. STAR wall time was 2:30.44, dominated by private reference
  loading.
- The full CRC CBQ execution completed at
  `/storage/star-spatial-ssd-benchmark-20260729/optimization/crc_hard_packed_cbq_v3`.
  It decoded 212,554,625 reads with exact lane totals, produced 79,236,462 read
  cliques and 77,819,072 hard molecules, and its complete hard MEX is
  byte-identical to the accepted FASTQ output. STAR wall time was 13:01.85 and
  maximum RSS was 48,089,616 kB.
- The comparable accepted native FASTQ production run at
  `/storage/star-spatial-ssd-benchmark-20260729/optimization/crc_hard_packed_fastq_pipeline_v2`
  took 16:57.60. Reusable CBQ input therefore reduced STAR wall time by 235.75
  seconds, or 23.17% (1.302x speedup).
- The three immutable level-0 CBQs total 9,789,019,419 logical bytes. Their
  encoders ran in parallel; the slowest lane took 2:47.68. Including that
  one-time conversion gives 15:49.53 for encode plus STAR, still 6.69% below
  direct FASTQ wall time. Subsequent analyses reuse the CBQs and pay only the
  13:01.85 STAR stage.

These successful definitions are closed to repeats without explicit user
authorization. The earlier incomplete attempts and their rejection reasons
remain in their run roots.
