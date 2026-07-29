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
