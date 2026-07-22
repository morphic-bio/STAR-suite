# Runbook: Gzip / Pipeline Iteration Strategy (2026-03-26)

## Goal

Establish a disciplined way to iterate on FASTQ decompression and packetization
strategy without conflating gzip costs with unrelated mapping or Solo changes.

This runbook is for two related but distinct tracks:

1. **Flex first**: use the existing multi-stage / fused pipeline work as the
   fast iteration surface for decompression experiments.
2. **GEX second**: once a strategy looks good on Flex, port the idea into a
   separate GEX-only harness and validate it there.

## Current Conclusions

### UCSF GEX-only

On the current UCSF GEX-only benchmark surface, external `zcat` is the fastest
validated option so far:

- `zcat` run:
  `/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_optimized_current_zcat_20260326/`
- native gzip run:
  `/storage/solo_overnight_20260326/ucsf_gexonly_no_bam/star_optimized_current_retry2/`

Working assumption for new work:

- treat **external `zcat` as the baseline**
- treat internal gzip as the optimization target

### Flex

Flex is the better place to iterate on decompression strategy because we already
have queue-aware pipeline modes and prior benchmark surfaces showing where the
decompression / routing stages saturated.

Best historical references:

- pipeline sweep artifact root:
  `/storage/flex_pipeline_sweep_20260320_230900/`
- fused pipeline summary:
  `docs/HANDOFF_FLEX_PIPELINE_BENCHMARK_SUMMARY_20260323.md`
- pipeline v1 implementation summary:
  `docs/HANDOFF_FLEX_PIPELINE_V1_STATUS_20260320.md`

## Principle

Do **not** mix these experiments with broad algorithmic changes.

For this runbook, only change one of the following at a time:

- how gzip bytes are read
- how FASTQ records are delimited
- how many records are grouped into a packet
- how decompressed packets are handed to fused or downstream consumers
- queue depth / batching policy

Hold constant:

- dataset
- thread count
- STAR/Flex/Solo counting semantics
- hash-screen / filtering logic
- output surface

## Success Metric

Primary signal:

- **downstream queue occupancy**

Interpretation:

- if the downstream queue is often near-empty, decompression / parsing is not
  feeding consumers fast enough
- if the downstream queue is pinned full, decompression is outrunning the next
  stage and another stage is the bottleneck
- if the queue stays active but not saturated, decompression and downstream work
  are closer to balance

For this track, the preferred state is:

- queue occupancy is **below hard saturation**
- queue starvation is **rare**
- throughput improves or at least does not regress

Secondary metrics:

- `delta10s` / reads-per-second
- queue occupancy over time
- CPU utilization by stage
- wall time for fixed-duration pilot runs
- full-run wall time only after a pilot shows improvement

## Existing Signals To Reuse

From the historical Flex pipeline work:

- `PIPELINE_STATS`
- `readerQ`
- `soloQ`
- `alignQ`
- `delta10s`

Examples:

- v1 sweep summary:
  `/storage/flex_pipeline_sweep_20260320_230900/summary.tsv`
- per-run queue telemetry:
  `/storage/flex_pipeline_sweep_20260320_230900/t2s2/pipeline_stats.txt`
  `/storage/flex_pipeline_sweep_20260320_230900/t4s2/pipeline_stats.txt`

Important historical lesson:

- in v1, `readerQ` pinned at `255/256` showed triage was the bottleneck, not
  decompression
- in fused runs, `alignQ` saturation showed the upstream fused reader path was
  keeping up and alignment had become the bottleneck

## Flex Iteration Surface

Use Flex first because the pipeline already supports staged and fused forms.

Relevant modes:

1. **Separate stages**
   - `--flexPipeline yes`
   - `--flexPipelineNTriage > 0`
   - `--flexPipelineNSolo > 0`
2. **Fused reader+router**
   - `--flexPipeline yes`
   - `--flexPipelineNTriage 0`
   - `--flexPipelineNSolo > 0`
3. **Fully fused**
   - `--flexPipeline yes`
   - `--flexPipelineNTriage 0`
   - `--flexPipelineNSolo 0`

Recommended baseline order:

1. external `zcat`
2. current internal gzip line-oriented reader
3. candidate internal gzip strategies below

## Candidate Strategy Matrix

Start with these strategies, in this order:

1. **Baseline external helper**
   - keep `--readFilesCommand zcat`
   - no semantic changes
   - this is the control

2. **Current internal gzip, line-oriented**
   - current `gzgets()`-style record assembly
   - establishes the present internal-gzip gap

3. **Buffered gzip block reads + delimiter scan**
   - read larger gzip chunks into a buffer
   - scan for `\n` delimiters directly
   - avoid repeated line-by-line function overhead

4. **Buffered gzip block reads + FASTQ field scanning**
   - parse 4 FASTQ lines from a larger block
   - keep packet assembly local to the lane thread
   - still emit one read-record at a time

5. **Small read-packet batching**
   - bundle several reads before handing off to downstream consumers
   - goal: reduce queue / mutex / condition-variable overhead

6. **Queue of several lines or several reads for fused consumers**
   - only after the buffered parser is stable
   - batch size must be varied explicitly, not hidden behind other changes

Do **not** start with:

- cross-stage semantic refactors
- alignment changes
- Solo hash or collapse changes
- mixed gzip + packetization + scheduling rewrites in one patch

## Recommended Flex Harness

Create a dedicated Flex decompression harness rather than reusing paper scripts.

Target properties:

- fixed dataset: full JAX SC2300771 or a stable 2-minute pilot subset
- fixed thread count
- fixed output mode: no BAM
- parameterized decompression mode
- parameterized packet/batch size
- authoritative telemetry dump every 10s

Suggested new harness:

- `tests/run_flex_gzip_strategy_sweep.sh`

Minimum knobs:

- `MODE=zcat|internal_line|internal_buffered|internal_buffered_batch`
- `THREADS`
- `BATCH_READS`
- `QUEUE_CAPACITY`
- `RUN_SECONDS` for pilot mode
- `OUT_ROOT`

Minimum outputs:

- `RUN_COMMAND.sh`
- `pipeline_stats.txt`
- `stderr.log`
- `Log.out`
- `Log.final.out`
- one-row summary TSV with throughput and queue occupancy

## Recommended GEX Harness

Do **not** reuse the Flex harness for GEX.

GEX should get its own benchmark wrapper because the downstream path and queue
shape are different.

Suggested new harness:

- `tests/run_gex_gzip_strategy_sweep.sh`

Start only after a Flex strategy looks promising.

GEX baseline:

- current optimized UCSF GEX-only `zcat` surface
- no BAM
- `GeneFull`
- same reference, same whitelist, same threads, same sample

Key question for GEX:

- does the candidate strategy keep the downstream processing queue from
  starvation while preserving the current `zcat` wall-time advantage?

## Experiment Loop

For each candidate:

1. Run the Flex pilot harness against `zcat` baseline.
2. Compare throughput and queue occupancy against current internal gzip.
3. Reject immediately if throughput regresses materially and queue behavior is
   not better.
4. If the pilot looks better, run a longer Flex validation.
5. Only then port the same idea into the GEX harness.
6. Compare against the UCSF GEX-only `zcat` control.

## Acceptance Criteria

A candidate is worth carrying forward only if it meets all of:

- better or equal throughput on Flex pilot
- queue occupancy is healthier than current internal gzip
- no obvious CPU waste or lock contention spike
- no output / parity regression
- reproducible behavior across repeated runs

For GEX promotion:

- candidate must beat or closely match the current GEX `zcat` baseline
- if it is still slower than `zcat`, do not replace the benchmark default

## Benchmark Hygiene

Follow repository benchmark rules:

- serialize benchmark jobs on the machine
- use fresh outdirs
- use wrapper-written summaries where present
- do not infer completion from partial logs

Because the point here is queue behavior, keep all competing benchmarks off the
host during timed runs.

## First Concrete Actions

1. Keep `zcat` as the published benchmark baseline for now.
2. Build a dedicated **Flex gzip strategy sweep harness** with queue telemetry.
3. Reproduce three control points before changing code:
   - v1 staged queueing
   - fused reader+router
   - fully fused
4. Add one new strategy only:
   - internal gzip block-read + delimiter scan
5. Compare queue occupancy and throughput against `zcat`.
6. Only if that wins, build the separate GEX harness and port the idea.
