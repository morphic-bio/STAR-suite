# Runbook: Broaden CBQ Parallel Readers in STAR Core and Flex

Date: 2026-05-31

## Goal

Indexed CBQ input is now mature enough to use independent range readers wherever
input order is not part of the output contract. This runbook broadens the first
STAR-core implementation and removes the most visible Flex range-reader fallback.

The guiding rule is:

- Use indexed CBQ range readers by default for CBQ input with `runThreadN > 1`.
- Keep the shared/whole-lane reader only for true semantic blockers or missing
  indexes.

## Scope

### STAR core

Enable `--readFilesCbqRangeMode auto` for the order-independent alignment
outputs we can validate without changing byte-order expectations:

- `--outSAMtype None`
- BAM `SortedByCoordinate`

Also support `--readMapNumber` by truncating the range plan to the requested
logical record count.

Keep these gates:

- `--outSAMorder PairedKeepInputOrder`
- direct SAM output
- BAM `Unsorted` side output
- `--emitYNoY*`
- batch and SLAM per-file/autotrim modes
- two-pass/runtime SJ insertion
- two-stage SJ filtering
- missing CBQ indexes

Rationale: `PairedKeepInputOrder` and Y/noY sidecars require deterministic
logical chunk ordering. BAM coordinate sorting is order-independent. Direct SAM
and unsorted BAM are not semantically sorted, but existing smokes and downstream
checks compare their body order exactly, so keep them on the shared reader until
they get ordered range flush/concatenation semantics. Flex unsorted BAM is
especially brittle because CB/UB tag injection depends on the unsorted BAM
output contract; revisit that path separately before allowing range tasks there.

### Flex

Use CBQ range tasks in the regular fully-fused Flex pipeline, not only the
no-genome/no-align path.

Keep these gates:

- non-fully-fused CBQ Flex pipeline modes
- missing CBQ indexes
- genome-backed auxiliary modes already rejected by the Flex no-genome guard

Rationale: `FlexCbqRangeTask` already carries lane id and logical first record,
and `processCbqModuleRecords()` already supports deterministic read ids. The
same worker-local reader can feed either count-only processing or the alignment
queue.

## Implementation Plan

1. Loosen `cbqCoreRangeGateReject()`.
   - Accept `None` and coordinate-sorted `BAM`.
   - Keep direct `SAM` and `BAM Unsorted` on the shared reader for now.
   - Remove the initial `--readMapNumber` rejection.

2. Truncate STAR-core range planning for `--readMapNumber`.
   - Inspect full lane counts.
   - Plan only `min(total_records, readMapNumber)` records.
   - Preserve lane-global ordinals for read ids.

3. Flush residual direct SAM buffers at the end of each CBQ range worker.
   - Normal shared-reader chunks mark the final chunk with `noReadsLeft`.
   - Range workers can finish with a partial per-thread SAM buffer and no final
     input chunk, so explicitly flush it after all claimed range tasks.

4. Prepare Flex CBQ range tasks in the regular fully-fused pipeline.
   - Reuse `flexPrepareCbqRangeTasks()`.
   - Log active/fallback state consistently with the no-genome path.

5. Remove the `noAlign` condition from Flex range task consumption.
   - Range-fed reads that miss the hash screen already enqueue
     `EnrichedPacket`s when alignment is enabled.

## Validation

Minimum local validation:

```bash
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR cbq-ordered-encoder cbq-star-adapter-harness

OUT_ROOT=/tmp/star_suite_cbq_solo_e2e_smoke_core_flex_range_$(date +%Y%m%dT%H%M%SZ) \
  tests/run_cbq_solo_e2e_smoke.sh

OUT_ROOT=/tmp/star_suite_cbq_star_input_smoke_core_flex_range_$(date +%Y%m%dT%H%M%SZ) \
  tests/run_cbq_star_input_smoke.sh

OUT_ROOT=/tmp/star_suite_cbq_ynoy_smoke_core_flex_range_$(date +%Y%m%dT%H%M%SZ) \
  tests/run_cbq_ynoy_smoke.sh
```

Additional validation when the public tiny Flex fixture dependencies are
available:

```bash
OUT_ROOT=/tmp/star_suite_cbq_flex_tiny_public_core_flex_range_$(date +%Y%m%dT%H%M%SZ) \
  tests/run_cbq_flex_tiny_public_smoke.sh
```

This public smoke uses unsorted BAM and should remain on the standard/shared
CBQ path until Flex unsorted BAM ordering is made explicit. A separate targeted
fully-fused Flex smoke should use the same tiny fixture with a prebuilt H0/H1
hash cache, `--outSAMtype None`, `--flexPipeline yes`,
`--flexPipelineNTriage 0`, `--flexPipelineNSolo 0`, and default
`--flexNoAlign 0`.

Expected logs:

- STAR-core CBQ range output modes log
  `CBQ indexed range reader: active ...`.
- Y/noY still logs fallback in `auto` mode.
- Public Flex unsorted-BAM CBQ logs shared-reader fallback.
- Fully-fused Flex CBQ with `--outSAMtype None` logs
  `Flex pipeline: runThreadN=...` and `Flex CBQ range: active ...`.

## Follow-Ups

- Add ordered CBQ Y/noY chunk flushing, then remove the Y/noY range gate.
- Add deterministic chunk ids for `PairedKeepInputOrder`, then remove that gate.
- Revisit Flex unsorted BAM tag injection and define an explicit ordered range
  flush contract before enabling CBQ range tasks for unsorted BAM outputs.
- Decide whether missing CBQ indexes should remain a fallback or become fatal in
  Flex production presets.
