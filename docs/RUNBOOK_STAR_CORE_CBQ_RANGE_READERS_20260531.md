# Runbook: STAR Core CBQ Indexed Range Readers

Date: 2026-05-31

## Goal

Add an opt-in/auto-gated STAR core path that reads indexed CBQ inputs through
independent logical record ranges instead of serializing all worker threads
behind the shared input mutex. This targets STAR/STARsolo runs where FASTQ/CBQ
decode is a bottleneck and the output path does not require input-order chunk
submission.

This complements the existing PF and FLEX range-reader work. It does not change
FASTQ input handling.

## Non-goals

- No parallel FASTQ reader.
- No range mode for CBQ files without a `CBQINDEX` footer.
- No range mode for order-dependent outputs in the first implementation.
- No changes to read mapping, STARsolo counting, barcode correction, or BAM
  sorting semantics.
- No SLAM per-file, batch-resume, or two-pass special-case support in the first
  implementation.

## User-Facing Gate

Add:

```text
--readFilesCbqRangeMode auto|off|range
```

Default: `auto`.

Behavior:

- `off`: always use the existing shared CBQ reader.
- `auto`: use indexed range readers only when the run is known to be
  order-independent and all input files have CBQ indexes; otherwise fall back to
  the existing shared reader with a log message.
- `range`: require indexed range readers; unsupported configurations or missing
  indexes are fatal.

Initial activation requirements:

- `--readFilesType Binseq SE|PE`
- `--runThreadN > 1`
- indexed CBQ files
- `--outSAMtype None`
- no `--emitYNoY`, `--emitYNoYFastq`, or `--emitYNoYCbq`
- no `--outSAMorder PairedKeepInputOrder`
- no `--readMapNumber` cap
- no batch mode, SLAM per-file skip, or SLAM auto-trim detection pass

The intentionally narrow initial surface is the STARsolo/MEX-only case. Sorted
BAM can be enabled later after explicit parity checks because coordinate sorting
should make input order irrelevant, but it is not part of this first gate.

## Implementation Plan

1. Add a copyable runtime range state to `Parameters`.
   - `readFilesCbqRangeMode`
   - `cbqRangeActive`
   - `cbqRangeFallbackReason`
   - vector of `{laneIndex, firstRecord, recordCount, globalFirst, taskOrdinal}`
   - shared atomic task cursor

2. Prepare range tasks in `Parameters::openReadsFiles()`.
   - Configure the normal `CbqInputModule`.
   - If the gate accepts range mode, open each lane with
     `CbqInputModule::open_range(lane, 0, UINT64_MAX)` to validate the index and
     obtain lane record counts.
   - Partition the logical concatenation of all lanes into roughly
     `runThreadN` global ranges. A global range may split at lane boundaries, so
     the task count can exceed `runThreadN`.
   - In fallback mode, open the existing shared reader exactly as before.

3. Add a `ReadAlignChunk::processChunks()` CBQ range branch.
   - Each worker atomically claims the next prepared range task.
   - Each worker constructs and opens its own `CbqInputModule` for that lane
     range.
   - Fill `cbqStarChunk` from the private reader without taking
     `g_threadChunks.mutexInRead`.
   - After appending each batch into the owned chunk, rebase copied read
     ordinals to the task's `globalFirst + localOffset + 1`.
   - Assign unique range-path chunk IDs from a shared atomic counter. These are
     for tracing only because range mode is gated away from order-dependent
     output surfaces.
   - Call the existing `mapCbqChunk()` path.

4. Preserve the existing path for all unsupported modes.
   - Shared-reader behavior remains unchanged when the gate is off or falls
     back.
   - Forced `range` mode exits early with an actionable message when the gate
     rejects the run.

## Validation

Clean rebuild first:

```bash
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR cbq-ordered-encoder cbq-star-adapter-harness
```

Focused smokes:

```bash
OUT_ROOT=/tmp/star_suite_cbq_star_input_smoke_core_range_$(date -u +%Y%m%dT%H%M%SZ) \
  tests/run_cbq_star_input_smoke.sh

OUT_ROOT=/tmp/star_suite_cbq_solo_e2e_smoke_core_range_$(date -u +%Y%m%dT%H%M%SZ) \
  tests/run_cbq_solo_e2e_smoke.sh

OUT_ROOT=/tmp/star_suite_cbq_ynoy_smoke_core_range_guard_$(date -u +%Y%m%dT%H%M%SZ) \
  tests/run_cbq_ynoy_smoke.sh
```

Expected checks:

- Existing CBQ STAR/STARsolo parity still passes.
- A new `--outSAMtype None --readFilesCbqRangeMode range` smoke run logs range
  activation and completes.
- Y/noY smoke remains on the shared reader or rejects forced range mode.
- No generated outputs are committed.

## Benchmark Follow-up

After smokes pass, benchmark the same dataset with:

- shared reader: `--readFilesCbqRangeMode off`
- indexed range readers: `--readFilesCbqRangeMode range`

Use fresh output directories and do not run benchmark jobs in parallel from the
same checkout. Document wall time, max RSS, thread count, input type, CBQ file
locations, and whether output is MEX-only or includes BAM.
