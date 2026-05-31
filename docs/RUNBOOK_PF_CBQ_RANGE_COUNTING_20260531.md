# Runbook: PF CBQ Range Counting

Date: 2026-05-31

Status: implemented in the STAR-suite worktree. The prerequisite indexed CBQ
range reader, direct decode benchmark gate, direct PF range-counting path, and
production `pf-multi` mode gate are in place.

Related docs:

- `docs/RUNBOOK_CBQ_PARALLEL_PF_FLEX_OPTIMIZATIONS_20260530.md`
- `docs/RUNBOOK_PROCESS_FEATURES_CBQ_NATIVE.md`
- `docs/RUNBOOK_CBQ_ORDERED_ENCODER.md`

## Goal

Wire CBQ indexed range readers directly into `process_features` counting so
each worker thread owns:

```text
CBQ logical record range -> per-thread PF counting state
```

This removes the current conservative path:

```text
single CBQ reader -> pf_process_record_views() -> shared queue -> PF consumers
```

The full-lane direct reader sweep showed that the CBQ range reader itself is
fast enough:

| Mode | 1 thread | 4 threads | 8 threads | 16 threads | 32 threads |
| --- | ---: | ---: | ---: | ---: | ---: |
| count-only | 5.72 s | 1.39 s | 0.86 s | 0.52 s | 0.36 s |
| sequence materialized | 6.56 s | 1.76 s | 0.93 s | 0.64 s | 0.41 s |

The next bottleneck to remove is the queue/stream handoff, not CBQ decoding.

## Implemented Priority Order

The eight PF suggestions were prioritized into this implementation order:

1. Add the worker-owned PF direct range API.
2. Keep the queue-backed CBQ path as the production fallback.
3. Add indexed CBQ range reads with one reader per worker.
4. Partition multi-lane CBQ inputs by logical record order.
5. Add `cbq-direct-pf` harness coverage and canonical large-output parity.
6. Promote the direct path into `pf-multi` behind `--crAssignCbqMode`.
7. Gate unsafe settings before direct work starts.
8. Record mode/fallback metadata and benchmark outputs in docs/artifacts.

The larger copy-elimination and shared consumer refactor remains a follow-on
because it changes the hottest PF record-processing code without changing the
current direct-range behavior.

## Non-Goals

- Do not add a second CBQ parser inside `process_features`.
- Do not materialize temporary FASTQ files.
- Do not change feature matching, barcode correction, UMI deduplication, or
  output semantics.
- Do not optimize Chromap-suite/libchromap in this runbook.
- Do not enable order-dependent legacy rescue on the range path until it has a
  deterministic design.

## Current Code Shape

Important files:

- `core/legacy/source/input/CbqInputModule.{h,cpp}`
  - now exposes `open_range(lane, first_record, record_count)`;
  - parses the tail `CBQINDEX`;
  - returns trimmed `CbqReadBatchView` slices.
- `core/legacy/source/input/cbq_pf_adapter_harness.cpp`
  - has `--mode cbq-direct-decode`;
  - proves independent CBQ readers can process logical ranges exactly once.
- `core/features/process_features/src/assignBarcodes.c`
  - `consume_reads()` owns almost all per-record PF counting logic;
  - it currently pulls records from `fastq_reader_set` ring buffers.
- `core/features/process_features/src/pf_api.c`
  - `pf_process_record_views()` copies borrowed records into the same ring
    buffer and starts `consume_reads()` threads;
  - finalization already merges per-thread hashes/stats and writes outputs.
- `core/legacy/source/PfMultiAssign.{h,cpp}`
  - resolves explicit `.cbq` sources or directories from feature-library
    `fastqs`;
  - selects `cbq_stream` or `cbq_range`;
  - writes `cbq_mode_requested`, `cbq_mode_effective`, and
    `cbq_mode_fallback_reason` to `assignBarcodes.api_run.txt`.
- `core/legacy/source/PfMultiProcess.cpp`
  - exposes `--crAssignCbqMode auto|stream|range`;
  - logs the effective input format and CBQ mode after each feature library.

The implementation should extract the queue-independent record-processing body
from `consume_reads()` and call it from both the existing queue consumer and the
new direct CBQ range API.

## Design

Add a small PF direct-worker API. STAR-suite C++ keeps responsibility for CBQ
reading and range partitioning. PF keeps responsibility for counting state,
merge, and final output.

Suggested public C API:

```c
typedef struct pf_direct_range_job pf_direct_range_job;

pf_error pf_direct_range_begin(pf_context *ctx,
                               const char *output_dir,
                               const char *sample_name,
                               int nworkers,
                               int nreaders,
                               pf_direct_range_job **job_out);

pf_error pf_direct_range_process_record_views(pf_direct_range_job *job,
                                              int worker_id,
                                              const pf_read_record_view *records,
                                              size_t n_records);

pf_error pf_direct_range_end(pf_direct_range_job *job,
                             pf_stats *stats_out);

void pf_direct_range_abort(pf_direct_range_job *job);
```

Rules:

- `worker_id` is stable and unique per caller thread: `0..nworkers-1`.
- Each worker has its own `statistics`, `data_structures`, memory pools,
  line scratch, bootstrap replay buffer, and hot hash state.
- The job holds `g_pf_runtime_mutex` from `begin` until `end`/`abort`, matching
  the current stream API and preventing process-global PF state bleed.
- `pf_direct_range_process_record_views()` is thread-safe only across distinct
  `worker_id` values in the same job.
- The first implementation supports paired feature mode (`nreaders=2`).
  Add `nreaders=3` only after paired parity is stable.

## Implementation Steps

### 1. Extract PF Per-Record Consumer Logic

Refactor `consume_reads()` in `assignBarcodes.c` into three internal helpers.

Suggested structs/functions:

```c
typedef struct pf_consumer_state {
    int thread_id;
    const sample_args *sample_args;
    statistics *stats;
    data_structures *hashes;
    memory_pool_collection *pools;
    char *lines_buffer;
    char *lines[6];
    char **barcode_lines;
    char **forward_lines;
    char **reverse_lines;
    seq_hash_t thread_hot_d0;
    int thread_hot_d0_active;
    unsigned int hot_check_counter;
    bootstrap_replay_buf_t replay_buf;
    int replay_buf_drained;
    int permit_batch_count;
    uint64_t permit_batch_wait_ns;
    double permit_batch_start_sec;
    uint64_t permit_batch_work_bytes;
} pf_consumer_state;

int pf_consumer_state_init(pf_consumer_state *state,
                           const fastq_processor *processor,
                           int nreaders);

int pf_consumer_process_lines(pf_consumer_state *state,
                              const char *barcode_seq,
                              const char *barcode_qual,
                              const char *feature_seq,
                              const char *feature_qual,
                              const char *feature_seq2,
                              const char *feature_qual2);

void pf_consumer_state_finish(pf_consumer_state *state);
```

Keep `consume_reads()` behavior unchanged by making its queue loop:

1. Initialize `pf_consumer_state`.
2. Copy one record from the ring buffer into local line buffers.
3. Call `pf_consumer_process_lines()`.
4. Finish state on exit.

Acceptance for this refactor alone:

```bash
make -C core/features/process_features lib
make -C core/legacy/source -j8 cbq-pf-adapter-harness

OUT_ROOT=/tmp/star_suite_cbq_pf_adapter_smoke_extract_$(date -u +%Y%m%dT%H%M%SZ) \
  tests/run_cbq_pf_adapter_smoke.sh
```

No output should change.

### 2. Add Direct Range Job Runtime In PF

Add `pf_direct_range_job` to `pf_api.c`. It can mirror
`pf_record_stream` but omit queue state and consumer pthreads:

```c
struct pf_direct_range_job {
    pf_context *ctx;
    sample_args sample_args;
    statistics *stats;
    data_structures *hashes;
    memory_pool_collection **pools;
    fastq_processor *processors;
    pf_consumer_state *consumer_states;
    char sample_directory[FILENAME_LENGTH];
    double start_time;
    int nworkers;
    int nreaders;
    int closed;
    int failed;
};
```

`pf_direct_range_begin()` should:

1. Validate `ctx`, `output_dir`, `sample_name`, `nworkers`, `nreaders`.
2. Lock `g_pf_runtime_mutex`.
3. Call `pf_apply_context_globals(ctx)`.
4. Create the sample output directory like `pf_process_records_begin()`.
5. Allocate `stats`, `hashes`, `pools`, and `processors` for `nworkers`.
6. Populate `sample_args` exactly as the stream API does.
7. Initialize one `pf_consumer_state` per worker.

`pf_direct_range_process_record_views()` should:

1. Validate `worker_id` and record spans with the same helpers used by
   `pf_process_record_views()`.
2. Copy sequence/quality spans into the worker state's fixed line buffers.
3. Default missing qualities with `pf_default_quality()`.
4. Call `pf_consumer_process_lines()`.

`pf_direct_range_end()` should:

1. Finish each `pf_consumer_state`.
2. Merge workers with `merge_process_feature_thread_data()`.
3. Call `finalize_processing()` unless `probe_only`.
4. Fill `pf_stats`.
5. Free job runtime state.
6. Unlock `g_pf_runtime_mutex`.

`pf_direct_range_abort()` must free state and unlock the runtime mutex.

### 3. Add Harness Mode `cbq-direct-pf`

Extend `core/legacy/source/input/cbq_pf_adapter_harness.cpp`:

```text
--mode cbq-direct-pf
--readFilesIn paired.cbq
--consumerThreads N
--featureOffset N
```

Implementation shape:

1. Configure PF as today.
2. Inspect total records with `CbqInputModule::open_range(..., UINT64_MAX)`.
3. Call `pf_direct_range_begin(ctx, outputDir, sampleName, N, 2, &job)`.
4. Partition logical records into `N` ranges.
5. Spawn `N` C++ workers:
   - each worker opens its own `CbqInputModule`;
   - each worker calls `open_range(0, first, count)`;
   - each worker loops `next_batch()`;
   - each worker materializes barcode and feature sequences into reused local
     scratch;
   - each worker passes `pf_read_record_view` batches to
     `pf_direct_range_process_record_views(job, worker_id, records, n)`.
6. Join workers.
7. Call `pf_direct_range_end(job, &stats)`.

Keep `--mode cbq` as the existing queue-backed control path.

### 4. Guardrails

Initial direct-PF range mode rejects or requires explicit fallback for:

- `legacy_cb_rescue=1`;
- no explicit feature offset when feature-mode bootstrap is enabled;
- `nreaders=3`;
- missing CBQ index;
- sequence length exceeding `LINE_LENGTH - 1`;
- quality length mismatch.

Multiple CBQ lanes are supported as an ordered file list or sorted `.cbq`
directory. Ranges are planned over the logical concatenation of those lanes.

For the current perturb benchmark, always use:

```text
--featureOffset 31
--runEmptyDrops omitted
--runHeatmaps omitted
```

The exact feature offset avoids feature-mode bootstrap and keeps the first
direct-PF implementation focused on reader/consumer mechanics.

### 5. Parity Gates

Synthetic smoke:

```bash
make -C core/legacy/source -j8 cbq-pf-adapter-harness

OUT_ROOT=/tmp/star_suite_cbq_pf_adapter_smoke_direct_pf_$(date -u +%Y%m%dT%H%M%SZ) \
  tests/run_cbq_pf_adapter_smoke.sh
```

Extend the smoke to run:

```bash
"$PF_BIN" \
  --mode cbq-direct-pf \
  --readFilesIn "$CBQ" \
  --whitelist "$WHITELIST" \
  --featureRef "$FEATURES" \
  --outputDir "$OUT_ROOT/cbq_direct_pf" \
  --sampleName sample \
  --barcodeLength 16 \
  --umiLength 12 \
  --consumerThreads 2 \
  --readBufferLines 128
```

Compare direct-PF output against FASTQ and queue-backed CBQ:

- `barcodes.txt`;
- `features.txt`;
- `feature_sequences.txt`;
- `matrix.mtx`;
- `feature_per_cell.csv`;
- `deduped_counts_histograms.txt`.

100K perturb parity:

```bash
ROOT=/tmp/star_suite_cbq_pf_direct_pf_100k_$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$ROOT"

PF_BIN=core/legacy/source/cbq_pf_adapter_harness
CBQ=/tmp/star_suite_cbq_pf_bench_perturb_20260530T165325Z/cbq/perturb100k.cbq
R1=/tmp/star_suite_cbq_pf_bench_perturb_20260530T165325Z/inputs/perturb100k/perturb100k_R1_001.fastq.gz
R2=/tmp/star_suite_cbq_pf_bench_perturb_20260530T165325Z/inputs/perturb100k/perturb100k_R2_001.fastq.gz
WHITELIST=/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt
FEATURE_REF=/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv

"$PF_BIN" --mode fastq \
  --barcodeFastq "$R1" \
  --featureFastq "$R2" \
  --whitelist "$WHITELIST" \
  --featureRef "$FEATURE_REF" \
  --outputDir "$ROOT/fastq" \
  --sampleName perturb100k \
  --barcodeLength 16 \
  --umiLength 12 \
  --featureOffset 31 \
  --consumerThreads 16 \
  --readBufferLines 65536

"$PF_BIN" --mode cbq-direct-pf \
  --readFilesIn "$CBQ" \
  --whitelist "$WHITELIST" \
  --featureRef "$FEATURE_REF" \
  --outputDir "$ROOT/cbq_direct_pf" \
  --sampleName perturb100k \
  --barcodeLength 16 \
  --umiLength 12 \
  --featureOffset 31 \
  --consumerThreads 16
```

Use canonical output comparison, not raw directory diff when logs/timestamps are
present. The helper for PF assign outputs is:

```bash
python3 tests/compare_pf_assign_outputs.py LEFT_SAMPLE_DIR RIGHT_SAMPLE_DIR
```

The existing full-lane parity checklist is:

- barcode set;
- feature list;
- matrix aggregated by `(feature, barcode)`;
- sorted `feature_per_cell.csv`;
- `deduped_counts_histograms.txt`;
- `stats.txt`.

### 6. Timing Gates

After synthetic and 100K parity pass, run full L001:

```bash
ROOT=/tmp/star_suite_cbq_pf_direct_pf_full_l001_$(date -u +%Y%m%dT%H%M%SZ)
CBQ=/tmp/star_suite_cbq_pf_bench_perturb_20260530T165325Z/cbq/perturb_l001_full.cbq

for c in 4 8 16 32; do
  /usr/bin/time -v "$PF_BIN" \
    --mode cbq-direct-pf \
    --readFilesIn "$CBQ" \
    --whitelist "$WHITELIST" \
    --featureRef "$FEATURE_REF" \
    --outputDir "$ROOT/cbq_direct_pf_c${c}" \
    --sampleName "cbq_direct_pf_c${c}" \
    --barcodeLength 16 \
    --umiLength 12 \
    --featureOffset 31 \
    --consumerThreads "$c" \
    > "$ROOT/cbq_direct_pf_c${c}.stdout" \
    2> "$ROOT/cbq_direct_pf_c${c}.time"
done
```

Record:

- CBQ path and size;
- record count;
- thread count;
- wall time;
- assignment/processing time from PF stats;
- max RSS;
- parity result;
- output root;
- storage tier (`/tmp`, SSD, or `/storage`).

Do not run full FASTQ and direct-PF benchmarks concurrently.

## Validation Results

Production `pf-multi` surface:

```text
--crAssignCbqMode auto
    Use indexed direct range counting when the CBQ has an index and the PF
    settings are safe for range processing. Fall back to queue-backed CBQ before
    direct output work starts.

--crAssignCbqMode stream
    Force the existing queue-backed CBQ stream.

--crAssignCbqMode range
    Require indexed direct range counting. Missing index or unsupported PF
    settings are fatal.
```

`auto` falls back only for preflight conditions such as missing CBQ index,
feature-mode bootstrap, or legacy CB rescue. Once direct range processing has
started, worker/runtime failures are reported rather than silently replayed
through the stream path.

Synthetic smoke:

```text
OUT_ROOT=/tmp/star_suite_cbq_pf_adapter_smoke_direct_pf_20260531T072131Z
tests/run_cbq_pf_adapter_smoke.sh
```

Result: PASS for FASTQ, queue-backed CBQ, direct-PF CBQ, direct indexed
count-only, and direct indexed sequence-materialization gates.

100K perturb validation:

- Output root:
  `/tmp/star_suite_cbq_pf_direct_pf_100k_20260531T072200Z`.
- Records: 100,000.
- Feature offset: `--featureOffset 31`.
- Consumers/workers: 16.
- Storage: `/tmp` SSD.

| Input | Wall Time | Assignment Time | Assignment Rate | Max RSS | Parity |
| --- | ---: | ---: | ---: | ---: | --- |
| FASTQ.gz | 1.88 s | 0.70 s | 143.5k reads/s | 241 MB | baseline |
| queue CBQ | 1.56 s | 0.68 s | 146.4k reads/s | 228 MB | PASS vs FASTQ |
| direct-PF CBQ | 1.47 s | 0.65 s | 155.0k reads/s | 228 MB | PASS vs FASTQ |

Canonical parity checked barcode set, feature list, matrix counts by feature
row and barcode, sorted `feature_per_cell.csv`, `deduped_counts_histograms.txt`,
and `stats.txt`. Raw row order differs under multithreaded PF and is not a
valid large-run parity signal.

Full L001 validation:

- Output root:
  `/tmp/star_suite_cbq_pf_direct_pf_full_l001_20260531T072623Z`.
- Input:
  `/tmp/star_suite_cbq_pf_bench_perturb_20260530T165325Z/cbq/perturb_l001_full.cbq`.
- CBQ size: 763,002,749 bytes.
- Records: 26,874,904.
- Feature offset: `--featureOffset 31`.
- Storage: `/tmp` SSD.

| Mode | Workers | Wall Time | Assignment Time | Assignment Rate | Max RSS | Parity vs queue c16 |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| queue CBQ | 16 | 3:19.94 | 189.83 s | 141.6k reads/s | 4.57 GB | baseline |
| direct-PF CBQ | 4 | 7:12.65 | 421.29 s | 63.8k reads/s | 3.26 GB | PASS |
| direct-PF CBQ | 8 | 4:20.96 | 248.73 s | 108.0k reads/s | 4.05 GB | PASS |
| direct-PF CBQ | 16 | 3:09.50 | 176.63 s | 152.1k reads/s | 4.43 GB | PASS |
| direct-PF CBQ | 32 | 2:21.30 | 128.35 s | 209.4k reads/s | 4.69 GB | PASS |

Direct-PF range counting improves over the queue-backed CBQ control once enough
workers are used. At 32 workers, full-lane wall time matches the previous
FASTQ/CBQ 32-consumer topline while avoiding the single CBQ stream handoff.

## Expected Result

The direct-reader gate proved the read/materialization side can process the
full L001 lane in under 1 second at high thread counts. Full PF direct range
counting will still spend time in feature matching, barcode correction, UMI
deduplication, hash merges, and output finalization, so do not expect a
sub-second end-to-end PF run.

The expected win is removing producer/queue pressure and letting PF consumers
stay fed directly from independent CBQ ranges. The benchmark should show:

- direct-PF CBQ output matches FASTQ and queue-backed CBQ;
- direct-PF CBQ is no slower than queue-backed CBQ;
- direct-PF scales more cleanly with `--consumerThreads` than the current
  single-reader CBQ path;
- RSS remains bounded by `threads * batch_scratch`, not by total reads.

## Rollback

If direct range counting regresses parity:

1. Keep `CbqInputModule::open_range()` and `--mode cbq-direct-decode`; these
   are validated and useful as diagnostics.
2. Disable or hide `--mode cbq-direct-pf`.
3. Keep production CBQ on the existing queue-backed `--mode cbq` path.
4. Compare the direct-PF worker's per-record line buffers against the
   queue-backed path for the first mismatching barcode/feature pair.

## Done Criteria

- `tests/run_cbq_pf_adapter_smoke.sh` passes FASTQ, queue-backed CBQ, and
  direct-PF CBQ parity: done.
- 100K perturb direct-PF parity passes: done.
- Full L001 direct-PF parity passes: done.
- STAR `pf-multi` exposes `--crAssignCbqMode auto|stream|range`: done.
- Direct mode supports ordered multi-lane CBQ planning: done.
- Timing results are added to
  `docs/RUNBOOK_CBQ_PARALLEL_PF_FLEX_OPTIMIZATIONS_20260530.md`: done.
- `git diff --check` passes: run before handoff.
