# Runbook: CBQ Parallel PF And FLEX Optimizations

Date: 2026-05-30

Status: PF easy wins, direct range-counting, and the production `pf-multi`
CBQ mode gate are implemented in the STAR-suite worktree. Start with
`process_features` because the producer/consumer structure and CBQ parity
harness are already in place.

Related docs:

- `docs/RUNBOOK_PROCESS_FEATURES_CBQ_NATIVE.md`
- `docs/RUNBOOK_CBQ_ORDERED_ENCODER.md`
- `docs/RUNBOOK_FLEX_CBQ_INPUT.md`
- `docs/RUNBOOK_FLEX_NO_GENOME_COUNT_ONLY_20260529.md`
- `docs/RUNBOOK_CHROMAP_ATAC_CBQ_IN_MEMORY.md`

## Goal

Make CBQ a production input path that is faster than gzipped FASTQ and still
easy to validate:

```text
CBQ block/index reader -> bounded producer queue -> existing consumers
```

The first target is `process_features` (`pf`) because it already has:

- a FASTQ producer and consumer ring buffer in
  `core/features/process_features/src/assignBarcodes.c`;
- a non-FASTQ streaming API in
  `core/features/process_features/include/pf_api.h`;
- a CBQ adapter harness in
  `core/legacy/source/input/cbq_pf_adapter_harness.cpp`;
- synthetic FASTQ-vs-CBQ parity coverage in
  `tests/run_cbq_pf_adapter_smoke.sh`;
- aggregate CBQ module regression coverage in
  `tests/run_cbq_e2e_module_regression.sh`.

After PF is stable, apply the same block/range producer shape to FLEX/STAR and
to Chromap-suite/libchromap.

## Non-Goals

- Do not add a second CBQ format parser in `process_features`.
- Do not materialize production CBQ input as temporary FASTQ.
- Do not change barcode, UMI, feature-search, Flex, or alignment semantics.
- Do not use benchmark output directories from previous runs.

## Current PF Surfaces

FASTQ path:

- `fastq_reader_set` in `core/features/process_features/include/common.h`
  owns the bounded buffer and condition variables.
- `read_fastqs_by_set()` in
  `core/features/process_features/src/assignBarcodes.c` reads FASTQ/gzip with
  `gzgets()` into that buffer.
- `consume_reads()` in the same file consumes the buffer and runs the existing
  feature assignment/counting path.

CBQ path:

- `CbqInputModule` is the shared reader used by STAR-suite input harnesses.
- `pf_process_records_begin()`, `pf_process_record_views()`, and
  `pf_process_records_end()` are the streaming PF adapter API.
- `cbq_pf_adapter_harness.cpp` currently reads a paired CBQ, materializes
  sequence spans into bounded scratch storage, and passes borrowed record views
  into PF.

The current CBQ adapter is correct but conservative. It is already an in-memory
path, but the production optimization is to make it a first-class PF producer
and then exploit CBQ block/index range reads.

## Implementation Order

### Phase 0: Baseline And Guardrails

Establish a clean baseline before editing shared PF or CBQ code.

```bash
make -C core/legacy/source clean
make -C core/legacy/source -j8 cbq-pf-adapter-harness cbq-ordered-encoder

OUT_ROOT=/tmp/star_suite_cbq_pf_adapter_smoke_baseline_$(date -u +%Y%m%dT%H%M%SZ) \
  tests/run_cbq_pf_adapter_smoke.sh

OUT_ROOT=/tmp/star_suite_cbq_e2e_module_regression_baseline_$(date -u +%Y%m%dT%H%M%SZ) \
RUN_CHROMAP_MAPPING_SMOKE=0 \
  tests/run_cbq_e2e_module_regression.sh
```

Record:

- git commit and branch;
- exact command;
- output root;
- wall time and max RSS when timed;
- whether the run used `/tmp`, SSD, or `/storage`.

Use wrapper-written summaries when available. For these smoke scripts, a
successful exit and the final `PASS:` line are the completion signal.

### Phase 1: PF Production CBQ Gate

Goal: promote the existing harness path into a production PF input without
changing PF counting behavior.

Status: STAR `pf-multi` accepts explicit CBQ paths/directories through
`PfMultiAssign::runAssignBarcodes()` by resolving `.cbq` sources. It can feed
`CbqInputModule` through the queue-backed stream path or through the indexed
direct range-counting path. The implementation records `input_format`,
`cbq_mode_requested`, `cbq_mode_effective`, and
`cbq_mode_fallback_reason` in `assignBarcodes.api_run.txt`.

Tasks:

- Production input is explicit CBQ source resolution. Keep the first surface
  narrow:

```text
pf-multi feature-library fastqs path -> paired.cbq
pf-multi feature-library fastqs path -> cbq_directory/
```

- For multi-lane samples, accept an explicit ordered list or manifest. Do not
  infer sample grouping from CBQ filenames until that convention is documented.
- Reuse `CbqInputModule` and `pf_process_record_views()`.
- Preserve mate roles from the existing harness:
  - mate 0: barcode + UMI read;
  - mate 1: feature/protospacer read;
  - optional mate 2: future reverse/alternate feature read.
- Keep default-quality synthesis and missing-quality handling identical to the
  current streaming API.
- Leave feature matching, barcode correction, UMI deduplication, and
  EmptyDrops behavior untouched.

Acceptance:

- `tests/run_cbq_pf_adapter_smoke.sh` passes with byte-identical `barcodes.txt`,
  `features.txt`, `feature_sequences.txt`, `matrix.mtx`,
  `feature_per_cell.csv`, and `deduped_counts_histograms.txt`.
- The new production CBQ case is registered in
  `tests/run_cbq_e2e_module_regression.sh`.
- Production CBQ processing creates no temporary FASTQ files.

Mode control:

```text
--crAssignCbqMode auto
    Use indexed direct range counting when supported, otherwise fall back to
    queue-backed CBQ before direct output work starts.

--crAssignCbqMode stream
    Force queue-backed CBQ streaming.

--crAssignCbqMode range
    Require indexed direct range counting and fail if unsupported.
```

### Phase 2: PF Easy Wins

Goal: reduce obvious read-side overhead before adding block-parallel reads.

Status: implemented for the first pass:

- PF FASTQ readers call `gzbuffer()` with a 1 MiB zlib buffer.
- PF API exposes `pf_config_set_read_buffer_lines()`.
- STAR pf-multi exposes this as `--crAssignReadBufferLines`.
- The CBQ PF harness exposes this as `--readBufferLines`.
- STAR pf-multi and the CBQ PF harness reuse CBQ record-view vectors and
  sequence scratch across batches.
- Borrowed quality spans are kept from CBQ batches.

Tasks:

- Reuse the `pf_read_record_view` vector and sequence scratch across CBQ
  batches instead of reallocating them every batch.
- Keep qualities as borrowed spans from the CBQ batch wherever possible.
- Minimize NUL-terminated sequence materialization to the fields that still
  need it inside PF.
- Add `gzbuffer()` to the FASTQ baseline reader so FASTQ-vs-CBQ benchmarks are
  fair and the current FASTQ path is not penalized by tiny zlib buffers.
- Keep `read_buffer_lines` configurable for PF timing runs; record the value in
  timing summaries.
- STAR pf-multi exposes this as `--crAssignReadBufferLines`; the standalone
  CBQ PF harness exposes `--readBufferLines`.

Acceptance:

- Synthetic FASTQ-vs-CBQ parity remains byte-identical.
- The FASTQ baseline does not regress on the synthetic smoke.
- `git diff --check` passes.

### Phase 3: PF CBQ Producer

Goal: replace the FASTQ/gzip producer with a CBQ producer feeding the existing
PF consumer queue.

Initial shape:

```text
single CBQ reader thread -> decoded PF record views -> existing PF consumers
```

Then:

```text
CBQ block/range producers -> ordered merge or record-range queue -> PF consumers
```

Tasks:

- Add a PF-local `CbqPfProducer` that owns `CbqInputModule` and emits batches
  into the same bounded queue concept used by `fastq_reader_set`.
- Keep the queue bounded. Do not allow a full sample to accumulate in memory.
- Preserve input record order for PF outputs that depend on stable streaming
  behavior.
- Add clear errors for:
  - single-end CBQ passed to paired PF mode;
  - fewer records in one mate stream than another;
  - missing mate role;
  - quality length mismatch;
  - sequence length over PF's current line capacity.
- Make the old harness a thin test wrapper over the production producer.

Acceptance:

- Synthetic smoke passes through the production producer.
- Downsampled PF fixture passes FASTQ-vs-CBQ count parity.
- Peak memory scales with queue depth and batch size, not total reads.

### Phase 4: CBQ Indexed Range Reader

Goal: exploit CBQ's parallel-friendly file layout.

Design rules:

- Partition work by logical record ranges, not by physical block boundaries.
- For split streams, read/barcode CBQ files may have different block
  boundaries. A range planner must seek each file to the blocks covering the
  same logical record interval.
- Each range producer emits records tagged with range ordinal and local record
  ordinal.
- PF can consume ranges directly if final output is order-independent. If a
  future path becomes order-sensitive, reorder by `(range_ordinal,
  local_record_ordinal)`.

Suggested shared API:

```text
CbqBlockPlan open_cbq_plan(paths)
CbqRangeReader open_range(plan, first_record, record_count)
next_batch(range_reader) -> CbqReadBatchView
```

Acceptance:

- Range 0..N output is byte-identical to sequential CBQ output.
- Parallel range output is byte-identical to sequential CBQ output after any
  required order merge.
- Default-compressed and level-0 CBQ both pass.

## Small PF Validation Artifacts

Use the synthetic fixture first:

- Input fixture:
  `core/features/process_features/tests/fixtures/assignbarcodes_baseline/input`
- Smoke output root:
  `/tmp/star_suite_cbq_pf_adapter_smoke*`
- Aggregate output root:
  `/tmp/star_suite_cbq_e2e_module_regression_*`

For host-local downsampled validation, preflight whichever fixture is present:

```bash
test -d /storage/ucsf-2M/fixtures/ucsf2m_iPSC2_AALG2_100k_pfconfig && \
  find /storage/ucsf-2M/fixtures/ucsf2m_iPSC2_AALG2_100k_pfconfig -maxdepth 3 -type f | head

test -d /mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2 && \
  find /mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2 -maxdepth 2 -type f,l | head
```

The preferred 100K PF fixture, when present, is:

```text
/storage/ucsf-2M/fixtures/ucsf2m_iPSC2_AALG2_100k_pfconfig
```

The corrected UCSF symlink fixture currently present on this host includes
`/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2/pf_multi_config.csv` and can be
used as a staging source if the 100K fixture has to be regenerated.

For a lane-level PF CBQ validation before multi-lane production support exists:

```bash
OUT=/tmp/star_suite_cbq_pf_100k_lane_$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$OUT"

BQTOOLS=${BQTOOLS:-/tmp/star_suite_bqtools/bin/bqtools}
PF_BIN=core/legacy/source/cbq_pf_adapter_harness
WHITELIST=/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt
FEATURE_REF=/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv

R1=/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2/guides/EBs2_2_AALG2_S26_L001_R1_001.fastq.gz
R2=/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2/guides/EBs2_2_AALG2_S26_L001_R2_001.fastq.gz

"$BQTOOLS" encode "$R1" "$R2" --mode cbq -o "$OUT/lane.cbq" -T 4

"$PF_BIN" --mode fastq \
  --barcodeFastq "$R1" \
  --featureFastq "$R2" \
  --whitelist "$WHITELIST" \
  --featureRef "$FEATURE_REF" \
  --outputDir "$OUT/fastq" \
  --sampleName EBs2_2_AALG2_L001 \
  --barcodeLength 16 \
  --umiLength 12

"$PF_BIN" --mode cbq \
  --readFilesIn "$OUT/lane.cbq" \
  --whitelist "$WHITELIST" \
  --featureRef "$FEATURE_REF" \
  --outputDir "$OUT/cbq" \
  --sampleName EBs2_2_AALG2_L001 \
  --barcodeLength 16 \
  --umiLength 12

diff -ru "$OUT/fastq/EBs2_2_AALG2_L001" "$OUT/cbq/EBs2_2_AALG2_L001"
```

For multi-lane PF validation, add or use the production CBQ manifest/list mode
first; do not concatenate unrelated lane records with ad hoc shell pipelines.

## Timing Protocol

Use fresh output directories and serialize benchmark jobs.

For PF:

- run synthetic smoke first, mostly for correctness;
- run one downsampled lane for early timing;
- run the full downsampled PF sample after multi-lane CBQ input is wired;
- record `read_buffer_lines`, CBQ compression level, thread counts, input
  device, output device, wall time, and max RSS.

Example wrapper:

```bash
/usr/bin/time -v bash -lc 'OUT_ROOT=/tmp/star_suite_cbq_pf_adapter_smoke_timed_$(date -u +%Y%m%dT%H%M%SZ) tests/run_cbq_pf_adapter_smoke.sh'
```

For apples-to-apples production comparisons, keep FASTQ and CBQ inputs on the
same storage tier. Use `/storage` for cache-resistant timing and SSD only when
the production comparison explicitly requires SSD.

## FLEX/STAR Follow-On

The FLEX bottleneck is FASTQ ingestion. Once PF has a reusable range reader,
apply the same shape to STAR/FLEX:

```text
CBQ range producers -> STAR input chunks -> ReadAlignChunk/Flex consumers
```

Rules:

- Preserve chunk IDs and record order for order-dependent outputs, especially
  unsorted BAM, Y/noY BAM/SAM splitting, and any CB/UB tag emission that uses
  input order.
- Coordinate-sorted BAM can tolerate more producer freedom because sorting
  becomes the order boundary.
- Keep no-genome/count-only FLEX on the same input abstraction as mapped FLEX.
- Validate against the topline FLEX count-only results in
  `docs/RUNBOOK_FLEX_NO_GENOME_COUNT_ONLY_20260529.md`.

Minimum gates:

```bash
OUT_ROOT=/tmp/star_suite_cbq_star_input_smoke_range_$(date -u +%Y%m%dT%H%M%SZ) \
  tests/run_cbq_star_input_smoke.sh

OUT_ROOT=/tmp/star_suite_cbq_ynoy_smoke_range_$(date -u +%Y%m%dT%H%M%SZ) \
  tests/run_cbq_ynoy_smoke.sh

OUT_ROOT=/tmp/star_suite_cbq_e2e_module_regression_range_$(date -u +%Y%m%dT%H%M%SZ) \
RUN_CHROMAP_MAPPING_SMOKE=0 \
  tests/run_cbq_e2e_module_regression.sh
```

For FLEX 100K and production timing, use the existing host-local artifacts
recorded in `tests/ARTIFACTS.md`:

- `/tmp/star_suite_cbq_flex_100k_*`
- `/storage/downsampled_100K/SC2300771`

## Chromap-suite/libchromap Follow-On

Chromap-suite already has a CBQ batch-producer direction in flight. The next
optimization is the same range-planned read:

```text
CBQ paired-read range + CBQ barcode range -> libchromap batch provider
```

Rules:

- For ATAC multiome, partition by logical record range across read-pair and
  barcode CBQ files.
- Do not assume the read-pair CBQ and barcode CBQ share physical block
  boundaries.
- Keep Y/noY output order deterministic.
- Keep the STAR-suite libchromap contract smoke as the cross-repo gate.

Minimum gates:

```bash
OUT_ROOT=/tmp/star_suite_libchromap_cbq_contract_smoke_range_$(date -u +%Y%m%dT%H%M%SZ) \
  tests/run_star_libchromap_cbq_contract_smoke.sh
```

In Chromap-suite, continue to run:

```bash
OUT_ROOT=/tmp/chromap_cbq_atac_smoke_range_$(date -u +%Y%m%dT%H%M%SZ) \
  tests/run_cbq_atac_smoke.sh

OUT_ROOT=/tmp/chromap_input_format_smoke_range_$(date -u +%Y%m%dT%H%M%SZ) \
  tests/run_input_format_smoke.sh
```

## Benchmark Table Template

Fill this table as each phase lands.

| Phase | Input | Dataset | Threads | Storage | Wall time | Max RSS | Parity gate | Output root |
| --- | --- | --- | ---: | --- | ---: | ---: | --- | --- |
| PF baseline | FASTQ.gz | synthetic | 1 | `/tmp` | TBD | TBD | PASS | TBD |
| PF baseline | CBQ | synthetic | 1 | `/tmp` | TBD | TBD | PASS | TBD |
| PF production | FASTQ.gz | 100K lane | TBD | TBD | TBD | TBD | TBD | TBD |
| PF production | CBQ | 100K lane | TBD | TBD | TBD | TBD | TBD | TBD |
| FLEX no-genome | FASTQ.gz | full | TBD | SSD | see Flex runbook | see Flex runbook | PASS | see Flex runbook |
| FLEX no-genome | level-0 CBQ | full | TBD | SSD | see Flex runbook | see Flex runbook | PASS | see Flex runbook |

## 2026-05-30 PF Perturb Benchmark Notes

Benchmark target:

- Dataset: UCSF `EBs2_2` AALG2 guide lane L001.
- Records: 26,874,904 read pairs.
- Feature reference:
  `/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv`.
- Whitelist:
  `/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt`.
- Harness: `core/legacy/source/cbq_pf_adapter_harness`.
- Correct feature offset for this benchmark: `--featureOffset 31`.

Diagnostic finding:

- A first run allowed the feature-ref pattern to auto-detect offset 0 from
  `(BC)GTTTNAGAGCTAAGC`.
- The observed read matches are around positions 30-32, so offset 0 made every
  constant-offset hash probe miss and forced full feature search.
- That wrong-offset run processed only about 26 thousand reads/s at 8 consumers
  and is not a valid CBQ benchmark.
- The tiered-search debug trace is now gated by `PF_TRACE_TIERED`.

100K thread sweep with `--featureOffset 31`, `--readBufferLines 65536`:

| Input | Consumer Threads | Wall Time | Assignment Rate |
| --- | ---: | ---: | ---: |
| FASTQ.gz | 4 | 2.46 s | 59.8k reads/s |
| CBQ | 4 | 2.54 s | 57.0k reads/s |
| FASTQ.gz | 8 | 1.72 s | 107.9k reads/s |
| CBQ | 8 | 1.73 s | 107.1k reads/s |
| FASTQ.gz | 16 | 1.43 s | 151.8k reads/s |
| CBQ | 16 | 1.42 s | 157.3k reads/s |
| FASTQ.gz | 32 | 1.24 s | 216.3k reads/s |
| CBQ | 32 | 1.24 s | 222.0k reads/s |

Thread sweep output root:

```text
/tmp/star_suite_cbq_pf_thread_sweep_100k_20260530T171730Z
```

Full L001 lane with `--featureOffset 31`, `--consumerThreads 32`,
`--readBufferLines 65536`, and both FASTQ.gz copies plus CBQ on `/tmp`:

| Input | Records | Wall Time | Assignment Time | Assignment Rate | Max RSS | Canonical Parity |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| FASTQ.gz | 26,874,904 | 135.70 s | 127.95 s | 210.0k reads/s | 4.0 GB | PASS |
| CBQ | 26,874,904 | 138.86 s | 128.61 s | 209.0k reads/s | 4.8 GB | PASS |

Full-lane output root:

```text
/tmp/star_suite_cbq_pf_bench_perturb_full_offset31_c32_20260530T171837Z
```

Canonical parity checks:

- barcode set;
- feature list;
- matrix aggregated by `(feature, barcode)`;
- sorted `feature_per_cell.csv`;
- `deduped_counts_histograms.txt`;
- `stats.txt`.

Interpretation:

- Once the feature offset is fixed and consumers use all 32 CPUs, feature-only PF
  is compute/finalization dominated on this guide library rather than FASTQ
  ingestion dominated.
- The current CBQ harness path is essentially tied with gzipped FASTQ for this
  full-lane perturb case. This is expected until the PF CBQ producer/range
  reader phase exploits CBQ's parallel block layout.

## 2026-05-31 PF Indexed Range Prototype

Implemented the first optimization gate for PF:

- `CbqInputModule::open_range(lane, first_record, record_count)` parses the
  existing tail `CBQINDEX`, seeks to the indexed block covering a logical record
  range, and returns trimmed `CbqReadBatchView` slices.
- `cbq_pf_adapter_harness --mode cbq-direct-decode` lets each consumer thread
  own an independent CBQ reader/range.
- `--materializeMode none` measures indexed decode/count only.
- `--materializeMode sequence` also materializes barcode/feature sequences into
  the same bounded scratch shape used by the PF adapter.

This is a tuning gate, not the final PF counting integration. Full PF parity is
still checked through the existing FASTQ-vs-CBQ adapter path.

Validation:

```bash
make -C core/legacy/source -j8 cbq-pf-adapter-harness

OUT_ROOT=/tmp/star_suite_cbq_pf_adapter_smoke_direct_threadsafe_20260531T063247Z \
  tests/run_cbq_pf_adapter_smoke.sh
```

Smoke result:

- FASTQ-vs-CBQ PF output parity: PASS.
- Direct indexed count gate: PASS (`processed_reads=7`, `indexed_reads=7`).
- Direct indexed sequence-materialization gate: PASS
  (`processed_reads=7`, `indexed_reads=7`, `lane_ordinal_sum=28`).

100K direct-reader sweep:

- Input:
  `/tmp/star_suite_cbq_pf_bench_perturb_20260530T165325Z/cbq/perturb100k.cbq`.
- Output root:
  `/tmp/star_suite_cbq_pf_direct_range_100k_20260531T063312Z`.
- Storage: `/tmp` cached/local.

| Materialize Mode | Reader Threads | Records | Batches | Bases Materialized | Ordinal Sum | Harness Elapsed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| none | 1 | 100,000 | 220 | 0 | 0 | 0.0339 s |
| none | 4 | 100,000 | 223 | 0 | 0 | 0.0131 s |
| none | 8 | 100,000 | 227 | 0 | 0 | 0.0071 s |
| none | 16 | 100,000 | 235 | 0 | 0 | 0.0053 s |
| none | 32 | 100,000 | 251 | 0 | 0 | 0.0072 s |
| sequence | 1 | 100,000 | 220 | 11,800,000 | 5,000,050,000 | 0.0334 s |
| sequence | 4 | 100,000 | 223 | 11,800,000 | 5,000,050,000 | 0.0109 s |
| sequence | 8 | 100,000 | 227 | 11,800,000 | 5,000,050,000 | 0.0069 s |
| sequence | 16 | 100,000 | 235 | 11,800,000 | 5,000,050,000 | 0.0046 s |
| sequence | 32 | 100,000 | 251 | 11,800,000 | 5,000,050,000 | 0.0073 s |

Interpretation:

- The indexed range reader now exposes the optimization path the full PF
  integration needs: direct consumers can read independent logical record ranges
  from one CBQ file.
- On this small cached 100K input, useful scaling appears by 4-8 reader threads
  and then flattens; the next meaningful gate is to attach these range readers
  to PF's per-thread counting state rather than feeding `pf_process_record_views`
  through the existing stream queue.

Full L001 direct-reader sweep:

- Input:
  `/tmp/star_suite_cbq_pf_bench_perturb_20260530T165325Z/cbq/perturb_l001_full.cbq`.
- CBQ size: 763,002,749 bytes.
- Records: 26,874,904.
- Output root:
  `/tmp/star_suite_cbq_pf_direct_range_full_l001_20260531T063530Z`.
- Storage: `/tmp` cached/local.

| Materialize Mode | Reader Threads | Records | Bases Materialized | Ordinal Sum | Harness Elapsed | Max RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| none | 1 | 26,874,904 | 0 | 0 | 5.7229 s | 7.1 MB |
| none | 2 | 26,874,904 | 0 | 0 | 2.7388 s | 9.4 MB |
| none | 4 | 26,874,904 | 0 | 0 | 1.3916 s | 13.5 MB |
| none | 8 | 26,874,904 | 0 | 0 | 0.8598 s | 21.7 MB |
| none | 16 | 26,874,904 | 0 | 0 | 0.5235 s | 38.0 MB |
| none | 32 | 26,874,904 | 0 | 0 | 0.3580 s | 70.3 MB |
| sequence | 1 | 26,874,904 | 3,171,238,672 | 361,130,245,942,060 | 6.5649 s | 7.9 MB |
| sequence | 2 | 26,874,904 | 3,171,238,672 | 361,130,245,942,060 | 3.4942 s | 11.7 MB |
| sequence | 4 | 26,874,904 | 3,171,238,672 | 361,130,245,942,060 | 1.7554 s | 17.3 MB |
| sequence | 8 | 26,874,904 | 3,171,238,672 | 361,130,245,942,060 | 0.9250 s | 30.7 MB |
| sequence | 16 | 26,874,904 | 3,171,238,672 | 361,130,245,942,060 | 0.6391 s | 59.3 MB |
| sequence | 32 | 26,874,904 | 3,171,238,672 | 361,130,245,942,060 | 0.4051 s | 114.3 MB |

Interpretation:

- On the full lane, the direct range reader scales nearly linearly through
  4 threads and continues improving through 32 threads in this cached local
  test.
- The constant ordinal checksum in sequence mode confirms that logical record
  ranges cover the lane exactly once despite extra boundary blocks being
  decoded internally.
- The optimization target is now clearly downstream of the reader: attach these
  range readers to PF's per-thread counting state and keep PF's output parity
  gate as the acceptance test.

## 2026-05-31 PF Direct Range Counting Results

Implemented the direct PF range-counting gate from
`docs/RUNBOOK_PF_CBQ_RANGE_COUNTING_20260531.md`:

- `pf_direct_range_begin/process_record_views/end/abort` gives externally
  partitioned CBQ readers one worker-local PF counting state each.
- `cbq_pf_adapter_harness --mode cbq-direct-pf` partitions one indexed paired
  CBQ file by logical record range and runs one `CbqInputModule` per worker.
- STAR `pf-multi` exposes the same path through `--crAssignCbqMode range`, with
  `auto` fallback to queue-backed CBQ for missing indexes or unsafe settings.
- Ordered multi-lane CBQ inputs are partitioned by logical record order before
  being assigned to direct workers.
- `tests/compare_pf_assign_outputs.py` compares large multithreaded PF outputs
  canonically, independent of barcode row ordering.

Validation commands:

```bash
make -C core/legacy/source -j8 cbq-pf-adapter-harness

OUT_ROOT=/tmp/star_suite_cbq_pf_adapter_smoke_direct_pf_20260531T072131Z \
  tests/run_cbq_pf_adapter_smoke.sh
```

Synthetic smoke result: PASS for FASTQ, queue-backed CBQ, direct-PF CBQ, direct
indexed count-only, and direct indexed sequence-materialization gates.

100K perturb results:

- Output root:
  `/tmp/star_suite_cbq_pf_direct_pf_100k_20260531T072200Z`.
- Records: 100,000.
- Feature offset: `--featureOffset 31`.
- Consumers/workers: 16.
- Storage: `/tmp` SSD.

| Input | Wall Time | Assignment Time | Assignment Rate | Max RSS | Canonical Parity |
| --- | ---: | ---: | ---: | ---: | --- |
| FASTQ.gz | 1.88 s | 0.70 s | 143.5k reads/s | 241 MB | baseline |
| queue CBQ | 1.56 s | 0.68 s | 146.4k reads/s | 228 MB | PASS vs FASTQ |
| direct-PF CBQ | 1.47 s | 0.65 s | 155.0k reads/s | 228 MB | PASS vs FASTQ |

Full L001 results:

- Output root:
  `/tmp/star_suite_cbq_pf_direct_pf_full_l001_20260531T072623Z`.
- Input:
  `/tmp/star_suite_cbq_pf_bench_perturb_20260530T165325Z/cbq/perturb_l001_full.cbq`.
- CBQ size: 763,002,749 bytes.
- Records: 26,874,904.
- Feature offset: `--featureOffset 31`.
- Storage: `/tmp` SSD.

| Mode | Workers | Wall Time | Assignment Time | Assignment Rate | Max RSS | Canonical Parity |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| queue CBQ | 16 | 3:19.94 | 189.83 s | 141.6k reads/s | 4.57 GB | baseline |
| direct-PF CBQ | 4 | 7:12.65 | 421.29 s | 63.8k reads/s | 3.26 GB | PASS vs queue c16 |
| direct-PF CBQ | 8 | 4:20.96 | 248.73 s | 108.0k reads/s | 4.05 GB | PASS vs queue c16 |
| direct-PF CBQ | 16 | 3:09.50 | 176.63 s | 152.1k reads/s | 4.43 GB | PASS vs queue c16 |
| direct-PF CBQ | 32 | 2:21.30 | 128.35 s | 209.4k reads/s | 4.69 GB | PASS vs queue c16 |

Interpretation:

- Direct range counting removes the single CBQ stream handoff and wins once
  enough workers are used.
- c32 is the useful full-lane point on this host: it matches the earlier
  FASTQ/queue-CBQ c32 topline while keeping CBQ input native.
- c4 and c8 are dominated by PF feature/counting work, not CBQ reading.
- Raw `barcodes.txt` and `feature_sequences.txt` row order can differ across
  multithreaded runs; use the canonical comparator for large-run parity.

Production guardrails:

- `legacy_cb_rescue=1` uses the stream fallback in `auto` and is fatal in
  forced `range` mode.
- Feature-mode bootstrap uses the stream fallback unless an explicit feature
  offset disables bootstrap. For the perturb benchmark, use
  `--crAssignFeatureOffset 31`.
- Direct range currently supports paired PF input (`nreaders=2`).
- Missing CBQ index uses stream fallback in `auto` and is fatal in forced
  `range` mode.

Remaining optimization follow-ons:

- Remove one more copy layer by teaching the PF direct consumer to accept
  borrowed spans directly instead of NUL-terminated line buffers.
- Reconcile the stream consumer and direct consumer into one shared record core
  after another parity gate, so future PF behavior changes cannot drift.
- Run the same direct-range timing on `/storage` when a cache-resistant
  production benchmark is needed.

## Done Criteria

PF is done when:

- production CBQ input is available without temporary FASTQ;
- synthetic smoke and aggregate CBQ regression pass;
- downsampled PF FASTQ-vs-CBQ count parity passes;
- timings are recorded in this runbook or a dated handoff;
- new test output locations are added to `tests/ARTIFACTS.md` if new scripts or
  stable artifact roots are introduced.

FLEX/Chromap follow-on is done when:

- range-planned CBQ readers are shared or have a clearly documented common
  contract;
- order-dependent BAM/Y/noY outputs pass smoke and parity checks;
- production timing is captured with FASTQ and CBQ on the same storage tier.
