# Dynamic Thread Interface Tiny Fixture Runbook

Date: 2026-02-17  
Branch: `core-dynamic-threads`  
Status: reproducible smoke + mock consumer harness

## Purpose
Provide a fast, reproducible way to validate the Stage 1 dynamic-thread
interface (constant map permits + telemetry) and a mock consumer that reads the
interface output, including parity checks and retune-sequence assertions.

## Fixture Source
Use the built-in tiny Cell Ranger fixture:

1. Reference: `/home/lhhung/cellranger-9.0.1/external/cellranger_tiny_ref`
2. Prebuilt STAR index (for inventory only):
   `/home/lhhung/cellranger-9.0.1/external/cellranger_tiny_ref/star`
3. FASTQs: `/home/lhhung/cellranger-9.0.1/external/cellranger_tiny_fastq`
4. Read inputs used by this smoke:
   `tinygex_S1_L001_R2_001.fastq.gz,tinygex_S1_L002_R2_001.fastq.gz`

Quick fixture sanity:

```bash
head -n 1 /home/lhhung/cellranger-9.0.1/external/cellranger_tiny_ref/fasta/genome.fa
# expected: >chr21.part
```

## Harness Components
1. Integration smoke script:
   `tests/run_dynamic_threads_tiny_fixture.sh`
2. Mock consumer report generator:
   `tests/dynamic_threads/mock_consumer_report.py`
3. Mock consumer unit test:
   `tests/dynamic_threads/test_mock_consumer_report.py`

## Reproduce
1. Build STAR:

```bash
make core
```

2. Run mock consumer unit test:

```bash
python3 tests/dynamic_threads/test_mock_consumer_report.py
```

3. Run tiny fixture dynamic-thread smoke:

```bash
tests/run_dynamic_threads_tiny_fixture.sh
```

Optional overrides:

```bash
RUN_THREADS=6 MAP_PERMITS=3 READ_MAP_NUMBER=200000 \
VARIABLE_THREADS=1 \
VARIABLE_THREADS_RETUNE_EVERY_ACQUIRES=1 \
VARIABLE_THREADS_PERMIT_SEQUENCE="2 4" \
MIN_RETUNES=1 \
OUT_SAMTYPE="BAM SortedByCoordinate" \
CHECK_BAM_PARITY=1 \
OUT_BASE=/tmp/dynamic_threads_tiny_manual \
tests/run_dynamic_threads_tiny_fixture.sh
```

Requested sequence smoke (both scenarios):

```bash
tests/run_dynamic_threads_variable_sequences_smoke.sh
```

## What The Smoke Does
1. Runs baseline mode (`--dynamicThreadInterface 0`) on the tiny fixture.
2. Builds a fresh tiny STAR index for the current STAR binary from fixture
   FASTA/GTF (avoids version mismatch with old prebuilt indexes).
3. Runs dynamic mode with telemetry:
   `--dynamicThreadInterface 1 --dynamicThreadConstMapPermits <N> --dynamicThreadTelemetry 1`.
4. Parses `on/Log.out` with the mock consumer parser.
5. Enforces telemetry exists with `acquires >= 1` and `work_units >= 1`.
6. Enforces `Log.final.out` parity by canonical metric diff (default on).
7. Confirms baseline and dynamic runs processed the same
   `Number of input reads`.
8. When variable retuning is enabled, enforces exact cyclic retune trace.
9. Optional: when `OUT_SAMTYPE` is BAM and `CHECK_BAM_PARITY=1`, validates BAM
   integrity and parity (`samtools quickcheck`, record counts, body MD5).

When enabled, variable retune controls are:

1. `--variableThreadsRetuneEveryAcquires <N>`
2. `--variableThreadsPermitSequence <p1> <p2> ...`
3. Permit targets are clamped to `1..runThreadN` (fixed worker pool bound).

## Output Artifacts
Default output root:
`/tmp/dynamic_threads_tiny_<YYYYMMDD_HHMMSS>_<pid>/`

Produced files:
1. Baseline logs: `off/Log.out`, `off/Log.final.out`
2. Dynamic logs: `on/Log.out`, `on/Log.final.out`
3. Mock consumer JSON report: `dynamic_thread_report.json`
4. Mock consumer text summary: `dynamic_thread_report.txt`
5. Canonical `Log.final.out` parity diff (expected empty): `log_final.diff`
6. Optional BAM parity summary: `bam_parity_summary.txt`

## Interface Contract Captured By Mock Consumer
The current mock consumer parses and validates these interface outputs:

1. Configuration line:
   `Dynamic thread interface enabled: map permits=... (runThreadN=..., telemetry=on|off, variableThreads=..., retuneEveryAcquires=..., retuneSequenceLength=...)`
2. Telemetry line:
   `Dynamic thread telemetry: acquires=..., retunes=..., retuneEveryAcquires=..., retuneSequenceLength=..., targetPermits=..., configuredPermits=..., retuneTrace=..., retuneTraceDropped=..., workUnits=..., workBytes=..., waitMs(...), workMs(...)`

This keeps the consumer-facing contract stable while STAR internals evolve
toward later producer/consumer mode.
