# Handoff: Feature Assignment Consumer Optimizations (2026-03-17)

## Summary

Five commits optimizing the `assignBarcodes` consumer hot path in
`core/features/process_features/`. Total diff: 2 files, +99/-47 lines.
All changes are in the producer-consumer pipeline that connects FASTQ reader
threads (producers) to barcode/feature matching threads (consumers).

Smoke-tested after each commit against the corrected UCSF `EBs2_2` dataset
(100K GEX reads, 10 guide lanes, 548 CRISPR features). Every smoke produced
identical results: **708,057 feature counts, 286,905 unmatched reads**.

## Commits

```
4ef7359 pf: fix consumer spin, probe-exit deadlock, and trylock error handling
632eff9 pf: fix head-of-line blocking in consumer with non-blocking round-robin
2481b57 pf: batch permit acquire/release per 64 reads
6e17064 pf: eliminate redundant strlen calls in consumer hot path
919856a pf: increase ring buffer from 1024 to 8192 lines
```

Base: `6452600 pf: fix autodetect regressions and add cpu-aware controller`

## Changes in Detail

### 1. Increase ring buffer: 1024 → 8192 lines (`common.h`)

**File**: `core/features/process_features/include/common.h`

One-line change: `#define READ_BUFFER_LINES 1024` → `8192`.

With `lines_per_block = 4` (2 readers × 2 lines each), the old 1024-line
buffer held ~256 reads. With 31 consumer threads competing for data, the
buffer emptied after ~8 reads per consumer, forcing producers to stall
frequently on `can_produce` waits. The 8192-line buffer holds ~2048 reads,
giving producers 8× more runway before blocking.

Smoke: 9.96s → baseline was 10.68s (**~7% faster**).

### 2. Eliminate redundant `strlen` calls (`assignBarcodes.c`)

**Before**: After `strcpy`-ing 4–6 lines from the ring buffer, the consumer
called `strlen()` on each copied string to compute `work_bytes` for the
permit system. This re-scanned the same bytes that were just copied.

**After**: `strcpy` → `stpcpy`, which returns a pointer to the destination's
null terminator. String length is computed as `end - dst` with zero
additional scanning.

Eliminates 4–6 `strlen` calls per read (at 278M reads, that's ~1.4 billion
redundant byte scans removed).

### 3. Batch permit acquire/release per 64 reads (`assignBarcodes.c`)

**Before**: Every single read did:
1. `permit_acquire_hook()` — atomic acquire + potential condvar wait
2. (process read)
3. `permit_release_hook()` — atomic release + accounting

At 278M reads, that's 556M atomic acquire/release round-trips.

**After**: Acquire once, process 64 reads, release once with aggregated
`items=64` and `work_bytes` totals. A final flush at consumer exit handles
any remainder < 64.

Eliminates ~274M atomic operations. The permit controller still sees
accurate throughput because it receives batched `(items, bytes, work_ns)`
tuples.

### 4. Fix head-of-line blocking with non-blocking round-robin (`assignBarcodes.c`)

**This is the structural fix.** The other three are micro-optimizations.

**Root cause**: The consumer loop iterated sets 0..N-1 in order and
**blocked** on `pthread_cond_wait(&set->can_consume, &set->mutex)` if
set `i` had no data. With 10 reader sets (one per guide lane), a consumer
blocked on set 0 could never check sets 1–9 even if they had data ready.
This serialized lane draining: set 0 drained first, then set 1, etc.

**Fix**:
- Replace `pthread_mutex_lock` + `pthread_cond_wait` with
  `pthread_mutex_trylock` + non-blocking data check
- If trylock fails (another consumer holds the lock) or no data available,
  skip to the next set immediately
- Round-robin starting position (`rr_start`) advances after each successful
  drain, preventing all consumers from piling on set 0
- If a full sweep finds no data from any set, increment `empty_sweeps`
- After 2+ empty sweeps, `usleep(100 * min(empty_sweeps, 8))` before the
  next sweep (200us..800us), avoiding CPU burn when producers are I/O-bound
- `pthread_mutex_trylock` errors distinguished: EBUSY = normal contention
  (skip), any other return = fatal (abort with diagnostic)

**Before** (blocking):
```c
for (int i = 0; i < nsets; i++) {
    pthread_mutex_lock(&set->mutex);
    while (set->filled < lines_per_block && !set->done)
        pthread_cond_wait(&set->can_consume, &set->mutex);
    // blocked here even if sets i+1..N have data
```

**After** (non-blocking round-robin with sleep backoff):
```c
if (empty_sweeps >= 2)
    usleep(100 * (empty_sweeps > 8 ? 8 : empty_sweeps));

for (int sweep = 0; sweep < nsets; sweep++) {
    int i = (rr_start + sweep) % nsets;
    int rc = pthread_mutex_trylock(&set->mutex);
    if (rc == EBUSY) continue;
    if (rc != 0) { /* fatal */ }
    if (set->filled < lines_per_block) { unlock; continue; }
    rr_start = (i + 1) % nsets;  // advance for next iteration
```

**Probe-exit safety**: consumers check `probe_only && chem_detect->done`
at the top of the while loop and signal `can_produce` on all sets before
breaking. Producers also check `chem_detect->done` in the `can_produce`
wait condition and use a `goto producer_probe_exit` path to close handles
and set `set->done`. Without this, the 8K ring buffer can deadlock: a
producer blocks on buffer-full, all consumers have exited, nobody drains.

At smoke scale (100K reads), throughput difference is within noise
(9.99s vs 10.24s). The real impact is at full scale with 10 lanes × 28M
reads where the old code forced sequential lane draining.

## Attempted but Deferred

### Batch consumption (drain N reads per lock acquisition)

Attempted draining up to 16 reads from the ring buffer under a single lock
instead of 1 read per lock. The implementation redirected `barcode_lines` /
`forward_lines` / `reverse_lines` char** pointers into a batch buffer.

**Result**: Caused a process hang during smoke testing. The processing code
is ~200 lines of inline logic with multiple function calls that expect
stable `char**` pointers; the redirection interacted poorly with something
in the call chain. Reverted and deferred.

This optimization is still worth pursuing but needs a different approach:
either extract the processing into a function that takes explicit line
buffers, or use a dedicated batch struct that doesn't alias the `lines[]`
array.

### Skip quality line for feature reads

The barcode quality line (`barcode_lines[1]`) IS used in
`checkAndCorrectBarcode` for quality-score-based error correction.
Feature quality lines (`forward_lines[1]`, `reverse_lines[1]`) are NOT
used, but skipping only those saves ~30 bytes per read — marginal.
Cancelled.

### Replace `strcpy` with `memcpy` + known length

Already covered by opt 2's `stpcpy` change, which avoids the double-scan
without requiring length metadata in the ring buffer. Cancelled.

## Verification

All four smoke tests used identical parameters and produced identical
feature counts:

| Optimization | Feature Counts | Unmatched | Time (1M reads) |
|-------------|---------------|-----------|-----------------|
| Baseline (6452600) | 708,057 | 286,905 | 10.68s |
| + ring buffer 8K | 708,057 | 286,905 | 9.96s |
| + stpcpy | 708,057 | 286,905 | 10.24s |
| + batch permits | 708,057 | 286,905 | 10.29s |
| + HOL fix | 708,057 | 286,905 | 9.84s |
| + spin/deadlock/error fixes | 708,057 | 286,905 | 9.99s |

Note: smoke-scale timings are noisy (±5%) and don't reflect the HOL fix's
impact which manifests at full scale with multiple lanes running in parallel.

## Files Changed

| File | Lines | What |
|------|-------|------|
| `core/features/process_features/include/common.h` | +1/-1 | Buffer size constant |
| `core/features/process_features/src/assignBarcodes.c` | +98/-46 | Consumer hot path in `consume_reads()`, producer probe-exit in `read_fastqs_by_set()` |

## Next Steps

1. **Full-scale benchmark** on `EBs2_2` corrected (278M guide reads, 10
   lanes) to measure actual wall-clock improvement from the HOL fix
2. **Parity comparison** with CellRanger on the common filtered barcode set
3. **Batch consumption** (revisit with a safer implementation approach)
4. Update `README.md` benchmark table and `docs/feature_barcodes.md`

## Build

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
```
