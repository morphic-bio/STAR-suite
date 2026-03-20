# Runbook: Hash Screen Replay Harness + Tiered Lookup + H2

**Date**: 2026-03-19
**Branch context**: `benchmark-flex`
**Parent**: `docs/HANDOFF_FLEX_MAPPING_PHASE_OPTIMIZATION_20260319.md`
**Status**: IN-PROGRESS — Phase 1-4 validated (0 mismatches), Phase 5+ pending

## Executive Summary

Build a standalone replay harness for the hash screen lookup (following
the clique replay pattern), then use it to develop and prove a tiered
lookup optimization and test H2 cache variants — all without running STAR
for each iteration.

**Approach**:
1. Dump hash screen inputs/outputs from a single STAR 100K run
2. Build a standalone replay tool that can run flat and tiered lookups
3. Prove tiered == flat on the dump (fast, seconds per run)
4. Iterate on the tiered implementation and H2 cache in the harness
5. Port the validated code back into `FlexHashScreen.cpp`
6. Final STAR 100K validation

## Why a Harness First

The hash screen is a pure function:

```cpp
FlexHashScreenDecision classifyRead(const char* readSeq, uint32_t readLen, uint16_t sampleIdx)
```

It depends only on the loaded cache records and three scalar inputs. No
STAR state, no genome, no alignment context. This makes it ideal for
extraction into a standalone replay tool — exactly the pattern that worked
for clique correction (P0 in the Solo runbook, validated in 152 ms on
623K groups).

A full STAR 100K run takes ~37s and requires genome load, alignment, and
Solo processing. The hash screen replay will process 800K reads in <1s,
giving a ~100x faster iteration cycle.

## The Problem: One Flat Array

Current `FlexHashScreen.cpp` stores all cache records in a single sorted
`std::vector<Record>` (8M records, 185 MB). For each read, `classifyRead`
does up to 3 offsets × 2 binary searches = 6 binary searches through the
entire 8M array. Each search is 23 steps; with the array 5× larger than
L3 (36 MB), the later steps hit main memory at ~100 ns each.

The fix: split into per-class tiers (H0 exact, H1 near-match, deny) at
load time. The fast path — H0 exact match at offset 0 — resolves ~70-80%
of reads with a single search through a small, L2/L3-resident array.

## Phase 1: Add Dump Mode to STAR

### 1a: Dump format

Add an env-gated dump in `ReadAlign_oneRead.cpp`, triggered by
`STAR_DUMP_HASH_SCREEN=<path>`. For each call to `classifyRead`, write:

```
Per-read record (binary):
  uint32_t  readLen        (4 bytes)
  uint16_t  sampleIdx      (2 bytes)
  uint8_t   action         (1 byte: 0=Disabled, 1=Pass, 2=Keep, 3=Deny)
  uint8_t   cacheClass     (1 byte)
  uint16_t  geneIdx15      (2 bytes)
  int8_t    offset         (1 byte)
  uint8_t   _pad           (1 byte)
  char      readSeq[readLen]  (variable, ASCII A/C/G/T/N)
```

File header:
```
char     magic[8]    = "HSCRN001"
uint64_t nReads      (filled on close, or use append mode)
```

The dump captures both inputs (readSeq, readLen, sampleIdx) and the
decision (action, geneIdx15, cacheClass, offset) so the replay tool can
verify against ground truth.

### 1b: Thread safety

Multiple threads call `classifyRead` concurrently. Options:
- **Per-thread dump files** (simplest): each thread writes its own file,
  merge after. Follow the same pattern as the clique dump which runs
  after the per-thread merge.
- **Single file with mutex**: simpler to consume but adds contention.
  Acceptable for 100K harness (not a timing run).

For the 100K harness (~800K reads), a single mutex-protected file is fine.

### 1c: Dump site

In `ReadAlign_oneRead.cpp`, after line 314:

```cpp
hashScreenDecision_ = FlexHashScreenCache::instance().classifyRead(
    Read0[0], readLengthOriginal[0], hashScreenSampleIdx);
// --- dump here if env set ---
```

Two dump sites:
1. **Active path** (inside the hash screen block): dump after `classifyRead`
   with the full decision (Keep/Deny/Pass) and the runtime sampleIdx. This
   covers short reads (Pass), N-base reads (Pass), cache misses (Pass), and
   all cache hits.
2. **Disabled path** (else branch): dump with `action=Disabled, sampleIdx=0`.
   These reads never called `classifyRead` (STAR routing skipped them).
   The replay tool counts them but does not compare against cache lookups.

### 1d: Run the dump

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
STAR_DUMP_HASH_SCREEN=/storage/downsampled_100K/SC2300771/results/hash_screen_dump.bin \
  tests/run_flex_hash_screen_internal_100k.sh
```

Expected: ~800K read records, ~50 MB dump file.

## Phase 2: Build the Replay Tool

### 2a: Tool location

`flex/tools/hash_screen_replay/` (parallel to `flex/tools/clique_replay/`):

```
flex/tools/hash_screen_replay/
  Makefile
  hash_screen_replay.cpp    — replay tool
  FlexHashScreenFlat.cpp    — copy of current flat lookup (standalone)
  FlexHashScreenTiered.cpp  — tiered lookup implementation
  FlexHashScreenCommon.h    — shared Record type, cache loader, encodeWindow
  test_hash_screen.cpp      — unit tests
```

### 2b: Standalone cache loader

Extract the cache-loading logic from `FlexHashScreen.cpp::loadFile` into
`FlexHashScreenCommon.h`. This reads the FH01SEQ1 binary format and
returns a `vector<Record>`. No dependency on `ParametersSolo`.

```cpp
struct Record {
    uint64_t seqLo, seqHi;
    uint32_t resolvedGeneIdx15;
    uint8_t  cacheClass;
    uint8_t  negativeCode;
    uint16_t sampleIdx;
};

bool loadCacheRecords(const char* path, std::vector<Record>& out);
```

### 2c: Flat lookup (reference implementation)

Copy the current `classifyRead` / `findRecord` / `classifyHits` /
`encodeWindow` logic into `FlexHashScreenFlat.cpp`. This is the ground
truth. It must produce identical output to STAR's in-process decisions.

### 2d: Tiered lookup (new implementation)

`FlexHashScreenTiered.cpp` implements the same `classifyRead` signature
but with records partitioned at init time:

```cpp
class TieredCache {
    std::vector<Record> h0_;     // cacheClass 0, gene > 0
    std::vector<Record> h1_;     // cacheClass 1, gene > 0
    std::vector<Record> deny_;   // cacheClass 2, ProbeAmbig
    // dead-weight records not loaded
public:
    void init(const std::vector<Record>& allRecords);
    FlexHashScreenDecision classifyRead(const char* seq, uint32_t len, uint16_t sampleIdx) const;
};
```

**First iteration (correctness)**: search all tiers at all offsets, same
total logic as flat. This validates the partition.

**Second iteration (fast path)**: add the H0-offset-0 short-circuit.
Must still produce identical decisions for all reads.

### 2e: Replay tool

```
Usage: hash_screen_replay <cache.bin> <dump.bin> [--mode flat|tiered|both]
                          [--summary] [--diff <out.tsv>]
```

For each read in the dump:
1. Call the selected lookup (flat, tiered, or both)
2. Compare against the ground-truth decision stored in the dump
3. Report mismatches

Output:
- Summary: total reads, KEEP/DENY/PASS counts per implementation,
  mismatches between flat and tiered, mismatches against ground truth
- Optional diff TSV with per-read details for any mismatch
- Elapsed time for the lookup-only portion (no I/O in the hot loop)

### 2f: Makefile

Follow the clique replay pattern — standalone build, no STAR dependency:

```makefile
CXX = g++
CXXFLAGS = -O3 -std=c++11 -Wall -Wextra

REPLAY_TARGET  = hash_screen_replay
TEST_TARGET    = test_hash_screen

all: $(REPLAY_TARGET) $(TEST_TARGET)

$(REPLAY_TARGET): hash_screen_replay.cpp FlexHashScreenFlat.cpp \
                  FlexHashScreenTiered.cpp FlexHashScreenCommon.h
	$(CXX) $(CXXFLAGS) -o $@ hash_screen_replay.cpp \
	       FlexHashScreenFlat.cpp FlexHashScreenTiered.cpp

$(TEST_TARGET): test_hash_screen.cpp FlexHashScreenFlat.cpp \
                FlexHashScreenTiered.cpp FlexHashScreenCommon.h
	$(CXX) $(CXXFLAGS) -o $@ test_hash_screen.cpp \
	       FlexHashScreenFlat.cpp FlexHashScreenTiered.cpp

test: $(TEST_TARGET)
	./$(TEST_TARGET)

clean:
	rm -f $(REPLAY_TARGET) $(TEST_TARGET)
```

## Phase 3: Unit Tests

`test_hash_screen.cpp` — targeted tests for both flat and tiered lookups:

| Test | What it covers |
|------|----------------|
| H0 exact match, offset 0, sample matched | Fast path KEEP |
| H0 exact match, offset +1 | Off-by-one KEEP |
| H1 near match, single gene | 1MM KEEP |
| H1 near match, gene conflict (two offsets) | Gene-conflict DENY |
| Deny record (ProbeAmbig) | Hard DENY |
| Sample-specific record, wrong runtime sample | PASS (invisible to `findRecord`†) |
| Cache miss (no record for sequence) | PASS |
| Read too short (< 50 bp) | PASS |
| Read with N base | PASS (encodeWindow fails) |
| H0 + deny at same offset | Deny overrides keep |
| H0 at offset 0 + H1 at offset +1, different gene | Gene conflict across tiers |
| Empty cache | All reads PASS |

†`findRecord` returns only an exact-sampleIdx match or a sampleIdx=0
fallback. A record at sampleIdx=5 is invisible when querying with
sampleIdx=7 — `findRecord` returns false (not found), so the read gets
PASS, not DENY. The `sampleSpecifiedMismatch` path in `classifyHits` is
unreachable with the current `findRecord` architecture; it would require
`findRecord` to return a record whose sampleIdx is nonzero and differs
from the runtime sample.

For each test, construct a small cache in memory, run both flat and tiered,
assert identical decisions.

## Phase 4: Validate Replay Against STAR Dump

### 4a: Flat replay vs STAR ground truth

```bash
cd flex/tools/hash_screen_replay && make
./hash_screen_replay \
  /storage/.../sequence_cache.bin \
  /storage/.../hash_screen_dump.bin \
  --mode flat --summary
```

**Expected**: zero mismatches. Every read's flat-replay decision must match
the STAR dump exactly. If not, the standalone flat implementation has a bug.

### 4b: Tiered replay vs flat replay

```bash
./hash_screen_replay \
  /storage/.../sequence_cache.bin \
  /storage/.../hash_screen_dump.bin \
  --mode both --summary --diff mismatches.tsv
```

**Expected**: zero mismatches between flat and tiered. The tier partition and
lookup order must produce the exact same KEEP/DENY/PASS for every read.

### 4c: Timing

The replay tool reports lookup-only elapsed time:

| Implementation | 800K reads | Per-read |
|---------------|-----------|----------|
| Flat (8M records) | ? ms | ? ns |
| Tiered (H0 + H1 + deny) | ? ms | ? ns |

Even on the 100K dump, the timing ratio (flat vs tiered) is meaningful
because it reflects the binary-search depth difference. At full scale
(2B reads), the per-read savings multiply.

### Phase 4 Results (2026-03-20)

**Dump**: 800K reads, 78 MB → `/storage/downsampled_100K/SC2300771/results/flex_hash_screen_dump_20260320_005618/hash_screen_dump_v2.bin`

**Flat vs STAR ground truth**: 0 mismatches (670,450 KEEP / 5,396 DENY / 124,154 PASS)
**Tiered vs flat**: 0 mismatches
**Tiered vs truth**: 0 mismatches

**Tier composition from loaded cache** (8,070,044 total records):

| Tier | Records | Size | Search depth (log2) |
|------|--------:|-----:|--------------------:|
| H0 (class 0) | 53,294 | 1.3 MB | 16 |
| H1 (class 1) | 7,360,125 | 176.6 MB | 23 |
| Deny (ambig) | 636,746 | 15.3 MB | 20 |
| Dropped (dead weight) | 19,879 | 0.5 MB | — |

**Timing** (800K reads, i9-13900KF):

| Implementation | Elapsed | Per-read |
|---------------|---------|----------|
| Flat (8M records) | 958 ms | 1,197 ns |
| Tiered (no short-circuit) | 1,567 ms | 1,959 ns |

Tiered is currently slower because it searches 3 tiers × 3 offsets = 9 binary searches
vs flat's 3. The H0-offset-0 short-circuit (Phase 5b) will reduce most reads to a
single 16-step binary search through the 53K-record H0 tier (1.3 MB, fits in L2).

**Note on KEEP count discrepancy**: The dump records raw `classifyRead` decisions
(670,450 KEEP). STAR's `Log.final.out` reports 659,863 KEEP because ~10,587
classifyRead-KEEP decisions are not "handled" by `record_flex_hash_screen_keep`
(these fall through to the aligner as PASS). The dump captures the hash-screen-level
decision, which is the correct surface for replay validation.

## Phase 5: Iterate on Tiered Implementation

With the harness validated, iterate rapidly:

### 5a: Measure tier composition

The replay tool (or a `--stats` flag) reports the actual tier sizes from
the loaded cache:

| Tier | Records | Size | Search depth |
|------|--------:|-----:|-------------:|
| H0 | ? | ? MB | log2(?) |
| H1 | ? | ? MB | log2(?) |
| Deny | ? | ? MB | log2(?) |
| Dropped (dead weight) | ? | ? MB | — |

### 5b: Add H0-offset-0 short-circuit

After the simple tiered partition is validated (all tiers searched at all
offsets), add the fast path:

```
1. Encode offset-0 window
2. Search H0 for (seq, sampleIdx)
3. If found + sampleMatched + no deny at this seq → immediate KEEP
4. Otherwise fall through to full multi-tier, multi-offset logic
```

**Correctness**: re-run `--mode both` (short-circuit tiered vs original
flat) — must be zero mismatches.

**Timing**: measure the speedup from the short-circuit.

### 5c: Test with H2 cache (if available)

Build a larger cache from 1M-read pilot (Phase 7), load it in the
replay tool, compare flat vs tiered on the same dump. The dump reads
are fixed (100K harness), but the cache changes — more reads will get
KEEP instead of PASS.

### 5d: Parameter sweeps

The harness makes it trivial to test:
- Different cache files (H1, H2, pruned variants)
- The effect of removing dead-weight records
- Lookup strategies (H0-first, deny-first, etc.)

Each test is a sub-second replay run.

## Phase 6: Port Back to FlexHashScreen.cpp

Once the tiered implementation is validated in the harness:

### 6a: Data structure change (FlexHashScreen.h)

```cpp
// Replace:
std::vector<Record> records_;

// With:
std::vector<Record> h0Records_;
std::vector<Record> h1Records_;
std::vector<Record> denyRecords_;
```

### 6b: Load-time partition (FlexHashScreen.cpp::loadFile)

After reading all records from disk, partition into tiers (same code as
`TieredCache::init` from the harness). Dead-weight records dropped.

### 6c: Tiered classifyRead

Port `TieredCache::classifyRead` into `FlexHashScreenCache::classifyRead`.
Same signature, same semantics, same decisions.

### 6d: Files changed

| File | Change |
|------|--------|
| `core/legacy/source/FlexHashScreen.h` | Replace `records_` with per-tier vectors |
| `core/legacy/source/FlexHashScreen.cpp` | Partition at load, tiered lookup |

No other files change — `ReadAlign_oneRead.cpp` calls `classifyRead()`
with the same signature.

## Phase 7: Validate in STAR

### 7a: 100K harness

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
tests/run_flex_hash_screen_internal_100k.sh
```

Compare against baseline:
- `Log.final.out` KEEP/DENY/PASS counts: **must be identical**
- Raw MEX sorted diff: **must be zero**
- Per-sample filtered MEX: **must be zero diff**

### 7b: Dump + replay round-trip

Run the dump again with the tiered STAR binary:

```bash
STAR_DUMP_HASH_SCREEN=dump_tiered.bin tests/run_flex_hash_screen_internal_100k.sh
# Compare tiered STAR dump vs original STAR dump
diff <(sort dump_flat.tsv) <(sort dump_tiered.tsv)
```

### 7c: Memory and timing

| Metric | Flat baseline | Tiered |
|--------|---------------|--------|
| Cache memory | 185 MB | ? (H0 + H1 + deny, no dead weight) |
| VmRSS after genome load | ? | ? |
| Hash screen timing (if logged) | ? | ? |

## Phase 8: H2 Cache (Larger Pilot)

### 8a: Build 1M-read downsample

Use the existing pilot pipeline with a larger input set.

### 8b: Build + prune H2 cache

Run `scripts/run_flex_h01_pilot.sh` with 1M reads. The tiered loader
automatically drops dead weight.

### 8c: Test H2 in replay harness first

```bash
./hash_screen_replay \
  H2_cache.bin \
  hash_screen_dump.bin \
  --mode tiered --summary
```

This shows how the H2 cache changes the KEEP/PASS split on the same 100K
reads. No STAR run needed.

### 8d: STAR 100K with H2 cache

If the harness looks good, run STAR with the H2 cache:

```bash
STAR_FLEX_HASH_SCREEN_CACHE=H2_cache.bin tests/run_flex_hash_screen_internal_100k.sh
```

Compare raw MEX against H2 baseline (not H1 baseline — the KEEP/PASS
split is intentionally different).

### 8e: Full-scale benchmark (only if 8d passes)

Full 2B-read benchmark with tiered lookup + H2 cache.

## Hard Constraints

1. **Always clean-rebuild before any STAR run.**
2. **Serialize benchmark runs.**
3. **Tiered lookup must produce identical decisions to flat** — zero
   mismatches in the replay tool on the 100K dump.
4. **Raw-MEX parity is the primary correctness surface** for STAR runs.
5. **Cache is read-only at runtime** — tiering is a load-time transform.
6. **Dump mode must be env-gated** — zero overhead when not dumping.
7. **No cross-tier key duplicates.** The flat sort key is `(seqHi, seqLo,
   sampleIdx)` — `cacheClass` is not included. If two records share the same
   key but differ in `cacheClass`, flat's `findRecord` returns whichever
   `std::sort` placed first (unstable), while the tiered model searches
   H0→H1→deny deterministically. The replay tool reports cross-tier
   duplicates via `--stats`; any nonzero count is a cache-builder bug that
   invalidates the flat==tiered proof. The current H1 cache has zero such
   duplicates (verified 2026-03-20).
8. **Dump covers active and disabled reads.** Reads that enter the hash
   screen block are dumped with their full decision. Reads that skip the
   block (disabled config, missing soloRead, etc.) are dumped with
   `action=Disabled, sampleIdx=0`. The replay tool counts but does not
   compare disabled reads against cache lookups (Disabled is a STAR
   routing decision, not a cache decision).

## Key Paths

| Item | Path |
|------|------|
| Current H1 cache | `/storage/downsampled_100K/SC2300771/results/flex_h01_full_cache_20260315_153914/reclassified/sequence_cache.bin` |
| 100K FASTQs | `/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_*` |
| 100K harness | `tests/run_flex_hash_screen_internal_100k.sh` |
| FlexHashScreen.h | `core/legacy/source/FlexHashScreen.h` |
| FlexHashScreen.cpp | `core/legacy/source/FlexHashScreen.cpp` |
| Call site | `core/legacy/source/ReadAlign_oneRead.cpp:314` |
| Clique replay (template) | `flex/tools/clique_replay/` |
| **New replay tool** | `flex/tools/hash_screen_replay/` |
| Pilot builder | `scripts/flex_h01_pilot.py` |

## Implementation Order

| Step | What | Effort | STAR build? |
|------|------|--------|-------------|
| Phase 1 | Add dump mode to STAR | ~1h | Yes (once) |
| Phase 2 | Build replay tool + flat reference | ~2h | No |
| Phase 3 | Unit tests (flat + tiered) | ~1h | No |
| Phase 4 | Validate flat replay vs STAR dump | ~10m | No |
| Phase 5a | Simple tiered partition | ~1h | No |
| Phase 5b | H0 short-circuit | ~1h | No |
| Phase 6 | Port to FlexHashScreen.cpp | ~1h | Yes (once) |
| Phase 7 | STAR 100K validation | ~30m | Yes (once) |
| Phase 8 | H2 cache | ~2h | Optionally |

Total estimated: ~10h of work, only 3 STAR builds needed.

## What This Runbook Does NOT Cover

- I/O pipeline optimization (M1 — separate track)
- Alignment parameter tuning (M3)
- Thread utilization profiling (M5)
- Solo post-map optimizations (already done in P1-P4)
- Changes to the cache-building pipeline (`flex_h01_pilot.py`)
