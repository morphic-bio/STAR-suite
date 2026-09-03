# Handoff: STAR-Flex performance work, JAX Flex benchmark (2026-09-03)

## The goal, in one line

Make STAR-Flex faster on the JAX Flex benchmark **without changing a single
output byte**, because our Cell Ranger parity is the thing cyto cannot claim and
it is not negotiable. Speed that costs parity is worthless here.

## Where we stand

Same box, same 8 CBQ lanes, same NVMe, cold page cache, 32 threads, hash-only
(`--flexNoAlign 1`), which is the mode matched to cyto since cyto does not align.

| Run | Wall | Mapping | Tail | CPU |
|---|---|---|---|---|
| cyto 0.4.6 (CBQ, SSD) | 3:46.76 | — | — | 1320% |
| STAR-Flex, before this work | 6:30.07 | 230 s | 139 s | 1910% |
| STAR-Flex, lock-free sample table | 5:07.81 | 150 s | 137 s | 1642% |
| STAR-Flex, plus full-budget EmptyDrops | 4:19.34 | 150 s | 89 s | 2002% |
| STAR-Flex, target | < 3:46 | | | |

The gap to cyto went from 1.72x to 1.14x. Every change so far produced
**byte-identical** per-sample matrices, barcodes and features.

## Read this before touching anything

**1. There are two copies of several source files. Only one is compiled.**

The Makefile sets `FLEX_SRC_DIR := /mnt/pikachu/STAR-suite/flex/source` and
`LIBFLEX_DIR := $(FLEX_SRC_DIR)/libflex`. So:

| Compiled (edit this) | Stale copy (do NOT edit) |
|---|---|
| `flex/source/libflex/FlexFilter.cpp` | `core/legacy/source/libflex/FlexFilter.cpp` |
| `flex/source/SoloFeature_flexfilter.cpp` | `core/legacy/source/SoloFeature_flexfilter.cpp` |
| `flex/source/SampleDetector.cpp` | — |

I lost twenty minutes reading the stale `FlexFilter.cpp`, which says it gives all
cores to the Monte Carlo, while the compiled one splits them evenly across tags.
The stale tree is from April; the live one is current. When a log line disagrees
with the source you are reading, you are reading the wrong file. Confirm with
`grep -a <string> core/legacy/source/STAR`.

**2. Never rebuild while a timed run is in flight.** `make` overwrites the
binary a running STAR is executing. Every diagnostic script here starts by
polling `pgrep -f "[s]ource/STAR --runThreadN"` and refusing if one is up. Keep
that guard.

**3. Everything in a timed run must be on the SSD, including `/tmp`.** Reads,
CBQ, the 35 GiB genome index, outputs, `TMPDIR`, and the bucket spill directory
all live on `/home` (nvme1n1p1). `/storage/flex_filtered_reference_2024` is a
**symlink onto the disk array**, so the raw path is not usable; the matrix script
stages the index and asserts the device. Array-sourced timings are worthless and
have been thrown away before.

**4. Verify output identity on every change.** The matrix script writes a
sha256 manifest of every `matrix.mtx`, `barcodes.tsv`, `features.tsv`,
`Summary.csv` and `Barcodes.stats`. Diff it against the recorded arm. A speedup
that changes one byte is a bug until proven otherwise.

**5. Single runs, no replicates, for now.** Deliberate: correctness and
optimization first. Two runs differing by ten seconds on a six-minute job is
noise, not a result. Do not build a story on a ten-second delta.

## What has been fixed (all byte-identical)

### A. Per-read mutex in the sample detector — the big one, 90 s

`SampleDetector::detectSampleFromPackedTag` called `registerSampleToken` on
every read carrying a sample tag. That helper took a process-wide `std::mutex`
to write an idempotent 32-entry table, a store that has already happened after
the first read of each tag. Thirty-two threads, 1.7 billion reads, one mutex.

A DWARF call graph put ~30% of all CPU samples in
`native_queued_spin_lock_slowpath`, every chain arriving through that function.
Fix: make the table `std::atomic<uint16_t>`, read it lock-free first, take the
mutex only on the genuine first registration. Mapping 240 s → 150 s.

Files: `flex/source/SampleDetector.{h,cpp}`.

### B. Per-thread pipeline counters — no measurable effect, kept

Four shared atomics plus a per-lane array were bumped three or four times per
read by every fused thread. Replaced with `FlexLocalCounters`, a per-thread
tally published every 65,536 reads and on scope exit (destructor flushes, so
early returns stay exact). Measured 240 s vs 230 s mapping, i.e. nothing outside
noise. Kept because it removes real contention in principle, but **do not claim
it as a win**.

File: `core/legacy/source/FlexPipeline.cpp`.

### C. FlexFilter thread allocation — 41 s, verified byte-identical

`flex/source/libflex/FlexFilter.cpp` split the thread budget evenly:
16 tag threads, `32 / 16 = 2` Monte Carlo threads per tag. But the work is not
even. Only 4 of 16 tags run the Monte Carlo at all; the other 12 report
`TailTested=0` and return almost immediately. So four tags ran their
simulations on two threads each while 24 cores idled.

Changed to one tag at a time with the whole budget on the Monte Carlo. Results
are provably unchanged: each simulation seeds from its own iteration index
(`seed ^ it`), and per-thread tallies are merged by **integer** addition, so the
p-values cannot depend on the thread split. Verify anyway.

The first attempt changed **two** things at once (tag threads 16 → 1 *and*
sampler threads 2 → 32) and the net 68 s → 27 s badly understated the sampler
effect. The 2x2 that decomposes it, flexfilter wall time, outputs identical in
every cell:

| tag threads | mc threads | flexfilter |
|---|---|---|
| 16 | 2 | 69 s (original) |
| 1 | 2 | 207 s |
| 1 | 32 | 27 s |
| **16** | **32** | **19 s** (new default) |

Two independent effects pulling opposite ways. Serializing tags costs 138 s;
giving the sampler the full budget saves 180 s. Solving the pair gives ~15 s of
serial per-tag work outside the sampler and ~384 core-seconds of Monte Carlo,
which fits both rows. So the sampler dominates, the tag work still wants its own
parallelism, and the right answer is to give both what they want.

Final: `numThreads = min(nTags, totalCores)`, `mcThreadsPerTag = totalCores`.
Flexfilter 68 s → 19 s, wall 4:19.34. Confirmed on the shipped default with no
overrides set (`run_confirm_default.sh`): flexfilter 20 s, mapping 150 s,
tail 89 s, wall 4:19.05, CPU 2002%, outputs byte-identical.

Run-to-run variance on this benchmark is negligible — the same configuration
measured 4:26.64 and 4:26.62 in separate runs — so single runs are sound here
and differences of a few seconds are real.

`STAR_FLEXFILTER_TAG_THREADS` and `STAR_FLEXFILTER_MC_THREADS` remain as
diagnostic overrides so the split can be re-measured without a rebuild.

**Before any high-plex assay**: this deliberately oversubscribes, which is
harmless at 16 tags because only 4 have an ambiguous tail, but hundreds of tags
would spawn hundreds of sampler threads. Bound the concurrent sampler threads
before running 384-plex.

## Where the remaining time goes

Tail is 137 s, measured from the `Log.out` timestamps of the fixed run:

| Stage | Seconds | Threading |
|---|---|---|
| `SoloFeature::sumThreads` | 26 | **serial** |
| Bucket-parallel collapse (1.63 B records → 224.9 M molecules) | 30 | 32 threads |
| MEX write | 6 | |
| flexfilter (EmptyDrops) | 68 → **19** | fixed, see C |
| Finish | 7 | |

After fix C the tail is ~89 s, and `sumThreads` at 26 s is now its largest
single serial component — the next thing to attack.

Mapping is 150 s. `perf` on the fixed binary, mapping window, top self time:
`matchCBtoWL` 9%, the hash-screen `find` 6%, `processCbqModuleRecords` 6%,
`record_flex_hash_screen_keep` 3%, CBQ batch build and sequence materialize 5%.
No single dominant symbol left, so mapping is now death by a thousand cuts.

### Next targets, in the order I would take them

1. **`sumThreads`, 26 s, serial.** Walks the 32 per-thread Solo structures one
   at a time on one core: merges each thread's inline hash and its
   pending-ambiguous table (2,627,928 entries in this run) into a master, adds
   stats and counts, then destroys the thread's khash tables and swaps out its
   vectors. The merge is an associative reduction, so a tree merge over threads
   would work. The frees are pure teardown that nothing reads afterwards and
   could be handed to a detached thread. Estimate 15 to 20 s recoverable.
   File: `core/legacy/source/SoloFeature_sumThreads.cpp`.

2. **Bucket collapse, 30 s, already on 32 threads.** Profile says ~24% is
   `__introsort_loop` over `PackedCbRecord` and ~11% `UMICorrector::findNeighbors`.
   A radix sort on the packed key is the obvious move since the key is an
   integer. Lower confidence, real work.

3. **Mapping, 150 s.** `matchCBtoWL` is the top symbol and takes `std::string`
   arguments by reference; look for avoidable string construction per read
   before touching the algorithm.

## Do not do these

- **Do not remove the EmptyDrops filter.** It rejects 2,130 of 22,454 tested
  barcodes, 4% of the 52,461 final cells. Small, but our headline is cell-calling
  concordance with Cell Ranger and 4% is not free. Make it fast, not absent.
- **Do not drop the H1 hash tier.** Measured: H0-only saves 20 s of mapping (9%)
  and loses 80.6 million keeps. Bad trade. The switch exists as
  `STAR_FLEX_HASH_H0_ONLY=1` for timing only and changes results.
- **Do not touch the alignment arm** (`--flexNoAlign 0`). Deferred until
  optimization is done, by explicit instruction.
- **Do not re-tune offsets.** The fused path already checks offset 0 only;
  `kProbeStartOffset` is 0. The multi-frame scan lives in the alignment fallback
  and does not run in this mode.

## How to run things

All scripts are in `docs/benchmarks/jax_matrix_20260903/`.

```bash
# One arm of the matrix. Stages inputs if needed, verifies the SSD locus,
# refuses on a dirty tree unless ALLOW_DIRTY=1, appends to results.tsv.
REPS=1 bash docs/benchmarks/jax_matrix_20260903/run_jax_matrix.sh cbq_noalign

# A diagnostic run whose logs must NOT enter the benchmark record:
ALLOW_DIRTY=1 REPS=1 ARTIFACT_ROOT=/home/lhhung/jax_matrix_20260903/diag/myrun \
  MATRIX_LOG_DIR=/home/lhhung/jax_matrix_20260903/diag/myrun/logs \
  bash docs/benchmarks/jax_matrix_20260903/run_jax_matrix.sh cbq_noalign

# Correctness sweep across arms (fails closed if too few assertions ran)
bash docs/benchmarks/jax_matrix_20260903/check_jax_matrix.sh
```

Arms: `{gzip,bgzf,cbq}_{noalign,align}`. Staged inputs live in
`/home/lhhung/jax_stage_20260903/` (141 GiB FASTQ, 100 GiB CBQ, 35 GiB index)
and are verified by size on every start, so leave them in place.

Profiling that worked, after several that did not:

```bash
# System-wide, started before the run, SIGINT to the real perf process
sudo -n perf record -F 99 -a -g -o perf.data &
# ... run ...
sudo -n kill -INT "$(pgrep -f '^/usr/lib/linux-tools/.*/perf record' | head -1)"
sudo -n perf report -i perf.data --force --no-children --stdio -g none \
     --comm STAR --sort sym --percent-limit 1

# Call graphs on a live process: attach to the STAR pid, not the wrapper
pid=$(for p in $(pgrep -f "[s]ource/STAR --runThreadN"); do \
        [ "$(ps -o comm= -p $p)" = STAR ] && echo $p && break; done)
sudo -n perf record -F 29 --call-graph dwarf,8192 -p "$pid" -o dwarf.data -- sleep 8
```

`--call-graph lbr` produced zero samples on this host; DWARF works. Kernel
symbols only resolve under `sudo` with `--force`.

## Traps that cost me time

- `strings BINARY | grep -q PATTERN` under `set -o pipefail` returns non-zero
  when grep closes the pipe early, so a present symbol reports as absent. Use
  `grep -a -c PATTERN BINARY`.
- `ps -eo ... -p PID` ignores `-p` because of `-e` and prints every process.
  Use `ps -o ... -p PID`.
- `pgrep -f` matches your own shell wrapper when the pattern appears in the
  command line. Bracket the first character **and** check `ps -o comm=`.
- `make STAR` does not refresh the compile timestamp in the version banner; only
  `make all` runs `cleanCompileInfo`. Release builds must use `make all`.
- A checker that skips every assertion must not print a pass. `check_jax_matrix.sh`
  counts assertions and reports `INCOMPLETE` with exit 2 if too few ran.

## State at handoff

Uncommitted in `/mnt/pikachu/STAR-suite` on `master` (last commit 47bfbfd,
pushed, the fused alignQ deadlock fix):

```
core/legacy/source/FlexPipeline.cpp      per-thread counters
flex/source/SampleDetector.{h,cpp}       lock-free token table
flex/source/SoloFeature_flexfilter.cpp   comment only
flex/source/libflex/FlexFilter.cpp       Monte Carlo thread budget
docs/benchmarks/jax_matrix_20260903/     protocol + diagnostics (untracked)
```

Nothing is committed, but the tree is in a committable state:

- **BGZF ingest gate re-run against these changes on 2026-09-03 19:05: all green.**
  Phase 0 and T1 through T9, plus the T7-off and T7-auto gzip/CBQ regressions.
  T9 (the fused alignQ deadlock test) passes in all three configurations.
- JAX `cbq_noalign` outputs byte-identical to the pre-optimization recorded arm.

Re-run the gate after any further change here; T9 exercises the same fused
pipeline these files touch.

Benchmark artifacts and diagnostic runs: `/home/lhhung/jax_matrix_20260903/`.
