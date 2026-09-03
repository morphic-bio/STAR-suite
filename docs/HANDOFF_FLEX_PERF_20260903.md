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
| STAR-Flex, plus parallel sumThreads | 4:16.50 | 150 s | 86 s | — |
| STAR-Flex, plus shared ambiguous store | 4:05.11 | 150 s | 89 s | 2002% |
| STAR-Flex, plus inlined accessors (pushed) | 3:59.94 | 150 s | 70 s | 2114% |
| STAR-Flex, plus parallel fan-in (uncommitted) | 3:58.40 | 150 s | 68 s | — |
| STAR-Flex, target | < 3:46 | | | |

The gap to cyto went from 1.72x to 1.13x. Every change so far produced
**byte-identical** per-sample matrices, barcodes and features.

## Read this before touching anything

**1. There are two copies of several source files. Only one is compiled.**

**It is decided per file, not per directory. Do not assume either tree wins.**
`make` runs with `cwd = core/legacy/source` and `-I.` first, so a bare name in
the source list resolves under `core/legacy/source`, while `flex/source` files
are pulled in by explicit path or through `LIBFLEX_DIR`. Verified against the
build log:

| Compiled (edit this) | Stale copy (do NOT edit) |
|---|---|
| `flex/source/libflex/FlexFilter.cpp` | `core/legacy/source/libflex/FlexFilter.cpp` |
| `flex/source/SampleDetector.{h,cpp}` | — |
| `flex/source/SoloFeature_collapseUMI_fromBuckets.cpp` | — |
| `core/legacy/source/SoloFeature_flexfilter.cpp` | `flex/source/SoloFeature_flexfilter.cpp` |
| `core/legacy/source/SoloReadFeature.cpp` | `flex/source/SoloReadFeature.cpp` |
| `core/legacy/source/flex/SoloReadFeature_record_flex.cpp` | `flex/source/SoloReadFeature_record_flex.cpp` |
| `core/legacy/source/SoloFeature_sumThreads.cpp` | — |

**The direction reverses between rows, including between two files with almost
the same name.** `FlexFilter.cpp` compiles from `flex/source`, but
`SoloFeature_flexfilter.cpp` compiles from `core/legacy/source`. I lost twenty
minutes reading the stale `FlexFilter.cpp`, which claims it gives all cores to
the Monte Carlo while the compiled one splits them evenly across tags — and
then, knowing all that, still edited the stale `SoloFeature_flexfilter.cpp`
later the same evening. That edit is comment-only so it is harmless, but it is
in commit f0f4fc8 and it is dead code. Treat this table as the only source of
truth and re-derive it from the build log if anything looks off.

**Check every time, with the build log, not with intuition:**
```bash
grep -o "[^ ]*YourFile\.cpp" <build log> | sort -u   # which path was compiled
grep -a -c "some string you added" core/legacy/source/STAR   # did it land
```
When a log line disagrees with the source you are reading, you are reading the
wrong file.

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

## THE TAIL IS NOW FULLY MEASURED — read this before optimising anything

Every phase is timed in the code (`Solo timing:` lines in `Log.out`). Do not
work from profile percentages: they are *sample* shares across 32 threads, and
reading them as wall-clock shares cost me four wrong diagnoses in a row.

| Phase of the ~68 s tail | Wall | Notes |
|---|---|---|
| **`resolvePendingAmbiguousToHash`** | **17.2 s** | **serial; see below** |
| Bucket loop (parallel, 32 threads) | 14.2 s | sort 7 s, UMI correct 4.5 s, load 2 s |
| Flexfilter | 20 s | already optimised this session, 68 s → 19 s |
| MEX write | 5 s | |
| After flexfilter | 8 s | not yet broken down |
| `sumThreads` | 2.5 s | was 26 s |
| Fan-in | 0.85 s | was ~2 s |
| Post-fan-in setup (triplet loop etc.) | 0.38 s | |

### 1. `resolvePendingAmbiguousToHash` — 17.2 s, and it resolves nothing

The single largest item left, and it sat untimed between two log timestamps all
session. The log line it prints is the striking part:

```
[AMBIG-CB-RESOLVE] pending=2627928 resolved=0 still_ambiguous=2627928 added_to_hash=0
```

**Seventeen seconds to walk 2.6 M ambiguous barcodes and resolve zero of them.**

**Do not "optimise" this by skipping it until someone establishes which of two
things is true**, because they need opposite fixes:

- **It is correctly a no-op on this data** (the resolver needs evidence this
  dataset does not provide), in which case the 17 s is dead work and the pass
  should exit early — a large, free win.
- **It is silently failing to resolve barcodes it should resolve**, in which
  case this is a *correctness* bug wearing a performance costume, and cutting
  the pass would lock in wrong answers. `resolved=0` out of 2.6 M is a
  suspicious number.

Start at `SoloFeature::resolvePendingAmbiguousToHash` in
`core/legacy/source/SoloFeature.cpp`; the early exits to check are
`g_disableAmbigResolve` and `!pSolo.cbCorrector`. Check what `cbCorrector` is on
this command line before assuming either branch.

There is an irony worth recording: this session made the *accumulation* of these
same 2.6 M entries nearly free (26 s → 2.5 s, commit 022b24d) while their
*resolution*, a few lines later, cost seven times more and nobody had timed it.

### Other targets, in the order I would take them

1. **`sumThreads`, 26 s, serial — now measured, sub-timers are in the code.**
   The loop's four steps on the JAX set:

   | step | time |
   |---|---|
   | merge (flush bucket segments + fold pending-ambiguous) | 14.56 s |
   | teardown of the thread's own tables | 11.27 s |
   | stats add | 0.000018 s |
   | count add | 0.0096 s |

   Both significant parts are **per-thread work**, which is why the right shape
   is each thread doing its own as its final act rather than one thread walking
   all 32 afterwards — the same pattern as the alignQ fix. Do NOT bother with a
   parallel reduction over the count arrays: that step costs 9 ms.

   The flush split out at 0.069 s, so it was never the cost. The read-identifier
   tracker is *not* involved either: it is gated on BAM output and these runs
   use none, so it is never allocated. The 14.56 s is entirely the
   pending-ambiguous fold, and reading `mergePendingAmbiguous` explains why —
   each of the 2.6 M entries carries several vectors and two maps, and a key
   present in more than one thread (the common case for an ambiguous barcode)
   walks and combines every field.

   **Done: the loop is now parallel** (`#pragma omp parallel for`, fold inside
   `critical(soloSumThreadsMerge)`, flush and teardown outside). Result: block
   wall 26 s → 20 s, tail 89 s → 86 s, wall 4:19 → 4:16.50, byte-identical.
   Less than the ~11 s hoped, and the reason is worth knowing: with 32 threads
   freeing large structures concurrently, allocator and memory-bandwidth
   contention slowed the serialized fold itself from 14.6 s to 19.1 s of
   thread-summed time. The teardown now hides behind the fold, so the fold is
   the floor.

   **Done, and it went further than sharding the fold.** If the accumulation is
   striped anyway, there is no reason to keep it per thread: both flex write
   sites now take a handle into a shared 256-way striped store, holding that
   stripe's lock for the entry update. A barcode seen by twenty threads is
   combined once, in one place, so the fold does not get faster — it stops
   existing. `sumThreads` then *adopts* the stripes by moving entries into the
   master.

   | step | before | after |
   |---|---|---|
   | fold | 14.6 s | 0.00006 s |
   | teardown | 11.3 s | 0.00003 s (thread maps are empty) |
   | adopt shared (new) | — | 1.64 s (2,627,928 entries) |
   | flush | 0.07 s | 0.87 s |
   | **block wall** | **26 s** | **~2.5 s** |

   Tail 86 s → 74 s, wall 4:16.50 → 4:05.11, CPU 2002% → 2123%,
   outputs byte-identical.

   **Mapping did not move: 150 s before and after.** That was the live risk —
   striping puts a lock on a path inside mapping, and the worry was that it
   would cost more in locality than it saved in the tail. It did not, because
   the lock is on the *ambiguous* path (a few million updates over 150 s across
   256 stripes), not the read path. Contrast the sample-token mutex in fix A:
   one mutex, 1.7 billion acquisitions. Rate and striping are what separate a
   safe shared structure from a catastrophic one.

   Enabled only for the fused Flex path (`sharedAmbigEnable` in
   `mapThreadsSpawn.cpp`); the staged pipeline keeps per-thread maps and its
   existing merge. The adopt step **exits** if a key is already present in the
   destination, which would mean both accumulation paths ran; silently keeping
   one side would corrupt counts.

   Three different risk levels, and they should not be treated alike:
   - **Teardown, 11.3 s, zero risk.** Freeing a thread's own tables touches
     nothing shared and has no ordering. Move it into the thread's exit.
   - **Bucket-segment flush, no shared-master contact.** Goes to the bucket
     store's per-bucket locks. Also movable. Timed separately now.
   - **Pending-ambiguous fold, order-sensitive.** Folding in thread-completion
     order is nondeterministic, and the merged table is later iterated by
     `resolvePendingAmbiguousToHash`. If any tie-break there depends on
     iteration order, outputs move. Either keep this in fixed thread order, or
     move it and *prove* byte-identity with the manifest check — do not assume.

   File: `core/legacy/source/SoloFeature_sumThreads.cpp`. Removing this block
   takes the tail from 89 s toward 65 s and the wall from 4:19 to about 3:55.

2. **Bucket collapse, now 67 s and the whole remaining tail.** Profile says
   ~24% `__introsort_loop` over `PackedCbRecord`, ~11% `UMICorrector::findNeighbors`.

   **A radix sort was tried and is SLOWER. Do not repeat it.** LSD radix on
   `groupSortKey`, 8 passes of 8 bits, skipping single-digit passes, one scratch
   buffer per thread: collapse 67.4 s → 74.2 s, wall 3:59.94 → 4:06.37, outputs
   byte-identical (so it was correct, just worse). The reason is memory, not
   comparisons: each pass moves all 6.4 M records of a bucket, so eight passes
   move ~800 MB per bucket, and the scatter writes 256 streams per thread which
   at 32 threads thrashes cache and TLB. `std::sort` partitions in place and
   touches far less memory. **The collapse is memory-bound, not
   comparison-bound**, which also means micro-optimising the comparator further
   is not where the win is.

   **Hash pre-aggregation was also tried and is a WASH. Do not repeat it
   either.** Collapsing exact-duplicate keys with the vendored open-addressing
   table *before* the sort, to shrink what gets ordered: collapse 67.41 s →
   67.59 s, wall 3:59.94 → 3:59.94, byte-identical. It cannot help, and the
   reason is in the log all along:

   ```
   [CB-BUCKET] streamed_records=1626999647 aggregated_counts=1626999647
   ```

   **Those two numbers are equal: there are no exact-duplicate keys to remove.**
   The compaction pass after the sort collapses nothing. I assumed a ~5x
   deduplication ratio and built against the guess; it is 1.0x. The 7x gap to
   `final_molecules=224941678` comes from UMI *correction* merging similar UMIs
   within a group, which happens after grouping and cannot be exploited before
   the sort. Check this ratio on any new dataset before assuming otherwise.

   **So the ordered volume is irreducible.** With both the algorithm (radix) and
   the volume (dedup) ruled out, the only remaining lever is the sort's cache
   behaviour. That points at the segment-merge design below.

   **Constraint any such design must respect.** The emission loop detects a new
   barcode with a `previousCb` / `previousCbTag` run check, so it needs globally
   *sorted* molecules, not merely grouped ones. Hash-grouped output fed straight
   in would push each barcode once per group and change the matrix. Group with a
   hash if you like, but sort before emitting.

3. **Mapping, 150 s.** `matchCBtoWL` is the top symbol and takes `std::string`
   arguments by reference; look for avoidable string construction per read
   before touching the algorithm.

## How I wasted most of an evening, so you do not repeat it

Five diagnoses of the tail, four wrong, all in the same way: reasoning from a
plausible mechanism instead of timing the thing itself.

| I believed | Reality |
|---|---|
| The bucket sort is the bottleneck (24% of profile) | 7 s of a 68 s tail. Three separate attacks on it all measured as noise. |
| A radix sort will beat `std::sort` | Slower, 67 → 74 s. Memory-bound, not comparison-bound. |
| Deduplicating before the sort cuts volume ~5x | 1.0x. `streamed_records == aggregated_counts` in the log said so all along. |
| The 21 s gap is the serial fan-in | Fan-in was ~2 s. |
| The rest is the unreserved count-matrix vector | 0.38 s. The reserve was still right, but bought nothing. |

The actual answer, a 17 s serial pass that resolves nothing, was invisible to
every one of those hypotheses because **no timer bracketed it**. Profile sample
shares pointed at the sort because the sort is what 32 threads were busy doing;
they cannot show you a serial block that runs while 31 cores idle.

**Method that worked, use it first next time.** Put a `steady_clock` timer
around every phase between two log lines, print thread-seconds *and* the
enclosing wall time, then divide. Total cost: one 4-minute run. It found the
answer immediately after four failed attempts costing an hour each.

**A trap in the timing itself**: my first bucket-loop timer started *after* four
untimed pre-loop steps, so the phase looked like it began at its log line. If
your timers do not sum to the wall time between the enclosing timestamps, the
difference is real work you have not bracketed — chase it before theorising.

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

**Uncommitted in `flex/source/SoloFeature_collapseUMI_fromBuckets.cpp`** (all
verified byte-identical, gate not yet re-run on them):

- Parallel fan-in: prefix offsets, size destinations once, buckets write into
  disjoint ranges. ~2 s → 0.85 s, wall 3:59.94 → 3:58.40. Worth keeping.
- Reserves on `countCellGeneUMI` (296.7 M entries) and the per-bucket molecule
  vectors. No measurable gain, but correct: these sizes are known exactly a few
  lines earlier and were being discovered by ~30 reallocations. An indexed write
  into a pre-sized buffer would also vectorise where `push_back` cannot.
- The phase timers themselves. **Keep them.** They are what finally located the
  17 s, they cost nothing measurable, and without them the next person starts
  from the same blind position I did.

Everything before that is committed and pushed (1016700).

The rest of the tree:

- **BGZF ingest gate re-run against these changes on 2026-09-03 19:05: all green.**
  Phase 0 and T1 through T9, plus the T7-off and T7-auto gzip/CBQ regressions.
  T9 (the fused alignQ deadlock test) passes in all three configurations.
- JAX `cbq_noalign` outputs byte-identical to the pre-optimization recorded arm.

Re-run the gate after any further change here; T9 exercises the same fused
pipeline these files touch.

Benchmark artifacts and diagnostic runs: `/home/lhhung/jax_matrix_20260903/`.
