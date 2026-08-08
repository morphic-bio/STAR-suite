# Standalone VB/EM finalizer: the pieces already exist

Date: 2026-08-07

## Finding

The eight-node finalizer runs full STAR purely to reach the VB/EM stage. It
loads the 29.7 GB genome index, maps zero reads, frees the index, and only then
does the work that matters. `final/Log.out` shows the sequence plainly:

```text
Loaded 199138 transcript sequences for TranscriptVB sequence-context models
RAM after mapping:
RAM after freeing genome index memory:
```

The index contributes nothing. A standalone VB/EM binary needs the transcriptome
and the small annotation tables, not `Genome`, `SA` or `SAindex`.

## Almost everything needed is already in the tree

| Piece | Where | State |
|---|---|---|
| Sidecar read / compatibility / **merge** | `TranscriptVBSidecar.{h,cpp}` — `read()`, `compatible()`, `merge()` | done, used by the current finalizer |
| Equivalence classes | `Sidecar::Evidence::ec_table` (`ECTable`) | already *in* the sidecar |
| GC and fragment-length state | `Evidence::gc_counts`, `Evidence::fld_state` | already in the sidecar |
| VB/EM engine | `core/features/vbem/source/libem` | done |
| Standalone VB/EM CLI | `core/features/vbem/tools/em_quant` | done — takes ECs + lengths, no genome |
| Effective lengths / GC bias | `tools/compute_gc_bias`, `tools/compute_expected_gc` | done |
| Fragment length distribution | `tools/sample_fld` | done |
| transcript→gene map | `tests/transcriptvb/make_gene_map_from_gtf.sh`, `tx2gene.tsv` | done |

The sidecar already carries an `ECTable`, so it *is* the equivalence-class
evidence `em_quant` consumes. `transcript_vb_sidecar::merge()` already performs
the N-way merge.

**The linkage is already proven to need no genome.** The existing
`transcriptvb-sidecar-unit-tests` target links only `TranscriptVBSidecar.o` and
`libem.a`:

```make
transcriptvb-sidecar-unit-tests : TranscriptVBSidecar.o libem $(TRANSCRIPT_VB_SIDECAR_TEST_SRC)
	$(CXX) -o transcriptvb_sidecar_unit_tests $(CXXFLAGS) $(TRANSCRIPT_VB_SIDECAR_TEST_SRC) \
	    TranscriptVBSidecar.o $(LIBEM_LIB) $(LDFLAGS)
```

No `Genome.o`, no index. That target builds and passes today.

## What a standalone finalizer would do

1. `read()` the N sidecars;
2. `compatible()` to validate revision, reference, sample and input identity,
   and that ordinals 0..N-1 each appear exactly once;
3. `merge()` into one `Sidecar`;
4. load the transcriptome FASTA for sequence-context and effective lengths;
5. run VB/EM via `libem`;
6. write `quant.sf`, and `quant.genes.sf` using a transcript→gene map.

Steps 1--3 and 5 are existing code called as-is.

## What it needs from genomeDir, and what it does not

Only the small annotation tables, for the transcript→gene mapping behind
`quant.genes.sf`:

```text
        46,073,063  exonGeTrInfo.tab
        12,910,609  transcriptInfo.tab
         1,308,277  geneInfo.tab
   -------------------------------------------
     3,211,173,485  Genome        not needed
    24,900,747,856  SA            not needed
     1,565,873,619  SAindex       not needed
```

60 MB of tables against 29.7 GB of index. Alternatively a plain `tx2gene` TSV,
for which tooling already exists.

## Expected benefit

**Time.** The genome load in the finalizer was 4 s in the measured runs, because
it landed on a node where a worker had just warmed the page cache. On a cold
node it is about 21 s with the parallel loader and striped index, and was 183 s
before those changes. So the saving is roughly 20 s in the general case, and the
finalizer stops being coupled to index load performance at all.

**Memory, which matters more.** The finalizer's peak RSS was 35.37 GB, almost
entirely the index it then frees. Without it the finalizer should need a small
multiple of the merged evidence rather than tens of gigabytes. That directly
addresses the constraint flagged for large runs in
`docs/RUNBOOK_8NODE_ASSAY_SHARDING_TESTS_20260807.md`, where finalizer memory is
the binding limit for the JAX 16-plex Flex arm.

**Scheduling.** A finalizer that needs neither 120 GB nor a full RM node can run
on `RM-shared`, which shortens queue time and reduces the cost of the serial
tail. It also decouples the finalizer from the STAR binary's genome-format
compatibility.

## Caveats to settle before building

- The sidecar embeds `source_revision` and refuses cross-revision evidence.
  A separate binary needs a defined versioning relationship with the STAR that
  wrote the sidecars, rather than inheriting the check by accident.
- `quant.genes.sf` gene-level aggregation must reproduce what STAR currently
  emits, including how it handles transcripts absent from the map.
- Effective-length and GC-bias handling must match the in-STAR path exactly, or
  results will differ for reasons unrelated to sharding. This is the most likely
  source of divergence and should be the first parity check.
- Parity is judged against the existing bands, not exact equality: Pearson
  `>= 0.99` generally, with the chr22 smoke's `0.99998` NumReads and `0.99995`
  TPM as the tighter reference. TranscriptVB is deliberately order-dependent.

## VB/EM parallelism and workload shape

The engine is already parallel, and already splits per transcript.

`vb_engine.cpp` (and `em_engine.cpp` identically):

- **E-step**: `#pragma omp parallel for schedule(dynamic)` over equivalence
  classes. Each thread accumulates into its own slice of a thread-local buffer,
  `expected_counts_tls`, so there are no atomics on the hot path.
- **Reduction / M-step**: `#pragma omp parallel for schedule(static)` over
  transcripts, summing the per-thread buffers. The comment states the fixed
  thread order is deliberate, for determinism.
- The `#pragma omp critical` blocks are inside debug-trace branches guarded by
  `params.debug_trace && !params.debug_file.empty() && !debug_indices.empty()`,
  so they are inert in production.

So the EM kernel is deterministic given identical equivalence classes. That is
consistent with the run-to-run variation originating upstream, in Salmon's
auto-detection during the first rounds, rather than in the reduction.

### Where the finalizer time actually goes

Measured on the eight-node run (`final/Log.out`):

```text
Merged 232633 equivalence classes from 32 threads
GC bias: collected 23335011 fragment observations
GC bias: dynamic update at iteration 11
Quantification converged: yes, iterations: 1863
```

The shape is **many iterations over a modest number of equivalence classes**:
232,633 ECs and 199,138 transcripts, iterated 1,863 times.

That has two consequences which explain why the finalizer phase showed only
about 5.6x on 32 threads while the kernel itself is well parallelized:

1. **Per-iteration parallel overhead is amortized over little work.** 232,633
   ECs across 32 threads is roughly 7,270 ECs per thread per iteration, with
   small per-EC cost. With two parallel regions per iteration, the run crosses
   about 3,700 fork/join barriers.
2. **The reduction is O(threads x transcripts) per iteration.** The thread-local
   buffer is `num_threads * n` doubles — 32 x 199,138 x 8 B is about 51 MB —
   which must be zeroed and summed every iteration. Over 1,863 iterations that
   is roughly 95 GB of memory traffic and 12 billion adds for the reduction
   alone, independent of how much real work the E-step does.

The dense thread-local buffer is the interesting target: each EC touches only a
handful of transcripts, so the overwhelming majority of those 51 MB are zeros
being written and re-read every iteration. A sparse or blocked accumulation
would cut that traffic sharply without changing the arithmetic, and it preserves
the deterministic fixed-order reduction.

Also note the finalizer's phase timing is measured on Bridges RM (AMD EPYC 7742,
Zen 2, 2.25 GHz base) while the serial baseline was measured on pikachu (Intel
i9-13900KF, up to 5.8 GHz). For a phase with a large serial and
memory-latency-bound component, that difference alone is a large factor, so
absolute cross-machine comparisons of the finalizer are not meaningful without
correcting for it.

### A parity trap for the standalone binary

STAR's in-process VB converged at **1,863 iterations**, while `em_quant`'s
defaults are `--max-iters 200` for VB with `--tolerance 1e-8`. A standalone
finalizer must adopt STAR's iteration limit and convergence criterion, not
`em_quant`'s defaults, or it will stop early and disagree for reasons that have
nothing to do with sharding. This belongs in the first parity check alongside
effective-length and GC handling.

## Parity milestone 1: the gather runs without the index

`tests/transcriptvb/sidecar_gather_parity.cpp`, built by
`make sidecar-gather-parity`, reads N sidecars, validates ordinals, and merges
them using only `TranscriptVBSidecar` and `libem`. `nm` reports zero
genome-index symbols in the resulting binary.

Run against the eight production sidecars from job `43119129`:

| Figure | STAR's `final/Log.out` | Standalone harness |
|---|---:|---:|
| Gathered pairs | 50,917,353 | **50,917,353** |
| Merged equivalence classes | 232,633 | **232,633** |
| GC fragment observations | 23,335,011 | **23,335,011** |

All three reproduce exactly, so the gather half of the finalizer is confirmed to
need no genome index.

Resource contrast for that half:

| | in-STAR finalizer | standalone gather |
|---|---:|---:|
| Peak RSS | 35.37 GB | **318 MB** |
| Wall | ~192 s for the whole phase | 0.83 s |
| Binary size | 10.8 MB | 125 KB |

Incidental observations worth keeping:

- Per-shard EC counts (136,487--137,111) sum to about 1.09 M but merge to
  232,633, so the merge genuinely combines equivalence classes across shards
  rather than concatenating them.
- `dropped_incompat` totals 749,346 fragments across the run.
- Every shard reports `first_pair=0`. The wrapper does not set
  `--quantVBSidecarFirstPair` per shard, so shard identity rests on the ordinal
  alone. Harmless today because ordinals are validated, but if `first_pair` is
  ever relied on for ordering it must actually be populated.

### Remaining parity work

The quantification half is untouched by this milestone. Next, in order:

1. effective lengths and GC bias reproduced from `Evidence::gc_counts` and
   `fld_state` — the most likely source of divergence;
2. VB/EM run with STAR's iteration limit and tolerance, not `em_quant`'s
   defaults, checked against the 1,863-iteration convergence;
3. `quant.sf` compared to `final/quant.sf` on the established Pearson bands;
4. `quant.genes.sf` via a transcript-to-gene map, including transcripts absent
   from the map.

## The two-pass gate: use non-zero support

The existing gate `applyGlobalPruning(ECTable&, double min_global_abundance)` in
`libem/extended_pruning.cpp` accumulates a per-transcript total and strips
anything below the threshold from every EC. It is deterministic — a serial loop
in fixed order with a plain threshold — but two properties make an abundance
threshold the wrong choice here:

- **It is off by design.** `enable_local_pruning` and `enable_global_pruning`
  both default `false`, commented "Default OFF for Salmon parity", and the
  function is only reachable from `ec_filter_test`, not from the production
  VB/EM path. Enabling it is a change against Salmon parity, not a performance
  switch.
- **It loses mass.** When every transcript in an EC is pruned the code sets
  `ec.count = 0.0`, discarding those reads. Where an EC keeps at least one
  transcript, mass is preserved and implicitly redistributed during EM, but
  whole-EC kills show up directly in total `NumReads`.

**The gate is therefore non-zero support: a transcript that appears in at least
one equivalence class.** This is exact rather than approximate. A transcript in
no EC contributes nothing to the E-step and receives nothing from it, so its
expected count is 0 and it is truncated to 0 on output regardless. No mass
moves, no EC is emptied, and the downstream dropping mechanics are untouched. It
is also trivially compatible with the deterministic fixed-thread-order
reduction, because it is only a compaction of the transcript index space.

### Measured on the merged production evidence

| | |
|---|---:|
| Transcripts total | 199,138 |
| With at least one EC | **129,516 (65.0%)** |
| With global weight > 0 | 129,476 |
| No support at all | **69,622 (35.0%)** |
| EC-transcript incidence | 1,043,422 |
| Mean transcripts per EC | 4.49 |

| Reduction buffer | dense | non-zero gated |
|---|---:|---:|
| 32 threads | 48.6 MiB | 31.6 MiB |
| 256 threads (8 x 32) | 388.9 MiB | 253.0 MiB |

Note the 40-transcript gap between "at least one EC" and "global weight > 0";
those sit in ECs contributing zero weight. Gating on EC membership is the more
inclusive and therefore safer of the two.

### Honest assessment of the benefit

**35%, which is worthwhile but not transformative.** Over 1,863 iterations it
removes roughly a third of the reduction traffic, but it does not change the
iteration count or the fork/join overhead, which are the other two components of
the finalizer's modest scaling.

It also does not address the *per-EC* sparsity, which is much larger: each EC
touches about 4.5 transcripts out of 129,516 supported, so even the gated buffer
is over 99.99% zeros for any single EC. Capturing that would need a genuinely
sparse per-thread accumulation rather than a dense buffer, and that changes the
reduction order — so it must be designed to preserve determinism explicitly, for
example by reducing over a sorted compacted index rather than in encounter
order. The non-zero gate needs no such care.

Recommended order: take the non-zero gate first since it is exact and free of
parity risk, and treat sparse accumulation and iteration count as separate
questions with their own evidence.

## Sparse-gate parity: results

`make sparse-gate-parity` runs the same VB three times over the same merged
sidecar evidence — dense, dense again, and gated to transcripts with non-zero EC
support — and compares them.

### A real bug found: initialization is not gate-invariant

`vb_engine.cpp` computed `numActive = state.n` for Salmon's
`uniformPrior = totalWeight / numActive`, with the comment "Use all transcripts
like Salmon". Compacting the index space changes `numActive` from 199,138 to
129,516, so the gated run starts from different initial abundances.

Effect, before the fix:

| | dense vs gated |
|---|---|
| iterations | 1360 vs **1361** |
| max relative delta | **1.0** (complete disagreement on some transcripts) |

Fixed by adding `EMParams::num_active_override`, default `0` meaning "use
`state.n`", which a gated caller sets to the original transcript count. After
the fix iterations agree and the delta falls to the 1e-10 range. This is
necessary for any sparse or two-pass scheme: without it the gate silently
changes the answer.

### The engine is not bit-reproducible against itself

This corrects an earlier claim in this document that the EM kernel is
deterministic given identical equivalence classes. It is not. The E-step uses
`schedule(dynamic)`, so thread-to-EC assignment varies between runs and the
per-thread partial sums differ. The fixed-order reduction preserves the order in
which partials are combined, but not the partials themselves.

Running dense twice, same binary, same input:

| dense vs dense | |
|---|---|
| iterations | 1360 vs 1360 |
| transcripts differing | 23,011 / 199,138 |
| max relative delta | 1e-10 to 1e-9, varying by run |

### Where the gate stands: exact within engine noise

The gate's *premise* is confirmed: across every run, **zero** unsupported
transcripts ever carried a non-zero dense count, and totals agree to about
1.6e-16 relative (23,351,021.000 both ways).

Because the engine is not bit-reproducible, the floor was characterised over
four dense runs, giving six pairwise comparisons:

```text
1.28e-10  2.97e-10  8.63e-10  9.91e-10  1.16e-09  1.29e-09
```

The gated-versus-dense delta was **6.10e-10**, which lies inside that
distribution rather than beyond it. Gating therefore perturbs the result no more
than running the engine twice does, which is the strongest statement available
given `schedule(dynamic)`.

| Axis | Result |
|---|---|
| 1. Gather | **exact** — pairs, EC count and GC observations all match STAR's log |
| 2. Sparse gate | **PASS** — 6.10e-10, inside the 1.28e-10 to 1.29e-09 engine floor |
| 3. VB/EM vs STAR `quant.sf` | **PASS** — Pearson 0.999999915 against the 0.99998 band |

### Axis 3 in detail, and what it does not yet prove

```text
transcripts differing 56,858 / 199,138
Pearson r             0.999999915069082
totals                23,351,021.000 vs 23,351,021.050
```

Pearson clears the chr22 TranscriptVB smoke's `0.99998` NumReads band by about
four orders of magnitude, and totals agree to 2e-9 relative.

**This validates the EM half, not the effective-length half.** The harness is
fed STAR's *final* effective lengths from `quant.sf` and holds them fixed, while
STAR updates them dynamically during VB — its log records a GC update at
iteration 11. So axis 3 shows that given the same effective lengths, the
standalone VB reproduces STAR's quantification within the established band. It
does not yet show that the effective lengths themselves can be reconstructed
from `Evidence::gc_counts` and `fld_state`.

That reconstruction is the remaining work, along with `quant.genes.sf`. Until it
is done the standalone finalizer is validated end to end only where its inputs
are borrowed from STAR.
### libem consumers checked

`run_vb`/`run_em` have two callers: `core/legacy/source/STAR.cpp:2911`, the
in-process finalizer, and `core/features/vbem/tools/em_quant/em_quant.cpp`.
`libem.a` is additionally linked by `compute_gc_bias` and `ec_filter_test`.
SLAM has its own separate `slam/source/libem`. The `num_active_override` field
is additive and defaults to preserving existing behaviour, so no caller changes;
STAR (1.6.1) and `em_quant` both rebuild cleanly.

## Axis 4/5: reproducing STAR's effective lengths and quantification

Earlier revisions of this section reported effective-length reconstruction as
unsolved and named `TranscriptQuantEC`'s `Parameters` dependency as a blocker.
**Both were wrong.** They came from approximating STAR's path instead of tracing
it. Tracing `STAR.cpp:2551-2915` found two concrete defects in the harness, and
with those fixed the reproduction is exact.

### Defect 1: the sidecar's EC weights are unnormalised

STAR writes the sidecar from `mergedEC.getECTable()` at line 2648, **before**
calling `mergedEC.finalize()` at line 2690. `TranscriptQuantEC.h` says so
explicitly: "deliberately before normalize/finalize and VB/EM". So every
consumer of a sidecar must normalise the EC weights itself before VB.

`finalize()` is pure and self-contained — per EC, rescale weights to sum to 1:

```cpp
for (EC& ec : ecTable_.ecs) {
    if (ec.weights.empty()) continue;
    double sum = 0.0;
    for (double w : ec.weights) if (std::isfinite(w) && w > 0.0) sum += w;
    if (sum > 0.0) { const double inv = 1.0/sum;
        for (double& w : ec.weights) w = (std::isfinite(w) && w > 0.0) ? w*inv : 0.0; }
}
```

Running VB on unnormalised weights conserves total mass but scrambles its
assignment, which is what produced the NumReads Pearson of 0.008.

### Defect 2: FASTA order is not the transcript index order

The transcript index space is `transcriptInfo.tab` order, which is also
`quant.sf` order and STAR's `transcriptomeMain->trID`. `transcriptome.fa` holds
the **same 199,138 transcripts in a different permutation** — verified by sorting
both name lists, which match as sets but not as sequences.

Indexing sequences by FASTA position therefore applies each transcript's GC
correction to the wrong transcript. `libem::Transcriptome::reorderByNames` exists
for exactly this and must be called with the index-order names.

### `TranscriptQuantEC` is not required

The retracted blocker claim was inferred from the constructor signature without
checking whether the class was needed. It is not:

- `finalize()` is the eight lines above, with no `Parameters` involvement;
- `importEvidence` only validates counts, restores `FLDAccumulator` from
  `Evidence::fld_state`, and restores GC from `Evidence::gc_counts` — all
  reachable directly from the merged sidecar.

A standalone finalizer needs neither the class nor STAR's `Parameters`.

### Result: exact reproduction outside STAR

With both defects fixed, running the traced sequence against the eight
production sidecars — no STAR, no genome index:

```text
GC applied at iteration 11, converged=yes, iterations=1863
effective lengths:  Pearson 0.999999999999959  mean abs 9.99e-05 nt
NumReads:           Pearson 1.000000000000007  mean abs 6.84e-05
                    totals 23,351,021.000 vs 23,351,021.050
```

The GC update fires at **iteration 11** and convergence is at **1863
iterations**, both matching STAR's log exactly. Agreement is at the 1e-13 level,
bounded by `quant.sf` printing three decimals.

### The traced sequence, for implementation

```text
1. transcript_vb_sidecar::read  x N
2. transcript_vb_sidecar::merge
3. normalise EC weights            <- finalize() equivalent; the sidecar is unnormalised
4. names + lengths in INDEX order  <- transcriptInfo.tab / quant.sf, not FASTA
5. libem::Transcriptome::loadFromFasta then reorderByNames(index_names)
6. FLDAccumulator::restore(fld_state) -> getPMF()
7. computeEffectiveLengthsFromPMFWrapper(pmf, raw_lengths)
8. install EMParams::effective_length_update calling
   computeDynamicGCBiasedEffectiveLengthsWrapper(txome, pmf, raw_lengths,
       alpha_counts, cb.eff_lengths, gc_counts) and swapping the result in
9. run_vb
10. write quant.sf / quant.genes.sf
```

Steps 1, 2, 5, 6, 7, 8 and 9 are existing library calls. Only step 3, the
eight-line normalisation, and step 10 need writing.


## Benchmark on Bridges-2 against the production artifacts

Job `43135160`, one exclusive RM node (`r270`, AMD EPYC 7742, 128 cores), against
the eight sidecars from run `43119129`. Packet `1002423268c0`, which still
carried the caller-side gate later superseded by the in-engine one.

| mode | threads | VB | total wall | peak RSS |
|---|---:|---:|---:|---:|
| dense | 8 | 249.8 s | 255.2 s | 1.88 GB |
| dense | 32 | 258.2 s | 261.6 s | 1.92 GB |
| dense | 64 | 241.8 s | 245.2 s | 1.97 GB |
| gated | 8 | 247.3 s | 250.8 s | 1.86 GB |
| gated | 32 | 240.4 s | 243.9 s | 1.88 GB |
| gated | 64 | 220.0 s | 223.5 s | 1.91 GB |

Parity on Bridges: `max abs NumReads delta 0` against STAR's `quant.sf` in both
modes, totals identical, GC at iteration 11 and 1863 iterations throughout.

### Memory: the result that matters

**Peak RSS is about 1.9 GB against the in-STAR finalizer's 35.37 GB, an 18x
reduction**, and it barely moves with thread count. That is the case for the
standalone binary: it removes the finalizer's dependence on a large-memory node
entirely, which is the constraint flagged for the JAX 16-plex Flex arm in
`docs/RUNBOOK_8NODE_ASSAY_SHARDING_TESTS_20260807.md`. A finalizer needing 2 GB
can run on `RM-shared` rather than holding a full node.

### RETRACTED: "VB does not scale with threads"

The thread sweep above ran every configuration **unbound** — a plain process on
a 128-core exclusive node with only `OMP_NUM_THREADS` set, so threads were free
to spread across both sockets of the EPYC 7742. Job `43142778` shows what that
costs:

| placement | VB | wall | %CPU |
|---|---:|---:|---:|
| unbound | 252.65 s | 258.49 s | **1121%** |
| `srun --cpus-per-task=32` | 184.95 s | 188.94 s | 567% |
| `OMP_PROC_BIND=close OMP_PLACES=cores` | 184.81 s | 188.90 s | 569% |

Binding cuts wall time by 27% and halves CPU consumption. The unbound run burns
about twice the CPU to do the same work more slowly, which is barrier spinning
across sockets, not useful parallelism.

So the sweep measured placement noise, not scaling. **The scaling curve is
unknown** and must be re-measured with binding before any statement about thread
behaviour is made. The related inference — that the workload is bound by the
serial dependency between iterations rather than available parallelism — was
drawn from those same unbound numbers and is equally unsupported.

### Standalone versus in-STAR timing: the gap was placement

With matched placement the standalone completes in **188.9 s** against STAR's
**188 s** for the quantification phase. The previously recorded ~30-39% deficit
was entirely an artifact of how the benchmark was launched.

**This is still not a clean comparison, and no precise claim should be made from
it.** Three things remain unmatched:

- **Cache state.** STAR's finalizer ran on a node where a worker had just
  finished, so its inputs were warm; every standalone arm here ran cold. The
  non-VB portion is only about 3.4 s of wall, so the effect is small, but the
  warm-condition timing of either arm has not been measured.
- **Output set.** STAR's 188 s includes writing `quant.genes.sf`; the standalone
  does not produce it at all.
- **Provenance.** The STAR figure is read from the logs of run `43119129`, hours
  earlier in a different job shape, not run fresh alongside.

A defensible comparison needs both arms fresh, back to back, in one job on one
node, at the same thread count and placement, with the output sets matched.
Until then the only supported statement is that placement, not algorithm,
accounted for the difference.

### Practical consequence

Any deployment of this binary must bind threads. Under Slurm, launching through
`srun --cpus-per-task=N` is sufficient; standalone invocations should set
`OMP_PROC_BIND=close` and `OMP_PLACES=cores`. Unbound on a large multi-socket
node it runs 27% slower while consuming twice the CPU.

## The VB decomposes: per-component solving beats distribution

Measured on the merged production evidence, the EC-transcript graph is not
coupled. Two equivalence classes interact only if they share a transcript, and
almost none do:

```text
components (with >=1 EC) = 23,995
largest component:        593 ECs (0.25% of 232,633), 0.06% of mass
largest by mass:          208 ECs, 0.37% of mass
outside the largest:      232,040 ECs (99.7%), 99.94% of mass
```

**There is no giant component.** Each of ~24k components is an independent VB
problem, so a distributed VB needs no per-iteration exchange at all: workers take
whole components and run to convergence locally.

What genuinely remains global is small: the initialization constants
(`totalWeight`, `fracObserved`, `uniformPrior`) from one pass over all ECs, the
dynamic GC update -- STAR's log records it using 90,070 background transcripts --
and the final TPM normalisation. One setup pass, one barrier, one gather, against
1863 iterations.

### Most components converge almost immediately

Running VB per component:

| iterations to converge | components |
|---|---:|
| <= 100 | 23,084 (96.2%) |
| <= 500 | 909 |
| <= 1000 | 2 |
| slowest | **611** |

| | EC-iterations |
|---|---:|
| per-component | 25,099,735 |
| global | 335,456,786 |
| **reduction** | **13.4x** |

Wall for all 23,995 components, single-threaded and sequential: **1.8 s**,
against roughly 56 s for the global VB on 16 bound threads.

Two effects compound. The global convergence check iterates everything until the
slowest transcript settles, so 96% of components keep iterating long after they
have converged. And a component spans tens of transcripts, so its working set
stays in cache, where the global run touches a 199,138-element dense buffer every
iteration.

### Consequence

The lever is not more workers and not a faster exchange. It is not iterating
sub-problems that have already converged, and that is a **single-node** change
which applies to STAR's in-process path as much as to the standalone binary.
Distribution becomes an option layered on top of an already much cheaper
computation, not the way to make it cheap.

### Not yet parity-validated

**Initialization is not matched.** Salmon derives `uniformPrior` from a global
`totalWeight / numActive`; running components in isolation changes both terms, so
the counts above are not comparable to `quant.sf`. The measurement establishes
the convergence structure and the work reduction, nothing about correctness.

Making it parity-correct means computing the initialization constants globally in
one pass and seeding each component with the values the global run would have
used -- the same class of fix as `num_active_override`, which was needed for
exactly this reason in the earlier gating work. That must be demonstrated against
`quant.sf` on the established Pearson bands before any of this is adopted.


## Per-component convergence: implemented, opt-in, and NOT equivalent

Implemented in `vb_engine` behind `EMParams::per_component_convergence`
(default `false`), exposed by the standalone as `--per-component`. The engine
builds the EC-transcript components once, then lets each stop as soon as it
meets `tolerance` instead of iterating everything until the slowest component in
the dataset converges.

Two properties make it safe to freeze a component mid-run: the E-step normalises
within each EC, so the global `logNorm` cancels and a frozen component cannot
perturb others; and `min_iters` (100) exceeds
`effective_length_update_target_iter` (10), so nothing freezes before the
effective-length callback has run. That last point resolves the GC rendezvous
ordering concern raised earlier -- it does not arise.

One implementation trap: `expected_counts` is `memset` to zero every iteration,
so a frozen component must be excluded from that zeroing or its alpha collapses
to zero rather than holding its fixed point.

### It buys speed by converging less

| tolerance | iterations | VB | NumReads Pearson vs STAR | max abs |
|---|---:|---:|---|---:|
| global (default) 0.01 | 1863 | 58.4 s | *(reference, byte-identical)* | 0 |
| per-component 0.01 | **614** | **39.8 s** | 0.999984785 | 26,782 |
| per-component 0.001 | 1589 | 46.8 s | 0.999991794 | 20,893 |
| per-component 0.0001 | 5759 | 74.0 s | 0.999999868 | 1,306 |

**This is a speed/accuracy tradeoff, not a free optimisation**, and earlier notes
in this document claiming a 13.4x work reduction were measuring the looser
stopping point rather than saved redundant work. Totals are conserved (23,351,021
either way) but individual transcripts move by up to 26,782 reads at the default
tolerance.

The reason is that neither rule reaches the fixed point; both stop at a
tolerance. The global rule keeps every component iterating until the slowest one
settles, which drives most components well past their own threshold. Per-component
honours the tolerance as written, which is arguably the more faithful reading of
the parameter, but it lands further from the fixed point than STAR does.

At the default tolerance it clears the chr22 smoke's `0.99998` NumReads band, but
only just, so it should not be enabled for work that is compared against existing
TranscriptVB results without re-establishing parity.

### Default behaviour is unchanged

With the flag off the standalone reproduces STAR byte-identically on both
`quant.sf` and `quant.genes.sf`, at the same 1863 iterations. STAR itself is
unaffected: the 100K SAM body hash is still
`27ccef4f72b243b00ce45d8e7f59adc9011b15fe1e27ebec3815832b861da87e` and the
sidecar unit tests pass.

## Component-partitioned E-step: exact, and the better lever

`EMParams::component_partition` (standalone: `--component-parallel`)
parallelises the E-step over connected components instead of over equivalence
classes. Because components share no transcripts, **a thread owning a component
is the only writer to those transcripts**, which removes both the
`num_threads x n_transcripts` accumulation buffer (51 MB at 32 threads) and the
reduction pass that zeroed and summed it every iteration.

Convergence stays global, so the answer is unchanged. Verified: `quant.sf` and
`quant.genes.sf` are **byte-identical to STAR's**, at the same 1863 iterations.

| threads | EC-parallel | component-parallel |
|---:|---:|---:|
| 4 | 63.07 s | 59.43 s |
| 8 | 58.04 s | **51.58 s** |
| 16 | 58.00 s | **46.47 s** |
| 32 | 55.56 s | 45.99 s |

At 16 threads that is 20% faster with an identical result, and
component-parallel on 8 threads beats EC-parallel on 32. EC-parallel is
essentially flat from 8 threads up, which is the reduction pass and the shared
buffer dominating; removing them restores some scaling.

Components vary enormously in cost, so the loop uses `schedule(dynamic)` -- each
thread takes one component, finishes it, and requests the next -- with
`comp_order` sorting components by descending EC count so the largest start
first. That is longest-processing-time-first, the standard heuristic when job
costs are uneven and known in advance.

### This supersedes `--per-component` for production

`--per-component` trades accuracy for speed by stopping components at their own
tolerance. `--component-parallel` gets a comparable saving from parallelism
while leaving the result bit-identical, so it is the right default lever.
`--per-component` remains useful only as an explicitly approximate fast mode.

Both are off by default, and STAR is unaffected: its 100K SAM body hash is still
`27ccef4f72b243b00ce45d8e7f59adc9011b15fe1e27ebec3815832b861da87e`, the sidecar
unit tests pass, and `em_quant` builds.

## Lustre benchmark: the component-parallel win is largely hardware-specific

Job `43148409`, one exclusive RM node (`r213`, EPYC 7742), against the eight
production sidecars on Lustre. Every arm bound via
`srun --exclusive --exact --cpus-per-task=N`, so this is not the unbound
measurement retracted earlier.

| threads | EC-parallel | component-parallel | gain |
|---:|---:|---:|---:|
| 8 | 185.49 s | 179.64 s | 3.2% |
| 16 | 174.90 s | 168.85 s | 3.5% |
| 32 | 171.85 s | 160.92 s | 6.4% |
| 64 | 174.75 s | 159.97 s | 8.5% |

Exactness holds on the production hardware: `quant.sf` and `quant.genes.sf` are
**byte-identical to STAR's** for both variants at 32 threads.

**On EPYC the component-parallel gain is 3--8%, against 20% measured on
pikachu.** The optimisation removes a 51 MB accumulation buffer and its
reduction pass, which is a memory-bandwidth effect; that mattered on a
single-die i9 with high per-core bandwidth and barely registers on a two-socket
Zen 2. Local numbers overstated it by roughly 3x, and the Lustre figures are the
ones to quote.

### Neither variant scales, and that is the real finding

From 8 to 64 threads, EC-parallel moves 185.5 s to 174.8 s and
component-parallel 179.6 s to 160.0 s -- roughly 6% and 11% for eight times the
cores. Removing the shared buffer and the reduction pass entirely did not change
the shape of the curve.

So the ceiling is not the parallel structure. What remains are the per-iteration
passes that stay O(n) over all 199,138 transcripts no matter how the E-step is
partitioned: the `alphaSum` reduction, the `expTheta` computation, the
`prev_alpha` copy, and the abundance normalisation. Across 1863 iterations that
is on the order of 1.5 billion transcript-touches that component partitioning
cannot remove.

Those are the same shape as the reduction already eliminated -- dense sweeps over
the full transcript space where only supported transcripts, or in the
partitioned case only a component's own transcripts, actually need touching.
That is where to look next, and until it is addressed adding workers will not
help the VB.

## Static-range architecture: measured on EPYC, and super-linear

Job `43150259`, exclusive EPYC 7742 node, single-threaded and bound, fixed at 300
iterations for every arm so this measures the architecture rather than a
different stopping rule. Components partitioned into W equal-cost ranges, each
range given the global effective lengths (so no GC rendezvous is needed) and
solved over **only its own transcripts**.

| | wall | vs global |
|---|---:|---:|
| global, all 199,138 transcripts | 6.57 s | 1.0x |
| 8 ranges, slowest range | 0.59 s | **11.1x** |
| 16 ranges, slowest range | 0.30 s | **21.9x** |
| total work across ranges | 4.70 s | **28% less work than global** |

Two things make this the result that justifies building the layout.

**It is super-linear** -- 11.1x on 8 ranges, 21.9x on 16 -- because splitting
does not merely divide the work, it *reduces* it. Total work falls from 6.57 s to
4.70 s. Each range's per-iteration passes sweep its own ~16k transcripts instead
of all 199,138, which is exactly the ceiling the Lustre benchmark identified and
that component partitioning could not break: removing the shared accumulation
buffer and its reduction left `alphaSum`, `expTheta`, the `prev_alpha` copy and
abundance normalisation still running over the full transcript space every
iteration.

**It replicates across hardware.** pikachu measured 11x and 22x; EPYC measures
11.1x and 21.9x. The component-parallel result did not replicate -- 20% locally
against 3--8% on EPYC -- because it was a memory-bandwidth effect. This one being
stable across very different machines is evidence it is structural.

Static assignment is also viable without any scheduler: 8 ranges balance to
within 0.5% of equal cost, 16 within 1.4%, and the largest single component is
0.54% of total work, so no unit can dominate a range.

### What to build

- component-ordered evidence file with a `(component -> byte offset, length)`
  index, so a worker reads one contiguous range;
- transcripts renumbered in component order, so each worker also *writes* one
  contiguous range of a pre-sized output file, with no shared blocks;
- two phases around the single GC rendezvous: iterations to the update point,
  one process computes the GC correction from the gathered alphas, then ranges
  resume to convergence;
- inverse permutation applied at assembly when writing `quant.sf`.

Two large reads and two large writes per worker for an entire run, with no
per-iteration coordination.
