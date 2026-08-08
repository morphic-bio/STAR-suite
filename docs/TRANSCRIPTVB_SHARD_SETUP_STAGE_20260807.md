# TranscriptVB shard setup stage: alpha seeding and effective-length correction

Date: 2026-08-07
Scope: `core/features/vbem/source/libem/`, `core/features/vbem/tools/transcriptvb_finalize/`

## Problem

The Lustre-friendly finalizer architecture splits the VB/EM into W independent
workers, each reading one self-contained `.tvbr` shard and writing its own
output file. For that to be correct, a worker must start from exactly the state
the in-process global run would have been in. Two global quantities stand in the
way:

1. **Salmon's alpha initialization.** `alpha[i] = projectedCounts[i]*fracObserved
   + uniformPrior*(1-fracObserved)`, where `fracObserved` depends on the *total*
   projected weight and `uniformPrior` is `totalWeight/numActive` over *all*
   transcripts. A shard holding 1/8 of the transcripts cannot derive either.
2. **The dynamic GC effective-length correction.** The engine applies it exactly
   once, at `iter > effective_length_update_target_iter` (default 10), from the
   alpha of that moment.

## What was wrong

The first cut of `--emit-ranges` produced `alpha0` by calling `run_vb` with
`min_iters = max_iters = 0` and reading `state.counts`. That was an unverified
assumption about engine behaviour, and it was false: `run_vb` writes results
into `EMResult`, not back into `state`, so **every shard emitted before this
change carried an all-zero alpha vector** (`sum(alpha0) = 0.000` vs the correct
`18999823.478`). The shard runner also read `alpha0` from the file and then
ignored it, letting `run_vb` recompute an initialization from the shard's own
slice — local `totalWeight`, local `numActive`. Both defects were silent.

## Fix

### Shared initialization

`compute_initial_alpha(const ECTable&, const TranscriptState&, const EMParams&)`
is now exported from `vb_engine.h`. `run_vb` calls it; the setup stage calls the
same function. There is one implementation, so a splitter cannot drift from the
engine.

`EMParams::initial_alpha` carries a pre-seeded vector. When its size matches
`state.n`, `compute_initial_alpha` returns it unchanged. The shard runner sets
it from the `.tvbr` file, so a worker uses the global seed rather than
recomputing a local one.

### Setup burn-in removes the GC rendezvous

The GC correction's only cross-transcript quantity is a 101-bin histogram
(~808 bytes) — but collecting it mid-run would force every worker through a
barrier at iteration 11, which is exactly the crosstalk this architecture exists
to avoid.

Instead the setup stage runs the burn-in once, globally, and lets the engine's
own callback fire (the firing rule is not restated in the splitter, so setup
cannot drift). It captures:

- `alpha` at the fire point — the global run's iteration-11 alpha;
- the GC-corrected effective lengths.

Both go into every shard. A worker's iteration 0 is therefore the global run's
iteration 11, and workers need no rendezvous at all. Cost: 12 of ~1870
iterations, paid once.

### The GC background pass was serial

That pass walks every fragment start of every transcript to accumulate the
101 expected-GC bins, and it dominated setup: **57.5 s of a 57.6 s split**. It is
a pure reduction and is now parallelised with per-thread bins combined in fixed
thread order.

`computeDynamicGCBiasedEffectiveLengthsWrapper` takes a new trailing
`background_threads` parameter, **defaulting to 1**. The default is deliberate:
serial summation order is what makes `quant.sf` byte-identical to an in-process
STAR run, and any parallel reduction shifts a few values by one unit in the last
printed digit. Existing call sites are unchanged and keep byte-identity.

- The shard setup stage passes `threads` — it already has its own parity band.
- `transcriptvb_finalize --parallel-gc` opts the in-process path in.

## Measurements

All numbers below are **pikachu (128 cores), local filesystem**, full
SRR4422207 sidecars: 199,138 transcripts, 232,633 ECs, 23,412,000 processed
fragments, 116.7 MB of `.stvb` input. They are not Lustre numbers and not
multi-node numbers -- the 8 "workers" are 8 processes on one node. The
architecture exists for Lustre, so these measure the algorithm, not the
deployment.

Setup burn-in + emit, W=8:

| threads | burn-in | emit total |
|---|---|---|
| 1 | 57.49 s | 57.63 s |
| 8 | 7.88 s | 8.03 s |
| 32 | 4.49 s | 4.63 s |

Thread-count sensitivity of the emitted shards (vs `threads=1`): max absolute
alpha delta 4.07e-10 (2.5e-12 relative), max effective-length delta 3.2e-06.

In-process finalizer, 32 threads:

| mode | VB wall | `quant.sf` |
|---|---|---|
| default (serial GC) | 55.40 s | IDENTICAL to STAR |
| `--parallel-gc` | 25.46 s | 19 of 199,139 lines differ by 1 ulp |

`quant.genes.sf` under `--parallel-gc`: 15 of 36,602 lines differ by 1 ulp.

End-to-end sharded pipeline. The honest comparison holds core count equal on
both sides -- 8 cores each:

| | wall |
|---|---|
| in-process, 8 threads, `--parallel-gc` | 31.08 s |
| sharded: setup 8 threads (9.86 s) + 8 single-thread workers (0.87 s) | **10.74 s** |

2.9x at equal cores. A 32-thread setup plus 8 single-thread workers reaches
7.62 s (setup 6.74 s + workers 0.88 s), but that consumes more cores than the
in-process run it would be compared against, so it is not a speedup figure.

Neither figure includes the gather: `tvb_assemble` does not exist yet, and the
final `quant.sf` / `quant.genes.sf` write is not measured here.

## Parity

Gather of the 8 shard outputs against STAR's own in-process `quant.sf`:

| quantity | Pearson | gate | result |
|---|---|---|---|
| NumReads | 0.9999990344 | 0.99998 | PASS |
| EffLength | 1.0000000000 | 0.99995 | PASS |

129,516 supported transcripts compared. The 69,622 transcripts absent from the
shards are exactly those STAR reports with `NumReads == 0` — the non-zero-support
gate is exact, not approximate.

### The residual is stopping-rule slack, not a decomposition error

Shards converge each component against its own criterion (401–827 iterations
here); the global run holds *every* component until the global maximum change
falls below tolerance (1,268 iterations with GC off). Driving both to
`--tolerance 1e-10` collapses the difference:

| | max abs NumReads delta | Pearson |
|---|---|---|
| default tolerance | 1046.66 | 0.9999999155 |
| `--tolerance 1e-10` | 0.000499998 | 1.000000000000 |

Six orders of magnitude. The component decomposition itself is exact; what
remains is only how far each side is driven toward the same fixed point.

## Regression status

`STAR`, `sidecar-gather-parity`, `sparse-gate-parity`, and
`transcriptvb-finalize` all build. `sidecar_gather_parity` passes (`PASS ordinals
0..7 each present exactly once`; 232,633 ECs; 129,516 transcripts with support).

## Not yet done

- `tvb_assemble`: inverse permutation from shard outputs to `quant.sf` /
  `quant.genes.sf`. The gather is currently done in the test harness.
- Setup is now 92% of sharded wall time at equal cores (9.86 s of 10.74 s). It
  is the next target, not the workers.
- Nothing here has been measured on Lustre or across nodes. The one Bridges data
  point that exists (job 43152701) predates this commit and used zero-seeded
  shards, so its values are void; only its timings stand.
- Six of the eight full-width `state.n` passes per VB iteration are still
  ungated; only the zeroing and reduction passes use `supported_idx`.
