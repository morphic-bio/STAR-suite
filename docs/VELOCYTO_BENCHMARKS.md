# Velocyto Benchmarks

This document preserves the detailed Velocyto bridge benchmark notes that were
previously embedded in the top-level `README.md`.

These runs isolate the exact/deterministic Velocyto bridge on corrected UCSF
`EBs2_2/GEX`, without perturb feature-assignment confounders.

Primary wrapper and implementation references:

- [publications/benchmarks/README.md](../publications/benchmarks/README.md)
- [docs/RUNBOOK_VELOCYTO_BRIDGE_IMPLEMENTATION_20260326.md](RUNBOOK_VELOCYTO_BRIDGE_IMPLEMENTATION_20260326.md)
- [docs/RUNBOOK_VELOCYTO_COUNT_RESOLUTION_20260326.md](RUNBOOK_VELOCYTO_COUNT_RESOLUTION_20260326.md)

## Correctness / parity ladder

| Surface | Threads | Checks | Result |
|---|---:|---|---|
| 100K downsampled GEX-only fixture | 1 vs 8 | frozen baseline vs `stream_t1`, `stream_t1` vs `det_t1`, `det_t1` vs `det_t8`, `Gene`/`GeneFull` parity | **PASS** |
| 10M downsampled GEX-only fixture | 8 | frozen baseline vs `stream_t8`, `stream_t8` vs `det_t8`, `det_t8` vs `hash_t8`, `Gene`/`GeneFull` parity | **PASS** |
| Full corrected `EBs2_2/GEX` | 32 | `stream_t32` vs integrated-hash `hash_t32`, full Velocyto packaged/raw parity, `Gene`/`GeneFull` parity | **PASS** |

## Timing / RSS

| Surface | Arm | Reads | Wall | Speed | Peak VmRSS |
|---|---|---:|---:|---:|---:|
| 10M GEX-only | `stream_t8` | 10.0M | 1m 53s | 620.69 M reads/hour | 11,183,816 kB |
| 10M GEX-only | `det_t8` | 10.0M | 1m 53s | 620.69 M reads/hour | 11,183,816 kB |
| 10M GEX-only | `hash_t8` | 10.0M | 1m 54s | 620.69 M reads/hour | 11,109,796 kB |
| Full corrected `EBs2_2/GEX` | `stream_t32` | 444.9M | 23m 05s | 1208.78 M reads/hour | 47,608,160 kB |
| Full corrected `EBs2_2/GEX` | `hash_t32` | 444.9M | 23m 03s | 1210.60 M reads/hour | 37,651,052 kB |

- On the full corrected GEX surface, the integrated-hash Velocyto path is a
  near timing tie with the stream path and reduces whole-run peak VmRSS by
  about 21% (`47.6 GB -> 37.7 GB`) in this run.
- This `stream_t32` vs `hash_t32` pair is the measured apples-to-apples
  benchmark for the full `Gene + GeneFull + Velocyto` surface. `Gene` and
  `GeneFull` stay on the normal path while `Velocyto` switches to the
  integrated-hash implementation.
- The current RSS helper uses whole-`Log.out` `VmRSS` fallback on these runs,
  so the numbers above are whole-run peak RSS rather than a Velocyto-only phase
  metric.

## Fused `GeneFull + Velocyto` experimental mixed path

| Surface | Arm | Reads | Wall | Speed | Peak VmRSS | Status |
|---|---|---:|---:|---:|---:|---|
| Full corrected `EBs2_2/GEX` | `GeneFull` bridge + Velocyto hash (32 threads) | 444.9M | 12m 55s | 2227.58 M reads/hour | 21,532,344 kB | Velocyto exact, `GeneFull` not yet exact |

- Relative to the earlier full `Gene + GeneFull + Velocyto` hash control
  (`23m 03s`, `37,651,052 kB`), the mixed path is materially faster and lower
  memory.
- We did not benchmark a fully hash-based `Gene + GeneFull + Velocyto`
  implementation, and we are not currently pursuing the engineering needed to
  make `Gene` hash-based. The practical fast Velocyto-capable path today is to
  omit `Gene` and run `GeneFull + Velocyto` on the bridge/hash path.
- `Velocyto` raw outputs matched exactly on the full set.
- `GeneFull` had small sparse deltas. One representative trace was
  `TTTCACAGTTGGATCT / corr=7652715 / gene=4324`, an ambiguous-CB read with
  four whitelist candidates that appeared in both legacy and bridge paths; the
  bridge did not retain that singleton under the target barcode.

## Bottom line

- Full corrected `EBs2_2/GEX`, `Gene + GeneFull + Velocyto`, normal stream:
  **23m 05s**
- Full corrected `EBs2_2/GEX`, `Gene + GeneFull + Velocyto`, hashed `Velocyto`
  only: **23m 03s**
- Full corrected `EBs2_2/GEX`, `GeneFull + Velocyto`, bridge/hash path:
  **12m 55s**

The integrated-hash Velocyto path is roughly timing-neutral on the full
`Gene + GeneFull + Velocyto` surface. The practical fast path for Velocyto
output today is to skip `Gene` and run `GeneFull + Velocyto` on the bridge/hash
implementation.
