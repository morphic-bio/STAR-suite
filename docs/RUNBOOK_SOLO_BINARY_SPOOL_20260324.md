# Solo Binary Spool Runbook

Date: 2026-03-24
Branch: `feature/solo-binary-spool-20260324`
Base: `fa061d0`

## Goal

Evaluate lower-risk alternatives to the non-Flex Solo bridge by keeping legacy
Solo counting semantics and only replacing the temp-spool representation used
between mapping and `countCBgeneUMI()`.

## Context

The more aggressive non-Flex inline-hash bridge on the main Solo optimization
branch removed temp files but introduced:
- larger implementation risk
- parity drift work
- full-scale memory failures before later fixes
- no clear end-to-end speed win on full UCSF after the memory fixes

This branch tests a conservative alternative:
- keep legacy Solo replay/counting
- replace text temp records with a lossless binary spool
- then test whether moving that binary spool in memory helps further

## Designs Tried

### 1. Legacy text spool

Characteristics:
- existing STARsolo path
- writes formatted text temp records during mapping
- reparses them later during Solo counting

Pros:
- known-correct

Cons:
- formatting/parsing overhead
- larger temp footprint

### 2. Binary spool on disk

Implementation:
- opt-in via `STAR_SOLO_BINARY_SPOOL=1`
- same logical record content as text
- same replay boundary
- only supported for gene-like Solo features
- restart intentionally blocked

Key files:
- [SoloBinarySpool.h](/mnt/pikachu/STAR-suite-binary-spool/core/legacy/source/SoloBinarySpool.h)
- [SoloReadFeature.cpp](/mnt/pikachu/STAR-suite-binary-spool/core/legacy/source/SoloReadFeature.cpp)
- [SoloReadFeature_record_base.cpp](/mnt/pikachu/STAR-suite-binary-spool/core/legacy/source/SoloReadFeature_record_base.cpp)
- [soloInputFeatureUMI.cpp](/mnt/pikachu/STAR-suite-binary-spool/core/legacy/source/soloInputFeatureUMI.cpp)
- [SoloReadInfoLoader.cpp](/mnt/pikachu/STAR-suite-binary-spool/flex/source/SoloReadInfoLoader.cpp)

Result on corrected UCSF `iPSC2_1/GEX` 2M:
- exact top-line parity
- fair comparison vs warmed baseline: about `1.21 s` faster
- RSS effectively unchanged

Interpretation:
- correct and low-risk
- small but real gain on this host

Artifacts:
- warmed baseline:
  - `/storage/100K/ucsf_solo_binary_spool_20260324/iPSC2_1_GEX_2M_baseline_warm/`
- disk binary spool:
  - `/storage/100K/ucsf_solo_binary_spool_20260324/iPSC2_1_GEX_2M_binaryspool/`

### 3. In-memory binary spool, monolithic growable buffer

Implementation:
- opt-in via `STAR_SOLO_BINARY_SPOOL_IN_MEMORY=1`
- same binary record format
- per-thread `std::vector<char>` backing store
- cleared after replay

Result on corrected UCSF `iPSC2_1/GEX` 2M:
- exact top-line parity
- `1:34.11` wall
- slower than disk-binary `1:33.65`
- slightly higher RSS

Interpretation:
- no evidence that naive in-memory storage helps on this host
- likely hurt by monolithic vector growth and reallocation/copy overhead
- disk-binary path already benefits from page cache, so this was not a true
  “disk vs RAM” test

Artifact:
- `/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_total/`
- note: current wrapper ignored the custom output root and reused this path

### 4. In-memory binary spool, chunked slab-backed pool

Current implementation direction:
- keep the same binary record format
- replace the monolithic growable vector with fixed-size append-only slabs
- preserve sequential replay with an absolute read cursor

Why:
- avoid repeated reallocation and whole-buffer copy
- better model a real RAM-backed spool
- keep semantics identical so any timing difference is attributable to storage
  shape, not algorithm changes

Current chunk shape:
- fixed `4 MiB` slabs
- append-only writes
- sequential replay across slabs

Key idea:
- this remains a spool optimization, not a replay-elimination strategy

Result on corrected UCSF `iPSC2_1/GEX` 2M:
- exact top-line parity
- `1:34.00` wall
- slightly better than naive in-memory `1:34.11`
- still slower than disk-binary `1:33.65`
- RSS `40366192 kB`, lower than naive in-memory but still above disk-binary

Interpretation:
- the chunked allocator removed the obvious monolithic-vector penalty
- but the host still favors page-cache-backed binary spool over RAM spool
- the remaining cost is replay/serialization work, not allocator shape alone

Artifact:
- `/storage/100K/ucsf_solo_binary_spool_20260324/iPSC2_1_GEX_2M_binaryspool_memory_chunked/`

### 5. Hybrid in-memory spool with disk spillover

Implementation:
- still enabled with `STAR_SOLO_BINARY_SPOOL_IN_MEMORY=1`
- new optional limit:
  - `STAR_SOLO_BINARY_SPOOL_IN_MEMORY_LIMIT_MB=<MiB>`
- starts in RAM using the chunked slab-backed buffer
- once a thread exceeds the configured limit, it flushes the existing binary
  buffer to the normal binary spool file and continues on disk
- legacy text behavior remains the default when neither binary env var is set

Why:
- keep the full-speed RAM-first path for moderate runs
- bound memory on larger runs
- preserve the same binary record format and same replay boundary

Validation on corrected UCSF `iPSC2_1/GEX` 2M:
- `32 MiB` cap:
  - `/storage/100K/ucsf_solo_binary_spool_20260324/iPSC2_1_GEX_2M_binaryspool_memory_spill_v1/`
  - exact top-line parity
  - no explicit spill event observed in the surviving log, so this run mainly
    validated the capped configuration
- `1 MiB` forced-spill run:
  - first attempt `v1` exposed a spill-boundary bug:
    - `/storage/100K/ucsf_solo_binary_spool_20260324/iPSC2_1_GEX_2M_binaryspool_memory_spill_1mb_v1/`
    - records at the first spill point were written back into the cleared
      in-memory buffer because `outputReadCB_base()` cached the pre-spill
      `binarySpoolInMemory` state
  - fixed run:
    - `/storage/100K/ucsf_solo_binary_spool_20260324/iPSC2_1_GEX_2M_binaryspool_memory_spill_1mb_v3/`
    - explicit spill events logged for all threads
    - exact summary parity vs the on-disk binary spool
    - raw MEX files also match exactly (`matrix.mtx`, `barcodes.tsv`,
      `features.tsv`)
    - wall `1:32.65`
    - RSS `40308360 kB`

Interpretation:
- the hybrid RAM-first spool now works as intended
- the spill transition is real and parity-safe at the summary level on the
  corrected 2M fixture
- performance is effectively the same as the pure binary spool path on this host

### 6. Full corrected UCSF run with chunked in-memory spool

Result on corrected UCSF `EBs2_2` full sample:
- artifact:
  - `/storage/paper_bench_solo_binary_spool_20260324/ucsf_ebs2_2_chunkedmem/`
- completed successfully
- wall `18:35.81`
- peak RSS `82564728 kB`
- filtered cells `13721`

Interpretation:
- this conservative spool optimization produced a meaningful end-to-end wall
  win versus the matched corrected baseline (`19:48.45`) on the full workload
- the tradeoff is higher memory
- a spill-capable RAM-first spool is therefore the practical next shape for
  this branch

## What We Learned So Far

- Replacing text with binary is low-risk and preserves parity.
- On this host, page-cache-backed binary temp files are already efficient.
- A naive in-memory buffer is not enough; allocator behavior can erase the
  expected gain.
- If in-memory spooling is worth testing further, the fair version is a
  chunked slab allocator, not a single growable byte vector.
- The larger potential gains still come from removing replay/collapse work
  entirely, but that is a different risk class than this branch.

## Benchmark Policy

- Keep benchmarks serialized; do not run benchmark jobs in parallel.
- Use the corrected UCSF inputs under:
  - `/mnt/pikachu/ucsf-perturb-seq-corrected/`
- For perturb-seq parity/full benchmarks include:
  - `--soloCrMultimapRescue yes`

## Next Steps

1. Clean rebuild after chunked-slab implementation.
2. Re-run the corrected UCSF 2M benchmark for:
   - warmed baseline
   - disk binary spool
   - chunked in-memory binary spool
3. Compare full corrected UCSF baseline vs:
   - disk binary spool
   - chunked in-memory spool
   - chunked in-memory spool with a configured spill threshold
4. If RAM-first spool still needs too much memory, keep the disk-binary mode as
   the safe production path and treat RAM-first as optional.
