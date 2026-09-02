# STAR Suite v1.8.1 Release Notes

Date: 2026-09-02

`v1.8.1` adds streaming cell-barcode bucketing with a bucket-parallel Solo
tail and an optional disk-spill mode that bounds per-read memory for supported
fully fused STAR-Flex runs. The post-map counting phase, previously
single-threaded, now runs bucket-parallel with byte-identical outputs, and
datasets whose record accumulation exceeds host memory complete under a
configurable budget.

The release artifact version is `v1.8.1`, Debian packages use `1.8.1-1`, and
`STAR --version` reports `1.8.1`. The upstream STAR base remains `2.7.11b`;
genome-index compatibility remains `2.7.4a`. Existing indexes do not need to
be rebuilt.

## Streaming CB buckets and bucket-parallel collapse

- Route each kept Flex record into a cell-barcode-hash bucket as it is
  produced during streaming. Producers append to per-worker segments; sealed
  segments are handed to the bucket store under an atomic claim, mirroring the
  BGZF reader's frontier-claim pattern.
- Replace the monolithic post-map gather and global barcode grouping with tail
  workers that claim whole buckets and run the existing per-barcode collapse
  in parallel. Bucket barcode ranges are disjoint, so final outputs are a
  concatenation; no merge step exists.
- Outputs are byte-identical to the prior path across thread counts, bucket
  counts, and storage modes, enforced by the `tests/bucket` equality gates.
- Measured on the full JAX Flex benchmark (2.011B pairs, 32 threads, BGZF
  input, cold cache): post-map tail 210.5 s to 98.2 s; total wall 7:41 to
  5:40; identical per-sample outputs.

## Disk spill mode

- In spill mode, sealed bucket segments append to one file per bucket under a
  scratch directory; nothing is read back until the tail, which overlaps
  read-back with per-bucket compute. Spill files are removed on completion.
- This bounds the packed-record accumulation, which otherwise grows with kept
  reads (~30 bytes per read). The public 320k scFFPE GEM-X Flex benchmark
  (7.3B read pairs) completes on a 126 GB host in 35:25 at a 97.6 GiB peak
  with all 16 per-sample outputs produced.

## Controls and compatibility

- `--soloBucketMode auto|ram|spill|off` selects storage: `auto` (default)
  begins in RAM and switches to spill when the budget is crossed; `ram` and
  `spill` force a backend; `off` restores the previous monolithic tail.
- `--soloBucketMemGB` sets the packed-record budget for `auto` (default 32).
- `--soloBucketSpillDir` sets the spill location (defaults to the run's
  temporary directory); `--soloBucketCount` sets the bucket count (default 256).
- The bucketed tail is restricted to supported fully fused STAR-Flex
  configurations. Classic STARsolo paths, plain gzip handling, and all
  non-Flex behavior are unchanged.

## Benchmark tooling

- `docs/benchmarks/cloud_320k/run_cloud_320k_matrix.sh` reproduces the 320k
  public benchmark matrix on a fresh AWS r6id.8xlarge from the released binary
  and 10x's public data; `docs/RUNBOOK_CLOUD_320K_SPILL_20260902.md` is its
  narrative companion.
- Benchmark preflight snapshots now record argument-free process columns.
