# Runbook: concurrent R1 and GeneFull sidecar producers for Visium HD GEX

Date: 2026-07-22

Branch: `feature/visium-hd-gex-concurrent-producers-20260722`

Base: STAR Suite `master` at `2fc05770d1b0528abf3709da66199b708ca17afb`

## Objective

Reduce the elapsed time before spatial feature joining by running the independent
raw-R1 coordinate decoder and STAR GeneFull sidecar producer concurrently. Keep
the existing source-only, global-ordinal, read-name-digest, coordinate, and
completion contracts unchanged.

This changes scheduling only. It must not change candidate sets, GeneFull
evidence, UMI policy, molecule resolution, or materialized matrices.

## Execution graph

```text
clean builds and source checks
             |
             +-- R1 decode -- H0 prior --+
             |                           |
             +-- STAR GeneFull sidecar --+-- validated ordinal join
                                             |
                                             +-- molecule resolver
                                             +-- materializer
```

The join remains a barrier. It may start only after the R1 candidate/H0 branch
and the finalized STAR sidecar branch both exit successfully. No process reads a
partial candidate table or incomplete sidecar.

## Existing contracts that make concurrency safe

- The fixture checksum manifest is validated before either producer starts.
- Both producers receive the same frozen `L001` through `L004` order.
- The R1 decoder writes a global ordinal for every candidate-bearing read.
- STAR writes one fixed-width sidecar slot per global ordinal.
- STAR finalization validates every ordinal, writes a complete header, fsyncs,
  and atomically renames the binary.
- STAR writes independent lane/block read-name digests.
- The join rescans the paired FASTQs and rejects any ordinal, read count, lane
  boundary, read-name digest, candidate-count, or coordinate-contract mismatch.
- Raw UMI and coordinate evidence remain R1-owned; GeneFull evidence remains
  STAR-owned. No GX, GN, UR, UB, CB, CR, SAM, or BAM carrier is introduced.

Execution timing therefore cannot change the join key or evidence ownership.

## Scheduler and failure contract

The wrapper uses a single-threaded `Popen` scheduler rather than concurrent
Python mutations of the command manifest:

1. Record and start `r1_decode`.
2. Record and start `star_sidecar`.
3. Poll both process groups.
4. Start `h0_prior` only after successful R1 completion; STAR may still run.
5. On any launch or non-zero-exit failure, terminate every live sibling process
   group, record final exit states, withhold the join, and leave no completion
   marker.
6. Start the join only after all three producer-stage records are complete.

`commands.json` is updated atomically by the scheduler and records each stage's
start/finish timestamps, state, exit code, and wall time. The wrapper-level `RUN_COMPLETE.json` records
producer mode, total thread budget, branch budgets, producer critical-path wall
time, and the successful join barrier.

## Thread-budget contract

`--threads` is the total concurrent producer cap and remains the build/downstream
thread setting.

- `--producer-mode concurrent` is the default.
- With no overrides, `--threads 16` assigns 8 R1 threads and 8 STAR threads.
- `--r1-threads` and `--star-threads` may override that split, but their sum may
  not exceed `--threads` in concurrent mode.
- `--producer-mode serial` is the parity/control path. Each producer defaults to
  the full `--threads` value because they do not overlap.
- Concurrent mode requires at least two total threads.

This prevents the former full-thread R1 process and full-thread STAR process from
silently oversubscribing the host when overlapped.

## Frozen 100K sources

Use only these source inputs:

```text
/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1
/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1/fastqs
/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1/checksums.sha256
/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1/summary.json
/mnt/pikachu/star-spatial/references/refdata-gex-GRCh38-2024-A_STAR-2.7.11a/star
/storage/star-spatial/runs/cleanroom_hd_mouse_brain/slide_oligos/bc1_full_oligos.txt
/storage/star-spatial/runs/cleanroom_hd_mouse_brain/slide_oligos/bc2_full_oligos.txt
/storage/star-spatial/runs/cleanroom_hd_mouse_brain/barcode_contract/
```

The barcode contract is an allowed declared source. Do not use any prior output
under `/mnt/pikachu/star-spatial/runs/`, any older sidecar/candidate ledger, or
any Space Ranger artifact as computational input.

## Fixture-free scheduler test

```bash
python3 tests/test_visium_hd_gex_sidecar_concurrency.py
```

The test proves:

- the R1 and STAR processes overlap;
- H0 starts only after R1 completion;
- serial mode preserves the historical order;
- thread-budget violations fail before launch;
- a producer failure terminates its sibling and suppresses H0/join work.

## Fresh serial control

```bash
python3 scripts/run_visium_hd_gex_sidecar_100k.py \
  --producer-mode serial \
  --threads 16 \
  --out-dir /mnt/pikachu/star-spatial/gex_sidecar_tests/20260722_ovarian_100k_serial_concurrency_control_v1
```

## Fresh concurrent run

```bash
python3 scripts/run_visium_hd_gex_sidecar_100k.py \
  --producer-mode concurrent \
  --threads 16 \
  --r1-threads 8 \
  --star-threads 8 \
  --out-dir /mnt/pikachu/star-spatial/gex_sidecar_tests/20260722_ovarian_100k_concurrent_v1
```

Both output roots must be absent or empty at launch. The wrapper performs a clean
STAR rebuild for each run and refuses dirty STAR Suite or companion worktrees.

## Parity comparison

Compare SHA-256 and byte size for all scheduling-independent products:

- `decoder/raw_r1_candidates.tsv`
- `decoder/h0_read_prior.tsv`
- `star/gex_features.bin`
- `star/gex_features.features.tsv`
- `star/gex_features.read_name_digests.tsv`
- `join/normalized_evidence.tsv`
- `resolver_a/{strict,hard,gated_hard}_molecules.tsv`
- every `materialized/*/*/{matrix.mtx,barcodes.tsv,features.tsv}`

Do not require equality for logs, `commands.json`, summaries containing elapsed
time/mode, or source paths.

## Acceptance gates

- Both fresh runs write `RUN_COMPLETE.json` with `status: complete`.
- Both reconcile exactly 100,000 paired reads and 100,000 sidecar records.
- Candidate, emitted-row, product, and scale counts match.
- Every scheduling-independent product listed above is byte-identical.
- Concurrent command intervals prove R1/STAR overlap and H0-after-R1 causality.
- Concurrent R1 + STAR threads do not exceed the declared total.
- The join reports matching name digests, complete candidate sets, and a complete
  coordinate contract.
- The concurrent producer wall time is reported against the fresh serial control;
  no speedup is assumed if shared FASTQ/index I/O becomes the limiting resource.

## Validation result

Pending implementation and fresh source-only 100K serial/concurrent runs.
