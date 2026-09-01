# Runbook: BGZF Parallel FASTQ Ingest (Flex)

Date: 2026-09-01

Status: planned. This runbook specifies the module for an implementing agent.
Tests are defined first (Phase 0) and gate every later phase.

Related docs:

- `docs/RUNBOOK_PF_CBQ_RANGE_COUNTING_20260531.md` (the indexed range-reader pattern this mirrors)
- `docs/RUNBOOK_CBQ_ORDERED_ENCODER.md`
- `docs/DYNAMIC_THREADS_TINY_FIXTURE_RUNBOOK_20260217.md` (fixture + smoke conventions)

## Goal

Parallel ingestion of BGZF-compressed FASTQ for the Flex pipeline, mirroring the
CBQ indexed range-reader design:

```text
BGZF block index -> per-worker block ranges -> raw-inflate + record parse -> fused Flex consumers
```

Illumina/BCL-Convert-delivered `fastq.gz` files are BGZF (blocked gzip: many
independent gzip members, each carrying its compressed size in a `BC` extra
subfield). They have been read through the single-threaded zlib stream like any
other gzip file, which makes ingest the bottleneck for Flex no-align runs.
BGZF blocks decompress independently, so ingest parallelizes the same way CBQ
ranges do — with **zero conversion** and files that remain standard `fastq.gz`
to every other tool.

Verified BGZF as delivered (header sniff): `/mnt/pikachu/JAX_sequences/JAX_scRNAseq01`
(all 8 lanes), 10x portal downloads (A375 fixture, 320k GEM-X set). Counter-examples
that must keep working through the existing zlib path: ENA/SRA mirrors and
tool-recompressed files (plain gzip).

Reference wall-clocks, JAX SC2300771 no-align, 32 threads (Supp Table S5 basis):

| Input | Wall time |
| --- | ---: |
| gzip FASTQ (current path) | 10m18s |
| CBQ | 7m22s |
| BGZF FASTQ (this module, target) | ~CBQ-class |

## Constraints

1. **No htslib changes.** The vendored htslib is not modified, not upgraded, and
   not used by this module. BGZF member parsing is a few hundred lines over the
   already-vendored zlib (`inflateInit2(windowBits=-15)` raw inflate per block).
   The BAM output path is untouched.
2. **Plain gzip behavior unchanged.** Non-BGZF `.gz` input continues through the
   existing zlib streaming path, byte-for-byte. Parallelizing plain gzip is out
   of scope for this module.
3. **Flex-first.** Flex has no legacy output contract (no prior STAR produced
   Flex output), so the fast path is simply the Flex behavior. Bulk/scRNA wiring
   is a non-goal here (those paths carry the byte-exact legacy guarantee and get
   their own runbook + legacy mode later).
4. No new external dependencies. No new file formats required at input; the
   optional index sidecar (below) is a cache, never a requirement.

## Detection

A file is BGZF iff its first gzip header has: magic `1f 8b`, `FLG.FEXTRA` set,
and an extra subfield with `SI1=66 ('B')`, `SI2=67 ('C')`, `SLEN=2`, whose
payload `BSIZE` is the total member length minus one. Detection reads only the
first 18 bytes. Files with FEXTRA but no `BC` subfield, or without FEXTRA, are
plain gzip → existing path.

The trailing 28-byte empty member is the BGZF EOF marker. Its absence on a
detected-BGZF file is a hard error (truncation) — a robustness gain over plain
gzip, which cannot detect truncation at all.

## Design

### 1. `input/BgzfIndex`

- Block scan by header hops: read the 18-byte header at offset, jump
  `BSIZE+1`, repeat. Produces `vector<{compressedOffset, compressedSize, isize}>`
  (`ISIZE` from the member trailer if wanted; may be deferred).
- **Record-count pass (required for mate pairing):** one parallel
  decompress-and-count pass over blocks records `recordsEndingPerBlock` (or
  cumulative record ordinals at block boundaries). Decompress-only throughput
  makes this cheap, and it is amortized by the sidecar.
- Optional sidecar cache `<file>.bgzi`: block offsets + record ordinals.
  Auto-written on first index build, auto-read (and validated against file
  size/mtime) on reuse. A cache only — absence or mismatch triggers a rebuild,
  never an error.

### 2. `input/BgzfRangeReader`

- A worker owns a contiguous block range. Raw-inflate each member (skip the
  header, `inflate` the deflate stream, optional CRC32 verify behind a flag,
  default on).
- **Record boundaries are exact, not heuristic:** with per-block record
  ordinals from the index, each worker's range is defined in *record* space
  (records k..m), and the worker knows which block its first record starts in
  and at what decompressed offset (stored during the count pass, or derived by
  parsing forward from the block in which record k begins). A record spanning
  past the range's last block is completed by inflating into the following
  block(s) — bounded overlap, deterministic ownership: a worker owns exactly
  the records that begin in its range.

### 3. Mate pairing (R1/R2)

R1 and R2 are separate files with equal record counts. Sharding is by **record
ordinal**: worker w owns records [k_w, k_{w+1}) in *both* files, using each
file's own index to locate the containing blocks. The record-count pass makes
this exact. A mismatch in total record counts between mates is a hard error.

### 4. Wiring

- New `input/BgzfStarAdapter` mirroring `input/CbqStarAdapter`: same producer
  interface feeding the fused Flex consumers (`FlexPipeline` /
  `PfMultiProcess` paths), per-worker feeds, order-independent (Flex carries no
  ordering contract; divergence from serial-order runs localizes to the
  documented order-dependent tie variance classes).
- Flags (parametersDefault):

```text
readFilesBgzfMode           auto
    string: auto|off|range - BGZF parallel range readers for gzipped FASTX input (Flex)
                            auto  ... sniff header; enable range readers for supported
                                      order-independent Flex runs; otherwise fall through
                                      to the standard zlib stream
                            off   ... always use the standard zlib stream
                            range ... require BGZF range readers; error if input is not BGZF

bgzfReaderThreads           0
    int: >=0: reader/inflate worker threads. 0 ... derive from --runThreadN.

bgzfCrcCheck                1
    int: 0/1 - verify per-block CRC32 during inflate.
```

- Scheduler: reader workers run outside the permit domains in v1 (like the CBQ
  readers). Accounting them under the FEATURE domain is future work, noted in
  the code with a TODO, not attempted here.

## Phase 0 — fixtures and gold tests (write these first)

Fixture prep: `tools/make_bgzf_fixture.sh` — re-blocks an existing gzip fixture
to BGZF (build the in-tree `bgzip` from vendored htslib sources for the tool
only, or a ~50-line python re-blocker; the tool is test-side only and must not
add a runtime dependency). Apply to the Flex smoke fixture and the tiny
Cell Ranger fixture used by `tests/run_dynamic_threads_tiny_fixture.sh`.

- **T1 detection unit tests:** BGZF, plain gzip, gzip+other-FEXTRA, empty file,
  bare EOF marker, missing EOF marker.
- **T2 index correctness:** block count and offsets equal an independent
  reference scanner (small python checker committed under `tests/`); sidecar
  write/read/invalidation (touch file → rebuild).
- **T3 record equality:** records parsed via BgzfRangeReader over N workers ==
  `zcat` reference: same count, same order-insensitive 64-bit checksum of
  (name, seq, qual) triples; run at 1, 3, 8 workers.
- **T4 end-to-end equivalence:** Flex smoke run, gzip input vs BGZF-re-blocked
  identical input → identical filtered matrices and per-sample outputs.
  Assert bit-identical if the fixture proves free of order-dependent ties;
  otherwise assert functional equality and record which variance class absorbs
  the difference.
- **T5 truncation:** strip the EOF marker → hard error; truncate mid-block →
  hard error with block offset in the message.
- **T6 mixed lanes:** some lanes BGZF, some plain gzip, in one run → correct
  outputs, BGZF lanes on range readers, gzip lanes on the zlib path.
- **T7 regression:** existing gzip/CBQ test suites pass unchanged with the
  module compiled in and `readFilesBgzfMode off` as well as `auto`.

## Phase 1 — `BgzfIndex` (+ sidecar, + record-count pass)
## Phase 2 — `BgzfRangeReader` (inflate, CRC, exact record ranges)
## Phase 3 — `BgzfStarAdapter` wiring, flags, docs in parametersDefault
## Phase 4 — benchmarks

Protocol (cold cache, 32 threads, no-align Flex, same invocations as the
Supp S5 rows):

1. JAX SC2300771, 8 lanes as delivered (BGZF):
   `/mnt/pikachu/JAX_sequences/JAX_scRNAseq01` — report the triple
   gzip-path (`readFilesBgzfMode off`) vs BGZF range readers vs CBQ.
2. 10x public 320k scFFPE GEM-X 16-plex (once extracted):
   `/mnt/pikachu/tenx_320k_scFFPE` — same triple; this is the public-facing
   benchmark row.

Record wall clock, `/usr/bin/time -v` peak RSS, and CPU%; store logs under
`docs/benchmarks/bgzf_ingest_20260901/`.

## Acceptance criteria

- T1–T7 green in CI (path-filtered like the other suites).
- JAX no-align with BGZF range readers within ~10% of the CBQ wall time.
- `readFilesBgzfMode off` output identical to current master on all existing
  tests (zero behavior change when disabled).
- No modifications under `core/legacy/source/htslib/`.
- parametersDefault documents the three flags; CHANGES entry added.

## Non-goals

- Parallelizing plain (non-BGZF) gzip input.
- BGZF/BAM output changes of any kind.
- Bulk/scRNA ingest wiring (legacy byte-exact contract; separate runbook).
- Writing BGZF FASTQ output.
