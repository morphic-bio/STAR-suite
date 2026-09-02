# STAR Suite v1.8.0 Release Notes

Date: 2026-09-02

`v1.8.0` adds on-demand parallel BGZF FASTQ ingestion to supported fully fused
STAR-Flex runs. BGZF inputs can now be consumed in place without a conversion,
pre-scan, external index, or sidecar cache while preserving exact FASTQ record
boundaries and R1/R2 pairing.

The release artifact version is `v1.8.0`, Debian packages use `1.8.0-1`, and
`STAR --version` reports `1.8.0`. The upstream STAR base remains `2.7.11b`;
genome-index compatibility remains `2.7.4a`. Existing indexes do not need to
be rebuilt.

## Parallel BGZF FASTQ ingest

- Detect BGZF from the inline gzip `BC` extra subfield and discover bounded
  work claims lazily from each member's `BSIZE`. Claimed member descriptors are
  cached only in memory; no `.bgzi` or other sidecar is read or written.
- Inflate independent members in parallel with the vendored zlib raw-deflate
  API, validate member sizes, and verify CRC32 by default.
- Reassemble completed work in compressed-stream order before FASTQ parsing.
  Records and lines may span BGZF members; correctness never depends on
  searching for heuristic `@` boundaries.
- Parse R1 and R2 as independent ordered streams and pair records by exact
  ordinal. Unequal mate record counts are a hard error even when the two files
  use different BGZF block layouts.
- Accept a structurally complete stream whose final declared member ends at
  physical EOF without requiring the optional canonical empty BGZF EOF member.
  Partial headers, members, deflate streams, and trailers remain hard errors.

## Controls and compatibility

- `--readFilesBgzfMode auto|off|range` selects header-detected fallback,
  unconditional standard zlib streaming, or required BGZF range ingest.
- `--bgzfReaderThreads 0` derives the reader pool from `--runThreadN`; a
  positive value sets it explicitly.
- `--bgzfCrcCheck 1|0` controls per-member CRC32 validation and defaults on.
- The new path is restricted to supported order-independent, fully fused
  STAR-Flex configurations. Plain gzip, unsupported Flex configurations, and
  mode `off` retain the existing zlib stream. Bulk and STARsolo ingest paths
  are unchanged.
- The implementation does not modify or use the vendored htslib BGZF reader
  and introduces no external runtime dependency.

## Dynamic permit control

- Extend the saturation-aware controller to coordinate MAP, FEATURE, and ATAC
  work domains. Optional domain work estimates guide probe ordering and
  remaining-work balancing while observed demands remain borrowable floors.
- Keep the controller opt-in through `--dynamicThreadAtacController 2` and its
  existing telemetry requirements; default scheduling behavior is unchanged.

## Validation and performance

- BGZF tests T1-T7 pass for detection, inline block traversal, record equality
  at 1/3/8 workers, end-to-end Flex equivalence, EOF-marker and truncation
  handling, mixed BGZF/plain-gzip lanes, and gzip/CBQ regression coverage.
- Existing Flex smoke and gzip/CBQ regression suites pass with the module
  compiled in. `--readFilesBgzfMode off` preserves the established path.
- Retained benchmark wrappers cover both the eight-lane JAX protocol and the
  public 320K scFFPE protocol with fresh-output and completion-marker guards.
- On the cold-cache, 32-thread, eight-lane JAX no-align benchmark, BGZF range
  ingest completed in 7:38.71 versus 7:33.91 for CBQ and 10:42.19 for the
  standard gzip path. All three runs produced identical Flex counters; BGZF
  was 1.06% from the CBQ wall time.

The implementation contract and retained benchmark evidence are in
`docs/RUNBOOK_BGZF_PARALLEL_INGEST_20260901.md` and
`docs/benchmarks/bgzf_ingest_20260901/`.
