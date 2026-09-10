# Flex stored khash and half-probe benchmark — 2026-09-10

The production loader now uses khash for H0 and H1/DENY, with a self-contained
file format that stores both completed tables. The new STAR binary passed
CBQ/BGZF count and cell-call regression checks and is installed locally.
The two-half alternative is implemented as a benchmark prototype and is not
the production default.

## Classifier results

Serial runs on cloud instance `i-06de289faa5d78117`, g++ 13.3.0,
`-O3 -std=c++11 -mavx2`, warm OS page cache. Query input was the first million
50-base R2 windows from each of L001–L004: 4,000,000 reads, without sample or
GEM filtering. Classification includes H0, H1X2, single-N handling and Hamming
rescue. FASTQ timing includes converting the shared packed query fixture to
ASCII before calling the ASCII classifier, equally for all implementations.
These are single-thread classifier microbenchmarks, not complete STAR runs
or compressed-input throughput measurements.

| Implementation | CBQ classification, 4M windows | FASTQ classification, 4M windows | CBQ initialization |
| --- | ---: | ---: | ---: |
| Previous unordered_map, load sorts removed | 0.596937 s | 1.196820 s | 110.858 s |
| Stored full-probe khash | 0.544574 s | 1.060840 s | 0.711236 s |
| Two 25-base half hashes | 0.187824 s | 0.566410 s | 1.038170 s |

The half implementation is 2.90× faster than the full khash for packed CBQ
classification and 1.87× faster for ASCII classification. The final ASCII
implementation uses a base-code lookup table; the earlier switch-based
prototype took 2.084960 s and is superseded for this comparison. Its changed
ASCII path was checked against all four million previous verdicts.

All implementations return exactly 3,859,246 KEEP and 140,754 DENY decisions.
There are zero differences in action, gene, cache class, negative reason,
region, Hamming distance, single-N metadata or offset, including between
ASCII and CBQ. A separate exhaustive comparison of all 316,072,780 source
cache records also found zero decision differences between the full khash
and the half-probe method. The source cache contains H0/H1X2 KEEP records;
the real-read comparison and synthetic tests additionally exercise DENY and
rescue decisions. This is not a claim about cell concordance on a new full run.

## Half-probe method and memory

Each probe contributes its exact 25-base half and 75 Hamming-1 neighbours
on each side. Two sets of 76 half entries represent the 5,775 non-exact
full-probe variants previously enumerated per parent. There are 54,580
parent probes, 4,148,062 left keys and 4,148,060 right keys.

Half-table values retain sorted lists of parent IDs. A single common parent
resolves a double match. Multiple common parents are denied even when they
share a gene. A unique one-sided match carries its parent directly into
full-probe Hamming scoring, without looking up the halves again. The ten-base
distance limit, ambiguous-half behavior, split-probe rejection and existing
single-N policy are preserved.

The four-million-read route counts are:

| Route | Reads |
| --- | ---: |
| Exact H0 match | 3,650,166 |
| Double-half match | 173,231 |
| HALF_MATCH → Hamming scoring | 55,439 |
| No half anchor | 107,969 |
| Different probes on the two halves | 13,194 |
| Ambiguous | 1 |

Of the half matches, 35,849 pass Hamming scoring and 19,590 fail.
The prototype's allocated lookup structures occupy 310,565,760 bytes
(296.18 MiB), including H0, parent records and candidate lists. This is
structure size, not whole-process peak RSS. The prototype currently builds
these structures from the mapped H0 parent records; it does not persist the
half tables yet. Reading those parents temporarily maps pages throughout the
original record array, so measured peak RSS is about 7.27 GiB despite the
small retained lookup structures. The stored full khash peaks around
19.11 GiB during real-read classification, versus 47.54 GiB previously.

## Stored full khash implementation

The new full cache is 20,615,036,928 bytes. It embeds the original
7,585,746,744-byte cache's records plus H0 indices, flags, keys and compact
values. It preserves original feature IDs and sample ordering. The loader
maps it read-only, without rebuilding either table or copying tier vectors.
Legacy `.bin` files still work, with an in-memory khash build.

On the full cache, conversion plus exhaustive stored-table verification took
86.48 seconds once. An isolated load of the stored file took 0.709 seconds,
and legacy-format loading with khash construction took 36.70 seconds. The
H1X2 seed index is still built from H0 parents. mmap is lazy: startup does
not imply every hash page has been populated. Query timings above include
demand page faults, with file contents already in the OS cache.

Generation and invocation: [FLEX_KHASH_CACHE.md](../FLEX_KHASH_CACHE.md).
New cloud resources, to be used together with the existing model gene list:

- Binary: `/scratch/flex_khash_star_fixture_20260910_v1/STAR`
- Cache: `/scratch/flex_khash_20260910_v1/model_h01x2_cache.khash`
- Gene list: `/scratch/full320k_star_deprecated_20260909_v1/config/model_gene_ids.txt`

Installed local binary: `/mnt/pikachu/STAR-suite/core/legacy/source/STAR`.
SHA256: `dca77afeb1d6e57c36d98b24806c2d0962ffb0304a855a7ef0addd9e649e8790`.
The previous binary is preserved under the local artifact directory.

## Validation and provenance

- Legacy versions 1/2/3, persisted payloads, full keys, sample-aware lookup,
  output non-overwrite, corrupt metadata and truncation tests pass.
- 3,522 actual-loader queries match the previous loader across legacy/stored
  files and both input modes; this includes H2 routing and gDNA metadata.
- Half-intersection, same-gene ambiguity, direct scoring, distance-limit,
  split-probe and single-N tests pass.
- A clean complete STAR build with default Chromap support succeeds.
- The prior CBQ and BGZF read regression fixture was each executed once with
  the new STAR and stored cache. Both have exact raw-MEX and cell-call parity
  with their prior runs, including eight called cells. No genome-index load,
  BAM, SAM or decision sidecar was used.

The first large comparison (`flex_pair_benchmark_20260910_v1`) was cancelled
and excluded: the old implementation's driver picked up the new header,
invalidating inline accessors. Version 2 places the driver outside source
directories and uses matching headers; its fixture parity and H1X2 readiness
checks pass. The initial STAR-fixture wrapper stopped while interpreting
fileless empty-sample directories; existing successful CBQ output was
validated post hoc, and BGZF was then run once. Neither correction required
repeating an identical successful dataset execution.

Local scripts, exact source snapshots, manifests, checks, completed logs and
reports: `/mnt/pikachu/star_suite_paper/analysis/flex_khash_20260910/`.
Durable evidence bucket:
`s3://star-suite-320k-benchmark-alt-171440768238-us-west-2-20260904/analysis-tools/`.
Prefixes: `flex_khash_20260910_v1`, `flex_pair_benchmark_20260910_v2`,
`flex_pair_ascii_20260910_v3`, and `flex_khash_star_fixture_20260910_v1`.
The packed 20.6-GB file remains on the running instance; the original source
cache and conversion code provide its reproducible inputs. No instance was
stopped, snapshotted or terminated. Changes are not committed or pushed.
