# STAR Suite v1.9.0 Release Notes

Date: 2026-09-10

STAR Suite 1.9.0 improves Flex probe assignment, grouped sample cell calling,
and cache/caller performance. It ships the changes since v1.8.4, including
the all-feature model and bounded parallel caller used in the 320K work.

`STAR --version` reports `1.9.0`. Debian source packaging uses `1.9.0-1`;
the Ubuntu 22.04 and 24.04 binary packages use
`1.9.0-1~ubuntu22.04.1` and `1.9.0-1~ubuntu24.04.1`. The upstream STAR version
remains `2.7.11b`, genome-index compatibility remains `2.7.4a`, and legacy
compatibility remains `2.7.1a`. Existing indexes do not need rebuilding.

## Flex cell calling

- The tag-aware default groups tags belonging to the same sample into one
  model while retaining complete CB16+TAG8 cell identities. The legacy caller
  remains selectable explicitly. The CR-config recipe selects tag-aware mode
  without injecting its former incompatible legacy expected-cell override.
- Bootstrap OrdMag uses deterministic quality ordering for equal total UMIs:
  non-mitochondrial counts when an annotation mask is supplied, then detected
  genes, then full barcode identity. Tie-breaking itself needs no bootstrap.
- Both primary and fallback selection enforce the configured inclusive UMI
  floor (500 by default), preventing low-count unused tags from generating
  spurious cell calls that distort occupancy.
- The model retains all probe features, including separate deprecated-feature
  IDs. `--soloFlexFilteredGeneList` restricts filtered exports without removing
  those features from OrdMag, ambient estimation, or EmptyDrops.
- Tag-aware defaults use a fixed ambient rank window, 100,000 Monte Carlo
  simulations, and BH-adjusted FDR 0.01. Joint occupancy runs after all sample
  groups finish, using observed distinct tag calls per GEM and a fitted
  zero-truncated Poisson mean with the ordinary Poisson 0.999 cutoff.
- Caller diagnostics expose ranks, candidates, ambient profiles, p-values,
  rescue decisions, and occupancy removals. See
  [model and export configuration](https://github.com/morphic-bio/STAR-suite/blob/v1.9.0/docs/FLEX_MODEL_FEATURES.md) and
  [quality ranking](https://github.com/morphic-bio/STAR-suite/blob/v1.9.0/docs/ORDMAG_QUALITY_RANKING.md).

## Cache and caller performance

- `flex_hash_cache_pack` writes completed, versioned khash lookup tables.
  STAR maps them read-only without sorting or rebuilding the stored hashes.
  Legacy record caches remain readable and must follow the generation-time
  sort contract.
- For complete H1X2 caches, `flex_hash_cache_pack --half` stores an exact H0
  table and two 25-base Hamming-1 tables. Intersecting parent lists recovers
  full-probe matches; one-sided candidates pass directly to full-probe
  Hamming scoring. Packing verifies the source decision universe before
  publishing an immutable file. Unsupported cache policies retain the
  general format. See [cache format and conversion](https://github.com/morphic-bio/STAR-suite/blob/v1.9.0/docs/FLEX_KHASH_CACHE.md).
- Each OrdMag bootstrap replicate sorts its array once and reuses it for all
  trial cell counts, preserving the previous estimator's arithmetic and draws.
- Independent sample groups run concurrently. A shared `runThreadN` permit
  budget covers group coordinators and bootstrap/EmptyDrops helpers. Active
  samplers can borrow returned permits as other groups finish, with stable
  per-iteration random streams and exact integer tally reduction.

## Probe assignment, input, and auditing

- Exact supplied H0 probe matches are authoritative. Parallel synthetic
  Hamming verification is deterministic and no longer mutates shared
  alignment parameters. Flex's default length gates accommodate probe-prefix
  alignment while preserving explicit user overrides.
- Sample-tag lookup uses the fixed supplied-table offset: authoritative H0
  first, followed by a separately constructed H1 table that rejects ambiguous
  ownership. Neighboring-offset search is unsupported. See
  [sample-tag policy](https://github.com/morphic-bio/STAR-suite/blob/v1.9.0/docs/FLEX_SAMPLE_TAG_POLICY.md).
- Conservative single-N probe resolution is shared by fused input and BAM
  paths. Experimental H1X2 remains opt-in: it permits at most one mismatch per
  25-base half and supports a unique-half seed extension with at least 40 of
  50 matching bases. It is not a blanket replacement for H1/H2 cache policy
  in other assays.
- Paired FASTQ lane reading uses record batches and separate mate-reader
  threads, and honors `--readFilesCommand`. CBQ and BGZF retain their native
  reader paths.
- An optional BAM-independent binary decision sidecar and native BAM/MEX audit
  tools support read-level diagnosis. The sidecar is disabled by default and
  is unnecessary for count-only benchmarking. Comparator cell keys retain
  the full barcode and tag; comparator corrections do not imply a Flex
  barcode-merging defect.

## Full 320K benchmark and validation

Both new benchmark arms processed all four lanes: **7,303,142,230 read pairs**,
eight samples and sixteen probe tags, on the same 48-thread `m6id.12xlarge`
instance. Both used the same frozen STAR 1.9.0 binary, cold OS page cache,
100,000 simulations, and the all-feature model. There was no alignment,
reference-index loading, BAM, or per-read decision sidecar. Input conversion
was outside the timed process. Each number is one completed execution, not a
repeated-trial estimate.

| Input | September 9 control | STAR 1.9.0 | Speedup | Peak RSS, control → new |
| --- | ---: | ---: | ---: | ---: |
| CBQ | 1,217.150 s (20m17s) | **628.033 s (10m28s)** | **1.94×** | 138.85 → 103.97 GiB |
| FASTQ-BGZF | 1,452.071 s (24m12s) | **867.680 s (14m28s)** | **1.67×** | 139.94 → 93.79 GiB |

For each arm, the global raw matrix and all eight filtered matrices matched
exactly, including all **333,439 called cell identities**. All 57 caller
diagnostic files per arm matched after excluding only the execution metadata
field `mc_threads`; ranks, ambient profiles, p-values, rescue decisions,
occupancy results, sample summaries, and read classifications were unchanged.
The previously measured pooled cell Jaccard against the saved reference results
therefore remains **0.973961**. See the
[full 320K benchmark report](https://github.com/morphic-bio/STAR-suite/blob/v1.9.0/docs/benchmarks/FULL320K_V190_20260910.md) and the
[September 9 caller integration report](https://github.com/morphic-bio/STAR-suite/blob/v1.9.0/docs/HANDOFF_FULL320K_STAR_DEPRECATED_20260909.md).

The benchmark binary was clean-built from commit
`a6c1c546cbb53403c542bb9819706631801c935b`, SHA-256
`54311642da84160e2dd75321d1ea6306ef0a5724b32727b402216748a5707d67`.
Subsequent release integration fixed the recipe wrapper and generated-default
build rule and updated documentation. The refreshed parameter header is
byte-identical to the header generated by the benchmark's clean build; no
counting or caller implementation changed after that benchmark source freeze.
The final tag preserves this commit in its ancestry. Packaged platform builds
have their own binary checksums and embed the exact final tagged source commit.

## Separate L004 optimization evidence

The earlier L004 optimization run processed **1,823,648,323 read pairs** with
CBQ and 48 threads. The optimized implementation took **201.901 seconds**,
versus **428.954 seconds** with the prior fixed-share parallel caller; caller
time fell from 244 to 18 seconds. All matrices, 249,194 called cell identities,
ranks, ambient profiles, and p-values matched across those optimization arms.
These are separate single-run, pre-version-bump implementation measurements;
they are not the full-set timings above. See the
[L004 report](https://github.com/morphic-bio/STAR-suite/blob/v1.9.0/docs/benchmarks/FLEX_ORDMAG_SORT_PERMITS_L004_20260910.md).

## Distribution and compatibility

Release packaging provides amd64 glibc234/glibc239 tarballs, an installer that
selects the compatible tarball, Ubuntu 22.04/24.04 Debian binaries, source
packaging, and checksums. Hosted tarball/Debian builds use the established
portable no-Chromap build; production multiome source builds retain the
Chromap-enabled default. Release checks cover packaged versions and source
commit provenance, runtime loading, installer selection, smoke tests, and
Debian installation/removal on the supported Ubuntu releases. Generated parameter headers are
content-checked and replaced atomically during builds, so source archives with
equal file timestamps cannot silently compile stale parameter defaults.

The compact probe cache and its gene-ID list are paired artifacts. Repack
after changing probe records, keep mapped cache files immutable, and retain
the original record cache when extending the model. The all-feature model
changes cell-calling behavior relative to v1.8.4; the sort and permit changes
preserve the validated model's statistical outputs.
