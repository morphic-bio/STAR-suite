# Stored Flex khash caches

`flex_hash_cache_pack input.bin output.khash` is the final cache-generation
step. It accepts sorted FH01SEQ1 version 1, 2 or 3 records and writes a
self-contained, versioned file with the completed H0 and H1/DENY tables.
The output is published atomically and existing files are never overwritten.
After writing, the tool reopens it and compares every distinct tier key and
payload against the embedded records. Any failed verification exits nonzero.

Point the usual `--soloHashScreenFile` option at this file. FASTQ, BGZF and
CBQ use the same artifact. No new STAR option or genome-index dependency is
introduced. The associated gene-ID list remains mandatory and unchanged.
Keep generated caches immutable while STAR processes have them mapped.

The loader uses the vendored klib khash with full two-word probe keys and
8-byte values. H0 is queried first. Within each tier, the first record in
sorted `(seqHi, seqLo, sampleIdx)` order supplies the sample-independent
value, preserving the previous `unordered_map::emplace` behavior. Exact
sample lookup and global fallback search the embedded sorted records;
H0-only sample lookup uses a small stored array of record indices. H2 records
are retained for the full classifier and do not enter the offset-zero tier.

The persistent format is `FH01KH01`, version 1, on 64-bit little-endian hosts:

- A 256-byte header identifies the format, source record version, key
  encoding, hash algorithm, counts and section offsets.
- Original 24-byte records and 64-bit H0 record indices retain file order.
- Each table stores khash flags, 16-byte keys and 8-byte values in separate
  page-aligned arrays. No pointers or allocator state are serialized.
- Lookup keys use CBQ's LSB-first order; original records retain their
  original MSB-first order. FASTQ lookup keys are converted to the same
  representation. A fixed hash mixer precedes khash's bucket mask.
- Loading validates header compatibility, lengths, section boundaries and
  the small H0 index, then binds read-only mapped arrays. It does not scan
  the large tables, reinsert entries, check record order or sort.

Legacy record files are also mapped directly, then indexed with khash in
memory. They avoid per-record stream decoding and tier copies but still pay
table-construction cost. For these two general formats, the H1X2 seed index is built from H0 parents.
The compact format below persists the half hashes as well.

Lazy mmap startup does not mean all table pages have been read. Benchmark
lookup time, page faults and full-process memory alongside startup. On the
316,072,780-record model cache, the stored file is 20,615,036,928 bytes versus
7,585,746,744 bytes for the record-only input; this trades disk space for
avoiding repeated index construction.

Storage regression: `bash tests/test_flex_khash_storage.sh [fresh-output-dir]`.
The actual-loader comparison driver is `tests/flex_khash_lookup_digest.cpp`;
compile it outside implementation source directories with each version's
matching include path. Mixing old implementations and new headers invalidates
the comparison because inline accessors depend on class layout.

## Compact H1X2 production cache

For a complete H0/H1X2 v3 cache, use the compact representation:

```bash
flex/tools/hash_screen_replay/flex_hash_cache_pack --half \
  model_h01x2_cache.bin model_h01x2_cache.half.khash
```

The input may also be an existing full `.khash` snapshot. Select the output
with `--soloHashScreenFile`; STAR detects its format automatically. The gene-ID
list remains the same paired artifact. The original record file remains the
input to tools that extend the probe universe.

`core/legacy/source/FlexProbePairKhash.h` implements the production matcher.
The small exact H0 hash is retained. Two 25-base Hamming-1 khashes store spans
into sorted parent-ID lists. Intersecting the two lists identifies full-probe
matches. More than one common parent is ambiguous even when the probes share
a gene. A unique one-sided `HalfMatch` carries its parent directly into the
existing full-probe Hamming scoring, without repeating the half lookup. The
ten-mismatch limit, split-probe rejection, single-N merge and region/distance
provenance are preserved. Fused CBQ and ASCII paths use the combined lookup;
the primary-only and seed-only APIs remain available separately.

The `FH25KH01` version-1 file contains a 256-byte header, unique probe parents,
original sample-specific H0 records, both parent-ID arrays, and all three
completed hashes. Sections are page-aligned and contain no pointers. Empty
hash slots are initialized before serialization. STAR maps the arrays read-only
and checks their section bounds, occupancy, parent IDs and spans; it builds
neither the large full-variant table nor the old half-seed index. Sample-aware
H0 queries use the retained records, and global non-exact records are recovered
from the half-table intersection.

Packing verifies every source record's primary decision and enumerates every
pair of half variants to ensure it exists in the source. Thus a truncated or
partial variant cache cannot acquire new matches silently. Unsupported H1/H2
or sample-specific non-exact policies fail conversion; keep the general cache
format for those inputs. After writing, the stored arrays are reopened and
compared byte-for-byte with the verified tables. Cache files must remain
immutable while in use.

The deprecated-complete 320K model's compact file is **331,268,096 bytes**
(315.92 MiB), versus 20,615,036,928 bytes for the full khash snapshot. This
includes the retained H0 sample mappings. Conversion of the already packed
source, full verification and writing took 45.80 seconds on the cloud host;
this is generation work, paid once. See
[the L004 production validation](benchmarks/FLEX_HALF_KHASH_L004_20260910.md)
for end-to-end measurements.
