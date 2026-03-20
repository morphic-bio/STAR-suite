# Runbook: Internal Hash Cache Generator

**Date:** 2026-03-20
**Status:** Implemented in STAR (`--runMode hashCacheGenerate`, `--hashCacheOutput`, `--hashCacheTiers`, `--hashCacheParentLimit`)

## Problem

The current H0/H1/H2 hash cache build is an external Python+STAR+MEX pipeline:

1. Python generates synthetic FASTQ (one read per variant sequence)
2. STAR aligns them as if they were real reads
3. Python reads the Solo MEX output to determine which variants map unambiguously
4. Python writes the FH01SEQ1 binary cache

This is slow because:
- **588M H2 variant reads** require 980 STAR invocations (one per shard of 600K reads, limited by CB whitelist size)
- Each Python invocation reloads the 8M-record binary cache (~30s)
- Total estimated wall time: **6-40 hours** depending on optimization
- The approach is fundamentally I/O bound: writes FASTQ, aligns, writes BAM, writes MEX, reads MEX

## Solution: Internal Cache Generator

Move cache generation inside STAR as a dedicated run mode. STAR already has everything needed at runtime:
- **Probe 50-mer sequences** are in the genome `G[]` array as probe pseudo-chromosomes (53,444 chromosomes, each exactly 50bp, identified by `ENSG` prefix in `chrName`)
- **`encodeWindow`** packs 50-mers into `(seqLo, seqHi)` pairs
- **The full alignment pipeline** (`mapOneRead`) can classify any synthetic read
- **Gene annotation** via `classifyAlign` / probe-chromosome detection
- **`FlexHashScreenCache`** can write FH01SEQ1 binary format

Instead of synthesizing FASTQ → aligning → reading MEX, the internal generator:
1. Extracts probe sequences from genome pseudo-chromosomes at startup
2. For each probe: generates H0, H1, H2 packed sequences in-memory
3. For H1/H2 variants: feeds them through the alignment pipeline internally
4. Checks if the aligner maps each variant unambiguously back to the correct gene
5. Writes qualifying records directly to the FH01SEQ1 binary cache file

### Why this is fast

- **No I/O overhead**: no FASTQ, no BAM, no MEX — all in-memory
- **No Python overhead**: variant generation in C++ is orders of magnitude faster
- **Single genome load**: genome is already in memory for the alignment
- **Alignment is cheap**: 2B real reads align in 39 minutes; 588M synthetic 90bp reads should take ~11 minutes of pure alignment time
- **No shard/CB constraint**: no need for unique cell barcodes — we check alignment directly

## Architecture

### Data Flow

```
┌─────────────────────────────────────────────────────────────────┐
│  STAR --runMode hashCacheGenerate                               │
│                                                                 │
│  1. Load genome (G[], chrName, chrStart, chrLength)             │
│  2. Load probe list (ProbeListIndex → geneIdx15 map)           │
│                                                                 │
│  3. Extract probe sequences from genome:                        │
│     for each chrName starting with "ENSG":                      │
│       probeSeq = G[chrStart[chr] .. chrStart[chr]+50]           │
│       geneId   = chrName split on '|' → first field             │
│       geneIdx15 = ProbeListIndex::geneIndex15(geneId)           │
│       encode(probeSeq) → (seqLo, seqHi)                        │
│       → H0 record: (seqLo, seqHi, geneIdx15, class=0)          │
│                                                                 │
│  4. Generate H1 variants (single substitution):                 │
│     for each probe (50bp):                                      │
│       for pos in 0..49:                                         │
│         for alt in {A,C,G,T} \ {ref[pos]}:                     │
│           mutSeq = probe with pos→alt                           │
│           encode(mutSeq) → (seqLo, seqHi)                      │
│           → candidate H1 record                                 │
│                                                                 │
│  5. Generate H2 variants (two substitutions):                   │
│     for each probe (50bp):                                      │
│       for i in 0..49, j in i+1..49:                             │
│         for alt_i, alt_j (3×3 = 9 combos):                     │
│           mutSeq = probe with i→alt_i, j→alt_j                 │
│           encode(mutSeq) → (seqLo, seqHi)                      │
│           → candidate H2 record                                 │
│                                                                 │
│  6. Validate candidates via internal alignment:                 │
│     Construct synthetic read (fill + variant@offset0 + sample)  │
│     Feed through mapOneRead → classifyAlign                     │
│     If maps uniquely to expected gene → KEEP (class 1 or 3)    │
│     If ambiguous / wrong gene / unmapped → DENY (class 2)      │
│                                                                 │
│  7. Deduplicate by (seqLo, seqHi):                              │
│     Same 50-mer from different parents → check gene agreement   │
│     All agree on same gene → KEEP                               │
│     Disagreement → DENY                                         │
│                                                                 │
│  8. Write FH01SEQ1 binary cache                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Threading Model

The variant generation and alignment validation are embarrassingly parallel by probe:

```
Main thread:
  - Load genome, probe list
  - Enumerate probe pseudo-chromosomes → probe table (shared, read-only)
  - Partition probes across N worker threads

Worker thread i:
  - For each assigned probe:
    - Generate H0 record (trivial — just encode)
    - Generate H1 candidates (150 per probe)
    - Generate H2 candidates (11,025 per probe)
    - For each H1/H2 candidate:
      - Build synthetic 90bp read in a thread-local buffer
      - Run alignment via ReadAlign::mapOneRead (thread has own ReadAlignChunk)
      - Check gene assignment
    - Append validated records to thread-local output vector

Main thread:
  - Merge per-thread vectors
  - Deduplicate across probes (same 50-mer from different parents)
  - Sort by (seqHi, seqLo, sampleIdx)
  - Write FH01SEQ1 binary
```

Each STAR alignment thread already owns a `ReadAlignChunk` with its own `ReadAlign` instance, genome pointers, and transcriptome reference. The existing STAR threading infrastructure handles this.

### Variant Counts

| Tier | Per probe | Total (53,294 unique H0) | Records |
|------|-----------|-------------------------|---------|
| H0   | 1         | 53,294                  | 53,294  |
| H1   | 150       | ~8.0M                   | ~8.0M   |
| H2   | 11,025    | ~588M                   | ~588M   |

H0 records need no alignment validation (they are the probe itself).
H1 records: 8M alignments — trivial (~seconds).
H2 records: 588M alignments — ~11 minutes at 2B reads / 39 min throughput.

### Synthetic reads vs Solo output

The external pipeline wrote MEX and inferred KEEP/DENY from Solo. The internal generator still builds **realistic Flex PE reads** (R2 = 90bp with variant + sample tag; R1 = CB+UMI from whitelist/fixed UMI) so `mapOneRead` → `alignedAnnotation` match production. **`soloRead->record` is skipped** (per-thread `ReadAlign::hashCacheSynthProbe_`) so the inline hash / MEX path is not polluted.

This removes the external shard/CB constraint for cache *generation*; CB/UMI are only needed to satisfy the same read layout as mapping.

### Alignment Validation Logic

For each H1/H2 candidate variant sequence:

1. Construct a synthetic 90bp read:
   - Positions 0-49: variant 50-mer
   - Positions 50-67: fill base (A)
   - Positions 68-75: sample tag (from `--soloSampleWhitelist`, first entry)
   - Positions 76-89: fill base (A)

2. Run the full PE path (`mapOneRead` → `outputAlignments` with synthetic CB/UMI) so annotation matches production.

3. **KEEP** only if the read would pass the same gates as `record_flex` → `outputReadCB_flex` for the configured Gene-type Solo feature, and **`flexResolveGeneIdx15_inlineResolver`** (shared with `outputReadCB_flex`) returns **`geneIdx15 == expectedGeneIdx15`**:
   - Valid CB (`cbMatch >= 0`), sample-tag policy matches production (whitelist / require-match: no tag → drop)
   - Gene set: `fSet` non-empty; multi-gene reads dropped unless `multiMap.yes.multi` (same as `record_flex`)
   - Resolver: same per-alignment probe vs genomic `CandidateView` build, filters, and **`resolveGeneFromCandidates`** as normal inline hash
   - If unmapped, resolver returns 0, wrong gene, or gates fail → **DENY**

4. The old union-of-`fAlign` / `fSet` heuristic is removed; validation is authoritative from the inline resolver only.

### Implementation Phases

**Phase 1: H0 only (no alignment needed)**
- Extract probe sequences from genome pseudo-chromosomes
- Encode each 50-mer → `(seqLo, seqHi)`
- Map `chrName` → `geneId` → `geneIdx15` via `ProbeListIndex`
- Write H0 records (cacheClass=0) to FH01SEQ1
- This is a few seconds of work and validates the probe extraction pipeline

**Phase 2: H1 with alignment validation**
- Generate all 150 single-substitution variants per probe
- Align each through `mapOneRead`
- Validate gene assignment
- Write KEEP records (cacheClass=1) and DENY records (cacheClass=2)
- ~8M alignments, should complete in seconds

**Phase 3: H2 with alignment validation**
- Generate all 11,025 two-substitution variants per probe
- Same alignment + validation pipeline as H1
- ~588M alignments, ~11 minutes
- Write KEEP (cacheClass=3) and DENY (cacheClass=2)

**Phase 4: Cross-probe deduplication**
- Merge all emitted rows by `(seqLo, seqHi)`; any DENY or disagreeing KEEP genes → final DENY (`cacheClass=2`, `FlexHashNegProbeAmbig`).
- **H0** rows use **`sampleIdx` = the hash-screen sample index** (`detectSampleIndex` → token → `sampleIndexForToken`, same as `ReadAlign_oneRead`). The generator emits **one H0 record per probe per whitelist row** that yields a non-zero index (same packed 50-mer, different `sampleIdx` keys). Dedup is keyed by `(seqLo, seqHi, sampleIdx)` so multiple samples do not collide. `FlexHashScreen::findRecord` + `classifyHits` then hit the **sample-matched exact-H0 fast path** for every sample represented in the cache.
- **H1/H2** rows stay **global** (`sampleIdx=0`); the sample tag in R2 is for detector/alignment realism.
- **KEEP validation** for H1/H2 uses the **production Flex gene resolver** (`flexResolveGeneIdx15_inlineResolver`); multimappers are accepted when that resolver returns the parent’s `geneIdx15` (same rules as a real alignment that would insert into the inline hash).

### CLI Interface

```bash
STAR \
  --runMode hashCacheGenerate \
  --genomeDir /path/to/star_index \
  --soloProbeList /path/to/probe_list.txt \
  --soloSampleWhitelist /path/to/sample_whitelist.tsv \
  --soloSampleProbes /path/to/probe-barcodes.txt \
  --soloSampleProbeOffset 68 \
  --hashCacheOutput /path/to/output_cache.bin \
  --hashCacheTiers H0,H1,H2 \
  --runThreadN 8
```

Parameters:
- `--hashCacheOutput`: path for the FH01SEQ1 binary output
- `--hashCacheTiers`: comma-separated list of tiers to generate (default: `H0,H1,H2`)
- Standard genome/probe/sample parameters are reused

### Files to Modify / Create

| File | Change |
|------|--------|
| `core/legacy/source/FlexHashScreen.h` | Add `writeFile` static method, `generateRecords` |
| `core/legacy/source/FlexHashScreen.cpp` | Implement `writeFile`, variant generation, alignment validation |
| `core/legacy/source/STAR.cpp` | Add `hashCacheGenerate` run mode dispatch |
| `core/legacy/source/Parameters.cpp` | Add `--hashCacheOutput`, `--hashCacheTiers` parameters |
| `core/legacy/source/ParametersSolo.h` | Add `hashCacheOutput`, `hashCacheTiers` fields |
| `core/legacy/source/Genome.cpp` | Add helper: `extractProbeSequences()` returning probe table |

### Output Format

Same FH01SEQ1 v2 binary format (24-byte header + 24-byte records), sorted by `(seqHi, seqLo, sampleIdx)`. Fully compatible with the existing `FlexHashScreenCache::loadFile` and the replay tools.

### Estimated Runtime

- Genome load: ~23 seconds
- Probe extraction + H0 encoding: <1 second
- H1 generation + alignment: ~5 seconds (8M alignments)
- H2 generation + alignment: ~11 minutes (588M alignments)
- Deduplication + sort + write: ~10 seconds
- **Total: ~12 minutes** (vs 6-40 hours for the external pipeline)

### Validation

Compare internal generator output against external pipeline output:
1. Run external pipeline on a small subset (e.g., 40 parents)
2. Run internal generator on the same subset
3. Diff the FH01SEQ1 binaries — should be identical (same KEEP/DENY decisions)
4. Load both in `hash_screen_replay --stats` and compare tier composition

### Risks and Mitigations

| Risk | Mitigation |
|------|-----------|
| Alignment results differ from full STAR run (different trimming, Solo processing) | The internal path skips CB/UMI/Solo; validate that gene assignment from `classifyAlign` matches MEX-based assignment on the pilot subset |
| Memory: 588M candidate records in flight | Process in batches (per-probe or per-thread); only accumulate KEEP records |
| Thread safety of `mapOneRead` | Each thread uses its own `ReadAlignChunk` (existing STAR pattern) |
| H2 deduplication across probes | Same 50-mer from different parents must agree on gene; use hash map keyed by `(seqLo, seqHi)` |
