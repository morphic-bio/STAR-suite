# Hash Cache Tiered Design

> **STAR Suite 1.9.4:** the default cache is `--hashCacheTiers H0,H1X2`, and
> `--flex yes` assigns reads from it without any genomic alignment (see
> [H1X2](#h1x2--half-probe-tier-default-from-194) below). The alignment-validated
> H1 and H2 tiers, and sending cache misses to STAR alignment, are **legacy**:
> they run only with `--flexLegacy yes`, to reproduce results from earlier
> releases. The H1/H2 sections, recipe policy and validation numbers below
> describe that legacy design. Current Flex speed and concordance are in the
> top-level [README Benchmarks](../README.md#benchmarks).

## Overview

The Flex hash screen cache accelerates single-cell Flex alignment by
pre-computing which probe-like sequences resolve to a known gene. This avoids
full STAR alignment for the majority of reads. The cache is organized into
three tiers (H0, H1, H2) corresponding to Hamming distance from the reference
probe sequences, plus the half-probe tier H1X2, which from 1.9.4 replaces
alignment altogether.

## Tier Definitions

### H0 — Exact Match (Hamming 0)

- **Contents**: One record per probe × sample (53,459 probes × 16 samples = 855,344 records).
- **Record type**: KEEP only. Every exact probe match is correct by definition.
- **Sample index**: Each H0 record carries the sample index derived from the
  sample tag, enabling sample-specific fast-accept.
- **Generation**: No alignment validation needed — records are stamped directly
  from the probe list and sample whitelist.
- **Resolution policy**: An unambiguous H0 hit is final. The included 10x probe
  set defines the assay target and is assumed not to contain an equivalent
  genomic target. Sequence after the 50-base probe is assay payload and is not
  allowed to overturn the exact probe identity.

### H1 — Single Mismatch (Hamming 1)

- **Contents**: Up to 150 variants per probe (3 substitutions × 50 positions).
- **Record types**: **KEEP + DENY**.
  - **KEEP** (`cacheClass=1`): Variant aligns and resolves to the same gene as
    the parent probe. Stored with the gene index.
  - **DENY** (`cacheClass=2`, `negativeCode=FlexHashNegProbeAmbig`): Variant
    aligns but resolves to a different gene or is ambiguous (multi-gene hit,
    wrong gene, resolver failure). Stored to enable fast-reject at runtime.
  - **DEAD** (not stored): Variant does not map at all. Pruned because a real
    read with the same sequence also won't map — the aligner reaches the same
    unmapped conclusion without a cache entry.
- **Generation**: Each variant is constructed as a synthetic PE read (R2=90bp
  with 50bp probe + padding + sample tag; R1=28bp CB+UMI). The full STAR
  alignment pipeline runs with relaxed per-read-length filters (the 50bp probe
  is embedded in a 90bp mate). The tri-state verdict (KEEP/DENY/DEAD) is
  returned by `flexHashCacheValidateSyntheticPair`.
- **Sample index**: H1 records use `sampleIdx=0` (global, not sample-specific).

### H2 — Double Mismatch (Hamming 2)

- **Contents**: Up to 11,025 variants per probe (C(50,2) × 9 = 1,225 position
  pairs × 9 substitution combos).
- **Record types**: **KEEP only**.
  - **KEEP** (`cacheClass=3`): Stored.
  - **DENY**: Not stored. At Hamming 2, the fraction of ambiguous variants is
    large and storing them provides diminishing returns — the reads that match
    H2 DENY entries are rare in practice and the aligner handles them.
  - **DEAD**: Not stored (same rationale as H1).
- **Pre-check optimization**: Before aligning an H2 variant, its sequence is
  looked up in the existing H0 + H1 cache. If it already exists (as H0 KEEP,
  H1 KEEP, or H1 DENY), the variant is **skipped** — it's already handled at a
  closer Hamming distance, so no alignment is needed. This saves significant
  compute time since a subset of H2 sequences coincide with H1 variants of
  neighboring probes.
- **Sample index**: `sampleIdx=0` (global).

### H1X2 — Half-probe tier (default from 1.9.4)

- **Contents**: every non-exact 50-mer with at most one substitution in each
  25-base half of a probe: 150 single substitutions plus 75 × 75 paired
  substitutions per probe.
- **Generation**: no alignment. Records are generated from the H0 parent
  probes (H1X2 requires H0 in `--hashCacheTiers`). A full 50-base key that can
  be generated from more than one probe is stored as DENY, even when the probes
  share a gene.
- **Runtime**: after an H0 miss, each 25-base half is looked up in the
  single-mismatch half tables. A half that identifies one probe is extended
  against that whole probe with fast Hamming scoring, and the read is kept when
  the full 50 bases have at most ten mismatches. Halves that resolve to
  different probes are denied. The decision is terminal: nothing is sent to
  genomic alignment.
- **Sample index**: `sampleIdx=0` (global).
- **Storage**: see [the khash cache notes](FLEX_KHASH_CACHE.md) for the compact
  `.half.khash` form.

## Recipe Policy and Assay Context

This section describes the legacy (pre-1.9.4) H0/H1 route. The 1.9.4 default
is `--hashCacheTiers H0,H1X2` with no alignment.

The legacy production policy is fixed rather than selected by a command-line option.
An unambiguous exact H0 match is final, generated H1 KEEP and certified H1 DENY
records retain their verified verdicts, and only an unencodable probe window or
an absent H0/H1 key passes to live STAR alignment. This makes the probe set,
rather than whole-read competition from post-probe assay sequence, the
authority for exact matches while keeping the cache and live resolver in
agreement for single substitutions.

For routine scRNA-seq Flex processing on the legacy route, recipes requested
`--hashCacheTiers H0,H1`. In the JAX scRNA-seq benchmark, H2 recovered
measurable read-level signal (about 604 additional KEEP reads per 100K reads),
but did not produce a material final count-level benefit. The full H0+H1+H2
cache contained about 489 million records and occupied 11.7 GB, so H2 was not
adopted for routine scRNA-seq production.

This is an assay-context result, not a general conclusion about H2. Spatial
assays may have a different coverage, sparsity, and signal-to-noise tradeoff.
Spatial recipe generation should preserve H2 as an explicit experimental
option and evaluate final count recovery and specificity on representative
data before selecting a default.

The runtime path also matters. The full `classifyRead()` implementation can
consume H2 KEEP records, but the fused production triage calls
`classifyReadH0H1Offset0()` and sends misses to alignment. Therefore, generating
an H2 cache does not by itself enable H2 in the fused path; an H2 experiment
must use an H2-aware runtime path or first add and validate H2 lookup in fused
triage.

## Runtime Lookup

Sample-tag eligibility is resolved first with the fixed-offset exact and
single-mismatch tag tables (see [the sample-tag policy](FLEX_SAMPLE_TAG_POLICY.md)).
A tag outside that accepted universe is denied before the probe cache is
consulted. On the legacy H0/H1 route, eligible reads then follow this routing:

```
1.  Encode the 50-base offset-0 probe window
       → contains N / cannot encode  → PASS

2.  H0 cache lookup
       → unambiguous included probe  → KEEP (final)

3.  H1 cache lookup
       → verified positive           → KEEP
       → certified negative/ambiguous → DENY
       → no record                   → PASS
```

With an H1X2 cache, step 3 is the half-probe lookup described under H1X2, and
a read it does not assign is discarded.

`--flexNoAlign` affects only `PASS`: `0` sends it to the normal Flex alignment
and resolver (legacy; needs `--flexLegacy yes`), while `1` discards it
(`--flex yes` sets `1` unless given explicitly). It never changes H0/H1 decisions
or sample tag rejection. H2 is not consulted by the fused production path.

## Cache File Format

All tiers are serialized into a single binary file with a fixed-size header
followed by `N` records of 24 bytes each:

```
Offset  Size   Field
 0       8     magic ("FH01SEQ1")
 8       2     version (2 = sample-aware, 3 = with probe-region metadata)
10       2     k-mer length (50)
12       4     record size (24)
16       8     N (record count)
24..     24×N  records
```

Each record:

```
Offset  Size   Field
 0       8     seqLo          — low 64 bits of 2-bit encoded 50bp probe
 8       8     seqHi          — high 64 bits
16       4     resolvedGeneIdx15
20       1     cacheClass     — 0=H0_KEEP, 1=H1_KEEP, 2=DENY, 3=H2_KEEP, 4=H1X2_KEEP
21       1     negativeCode   — 0=none, 1=FlexHashNegProbeAmbig
22       2     sampleIdx      — sample index (>0 for H0, 0 for H1/H2)
```

## Generation Summary

| Tier | Variants/probe | Total variants | Stored | Compute |
|------|---------------|----------------|--------|---------|
| H0   | 16 (samples)  | 855,344        | KEEP only | Instant (no alignment) |
| H1   | 150           | 8,018,850      | KEEP + DENY | ~40s / 24 threads |
| H2   | 11,025        | ~589M          | KEEP only | ~50 min / 24 threads (est.) |
| H1X2 | 5,775 (150 + 75 × 75) | - | KEEP + DENY (multi-probe keys) | No alignment; 2 min 18 s / 16 threads for the human probe set in the 1.9.4 release notes |

H2 generation time is reduced by the pre-check optimization (skip variants
already in H0/H1 cache).

## Validation Results (2024 Reference, 53,459 Probes; legacy H0/H1 cache, pre-1.9.4)

**H1 tier breakdown:**
- KEEP: 7,124,819 (88.9%)
- DENY: 798,904 (10.0%)
- DEAD: 95,127 (1.2%, pruned)

**100K alignment smoke test (H0+H1 cache):**
- Hash screen KEEP: 80.3%
- Hash screen DENY: 0.8%
- Hash screen PASS: 18.9%

## Key Implementation Files

- `core/legacy/source/FlexHashCacheGenerate.cpp` — tier generation, variant
  enumeration, `appendVariantRecord` with tri-state pruning.
- `core/legacy/source/ReadAlign_hashCacheSynth.cpp` — synthetic alignment
  validation (`flexHashCacheValidateSyntheticPair`), returns 1/0/-1.
- `core/legacy/source/FlexHashScreen.cpp` — runtime cache loading and lookup.
- `core/legacy/source/ReadAlign.h` — `flexHashCacheValidateSyntheticPair` declaration.
