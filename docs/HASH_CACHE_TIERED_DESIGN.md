# Hash Cache Tiered Design

## Overview

The Flex hash screen cache accelerates single-cell Flex alignment by
pre-computing which probe-like sequences resolve to a known gene. This avoids
full STAR alignment for the majority of reads. The cache is organized into
three tiers (H0, H1, H2) corresponding to Hamming distance from the reference
probe sequences.

## Tier Definitions

### H0 — Exact Match (Hamming 0)

- **Contents**: One record per probe × sample (53,459 probes × 16 samples = 855,344 records).
- **Record type**: KEEP only. Every exact probe match is correct by definition.
- **Sample index**: Each H0 record carries the sample index derived from the
  sample tag, enabling sample-specific fast-accept.
- **Generation**: No alignment validation needed — records are stamped directly
  from the probe list and sample whitelist.

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

## Runtime Lookup (Two-Tier)

At alignment time, each read's 50bp probe window is looked up in two stages:

```
1.  H0 cache lookup (exact match, sample-aware)
       → match found  →  KEEP  (immediate return, resolve to cached gene + sample)

2.  H1/H2 combined cache lookup (sequence match)
       → H1 DENY match      →  DENY  (early return, fast-reject)
       → H1 KEEP match      →  KEEP  (resolve to cached gene)
       → H2 KEEP match      →  KEEP  (resolve to cached gene)
       → no match            →  PASS  (fall through to full STAR alignment)
```

The H0 check is separate because H0 records carry a sample index; H1/H2
records are global (`sampleIdx=0`).

## Cache File Format

All tiers are serialized into a single binary file with a fixed-size header
followed by `N` records of 24 bytes each:

```
Offset  Size   Field
 0       8     magic (0x464C455848415348 = "FLEXHASH")
 8       4     version
12       4     reserved
16       8     N (record count)
24..     24×N  records
```

Each record:

```
Offset  Size   Field
 0       8     seqLo          — low 64 bits of 2-bit encoded 50bp probe
 8       8     seqHi          — high 64 bits
16       4     resolvedGeneIdx15
20       1     cacheClass     — 0=H0_KEEP, 1=H1_KEEP, 2=DENY, 3=H2_KEEP
21       1     negativeCode   — 0=none, 1=FlexHashNegProbeAmbig
22       2     sampleIdx      — sample index (>0 for H0, 0 for H1/H2)
```

## Generation Summary

| Tier | Variants/probe | Total variants | Stored | Compute |
|------|---------------|----------------|--------|---------|
| H0   | 16 (samples)  | 855,344        | KEEP only | Instant (no alignment) |
| H1   | 150           | 8,018,850      | KEEP + DENY | ~40s / 24 threads |
| H2   | 11,025        | ~589M          | KEEP only | ~50 min / 24 threads (est.) |

H2 generation time is reduced by the pre-check optimization (skip variants
already in H0/H1 cache).

## Validation Results (2024 Reference, 53,459 Probes)

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
