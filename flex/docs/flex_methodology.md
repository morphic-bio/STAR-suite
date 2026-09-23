# Flex Pipeline Technical Methodology

This document describes the technical implementation and data flow of the STAR-Flex inline hash pipeline.

## STAR Suite 1.9.4 default route

From STAR Suite 1.9.4, `--flex yes` assigns reads with the half-probe (H1X2)
cache and aligns nothing; a count-only command (`--outSAMtype None`,
`--outSJtype None`, `--chimSegmentMin 0`, `--soloFeatures Gene`) loads no genome
index. Several sections below
(Data Flow, Probe-Aware Gene Indexing, and the per-tag FlexFilter algorithms)
were written for the alignment route and the per-tag caller, which are now
legacy (`--flexLegacy yes`, `--soloFlexCellCaller legacy`). The default route is:

1. **Input.** Lanes are read in parallel: sequencer-delivered BGZF FASTQ with the
   in-process parallel reader (`--readFilesBgzfMode`), plain gzip through
   `--readFilesCommand` (for example rapidgzip), or CBQ (experimental).
2. **Triage, per read.** The cell barcode is extracted and corrected inline,
   the UMI is extracted, and the 8-base sample tag is read at `--soloSampleProbeOffset` and
   resolved through an exact table and then a single-mismatch table (see
   [the sample-tag policy](../../docs/FLEX_SAMPLE_TAG_POLICY.md)).
3. **Probe assignment.** The 50-base probe window is looked up in the exact H0
   table. Otherwise each 25-base half is looked up in exact and single-mismatch
   tables; the other half is scored by fast Hamming distance (2-bit XOR and
   popcount), and the read is accepted when the whole probe has at most 10
   mismatches. A read anchored by more than one probe, or otherwise ambiguous, is
   rejected; so is a read the cache does not assign. Cache entries that could
   come from more than one probe are marked ambiguous when the cache is built
   (see [the cache format notes](../../docs/FLEX_KHASH_CACHE.md)).
4. **Counting.** Each assigned read becomes a packed record in one of
   `--soloBucketCount` cell-barcode buckets, held in memory or spilled to disk
   (`--soloBucketMode`). Each bucket is sorted and its UMIs are collapsed within
   gene, and the raw matrix is written.
5. **Cell calling** with the tag-aware caller (default): one OrdMag and
   EmptyDrops model per sample label, with tags given the same label (in the
   `--soloFlexAllowedTags` file) called together as one group, then a joint GEM
   occupancy fit across all samples (see [Flex modeling features](../../docs/FLEX_MODEL_FEATURES.md)).
   A Flex cell is its 16-base barcode together with its sample tag (CB16+TAG8).
6. **Output.** Raw and per-sample filtered MEX.

Speed, memory and concordance for this route are in the top-level
[README Benchmarks](../../README.md#benchmarks).

## Architecture Overview

The flex pipeline replaces STAR's standard two-pass Solo workflow with a single-pass approach optimized for 10x Flex (Fixed RNA Profiling) samples. Key differences:

| Aspect | Standard Solo | Flex Pipeline |
|--------|---------------|---------------|
| Read storage | Temporary files | In-memory hash |
| Sample detection | Post-hoc | Per read, inline |
| CB correction | Two-pass | Inline |
| Cell filtering | CellRanger-style | Per-sample FlexFilter |
| Output | Single matrix | Raw + per-sample |

## Data Flow

This diagram shows the legacy alignment route (`--flexLegacy yes`). The default
route replaces alignment with the half-probe cache lookup described above.

```
ReadAlign::outputAlignments()
         │
         ▼
    SampleDetector::detect()  ──────────────────┐
         │                                      │
         ▼                                      │
    SoloReadFeature::record()                   │
         │                                      │
         ▼                                      │
    InlineHash (khash)  ◄───────────────────────┘
         │                                      sample_tag
         │
    ┌────┴────────────────────┐
    │  Hash Key (64-bit):     │
    │  ├─ CB index (20 bits)  │
    │  ├─ Sample tag (8 bits) │
    │  ├─ Gene index (16 bits)│
    │  └─ UMI packed (24 bits)│ (12bp, 2 bits/base)
    └─────────────────────────┘
         │
         ▼
    SoloFeature::collapseUMIall_fromHash()
         │
    ┌────┴────┐
    │         │
    ▼         ▼
CbCorrector  UMICorrector
    │         │
    └────┬────┘
         │
         ▼
    MexWriter::write()  ──► Solo.out/Gene/raw/
         │
         ▼
    FlexFilter::run()
         │
    ┌────┴────┐
    │         │
    ▼         ▼
OrdMag    EmptyDrops
(simple)  (multinomial)
    │         │
    └────┬────┘
         │
         ▼
    OccupancyGuard
    (partition filter)
         │
         ▼
    Per-sample MEX output
```

## Hash Schema

Records use a 64-bit key. In the bucket store (`core/legacy/source/CbBucketStore.h`)
the fields are, from most to least significant:

```
[CB whitelist index 20][UMI 24 (12bp, 2 bits/base)][probe gene index 15][sample tag index 5]
```

Values store the read count and flags (probe-region status). Spill files keep
this 12-byte record (8-byte key, 4-byte value); in memory a bucket's records are
compacted to 8 bytes by storing the barcode relative to its bucket.

## CB Correction

The `CbCorrector` module implements Cell Ranger-compatible CB correction:

1. **Whitelist matching**: Exact match against 737K-fixed-rna-profiling whitelist
2. **1MM correction**: Single mismatch correction with pseudocount-based resolution
3. **N-base handling**: Ns in CB treated as wildcards

Configuration: `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts`

## UMI Correction

The `UMICorrector` module implements clique-based UMI deduplication:

1. **Graph construction**: UMIs within 1 Hamming distance form edges
2. **Clique detection**: Connected components identify UMI families
3. **Representative selection**: Highest-count UMI represents the clique

Configuration: `--soloUMIdedup 1MM_CR`

UMI encoding uses 2 bits per base (A=0, C=1, G=2, T=3) packed into 24 bits for 12bp UMIs. See `UmiCodec.h`.

## Sample Tag Detection

The `SampleDetector` class identifies sample barcodes for each read:

1. **Probe sequence location**: exactly eight bases at `--soloSampleProbeOffset` (default 68); neighboring offsets are not searched (`--soloSampleSearchNearby yes` is rejected)
2. **Whitelist lookup**: exact match against `--soloSampleWhitelist`, then, after a miss, a single-mismatch table in which keys generated by more than one sample are rejected (`--soloSampleTagMismatch 1`, the default). See [the sample-tag policy](../../docs/FLEX_SAMPLE_TAG_POLICY.md)

Sample tags are stored in the hash key and used for per-sample filtering.

## Probe-Aware Gene Indexing

This section describes the legacy alignment route. On the default route the
probe gene comes from the half-probe cache entry (see
[STAR Suite 1.9.4 default route](#star-suite-194-default-route)).

Flex samples use a probe panel that maps to a subset of genes. The `GeneResolver` handles:

1. **Probe chromosome alignments**: Direct mapping from probe pseudo-chromosomes
2. **Genomic alignments**: Mapping genomic positions to probe genes
3. **Conflict resolution**: When alignments span both probe and genomic regions

The `--soloProbeList` defines the gene set; indices are 1-based probe list positions.

### Clarification: Genomic-Only Does Not Mean Probe-Less

BAM-visible genomic alignments are not, by themselves, evidence that a read was
counted. The Flex counting path only keeps a genomic alignment if its annotated
genomic gene maps back to an entry in `--soloProbeList`; otherwise that
alignment is skipped before gene resolution. In logs and debug counters,
`genomic-only` means no probe pseudo-chromosome alignment was available or won,
but the genomic evidence still resolved to a probe-list gene. Do not infer from
BAM records or these counters that Flex, or a Cell Ranger comparator, accepted
probe-less genomic assignments without checking the counted molecule or matrix
output.

## FlexFilter Cell Calling

`--soloFlexCellCaller` selects the caller. The default, `tag-aware`, builds one
model per sample label: tags given the same label are pooled, OrdMag estimates
the cell count and EmptyDrops tests candidates against a fixed ambient rank
window of 45,000-90,000 per tag, with 100,000 simulations and BH FDR 0.01 unless
`--soloFlexEdNiters` or `--soloFlexEdFdrThreshold` is set. After all samples are
called, the number of distinct tags per 16-base barcode is fitted with a
zero-truncated Poisson and barcodes above the 0.999 quantile are removed. Full
CB16+TAG8 identities are kept. Legacy expected-cell, ambient-rank and occupancy
options are rejected by this caller. Details are in
[Flex modeling features](../../docs/FLEX_MODEL_FEATURES.md).

The rest of this section describes the `legacy` caller, which runs per-sample
cell calling with two algorithms:

### OrdMag (Simple EmptyDrops)

Fast rank-based filtering for datasets with clear cell populations:

1. Sort barcodes by total UMI count descending
2. Identify knee point using robust percentile estimation
3. Call cells above the knee threshold

Used when: Low cell counts, clear separation, or as fallback.

### Full EmptyDrops

Multinomial-based statistical test for ambient vs. cell:

1. Estimate ambient RNA profile from low-UMI barcodes
2. Monte Carlo simulation of expected distribution
3. FDR-corrected p-value threshold for cell calling

Parameters (legacy caller):
- `--soloFlexEdNiters`: Simulation iterations (default: 10000)
- `--soloFlexEdFdrThreshold`: FDR threshold (default: 0.001)
- `--soloFlexEdLower`: start of the ambient rank window (default: 45000)

### Algorithm Selection

FlexFilter automatically selects the appropriate algorithm:

1. Attempt full EmptyDrops
2. If insufficient data (< minAmbient cells), fall back to OrdMag
3. If OrdMag also fails, report error

### Occupancy Guard (Post-EmptyDrops Filter)

After EmptyDrops/OrdMag cell calling, an **occupancy filter** removes cells with suspiciously high partition occupancy. This guards against GEM multiplets and barcode collisions.

**How it works:**

1. **Partition mapping**: Each cell barcode is hashed to one of 115,000 partitions (matching 10x GEM count)
2. **Occupancy estimation**: For each partition, estimate expected cell count based on:
   - Total called cells
   - Poisson lambda (cells per partition)
   - Recovery factor (1/1.65, CR default)
3. **Threshold computation**: Find 99.9th percentile of occupancy distribution
4. **Flagging**: Cells in partitions exceeding threshold are removed

**Monte Carlo simulation:**

The filter uses Poisson-based Monte Carlo simulation (default: 1M GEMs) to estimate the expected occupancy distribution. This matches Cell Ranger's production behavior and provides accurate partition statistics.

**Configuration (Cell Ranger defaults):**

```cpp
totalPartitions = 115000    // GEM partitions
recoveryFactor = 1.0/1.65   // Expected recovery
percentile = 0.999          // Stringent threshold
```

**Combined filter**: The occupancy filter runs **once after all sample tags** are processed, pooling all called cells for accurate partition statistics. This avoids per-tag bias from small sample sizes.

## Output Files

### Raw MEX (`Solo.out/Gene/raw/`)

Standard MEX format with all observed barcodes:

- `barcodes.tsv`: CB sequences (may include sample tag suffix)
- `features.tsv`: Gene IDs and names from probe list
- `matrix.mtx`: Sparse count matrix

### Per-Sample MEX (`<prefix>/<tag>/Gene/filtered/`)

One directory per sample label with filtered cells:

- Same MEX format
- Only cells passing FlexFilter
- Barcodes keep the full CB16+TAG8 identity with the tag-aware caller; the legacy caller strips the sample tag unless `--soloFlexKeepCBTag yes`

### Summary (`flexfilter_summary.tsv`)

Tab-separated statistics per sample:

```
sample_tag  total_barcodes  cells_called  median_umi  ...
BC001       15234           1245          3421        ...
BC002       12456           1102          2987        ...
```

## Configuration Reference

### Key Parameters

| Parameter | Effect |
|-----------|--------|
| `--flex yes` | Enable entire flex path (half-probe route, no alignment) |
| `--soloHashScreenFile` | Half-probe (H1X2) cache; required by `--flex yes` |
| `--flexLegacy yes` | Permit legacy alignment routes (reproducing earlier releases only) |
| `--soloSampleWhitelist` | Enable sample demultiplexing |
| `--soloFlexAllowedTags` | Sample labels used to group tags for cell calling |
| `--soloProbeList` | Enable probe-aware gene indexing; must match the cache |
| `--soloFlexExpectedCellsPerTag` | Legacy caller only: per-sample cell estimate (rejected by the tag-aware caller) |

### Recommended Settings for Flex

```bash
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
--soloUMIfiltering MultiGeneUMI_CR
--soloUMIdedup 1MM_CR
--soloMultiMappers Rescue
--soloStrand Unstranded
```

## Ambiguous CB Resolution

When a cell barcode matches multiple whitelist entries with equal Hamming distance, the `CbBayesianResolver` attempts to resolve the ambiguity using Bayesian inference:

1. **Accumulation**: Ambiguous CB observations are accumulated during read processing
2. **Merging**: Per-thread ambiguous maps are merged after mapping completes
3. **Resolution**: Bayesian inference weighs candidates by:
   - Quality scores at mismatch positions
   - UMI count evidence
   - Prior frequencies (if available)

Resolution occurs only when confidence thresholds are met:
- `MIN_POSTERIOR = 0.9` - Best candidate must have 90%+ posterior probability
- `MIN_RATIO = 2.0` - Best must be 2x more likely than runner-up

This strict threshold matches Cell Ranger's approach: prefer dropping ambiguous reads over incorrect assignment.

## Limitations

1. **BAM tag injection**: Not supported in inline flex path (use standard Solo for tagged BAMs)
2. **Memory usage**: In-memory records scale with unique CB/UMI/gene combinations; `--soloBucketMode auto` (default) spills buckets to disk once `--soloBucketMemGB` is crossed
