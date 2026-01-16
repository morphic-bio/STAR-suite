# Flex Pipeline Technical Methodology

This document describes the technical implementation and data flow of the STAR-Flex inline hash pipeline.

## Architecture Overview

The flex pipeline replaces STAR's standard two-pass Solo workflow with a single-pass approach optimized for 10x Flex (Fixed RNA Profiling) samples. Key differences:

| Aspect | Standard Solo | Flex Pipeline |
|--------|---------------|---------------|
| Read storage | Temporary files | In-memory hash |
| Sample detection | Post-hoc | During alignment |
| CB correction | Two-pass | Inline |
| Cell filtering | CellRanger-style | Per-sample FlexFilter |
| Output | Single matrix | Raw + per-sample |

## Data Flow

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

The inline hash uses `khash` with a 64-bit key encoding:

```cpp
// Key packing (source/SoloReadFeature.h)
uint64_t key = 0;
key |= (uint64_t)cbIndex;           // bits 0-19:  CB whitelist index
key |= (uint64_t)sampleTag << 20;   // bits 20-27: Sample tag index
key |= (uint64_t)geneIndex << 28;   // bits 28-43: Probe gene index
key |= (uint64_t)packedUmi << 44;   // bits 44-67: UMI (12bp packed)
```

Values store read count and auxiliary flags.

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

The `SampleDetector` class identifies sample barcodes during alignment:

1. **Probe sequence location**: Read position defined by `--soloSampleProbeOffset`
2. **Whitelist lookup**: Match against `--soloSampleWhitelist`
3. **Nearby search**: Optional fuzzy matching with `--soloSampleSearchNearby`

Sample tags are stored in the hash key and used for per-sample filtering.

## Probe-Aware Gene Indexing

Flex samples use a probe panel that maps to a subset of genes. The `GeneResolver` handles:

1. **Probe chromosome alignments**: Direct mapping from probe pseudo-chromosomes
2. **Genomic alignments**: Mapping genomic positions to probe genes
3. **Conflict resolution**: When alignments span both probe and genomic regions

The `--soloProbeList` defines the gene set; indices are 1-based probe list positions.

## FlexFilter Cell Calling

FlexFilter runs per-sample cell calling with two algorithms:

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

Parameters:
- `--soloFlexEdNiters`: Simulation iterations (default: 10000)
- `--soloFlexEdFdrThreshold`: FDR threshold (default: 0.001)
- `--soloFlexEdLower`: UMI threshold for ambient (default: 100)

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

One directory per sample tag with filtered cells:

- Same MEX format
- Only cells passing FlexFilter
- Barcodes without sample tag suffix

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
| `--flex yes` | Enable entire flex path |
| `--soloFlexExpectedCellsPerTag` | Per-sample cell estimate (affects filtering) |
| `--soloSampleWhitelist` | Enable sample demultiplexing |
| `--soloProbeList` | Enable probe-aware gene indexing |

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
2. **Memory usage**: In-memory hash scales with unique CB/UMI/gene combinations

