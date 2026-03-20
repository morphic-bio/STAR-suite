# Plan: Integrate Hash Cache Generation into genomeGenerate

**Status**: Planned (not yet implemented)  
**Date**: 2026-03-20  
**Branch**: To be implemented after hotfix/flex-index-deprecated-gtf merges

## Summary

Integrate `runFlexHashCacheGenerate` (currently `--runMode hashCacheGenerate` on
`benchmark-flex`) into the `--runMode genomeGenerate --flexGeneProbeSet` path so
a single command builds the full Flex index including the H0/H1/H2 hash cache.
Keep the standalone `--runMode hashCacheGenerate` for regeneration/customization.

## Current Architecture

Hash cache generation (`runFlexHashCacheGenerate`) requires:

1. A loaded genome index (Genome + SA + SAi, ~30GB)
2. A loaded Transcriptome (gene annotations)
3. ReadAlignChunk instances (one per thread, full ReadAlign aligner)
4. Solo/Flex parameters (CB whitelist, sample whitelist, probe list)
5. SampleDetector (maps sample tags to hash-screen sample indices for H0)

The core algorithm: for each probe pseudo-chromosome, enumerate H0/H1/H2
variant 50-mers, construct synthetic R1+R2 read pairs, run them through the
full STAR alignment + Flex gene resolver (`flexHashCacheValidateSyntheticPair`),
and record keep/deny. Results are deduped and written as binary FH01SEQ1.

## Key Design Challenge

`genomeGenerate` builds the index and exits — the index only exists on disk.
To run hash cache generation inline, we must load the freshly-built index back
after writing it. This adds ~30GB memory and ~2 min of load time but saves the
user from managing a separate step. Peak memory does NOT increase because SA
generation (~45GB) completes and frees memory before index load (~35GB).

## New Parameters (genomeGenerate mode)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--flexHashCache` | string | `no` | `yes` to auto-generate hash cache during genomeGenerate |
| `--hashCacheTiers` | string | `H0,H1,H2` | Comma-separated tiers (already on benchmark-flex) |
| `--hashCacheParentLimit` | uint32 | `0` | Cap probe count for testing (already exists) |
| `--hashCacheOutput` | string | auto | Defaults to `{genomeDir}/flex_hash_cache.bin` |

## Solo Parameters Required for Hash Cache

Must be supplied during genomeGenerate when `--flexHashCache yes`:

- `--soloCBwhitelist` (synthetic R1 construction)
- `--soloSampleWhitelist` (H0 sample index derivation)
- `--soloSampleProbes` (sample tag detection)
- `--soloCBstart`, `--soloCBlen`, `--soloUMIstart`, `--soloUMIlen` (R1 layout)
- `--soloFeatures Gene` (gene resolver surface)
- `--soloType CB_UMI_Simple` or equivalent

## Modified genomeGenerate Flow

```
Parse CLI params
  └─ --flexGeneProbeSet provided?
       ├─ no → Normal genomeGenerate
       └─ yes → FlexProbeIndex::run()
                  (filter probes, strip deprecated GTF, build hybrid FASTA+GTF)
                └─ Normal SA generation using hybrid files
                     └─ --flexHashCache yes?
                          ├─ no → Write index, exit
                          └─ yes → Validate Solo params present
                                   → genomeLoad() on freshly-built index
                                   → Create Transcriptome
                                   → runFlexHashCacheGenerate()
                                   → Write cache.bin to genomeDir/
                                   → Unload genome, exit
```

## Code Changes

### 1. STAR.cpp — genomeGenerate block

After genome generation completes and before exit, add conditional block:

```cpp
if (P.pGe.flexGeneProbe.enabled && P.pSolo.flexHashCacheEnabled) {
    Genome genomeForCache(P, P.pGe);
    genomeForCache.genomeLoad();
    Transcriptome* trForCache = new Transcriptome(P);
    runFlexHashCacheGenerate(P, genomeForCache, trForCache, nullptr);
    delete trForCache;
}
```

### 2. ParametersGenome.h — extend flexGeneProbe struct

Add `bool hashCacheEnabled` field.

### 3. Parameters.cpp — wire --flexHashCache

### 4. Genome_genomeGenerate.cpp — early validation

If hashCacheEnabled, validate required Solo params are present BEFORE SA
generation starts (fail early, not after 30 min of SA building).

### 5. FlexHashCacheGenerate.cpp — no changes needed

Already takes `Genome&` and `Transcriptome*` as inputs.

### 6. parametersDefault — document new parameter

## Standalone Mode (Preserved)

`--runMode hashCacheGenerate` continues to work for:
- Regenerating cache with different tiers without rebuilding full index
- Capping probe count for testing
- Running on a different machine with pre-built index

## Output Location

Cache goes to `{genomeDir}/flex_hash_cache.bin`. Runtime FlexHashScreen
auto-discovery already checks `{genomeDir}/` for cache files.

## Error Handling

If `--flexHashCache yes` but required Solo params missing:
- Fail early during parameter validation (before SA generation)
- Clear error message listing which params are needed
