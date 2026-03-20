# Hotfix: Flex Index — Deprecated GTF Filtering, Probe gene_id, and chrBinNbits

**Date**: 2026-03-20  
**Branch**: hotfix/flex-index-deprecated-gtf  
**Affects**: `--runMode genomeGenerate --flexGeneProbeSet`, `build_filtered_reference.sh`

## Bug 1: Deprecated genes not stripped from base GTF

### Problem

When building a Flex hybrid reference (base genome + probe pseudo-chromosomes),
the base GTF was copied verbatim into the hybrid GTF without removing annotations
for genes that only have `included=FALSE` probes. These "deprecated-only" genes
create alignment ambiguity when their genomic annotations overlap with included
probe genes (e.g., CHRNA5 and CHRNA3 on chr15).

For the 2024-A reference with v1.1.0 probe set: 939 deprecated-only genes
contributed 30,927 GTF annotation lines that competed with included probe genes
during alignment.

### Fix

Both `FlexProbeIndex.cpp` (C++ genomeGenerate path) and
`build_filtered_reference.sh` (external bash pipeline) now:

1. Collect gene_ids that only appear in `included=FALSE` rows of the probe CSV
2. Filter the base GTF line-by-line, dropping any line whose `gene_id` matches
   a deprecated-only gene
3. Report stats: `deprecatedOnlyGenes` and `gtfLinesDropped` in the manifest

### Verification

- All 397 non-prefixed deprecated ENSG IDs confirmed absent from hybrid GTF
- Included genes (e.g., CHRNA5, CHRNA3) confirmed present
- Line count exact: 3,293,161 − 30,927 + 53,459 = 3,315,693

## Bug 2: Probe GTF gene_id used gene_name instead of ENSG ID

### Problem

The probe-only GTF emitted `gene_id "TSPAN6"` (the gene name) instead of
`gene_id "ENSG00000000003"` (the Ensembl ID). This caused STAR's `geneInfo.tab`
to create duplicate gene entries — one from the base GTF ("ENSG00000000003")
and one from the probe GTF ("TSPAN6"). The Flex gene resolver, which matches
against the probe list (ENSG IDs), could not find the probe pseudo-gene
entries in geneInfo.tab.

Before fix: geneInfo.tab had 56,351 entries with `MissingGeneType` probe genes.
After fix: geneInfo.tab has 38,247 entries with probe genes correctly merged
under their ENSG IDs.

### Fix

- `FlexProbeIndex.cpp`: Changed probe GTF `gene_id` attribute from `geneName`
  to `probe.gene_id`
- `make_flex_reference.sh`: Changed `$gene_name` to `$gene_id` in GTF printf

## Bug 3: genomeChrBinNbits=18 wastes 82% of genome on probe padding

### Problem

STAR pads each chromosome to a `2^genomeChrBinNbits` boundary (default 2^18 =
256 KB). With 53,459 probe pseudo-chromosomes each only 50 bp, each probe
wastes 262,094 bytes of padding — a 5,243x inflation. Total waste: ~14 GB,
inflating the genome from 3.1 GB to 17.2 GB and the SA proportionally.

### Fix

During genomeGenerate, when `--flexGeneProbeSet` is provided and the user
hasn't explicitly set `genomeChrBinNbits`, auto-compute the optimal value
following STAR's own recommendation:
`min(18, floor(log2(GenomeLength / nContigs)))`.

For 53K probes: `log2(3.1G / 53659) = 15.8` → `genomeChrBinNbits=15`.

### Impact

| Metric | Before (bits=18) | After (bits=15) |
|--------|-------------------|------------------|
| Genome padded | 17.2 GB | 4.9 GB |
| SA on disk | 27.2 GB | 25.7 GB |
| VmPeak (alignment) | 48.5 GB | 34.8 GB |
| Build time | ~21 min | ~18 min |

## Smoke Test

100K-read Flex alignment against the 2024-A index confirms:
- Mapping rates match 2020-A index (16.4% unique, 81.3% multi)
- VmPeak reduced 28% (48.5 → 34.8 GB)
- Gene assignments correct (GX/GN tags populated with ENSG IDs)

## Files Changed

- `core/legacy/source/FlexProbeIndex.cpp` — deprecated GTF filter + gene_id fix
- `core/legacy/source/FlexProbeIndex.h` — new FilterStats fields
- `core/legacy/source/Genome_genomeGenerate.cpp` — chrBinNbits auto-tuning + logging
- `flex/scripts/build_filtered_reference.sh` — deprecated GTF filter (bash)
- `flex/scripts/make_flex_reference.sh` — gene_id fix (bash)
