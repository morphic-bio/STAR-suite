# External htslib pileup SNP-mask tool - Implementation Summary

## Overview

Implemented an external SNP mask caller using htslib pileup APIs, designed to match GEDI/STAR SNP detection semantics for fast parity testing and eventual integration into STAR.

## Deliverables

### 1. External Tool (`tools/pileup_snp/`)

**Files:**
- `pileup_snp.cpp` - Main tool implementation
- `binomial_math.cpp` - Binomial statistics (reused from STAR's libem)
- `Makefile` - Build configuration
- `README.md` - Usage documentation
- `test_parity.sh` - Parity test harness
- `IMPLEMENTATION_SUMMARY.md` - This file

**Features:**
- BED-driven pileup counting from BAM files
- Configurable counting modes (conv vs any mismatch)
- Alignment filtering (primary/secondary, MAPQ, BaseQ)
- GEDI/STAR-compatible binomial calling
- Outputs: bgzip+tabix BED mask + optional debug TSV

**Build:**
```bash
cd tools/pileup_snp
make
```

**Usage:**
```bash
./pileup_snp \
  --bam wt.bam \
  --bed candidates.bed \
  --ref reference.fa \
  --output mask.bed.gz \
  --kMode conv \
  --includeSecondary 0 \
  --pErr 0.001 \
  --pval 0.001 \
  --minTcRatio 0.3 \
  --minCov 6 \
  --minAlt 1
```

### 2. STAR Integration (Partial)

**Added Parameter:**
- `--slamSnpMaskBuildBam` - Path to BAM file for mask build (currently placeholder)

**Status:**
- Parameter parsing implemented
- Full integration deferred (external tool can be used directly)
- Future work: integrate pileup logic directly into STAR

**Current Workflow:**
1. Build mask using external tool: `tools/pileup_snp/pileup_snp ...`
2. Load mask in STAR: `--slamSnpMaskIn <mask.bed.gz>`

## Key Design Decisions

1. **BED-driven scanning** (v1): Tool scans only positions in candidate BED file for fast iteration/testing
2. **Primary-only default**: Excludes secondary/supplementary alignments by default (matches user preference)
3. **Reused binomial math**: Uses same `log_binom_tail_cdf` from STAR's libem for numeric consistency
4. **External tool first**: Enables fast parity testing without modifying STAR core

## Testing

**Parity Test Harness:**
```bash
tools/pileup_snp/test_parity.sh <bam> <ref_fa> <star_bed> <gedi_bed> <candidate_bed> [output_dir]
```

Compares:
- pileup_snp (STAR-compat mode) vs STAR mask
- pileup_snp (GEDI-compat mode) vs GEDI mask

## Next Steps

1. **Validate external tool** on real WT BAM + candidate loci
2. **Run parity tests** against GEDI and STAR outputs
3. **Optimize pileup counting** (currently manual CIGAR parsing; could use htslib pileup API more directly)
4. **Add genome-wide mode** to external tool (currently BED-only)
5. **Full STAR integration** (port proven logic into STAR, skip FASTQ alignment pre-pass)

## Files Modified

- `source/Parameters.h` - Added `buildBam` parameter
- `source/Parameters.cpp` - Added parameter parsing
- `source/parametersDefault` - Added parameter documentation
- `source/STAR.cpp` - Added placeholder for BAM-based build (currently exits with instructions)

## Files Created

- `tools/pileup_snp/pileup_snp.cpp`
- `tools/pileup_snp/binomial_math.cpp`
- `tools/pileup_snp/Makefile`
- `tools/pileup_snp/README.md`
- `tools/pileup_snp/test_parity.sh`
- `tools/pileup_snp/IMPLEMENTATION_SUMMARY.md`
