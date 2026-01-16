# Cutadapt 3.2 Parity Implementation Report

**Date**: January 3, 2026  
**Author**: STAR-Flex Development Team  
**Status**: ✅ Complete and Validated

---

## Executive Summary

We successfully implemented a faithful port of cutadapt 3.2's adapter matching algorithm in STAR-Flex, achieving **byte-identical output** with:
- Real cutadapt 3.2 (verified via Docker)
- The original JAX PE pre-trimmed files (Trim Galore 0.6.10 + cutadapt 3.2)

This implementation enables STAR-Flex to exactly reproduce trimming results from datasets processed with Trim Galore + cutadapt 3.2, closing previously observed 1-3 bp adapter-prefix residue differences.

---

## Background

### Problem Statement

When comparing STAR-Flex's integrated trimming against pre-trimmed JAX PE data (processed with Trim Galore + cutadapt 3.2), we observed small but consistent differences in adapter trimming positions. These differences were traced to algorithmic changes between cutadapt 3.2 and the current cutadapt 5.x.

### Solution

Implement an optional `--trimCutadaptCompat Cutadapt3` mode that uses cutadapt 3.2's exact algorithm for adapter matching.

---

## Technical Implementation

### Algorithm Differences: cutadapt 3.2 vs 5.x

| Aspect | Cutadapt 3.2 | Cutadapt 5.x |
|--------|--------------|--------------|
| **MATCH_SCORE** | 16384 | 1 |
| **Scoring** | match=+16384, mismatch=0, indel=-1 | match=+1, mismatch=-1, indel=-2 |
| **Best match selection** | Match count priority | Score-based with overlap checks |
| **Tie-breaking** | Lower error count wins | Overlap proximity rule |

### Key Code Changes

1. **`source/libtrim/adapter_trim.cpp`**
   - New function: `semiglobal_locate_cutadapt32()` - faithful 3.2 algorithm port
   - 3.2 scoring constants: `MATCH_SCORE_3_2 = 16384`
   - Match counting: `score_to_matches()` conversion
   - Case-insensitive comparison

2. **`source/libtrim/trim.h`**
   - Added `TrimCompatMode` enum (`TRIM_COMPAT_OFF`, `TRIM_COMPAT_CUTADAPT3`)

3. **CLI Integration**
   - New parameter: `--trimCutadaptCompat` (values: `-`, `Off`, `Cutadapt3`)

---

## Validation Results

### Test 1: Synthetic Fixtures

Generated 15 test cases using real cutadapt 3.2 (via Docker) covering:
- Exact matches, partial adapters (1-13 bp)
- Mismatches, insertions, deletions
- Min overlap filtering, case sensitivity
- Ambiguous match positions

**Result**: ✅ 15/15 tests pass

### Test 2: cutadapt 3.2 Direct Comparison

| Metric | cutadapt 3.2 | STAR-Flex (compat) | Match |
|--------|--------------|-------------------|-------|
| Pairs processed | 25,000 | 25,000 | ✅ |
| Pairs dropped | 54 | 54 | ✅ |
| R1 output | 8,331,105 bytes | 8,331,105 bytes | ✅ |
| R2 output | 8,325,923 bytes | 8,325,923 bytes | ✅ |
| **Diff lines** | - | **0** | ✅ |

**Result**: ✅ Byte-identical output

### Test 3: Trim Galore (cutadapt 5.1) Comparison

To understand remaining differences with Trim Galore, we compared cutadapt 3.2 vs Trim Galore (which uses cutadapt 5.1 internally):

| Comparison | R1 Diffs | R2 Diffs |
|------------|----------|----------|
| Trim Galore vs cutadapt 3.2 | 3 reads (0.012%) | 0 |

The 3 differing reads contain internal adapter-like sequences:

```
Read 1: Contains 'AGTCGGAAGAGC' (adapter AGATCGGAAGAGC with 1 deletion)
        cutadapt 3.2 keeps: 114 bp
        cutadapt 5.1 trims to: 58 bp

Read 2: Contains 'AGTCGGAAGAGC' (same pattern)
        cutadapt 3.2 keeps: 137 bp  
        cutadapt 5.1 trims to: 57 bp

Read 3: Contains 'AGATGGAAGAGC' (adapter with 1 substitution)
        cutadapt 3.2 keeps: 105 bp
        cutadapt 5.1 trims to: 39 bp
```

**Conclusion**: cutadapt 5.x is more aggressive at finding internal adapter-like matches. The Trim Galore wrapper itself introduces no additional differences.

### Test 4: Pre-trimmed JAX PE Data Verification

The ultimate validation: comparing against the actual pre-trimmed JAX PE files.

**Pre-trimmed file origin** (from trimming report):
- Trim Galore version: 0.6.10
- Cutadapt version: 3.2
- Adapter: AGATCGGAAGAGC (auto-detected)

| Comparison | R1 Diffs | R2 Diffs | Status |
|------------|----------|----------|--------|
| cutadapt 3.2 vs JAX PE Pre-trimmed | 0 | 0 | ✅ Identical |
| STAR-Flex (compat) vs JAX PE Pre-trimmed | 0 | 0 | ✅ Identical |

**Result**: ✅ **STAR-Flex with `--trimCutadaptCompat Cutadapt3` exactly reproduces the JAX PE pre-trimmed data**

---

## Usage

### STAR Command Line

```bash
STAR --genomeDir /path/to/genome \
     --readFilesIn R1.fq R2.fq \
     --trimCutadapt Yes \
     --trimCutadaptCompat Cutadapt3 \
     --trimCutadaptAdapter "AGATCGGAAGAGC AGATCGGAAGAGC"
```

### E2E Pipeline

```bash
TRIM_COMPAT_MODE=Cutadapt3 \
TRIM_CUTADAPT_ADAPTER="AGATCGGAAGAGC AGATCGGAAGAGC" \
./runE2E_PE_comp.sh
```

### Standalone Validation

```bash
./tools/trimvalidate/trimvalidate \
  -1 R1.fastq -2 R2.fastq \
  -o1 trimmed_R1.fastq -o2 trimmed_R2.fastq \
  --adapter-r1 AGATCGGAAGAGC --adapter-r2 AGATCGGAAGAGC \
  --quality 20 --length 20 \
  --compat Cutadapt3
```

---

## Files Modified

| File | Description |
|------|-------------|
| `source/libtrim/adapter_trim.cpp` | Full cutadapt 3.2 algorithm port |
| `source/libtrim/trim.h` | TrimCompatMode enum |
| `source/libtrim/trim.cpp` | Algorithm switching logic |
| `source/Parameters.h/.cpp` | CLI parameter handling |
| `source/parametersDefault` | Default parameter value |
| `source/ReadAlign_oneRead.cpp` | STAR integration |
| `tools/trimvalidate/trimvalidate.cpp` | Validation tool updates |
| `runE2E_PE_comp.sh` | E2E script support |

---

## New Files

| File | Description |
|------|-------------|
| `tests/fixtures/generate_cutadapt32_fixtures.py` | Docker-based fixture generation |
| `tests/fixtures/cutadapt32_adapter_fixtures.json` | 15 test fixtures |
| `tests/test_cutadapt32_parity.sh` | Automated fixture test runner |
| `plans/cutadapt_3.2_vs_5.x_delta.md` | Algorithm difference analysis |

---

## Conclusions

1. **✅ Implementation verified**: STAR-Flex cutadapt 3.2 compat mode produces byte-identical output to real cutadapt 3.2

2. **✅ Pre-trimmed data reproduced**: STAR-Flex with compat mode exactly matches the original JAX PE pre-trimmed files

3. **Root cause identified**: Differences between JAX PE pre-trimmed data and current tools are due to the cutadapt version upgrade (3.2 → 5.x), not Trim Galore wrapper behavior

4. **Behavioral difference**: cutadapt 5.x is more aggressive at finding internal adapter-like sequences (0.012% of reads affected in test dataset)

5. **Recommendation**: Use `--trimCutadaptCompat Cutadapt3` when exact reproduction of cutadapt 3.2-era trimming is required

---

## Appendix: Validation Commands

### Generate Fixtures (Docker required)

```bash
docker run --rm -v $(pwd):/work -w /work python:3.9-slim bash -c \
    "pip install -q cutadapt==3.2 && python3 tests/fixtures/generate_cutadapt32_fixtures.py"
```

### Run Fixture Tests

```bash
./tests/test_cutadapt32_parity.sh
```

### Direct cutadapt 3.2 Comparison (Docker required)

```bash
# Downsample
gzip -dc raw_R1.fq.gz | head -n 100000 > downsampled_R1.fq
gzip -dc raw_R2.fq.gz | head -n 100000 > downsampled_R2.fq

# Run cutadapt 3.2
docker run --rm -v $(pwd):/data python:3.9-slim bash -c \
    "pip install -q cutadapt==3.2 && cutadapt -j 1 -e 0.1 -q 20 -O 1 -m 20 \
     -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
     -o /data/cutadapt_R1.fq -p /data/cutadapt_R2.fq \
     /data/downsampled_R1.fq /data/downsampled_R2.fq"

# Run STAR-Flex
./tools/trimvalidate/trimvalidate \
    -1 downsampled_R1.fq -2 downsampled_R2.fq \
    -o1 starflex_R1.fq -o2 starflex_R2.fq \
    --quality 20 --length 20 \
    --adapter-r1 AGATCGGAAGAGC --adapter-r2 AGATCGGAAGAGC \
    --compat Cutadapt3

# Compare
diff cutadapt_R1.fq starflex_R1.fq && echo "R1: MATCH"
diff cutadapt_R2.fq starflex_R2.fq && echo "R2: MATCH"
```
