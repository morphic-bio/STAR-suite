# Count Gap Diagnostics Test Matrix

**Purpose:** Systematically identify root causes of the 3,132 read count gap between STAR and Salmon.

---

## Test Configuration

### Base Parameters
- **Dataset:** PPARG_R_WT_3 (100k PE reads)
- **Genome:** `/mnt/pikachu/patrick_bwb_mnt/processing/genome`
- **STAR:** `/mnt/pikachu/STAR-Flex/source/STAR`
- **Salmon:** `salmon` (v1.10.3)

### Common STAR Parameters
```bash
--runThreadN 1
--quantMode TranscriptVB TranscriptomeSAM GeneCounts
--quantVBem 1
--outFilterMultimapNmax 20
--alignSJDBoverhangMin 1
--twopassMode Basic
```

### Common Salmon Parameters
```bash
-p 1
--incompatPrior 0
```

---

## Test 1: Fixed libType (No Auto-Detect)

### Hypothesis
Auto-detection window differences may cause reads to be processed differently, affecting counts.

### Configuration
- **STAR:** `--quantVBLibType ISF` (no auto-detect)
- **Salmon:** `-l ISF` (no auto-detect)

### Expected Outputs
| Metric | Expected Change |
|--------|----------------|
| Read gap | Should be similar (±100 reads) |
| EC count | Should be similar (±50 ECs) |
| Dropped incompat | Should match baseline |

### Commands
```bash
# STAR
STAR --quantVBLibType ISF --outFileNamePrefix test1/star_

# Salmon  
salmon quant -l ISF -a test1/star_Aligned.toTranscriptome.out.bam -o test1/salmon
```

### Analysis
```python
# Compare with baseline
gap_fixed = salmon_total - star_total
gap_baseline = 3132
if abs(gap_fixed - gap_baseline) > 200:
    print("Auto-detect affects counts!")
```

---

## Test 2: Score Filtering Disabled

### Hypothesis
STAR's score filtering may drop alignments that Salmon keeps, reducing EC counts.

### Configuration
- **STAR:** `--quantVBMinScoreFraction 0` (disable score filtering)
- **Salmon:** Default (no score filtering in alignment mode)

### Expected Outputs
| Metric | Expected Change |
|--------|----------------|
| EC count | Should increase (more alignments kept) |
| STAR-only transcripts | Should decrease |
| Read gap | May decrease if filtering was dropping valid reads |

### Commands
```bash
# STAR
STAR --quantVBLibType ISF --quantVBMinScoreFraction 0 --outFileNamePrefix test2/star_

# Salmon
salmon quant -l ISF -a test2/star_Aligned.toTranscriptome.out.bam -o test2/salmon
```

### Analysis
```python
# Check EC count change
ec_count_no_filter = grep("equivalence classes", test2/star_Log.out)
ec_count_baseline = 20025
if ec_count_no_filter > ec_count_baseline + 100:
    print("Score filtering was reducing ECs!")

# Check STAR-only transcripts
star_only_no_filter = count_star_only(test2/star_quant.sf, test2/salmon/quant.sf)
star_only_baseline = 3007
if star_only_no_filter < star_only_baseline - 500:
    print("Score filtering was creating STAR-only transcripts!")
```

---

## Test 3: Trace STAR-Only Transcripts

### Hypothesis
STAR-only transcripts exist because:
1. Reads were dropped as incompatible
2. Reads were filtered by score
3. Reads had missing mate fields
4. EC building grouped reads differently

### Configuration
- **STAR:** `--quantVBTrace 1 --quantVBTraceLimit 10000`
- **Focus:** Top 100 STAR-only transcripts from Test 2

### Expected Outputs
| Drop Reason | Expected Count |
|-------------|----------------|
| Incompatible | ~70% of STAR-only reads |
| Score filter | ~20% (if Test 2 shows improvement) |
| Missing mate | ~5% |
| Other | ~5% |

### Commands
```bash
# Identify STAR-only transcripts
python3 <<PY
import pandas as pd
star = pd.read_csv("test2/star_quant.sf", sep="\t")
salmon = pd.read_csv("test2/salmon/quant.sf", sep="\t")
merged = star.merge(salmon, on="Name", suffixes=("_star", "_salmon"))
star_only = merged[(merged['NumReads_star'] > 0.1) & (merged['NumReads_salmon'] < 0.1)]
top_100 = star_only.nlargest(100, 'NumReads_star')
top_100['Name'].to_csv('star_only_top100.txt', index=False)
PY

# Get sample reads from these transcripts
while read txp; do
    samtools view test2/star_Aligned.toTranscriptome.out.bam | grep "$txp" | head -2 | cut -f1
done < star_only_top100.txt | sort -u > sample_reads.txt

# Run STAR with trace
STAR --quantVBLibType ISF --quantVBMinScoreFraction 0 \
     --quantVBTrace 1 --quantVBTraceLimit 10000 \
     --outFileNamePrefix test3/star_
```

### Analysis
```python
# Parse trace for sample reads
sample_reads = set(open('sample_reads.txt').readlines())

dropped_reasons = {
    'incompat': 0,
    'score': 0,
    'mate_fields': 0,
    'unknown_obs_fmt': 0,
    'compatible': 0
}

with open('test3/star_quant.ec_trace.tsv') as f:
    for line in f:
        qname = line.split('\t')[0]
        if qname in sample_reads:
            # Parse drop reasons from trace
            if 'droppedIncompat=1' in line:
                dropped_reasons['incompat'] += 1
            elif 'isCompat=0' in line and 'isCompat=1' not in line:
                dropped_reasons['incompat'] += 1
            elif 'isCompat=1' in line:
                dropped_reasons['compatible'] += 1

print("Drop reason distribution:")
for reason, count in dropped_reasons.items():
    print(f"  {reason}: {count}")
```

---

## Test 4: Incompatibility Prior Comparison

### Hypothesis
Salmon may handle incompatible alignments differently when `incompatPrior=0` vs non-zero.

### Configuration
- **STAR:** `incompat_prior = 0.0` (current)
- **Salmon:** Test with `--incompatPrior 0.0001` (very small but non-zero)

### Expected Outputs
| Metric | Expected Change |
|--------|----------------|
| Dropped incompat | Should decrease slightly |
| EC count | May increase slightly |

### Commands
```bash
# Salmon with tiny incompatPrior
salmon quant -l ISF -a test1/star_Aligned.toTranscriptome.out.bam \
    -o test4/salmon --incompatPrior 0.0001
```

---

## Test 5: EC Building Comparison

### Hypothesis
Different EC building logic may group reads differently.

### Configuration
- Compare EC composition for multi-mapper groups
- Check if transcript ID ordering affects EC formation

### Analysis
```python
# Load ECs from both tools
star_ecs = load_ecs('test1/star_ecs.txt')  # Need to dump ECs
salmon_ecs = load_ecs('test1/salmon/aux_info/eq_classes.txt.gz')

# Find ECs that differ
for txp_id in [13042, 13038, 13039]:  # Known problematic group
    star_ecs_with_txp = [ec for ec in star_ecs if txp_id in ec.transcripts]
    salmon_ecs_with_txp = [ec for ec in salmon_ecs if txp_id in ec.transcripts]
    
    print(f"Transcript {txp_id}:")
    print(f"  STAR ECs: {len(star_ecs_with_txp)}")
    print(f"  Salmon ECs: {len(salmon_ecs_with_txp)}")
    
    # Compare weights
    star_weight = sum(ec.weight for ec in star_ecs_with_txp)
    salmon_weight = sum(ec.weight for ec in salmon_ecs_with_txp)
    print(f"  STAR weight: {star_weight:.1f}")
    print(f"  Salmon weight: {salmon_weight:.1f}")
```

---

## Expected Results Summary

| Test | Metric | Baseline | Expected | Pass Criteria |
|------|--------|----------|----------|---------------|
| 1 | Read gap | 3,132 | 3,132 ± 200 | Gap unchanged |
| 1 | EC count | 20,025 | 20,025 ± 50 | ECs unchanged |
| 2 | EC count | 20,025 | >20,100 | ECs increase |
| 2 | STAR-only | 3,007 | <2,500 | STAR-only decrease |
| 3 | Incompat drops | ~70% | 60-80% | Most due to incompat |
| 4 | Dropped incompat | 2,226 | <2,500 | Slight decrease |

---

## Interpretation Guide

### If Test 1 shows gap unchanged:
- Auto-detect is NOT the issue
- Focus on filtering/EC building

### If Test 2 shows EC increase:
- Score filtering was dropping valid alignments
- May explain some STAR-only transcripts

### If Test 3 shows >70% incompat drops:
- Incompatibility detection is the main issue
- Need to match Salmon's filtering logic

### If Test 4 shows different counts:
- `incompatPrior` handling differs
- May need to adjust STAR's incompat handling

---

## Quick Run Script

See `/mnt/pikachu/STAR-Flex/test_count_gap_diagnostics.sh` for automated execution.

```bash
bash /mnt/pikachu/STAR-Flex/test_count_gap_diagnostics.sh
```

Results will be saved in `/storage/production/bulk_vb_diagnostics/`

