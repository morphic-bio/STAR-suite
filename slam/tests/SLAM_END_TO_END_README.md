# SLAM End-to-End Comparison (0h Mask + 6h Trim)

## Overview

This script (`run_slam_end_to_end.sh`) performs a comprehensive end-to-end comparison of STAR-SLAM and GRAND-SLAM (GEDI), demonstrating the full workflow:

1. **Build SNP mask** from 0h (time-zero/WT) sample
2. **Detect trim values** from 6h sample using auto-trim feature
3. **Run STAR-SLAM** on both 0h and 6h using the same mask and trims
4. **Run GRAND-SLAM (GEDI)** on the same BAMs with identical parameters
5. **Compare correlations** between STAR and GEDI outputs

## Files & Setup

### Input FASTQs
- **0h sample**: `/storage/SLAM-Seq-prod-compare-20260109/input/WDHD1-0h-3_S201_R1_001.fastq.gz`
- **6h sample**: `/storage/SLAM-Seq-prod-compare-20260109/input/ARID1A-6h-1_S43_R1_001.fastq.gz`

### Reference Data
- **STAR index**: `/mnt/pikachu/STAR-Flex/test/fixtures/slam/ref/star_index` (110-44 genome build)
- **GEDI genome**: `/home/lhhung/.gedi/genomic/homo_sapiens_110_44.oml`

### Binaries
- **STAR**: `/mnt/pikachu/STAR-Flex/source/STAR`
- **GEDI**: `/mnt/pikachu/STAR-Flex/gedi` (optional; comparison skipped if unavailable)

## Usage

```bash
# Basic usage (creates work directory with timestamped name)
bash /mnt/pikachu/STAR-Flex/tests/run_slam_end_to_end.sh

# Specify custom work directory
WORK_DIR=/storage/my_slam_comparison bash /mnt/pikachu/STAR-Flex/tests/run_slam_end_to_end.sh

# With custom work directory option
bash /mnt/pikachu/STAR-Flex/tests/run_slam_end_to_end.sh --work-dir /storage/my_slam_comparison
```

## Workflow Steps

### 1. Build SNP Mask (0h)
- Runs STAR-SLAM with `--slamSnpMaskBuildFastqs` and `--slamSnpMaskOnly`
- Outputs:
  - `wt0.mask.bed.gz` — SNP mask BED file (bgzip+tabix)
  - `wt0.mask.summary.tsv` — Summary statistics

### 2. Detect Trims (6h)
- Runs STAR-SLAM with `--autoTrim variance` and `--trimSource`
- Generates QC report with auto-detected trim values
- Parses `trim5p` and `trim3p` from JSON output

### 3. Apply Trims (0h)
- Re-runs STAR-SLAM on 0h using detected trim values
- Uses same SNP mask as 6h run
- Ensures both timepoints use identical trim and mask parameters

### 4. Run GEDI (Optional)
- Runs GRAND-SLAM on both BAM files with same trims
- Attempts Sense strandness first, falls back to Unspecific if needed
- Outputs gene-level SLAM quantification

### 5. Compare Results
- Runs `compare_fixture.py` to compute correlations between STAR and GEDI
- Generates Pearson and Spearman correlations for:
  - NTR (normalized transcription rate)
  - k/nT (conversion fraction)
- Multiple coverage thresholds (20, 50, 100 reads)

## Output Structure

```
$WORK_DIR/
├── mask/
│   ├── wt0.mask.bed.gz           # SNP mask (for both runs)
│   ├── wt0.mask.summary.tsv       # Mask build statistics
│   └── wt0.fofn                   # FOFN input
├── star/
│   ├── 0h_SlamQuant.out           # STAR 0h output
│   ├── 0h_Aligned.sortedByCoord.out.bam
│   ├── 6h_SlamQuant.out           # STAR 6h output
│   └── 6h_Aligned.sortedByCoord.out.bam
├── gedi/
│   ├── 0h.tsv.gz                  # GEDI 0h output
│   └── 6h.tsv.gz                  # GEDI 6h output
├── qc/
│   ├── trim_6h.slam_qc.json       # QC data (includes trim values)
│   └── trim_6h.slam_qc.html       # Interactive QC report
└── report/
    ├── run.log                    # Full execution log
    ├── build_mask.log             # Mask build stderr/stdout
    ├── star_0h.log                # 0h run log
    ├── star_6h_trim.log           # 6h trim detection log
    ├── gedi_0h.log                # GEDI 0h log
    ├── gedi_6h.log                # GEDI 6h log
    ├── compare_0h.txt             # STAR vs GEDI 0h correlation
    └── compare_6h.txt             # STAR vs GEDI 6h correlation
```

## Key Parameters

| Parameter | 0h Run | 6h Run | Purpose |
|-----------|--------|--------|---------|
| `--slamQuantMode` | 1 | 1 | Enable SLAM quantification |
| `--slamSnpMaskIn` | (build) | Shared mask | Apply SNP mask |
| `--autoTrim` | `-` (disabled) | variance | Trim detection mode |
| `--slamCompatTrim5p` | $TRIM5 | (auto) | 5' trim value |
| `--slamCompatTrim3p` | $TRIM3 | (auto) | 3' trim value |
| `--slamQcReport` | (none) | Prefix | QC output prefix |

## Expected Outputs

### SNP Mask Summary
- `wt0.mask.summary.tsv` shows:
  - Total sites observed
  - Sites passing minCov threshold
  - Sites passing posterior probability threshold
  - EM convergence statistics

### QC Report (trim_6h.slam_qc.json)
- Contains detected trim values
- Per-position statistics for all 4 plots
- Segmented regression coefficients for visualization

### Correlation Report
Example output format:

```
0h Comparison (STAR vs GEDI)
Threshold: 20 reads
  NTR Pearson r: 0.987
  NTR Spearman ρ: 0.991
  k/nT Pearson r: 0.945
  k/nT Spearman ρ: 0.964

6h Comparison (STAR vs GEDI)
Threshold: 20 reads
  NTR Pearson r: 0.989
  NTR Spearman ρ: 0.992
  k/nT Pearson r: 0.952
  k/nT Spearman ρ: 0.968
```

## Troubleshooting

### GEDI Not Found
- The script gracefully skips GEDI if binary is unavailable
- Comparison will also be skipped
- Check: `ls -l /mnt/pikachu/STAR-Flex/gedi`

### GEDI Strandness Mismatch
- Script automatically falls back from `Sense` to `Unspecific` if needed
- Both modes generate comparable outputs
- Check logs for which strandness was used

### SNP Mask Build Fails
- Verify 0h FASTQ exists and is readable
- Check disk space for temporary files
- See `report/build_mask.log` for details

### Compare Script Not Found
- Ensure working in STAR-Flex repository
- Check: `ls -l /mnt/pikachu/STAR-Flex/tests/slam/compare_fixture.py`

## Performance Notes

- **Total runtime**: ~30-60 minutes (depends on system load)
  - SNP mask build: ~10 min
  - STAR 0h run: ~10 min
  - STAR 6h run (with auto-trim): ~10 min
  - GEDI runs: ~5-10 min each
  - Comparison: <1 min
  
- **Disk space**: ~50-100 GB
  - BAM files: ~20-30 GB each
  - Temporary STAR files: ~30-50 GB

## Example Run

```bash
# Full workflow with timestamped directory
$ bash /mnt/pikachu/STAR-Flex/tests/run_slam_end_to_end.sh

# Expected output
========================================================================
SLAM/GRAND-SLAM End-to-End Comparison
========================================================================
Date: Mon Jan 12 10:30:45 UTC 2026
Working directory: /storage/slam_e2e_20260112
...

✓ SNP mask created: /storage/slam_e2e_20260112/mask/wt0.mask.bed.gz
✓ STAR 6h run complete
✓ Trim values detected:
  - trim5p: 10
  - trim3p: 8
✓ STAR 0h run complete
✓ GEDI runs complete
Comparing 0h STAR vs GEDI...
[correlation output]

========================================================================
SUMMARY
========================================================================
Working directory: /storage/slam_e2e_20260112
SNP mask: /storage/slam_e2e_20260112/mask/wt0.mask.bed.gz
Trim values (from 6h): trim5p=10, trim3p=8
...
```

## Files Generated by This Script

1. **run_slam_end_to_end.sh** — Main execution script
2. **SLAM_END_TO_END_README.md** — This documentation

## Integration with SLAM QC

The script leverages three key SLAM features:

1. **SNP Mask Building** (`--slamSnpMaskBuildFastqs`)
   - Detects static variants using 3-component EM model
   - Prevents T→C inflation from genetic variants

2. **Auto-Trim Detection** (`--autoTrim variance`)
   - Computes trim via segmented regression on T→C stdev curve
   - Generates QC report for visualization

3. **QC Report Generation** (`--slamQcReport`)
   - Interactive Plotly visualizations
   - 4-panel report: transition rates, stdev+fits, PHRED+counts, trim overlay

## References

- SLAM SNP Mask Build Plan: `.cursor/plans/slam_snp-mask_build_implementation_*.plan.md`
- SLAM QC Output Plan: (QC report generation)
- Compare Script: `tests/slam/compare_fixture.py`
- Runbook: `plans/STAR-Flex_SLAM-SNP-mask_runbook.md`

## License & Attribution

Part of STAR-Flex SLAM enhancement project, January 2026.
