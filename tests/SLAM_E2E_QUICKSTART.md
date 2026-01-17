# SLAM End-to-End Quick Start

## TL;DR

```bash
bash /mnt/pikachu/STAR-Flex/tests/run_slam_end_to_end.sh
```

This runs the complete SLAM workflow:
- Builds SNP mask from 0h sample
- Detects trims from 6h sample
- Runs STAR-SLAM on both with same parameters
- Compares vs GEDI (if available)

Output: `/storage/slam_e2e_YYYYMMDD/`

## What It Does

| Step | Input | Action | Output |
|------|-------|--------|--------|
| 1 | 0h FASTQ | Build SNP mask | `mask/wt0.mask.bed.gz` |
| 2 | 6h FASTQ | Auto-trim detection | `trim5p`, `trim3p` from QC JSON |
| 3 | 0h FASTQ + mask + trims | STAR-SLAM alignment | `star/0h_SlamQuant.out` |
| 4 | 6h FASTQ + mask + trims | STAR-SLAM alignment | `star/6h_SlamQuant.out` |
| 5 | BAMs | GEDI/GRAND-SLAM (optional) | `gedi/0h.tsv.gz`, `gedi/6h.tsv.gz` |
| 6 | STAR + GEDI outputs | Correlation analysis | `compare_0h.txt`, `compare_6h.txt` |

## Key Files

**Input:**
- 0h FASTQ: `/storage/SLAM-Seq-prod-compare-20260109/input/WDHD1-0h-3_S201_R1_001.fastq.gz`
- 6h FASTQ: `/storage/SLAM-Seq-prod-compare-20260109/input/ARID1A-6h-1_S43_R1_001.fastq.gz`

**Reference:**
- STAR Index: `/mnt/pikachu/STAR-Flex/test/fixtures/slam/ref/star_index` (v110-44)
- GEDI Genome: `/home/lhhung/.gedi/genomic/homo_sapiens_110_44.oml`

## Output Locations

```
/storage/slam_e2e_YYYYMMDD/
├── mask/wt0.mask.bed.gz          ← SNP mask
├── star/{0h,6h}_SlamQuant.out     ← STAR quantification
├── gedi/{0h,6h}.tsv.gz           ← GEDI quantification
├── qc/trim_6h.slam_qc.{json,html} ← QC report
└── report/
    ├── compare_{0h,6h}.txt       ← Correlations
    └── run.log                    ← Full log
```

## Expected Trim Values

From recent runs (reference):
- trim5p: 8-12 bp
- trim3p: 5-10 bp

These are specific to the sample and library prep protocol.

## Correlation Results

Expected Pearson correlation (STAR vs GEDI):
- **NTR**: 0.98-0.99
- **k/nT**: 0.94-0.96

(At 20-read threshold)

## Interpreting Results

### Good Results
- High correlation (>0.95) indicates agreement between STAR and GEDI
- Similar trim values across replicates
- Smooth QC curves without artifacts

### Potential Issues
- Low correlation (<0.90) suggests parameter mismatch or systematic error
- Trim values differ dramatically between timepoints (check library prep)
- QC curves show spikes or discontinuities (coverage artifacts)

## Troubleshooting

| Issue | Solution |
|-------|----------|
| GEDI binary not found | Install GEDI or ignore (comparison skipped) |
| Out of disk space | Use `--work-dir` to specify different location |
| Long runtime | Normal (~1 hour); check system load |
| FASTQ not found | Verify paths match current storage structure |

## Next Steps After Run

1. **Review QC Report** (interactive HTML)
   ```
   open /storage/slam_e2e_YYYYMMDD/qc/trim_6h.slam_qc.html
   ```

2. **Check Correlations**
   ```
   cat /storage/slam_e2e_YYYYMMDD/report/compare_0h.txt
   cat /storage/slam_e2e_YYYYMMDD/report/compare_6h.txt
   ```

3. **Inspect SNP Mask Summary**
   ```
   cat /storage/slam_e2e_YYYYMMDD/mask/wt0.mask.summary.tsv
   ```

## Script Location

- **Main script**: `/mnt/pikachu/STAR-Flex/tests/run_slam_end_to_end.sh`
- **Full documentation**: `/mnt/pikachu/STAR-Flex/tests/SLAM_END_TO_END_README.md`
- **Compare utility**: `/mnt/pikachu/STAR-Flex/tests/slam/compare_fixture.py`

## Performance

- **Runtime**: 30-60 minutes
- **Disk space**: 50-100 GB
- **CPU threads**: 8 (configurable via THREADS variable)

---

**Last updated**: 2026-01-12  
**STAR-Flex version**: January 2026 (SLAM SNP-Mask & QC features)
