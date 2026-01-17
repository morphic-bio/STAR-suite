# SLAM Fixture Trim Detection Task

**Handoff Date**: 2026-01-12  
**Prerequisites**: Review `SLAM_FIXTURE_REGENERATION_SUMMARY.md` for context

---

## Context

We now have **perfect reproduction** of the original fixture:
- GEDI on noSNP BAM with `-err 0.001` gives **r = 1.000** vs original fixture
- The regeneration pipeline is fully documented and reproducible

---

## Goals

1. **Run STAR with auto-trim detection** on the fixture to get `trim5p`/`trim3p`
   - Prefer `--autoTrim variance` with `--trimScope first`
   - Capture computed trim values from log or QC JSON

2. **Re-run STAR-SLAM** with those trims applied (or reuse same output if trims already applied)

3. **Run GEDI** on the prefiltered BAM with `-trim5p`/`-trim3p` set to STAR-detected values
   - Regenerate `fixture_ref_human.tsv.gz` to evaluate impact

4. **Compare STAR-SLAM vs GEDI** output using `compare_fixture.py`
   - Report correlations

5. **Summarize** if trim detection improves/worsens correlations and whether discrepancies are explained by trimming

---

## Constraints & Notes

- Use same fixture inputs: `slam_100000_reads_SRR32576116.fastq.gz`
- Keep same STAR alignment parameters (adapter clipping, NO EndToEnd):
  ```bash
  --clip3pAdapterSeq AGATCGGAAGAG
  --clip3pAdapterMMp 0.1
  --outSAMattributes NH HI AS nM MD
  ```
- Use existing SNP BED workflow and prefiltering step
- **Keep `-err 0.001`** for GEDI to preserve reproducibility
- Record commands/outputs in a note

---

## Useful Paths & Scripts

| Item | Path |
|------|------|
| Fixture regeneration script | `tests/regenerate_slam_fixture.sh` |
| Compare script | `tests/slam/compare_fixture.py` |
| Fixture reference | `test/fixtures/slam/expected/fixture_ref_human.tsv.gz` |
| SNP BED | `test/fixtures/slam/ref/snps.bed` |
| noSNP BAM (use this!) | `test/tmp_slam_fixture/fixture_human_Aligned.sortedByCoord.noSNP.bam` |
| STAR binary | `source/STAR` |
| GEDI binary | `gedi` |
| STAR index | `/storage/autoindex_110_44/bulk_index` |
| GEDI genome | `homo_sapiens_110_44` |
| Fixture FASTQ | `test/fixtures/slam/raw/slam_100000_reads_SRR32576116.fastq.gz` |

---

## Key Findings from Previous Work

1. **Original fixture was created from noSNP BAM** (not full BAM)
2. **GEDI requires `-err 0.001`** for reproducibility
3. **Parity test requires adapter clipping** (not EndToEnd alignment)
4. Current correlations:
   - STAR-Slam ↔ Original GEDI: **r = 0.999**
   - GEDI (noSNP BAM) ↔ Original fixture: **r = 1.000**

---

## Deliverables

1. ✅ Trim values detected (`trim5p`/`trim3p`)
2. ✅ Correlation table before/after applying trims to GEDI
3. ✅ Short explanation of any change in correlations

---

## Expected Commands (Template)

### Step 1: Run STAR with auto-trim detection
```bash
./source/STAR \
    --runThreadN 4 \
    --genomeDir /storage/autoindex_110_44/bulk_index \
    --readFilesIn test/fixtures/slam/raw/slam_100000_reads_SRR32576116.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix test/tmp_slam_fixture/trim_detect_ \
    --outSAMtype BAM SortedByCoordinate \
    --clip3pAdapterSeq AGATCGGAAGAG \
    --clip3pAdapterMMp 0.1 \
    --slamQuantMode 1 \
    --slamSnpMaskIn test/fixtures/slam/ref/snps.bed \
    --autoTrim variance \
    --trimScope first
```

### Step 2: Extract trim values
```bash
# From QC JSON or log file
# e.g., grep "trim5p\|trim3p" test/tmp_slam_fixture/trim_detect_*.json
```

### Step 3: Run GEDI with trims
```bash
./gedi -e Slam \
    -reads test/tmp_slam_fixture/fixture_human_Aligned.sortedByCoord.noSNP.bam \
    -genomic homo_sapiens_110_44 \
    -prefix test/tmp_slam_fixture/gedi_with_trims \
    -strandness AutoDetect \
    -nthreads 4 \
    -full \
    -err 0.001 \
    -trim5p <VALUE> \
    -trim3p <VALUE>
```

### Step 4: Compare
```bash
python3 tests/slam/compare_fixture.py \
    --reference test/tmp_slam_fixture/gedi_with_trims.tsv.gz \
    --test test/tmp_slam_fixture/trim_detect_SlamQuant.out
```
