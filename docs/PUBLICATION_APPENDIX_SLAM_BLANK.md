# Publication Appendix: SLAM‑seq Blank + GEDI Parity (2026‑02‑02)

This appendix captures the **key reproducible commands** and **parity metrics**
used to compare **STAR‑SLAM** against **GEDI** with **SNP‑masked BAMs**.

## Dataset & Inputs

- **FASTQs (1M subsets):** `/mnt/pikachu/NW-5-21/SLAM-Seq-1M/`
- **STAR index:** `/storage/autoindex_110_44/bulk_index`
- **SNP BED mask:** `/mnt/pikachu/slam_blank_artifacts_20260201/mask/snps_from_vcf.bed.gz`
- **GEDI genome:** `/home/lhhung/.gedi/genomic/homo_sapiens_110_44.oml`
- **Output base:** `/mnt/pikachu/slam_blank_artifacts_20260201/`

## Command Lines (Reproducible)

### 1) STAR‑SLAM 1M alignments (compat mode, no trim)
```bash
OUT_BASE=/mnt/pikachu/slam_blank_artifacts_20260201 \
FASTQ_DIR=/mnt/pikachu/NW-5-21/SLAM-Seq-1M \
THREADS=24 \
./scripts/run_slam_1M_alignments.sh
```

### 2) GEDI parity (compat mode, no trim; GEDI is reference)
```bash
TRIM5_OVERRIDE=0 TRIM3_OVERRIDE=0 \
OUT_BASE=/mnt/pikachu/slam_blank_artifacts_20260201 \
./scripts/run_slam_gedi_parity.sh
```

### 3) STAR‑SLAM trimmed run (blank‑derived constants)
```bash
OUT_BASE=/mnt/pikachu/slam_blank_artifacts_20260201 \
./scripts/run_slam_1M_with_trim.sh
```

**Notes**
- GEDI comparisons use **SNP‑masked BAMs** (bedtools intersect) so STAR vs GEDI
  parity is evaluated on identical masked positions.
- Compat parity uses **no trims**; blank effects apply only to trimmed mode.

## Key Parity Results (Compat Mode, No Trim)

From `./scripts/run_slam_gedi_parity.sh` on 2026‑02‑02:

| Sample | Pearson (NTR) | Spearman (NTR) | Status |
|--------|---------------|----------------|--------|
| ARID1A_0h_1 | 0.9777 | 0.9902 | PASS |
| ARID1A_24h_1 | 0.9667 | 0.9849 | PASS |
| ARID1A_6h_1 | 0.9722 | 0.9861 | PASS |
| ARID1A_no4su | N/A | N/A | GEDI output missing MAP column |

**Interpretation:** STAR‑SLAM matches GEDI with high correlation when both use
SNP‑masked BAMs and no trimming (compat mode). The **no4su GEDI output** lacked
MAP values and was excluded from correlation.

## Auto‑Trim Effects (Blank + Trimmed Parity)

### Auto‑trim constants
Derived from `qc/trim_6h.slam_qc.json`:
- **trim5p = 11**, **trim3p = 10**

### Blank (no4su) T→C reduction
From `scripts/run_slam_1M_alignments.sh` (untrimmed) and `scripts/run_slam_1M_with_trim.sh` (trimmed):

| Sample | Untrimmed T>C | Trimmed T>C | Δ |
|--------|----------------|-------------|---|
| ARID1A_no4su | 0.2174% | 0.0573% | −0.1601% |

### Overall T→C rates (STAR, untrimmed vs trimmed)
Computed from STAR `SlamQuant.out.transitions.tsv` outputs:

| Sample | Untrimmed T→C | Trimmed T→C | Δ |
|--------|---------------|-------------|---|
| ARID1A_0h_1 | 0.9790% | 0.8061% | −0.1729% |
| ARID1A_6h_1 | 0.9157% | 0.7441% | −0.1717% |
| ARID1A_24h_1 | 1.0108% | 0.8374% | −0.1734% |
| ARID1A_no4su | 0.2174% | 0.0573% | −0.1601% |

### Trimmed parity vs GEDI
Run in a dedicated output base with trimmed BAMs:
```bash
OUT_BASE=/mnt/pikachu/slam_blank_artifacts_20260201_trimmed_parity \
./scripts/run_slam_gedi_parity.sh
```

| Sample | Pearson (NTR) | Spearman (NTR) | Status |
|--------|---------------|----------------|--------|
| ARID1A_0h_1 | 0.9550 | 0.9549 | PASS |
| ARID1A_24h_1 | 0.9461 | 0.9572 | PASS |
| ARID1A_6h_1 | 0.9557 | 0.9608 | PASS |
| ARID1A_no4su | N/A | N/A | GEDI output missing MAP column |

### Trimmed vs untrimmed STAR (MAP correlation, readcount ≥20)
Computed from STAR `SlamQuant.out`:

| Sample | Pearson (MAP) |
|--------|----------------|
| ARID1A_0h_1 | 0.9497 |
| ARID1A_6h_1 | 0.9559 |
| ARID1A_24h_1 | 0.9570 |
| ARID1A_no4su | 0.8110 |

## Dump/Requant Parity (Binary Dump → `slam_requant`)

This validates that **binary dumps** produced by STAR‑SLAM can be re‑quantified
by the external `slam_requant` tool with **exact parity** to STAR outputs.

**Run (2026‑02‑04):**
- **Output base:** `/mnt/pikachu/slam_batch_dump_20260204_054426`
- **STAR flags:** `--slamDumpBinary 1 --slamDumpWeights 1` (with batch mode + auto prefix)
- **Requant tool:** `slam/tools/slam_requant/slam_requant`

**Requant invocation (template):**
```bash
slam/tools/slam_requant/slam_requant \
  --dumpIn alignments/<sample>_slam_dump.bin \
  --dumpWeightsIn alignments/<sample>_slam_weights.bin \
  --dumpOut requant/<sample>_SlamQuant.out \
  --dumpWeightsOut requant/<sample>_SlamQuant.out.weights \
  --threads 24
```

**Parity vs STAR `SlamQuant.out` (Gene‑level overlap):**

| Sample | Genes (overlap) | Pearson (ReadCount/Conversions/Coverage/NTR/Sigma) | Spearman (same) |
|--------|-----------------|-----------------------------------------------------|-----------------|
| ARID1A_0h_1 | 12,678 | 1.0000 | 1.0000 |
| ARID1A_6h_1 | 12,573 | 1.0000 | 1.0000 |
| ARID1A_24h_1 | 12,813 | 1.0000 | 1.0000 |
| ARID1A_no4su | 12,965 | 1.0000 | 1.0000 |

## Artifact Locations

- **STAR‑SLAM outputs (1M):**
  - `/mnt/pikachu/slam_blank_artifacts_20260201/fastq_1M_runs/`
- **GEDI outputs:**
  - `/mnt/pikachu/slam_blank_artifacts_20260201/fastq_1M_runs/gedi/`
- **SNP‑masked BAMs for GEDI:**
  - `/mnt/pikachu/slam_blank_artifacts_20260201/fastq_1M_runs/gedi/bam_masked/`
