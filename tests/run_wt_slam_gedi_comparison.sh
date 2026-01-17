#!/bin/bash
# Run STAR-SLAM + GEDI comparison on WT/no4sU sample with new p_err methodology
# This script builds SNP mask, runs STAR-SLAM with auto-trim detection,
# extracts p_err and trims from QC JSON, runs GEDI, and compares outputs.

set -euo pipefail

# Configuration
WORK_DIR="/storage/slam_e2e_arid1a_$(date +%Y%m%d)"
STAR="/mnt/pikachu/STAR-suite/core/legacy/source/STAR"
STAR_INDEX="/storage/autoindex_110_44/bulk_index"
FASTQ="/mnt/pikachu/NW-5-21/SLAM-Seq/ARID1A-no4su_S50_R1_001.fastq.gz"
GEDI="/mnt/pikachu/STAR-suite/gedi"
GEDI_GENOME="/home/lhhung/.gedi/genomic/homo_sapiens_110_44.oml"
COMPARE_SCRIPT="/mnt/pikachu/STAR-suite/tests/slam/compare_fixture.py"

# Create directories
mkdir -p "$WORK_DIR"/{mask,star,gedi,qc,report}

echo "=== WT/no4sU SLAM-GEDI Comparison Workflow ==="
echo "Work directory: $WORK_DIR"
echo "FASTQ: $FASTQ"
echo ""

# Check if FASTQ exists
if [ ! -f "$FASTQ" ]; then
    echo "ERROR: FASTQ not found: $FASTQ"
    exit 1
fi

# === TASK 1: Check for existing mask ===
echo "=== TASK 1: Checking for existing SNP mask ==="
MASK_BED="$WORK_DIR/mask/wt.mask.bed.gz"
MASK_SUMMARY="$WORK_DIR/mask/wt.mask.summary.tsv"
MASK_FOUND=false

# Search in common work directories
for search_dir in /storage/slam_e2e_* /storage/SLAM-Seq-prod-compare-*; do
    if [ -d "$search_dir" ]; then
        for mask_file in "$search_dir"/mask/*ARID1A-no4su*.bed.gz "$search_dir"/mask/*wt*.bed.gz "$search_dir"/*ARID1A-no4su*.bed.gz; do
            if [ -f "$mask_file" ] && [ -f "${mask_file%.bed.gz}.summary.tsv" ]; then
                echo "Found existing mask: $mask_file"
                cp "$mask_file" "$MASK_BED"
                cp "${mask_file%.bed.gz}.summary.tsv" "$MASK_SUMMARY"
                MASK_FOUND=true
                echo "Copied to: $MASK_BED"
                break 2
            fi
        done
    fi
done

# === TASK 2: Build mask if not found ===
if [ "$MASK_FOUND" = false ]; then
    echo ""
    echo "=== TASK 2: Building SNP mask ==="
    echo "$FASTQ" > "$WORK_DIR/mask/wt.fofn"
    
    echo "Running STAR SNP mask build..."
    "$STAR" \
        --runThreadN 8 \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$FASTQ" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$WORK_DIR/mask/wt_" \
        --outSAMtype None \
        --slamQuantMode 1 \
        --slamSnpMaskBuildFastqs "$WORK_DIR/mask/wt.fofn" \
        --slamSnpMaskBedOut "$MASK_BED" \
        --slamSnpMaskSummaryOut "$MASK_SUMMARY" \
        --slamSnpMaskOnly 1 \
        > "$WORK_DIR/report/build_mask.log" 2>&1
    
    if [ ! -f "$MASK_BED" ]; then
        echo "ERROR: Mask build failed. Check log: $WORK_DIR/report/build_mask.log"
        exit 1
    fi
    
    echo "Mask built successfully: $MASK_BED"
    echo "Summary: $MASK_SUMMARY"
else
    echo "Using existing mask: $MASK_BED"
fi

# Display mask summary
if [ -f "$MASK_SUMMARY" ]; then
    echo ""
    echo "Mask summary:"
    head -20 "$MASK_SUMMARY"
fi

# === TASK 3: Run STAR-SLAM with auto-trim detection ===
echo ""
echo "=== TASK 3: Running STAR-SLAM with auto-trim detection ==="
QC_PREFIX="$WORK_DIR/qc/trim_wt"

echo "Running STAR-SLAM (this will compute p_err and trims)..."
"$STAR" \
    --runThreadN 8 \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "$FASTQ" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$WORK_DIR/star/wt_" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS nM MD \
    --slamQuantMode 1 \
    --slamSnpMaskIn "$MASK_BED" \
    --autoTrim variance \
    --trimScope first \
    --trimSource "$FASTQ" \
    --autoTrimDetectionReads 100000 \
    --slamQcReport "$QC_PREFIX" \
    > "$WORK_DIR/report/star_wt.log" 2>&1

if [ ! -f "$QC_PREFIX.slam_qc.json" ]; then
    echo "ERROR: QC JSON not generated. Check log: $WORK_DIR/report/star_wt.log"
    exit 1
fi

echo "STAR-SLAM completed. QC JSON: $QC_PREFIX.slam_qc.json"

# === TASK 4: Parse trims and p_err from QC JSON ===
echo ""
echo "=== TASK 4: Parsing QC JSON ==="
QC_JSON="$QC_PREFIX.slam_qc.json"

TRIM5=$(python3 -c "
import json
try:
    with open('$QC_JSON', 'r') as f:
        d = json.load(f)
    print(int(d.get('trim5p', 0)))
except Exception as e:
    print(0)
")

TRIM3=$(python3 -c "
import json
try:
    with open('$QC_JSON', 'r') as f:
        d = json.load(f)
    print(int(d.get('trim3p', 0)))
except Exception as e:
    print(0)
")

SNP_ERR_EST=$(python3 -c "
import json
try:
    with open('$QC_JSON', 'r') as f:
        d = json.load(f)
    print(f\"{d.get('snp_err_est', 0.0):.6f}\")
except Exception as e:
    print('0.000000')
")

SNP_ERR_USED=$(python3 -c "
import json
try:
    with open('$QC_JSON', 'r') as f:
        d = json.load(f)
    print(f\"{d.get('snp_err_used', 0.001):.6f}\")
except Exception as e:
    print('0.001000')
")

SNP_ERR_FALLBACK=$(python3 -c "
import json
try:
    with open('$QC_JSON', 'r') as f:
        d = json.load(f)
    print(d.get('snp_err_fallback_reason', ''))
except Exception as e:
    print('')
")

echo "Detected parameters:"
echo "  trim5p: $TRIM5"
echo "  trim3p: $TRIM3"
echo "  snp_err_est: $SNP_ERR_EST"
echo "  snp_err_used: $SNP_ERR_USED"
if [ -n "$SNP_ERR_FALLBACK" ]; then
    echo "  snp_err_fallback_reason: $SNP_ERR_FALLBACK"
fi

# === TASK 5: Run GEDI ===
echo ""
echo "=== TASK 5: Running GEDI ==="

# Index BAM if needed
BAM="$WORK_DIR/star/wt_Aligned.sortedByCoord.out.bam"
if [ ! -f "${BAM}.bai" ]; then
    echo "Indexing BAM..."
    samtools index "$BAM"
fi

echo "Running GEDI with:"
echo "  -trim5p: $TRIM5"
echo "  -trim3p: $TRIM3"
echo "  -err: $SNP_ERR_USED"

"$GEDI" -e Slam \
    -reads "$BAM" \
    -genomic "$GEDI_GENOME" \
    -prefix "$WORK_DIR/gedi/wt" \
    -trim5p "$TRIM5" \
    -trim3p "$TRIM3" \
    -err "$SNP_ERR_USED" \
    -strandness Sense \
    -D \
    > "$WORK_DIR/report/gedi_wt.log" 2>&1

if [ ! -f "$WORK_DIR/gedi/wt.tsv.gz" ]; then
    echo "ERROR: GEDI output not generated. Check log: $WORK_DIR/report/gedi_wt.log"
    exit 1
fi

echo "GEDI completed. Output: $WORK_DIR/gedi/wt.tsv.gz"

# === TASK 6: Compare STAR vs GEDI ===
echo ""
echo "=== TASK 6: Comparing STAR vs GEDI ==="

if [ ! -f "$COMPARE_SCRIPT" ]; then
    echo "ERROR: Compare script not found: $COMPARE_SCRIPT"
    exit 1
fi

python3 "$COMPARE_SCRIPT" \
    --reference "$WORK_DIR/gedi/wt.tsv.gz" \
    --test "$WORK_DIR/star/wt_SlamQuant.out" \
    --thresholds 20,50,100 \
    | tee "$WORK_DIR/report/compare_wt.txt"

# === TASK 7: Generate summary ===
echo ""
echo "=== TASK 7: Summary ==="
SUMMARY_FILE="$WORK_DIR/report/summary.txt"

{
    echo "WT/no4sU SLAM-GEDI Comparison Summary"
    echo "======================================"
    echo ""
    echo "Date: $(date)"
    echo "Work directory: $WORK_DIR"
    echo ""
    echo "Inputs:"
    echo "  FASTQ: $FASTQ"
    echo "  STAR index: $STAR_INDEX"
    echo "  GEDI genome: $GEDI_GENOME"
    echo ""
    echo "SNP Mask:"
    echo "  BED: $MASK_BED"
    echo "  Summary: $MASK_SUMMARY"
    if [ "$MASK_FOUND" = true ]; then
        echo "  Source: Reused from previous run"
    else
        echo "  Source: Built in this run"
    fi
    echo ""
    echo "STAR-SLAM Parameters:"
    echo "  trim5p: $TRIM5"
    echo "  trim3p: $TRIM3"
    echo "  snp_err_est: $SNP_ERR_EST"
    echo "  snp_err_used: $SNP_ERR_USED"
    if [ -n "$SNP_ERR_FALLBACK" ]; then
        echo "  snp_err_fallback_reason: $SNP_ERR_FALLBACK"
    fi
    echo ""
    echo "GEDI Parameters:"
    echo "  -trim5p: $TRIM5"
    echo "  -trim3p: $TRIM3"
    echo "  -err: $SNP_ERR_USED"
    echo "  -strandness: Sense"
    echo ""
    echo "Outputs:"
    echo "  STAR-SLAM: $WORK_DIR/star/wt_SlamQuant.out"
    echo "  STAR BAM: $BAM"
    echo "  GEDI TSV: $WORK_DIR/gedi/wt.tsv.gz"
    echo "  QC JSON: $QC_PREFIX.slam_qc.json"
    echo "  QC HTML: $QC_PREFIX.slam_qc.html"
    echo ""
    echo "Correlation Results:"
    echo "  See: $WORK_DIR/report/compare_wt.txt"
    echo ""
} > "$SUMMARY_FILE"

cat "$SUMMARY_FILE"
echo ""
echo "Full summary written to: $SUMMARY_FILE"
echo ""
echo "=== Workflow completed successfully ==="
