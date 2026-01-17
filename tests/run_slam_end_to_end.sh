#!/bin/bash
#
# End-to-end SLAM/GRAND-SLAM comparison (0h mask + 6h trim)
#
# Goal:
# - Build SNP mask from 0h sample
# - Compute trims from 6h sample
# - Run STAR-SLAM on 0h and 6h using the same trims + mask
# - Run GRAND-SLAM (GEDI) on the same BAMs
# - Compare STAR vs GEDI correlation
#
# Usage:
#   bash run_slam_end_to_end.sh [options]
#
# Options:
#   --work-dir DIR        Working directory (default: /storage/slam_e2e_$(date +%Y%m%d))
#   --no-cleanup          Don't clean up intermediate files
#   --help                Show this help message

set -euo pipefail

# Defaults
WORK_DIR="${WORK_DIR:-/storage/slam_e2e_$(date +%Y%m%d)}"
NO_CLEANUP=0
OVERWRITE=0
THREADS=8

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --work-dir)
            WORK_DIR="$2"
            shift 2
            ;;
        --no-cleanup)
            NO_CLEANUP=1
            shift
            ;;
        --overwrite)
            OVERWRITE=1
            shift
            ;;
        --help)
            grep "^#" "$0" | tail -20
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Project paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
STAR_BIN="${PROJECT_ROOT}/core/legacy/source/STAR"
GEDI_BIN="${PROJECT_BIN:-/mnt/pikachu/STAR-suite/gedi}"

# Optional overrides (useful for A/B testing)
# Extra STAR args used ONLY during SNP mask build (step 1), e.g.:
#   STAR_MASK_BUILD_EXTRA_ARGS="--slamSnpMaskCompat gedi --slamSnpMaskKMode any --slamSnpMaskMinAlt 6"
STAR_MASK_BUILD_EXTRA_ARGS="${STAR_MASK_BUILD_EXTRA_ARGS:-}"

# Production STAR index (110-44)
STAR_INDEX="${STAR_INDEX:-/storage/autoindex_110_44/bulk_index}"

# FASTQ files (from /storage/SLAM-Seq-prod-compare-20260109/input/)
WT_FASTQ="/storage/SLAM-Seq-prod-compare-20260109/input/WDHD1-0h-3_S201_R1_001.fastq.gz"
H6_FASTQ="/storage/SLAM-Seq-prod-compare-20260109/input/ARID1A-6h-1_S43_R1_001.fastq.gz"

# GEDI genome
GEDI_GENOME="/home/lhhung/.gedi/genomic/homo_sapiens_110_44.oml"

# Verify binaries and inputs exist
if [[ ! -f "$STAR_BIN" ]]; then
    echo "ERROR: STAR binary not found: $STAR_BIN"
    exit 1
fi

if [[ ! -d "$STAR_INDEX" ]]; then
    echo "ERROR: STAR index not found: $STAR_INDEX"
    exit 1
fi

if [[ ! -f "$WT_FASTQ" ]]; then
    echo "ERROR: 0h FASTQ not found: $WT_FASTQ"
    exit 1
fi

if [[ ! -f "$H6_FASTQ" ]]; then
    echo "ERROR: 6h FASTQ not found: $H6_FASTQ"
    echo "       This script requires the SLAM-Seq 6h sample, not ATAC-Seq data."
    echo "       Please provide the correct FASTQ file path."
    exit 1
fi

if [[ ! -f "$GEDI_GENOME" ]]; then
    echo "WARNING: GEDI genome not found: $GEDI_GENOME"
    echo "         GEDI comparison will be skipped"
    RUN_GEDI=0
else
    if [[ ! -f "$GEDI_BIN" ]]; then
        echo "WARNING: GEDI binary not found: $GEDI_BIN"
        echo "         GEDI comparison will be skipped"
        RUN_GEDI=0
    else
        RUN_GEDI=1
    fi
fi

# Setup directories
mkdir -p "$WORK_DIR"/{mask,star,gedi,qc,report}
MASK_DIR="$WORK_DIR/mask"
STAR_DIR="$WORK_DIR/star"
GEDI_DIR="$WORK_DIR/gedi"
QC_DIR="$WORK_DIR/qc"
REPORT_DIR="$WORK_DIR/report"

MASK_BED="$MASK_DIR/wt0.mask.bed.gz"
MASK_SUMMARY="$MASK_DIR/wt0.mask.summary.tsv"
QC_PREFIX="$QC_DIR/trim_6h"

ensure_bam_index() {
    local bam="$1"
    if [[ -f "${bam}.bai" ]]; then
        return 0
    fi
    if ! command -v samtools >/dev/null 2>&1; then
        echo "ERROR: samtools not found; required to index BAM for GEDI"
        exit 1
    fi
    echo "Indexing BAM for GEDI: $bam"
    samtools index "$bam"
}

LOG_FILE="$REPORT_DIR/run.log"
exec > >(tee -a "$LOG_FILE")
exec 2>&1

echo "========================================================================"
echo "SLAM/GRAND-SLAM End-to-End Comparison"
echo "========================================================================"
echo "Date: $(date)"
echo "Working directory: $WORK_DIR"
echo "STAR binary: $STAR_BIN"
echo "GEDI binary: ${GEDI_BIN:-SKIPPED}"
echo "STAR index: $STAR_INDEX"
echo "0h FASTQ: $WT_FASTQ"
echo "6h FASTQ: $H6_FASTQ"
echo "========================================================================"
echo

# Step 1: Build SNP mask from 0h
echo "[1/6] Building SNP mask from 0h sample..."
echo "$WT_FASTQ" > "$MASK_DIR/wt0.fofn"

if [[ -f "$MASK_BED" && "$OVERWRITE" -eq 0 ]]; then
    echo "✓ SNP mask exists, skipping build (use --overwrite to rebuild): $MASK_BED"
else
    "$STAR_BIN" \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$WT_FASTQ" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$MASK_DIR/wt0_" \
        --outSAMtype None \
        --slamQuantMode 1 \
        --slamSnpMaskBuildFastqs "$MASK_DIR/wt0.fofn" \
        --slamSnpMaskBedOut "$MASK_BED" \
        --slamSnpMaskSummaryOut "$MASK_SUMMARY" \
        --slamSnpMaskOnly 1 \
        $STAR_MASK_BUILD_EXTRA_ARGS \
        > "$REPORT_DIR/build_mask.log" 2>&1
fi

if [[ ! -f "$MASK_BED" ]]; then
    echo "ERROR: SNP mask not created: $MASK_BED"
    exit 1
fi
echo "✓ SNP mask created: $MASK_BED"
echo

# Step 2: Run STAR-SLAM on 6h with auto-trim detection
echo "[2/6] Running STAR-SLAM on 6h with auto-trim detection..."
if [[ -f "$STAR_DIR/6h_SlamQuant.out" && -f "$STAR_DIR/6h_Aligned.sortedByCoord.out.bam" \
      && -f "$QC_PREFIX.slam_qc.json" && "$OVERWRITE" -eq 0 ]]; then
    echo "✓ STAR 6h outputs exist, skipping (use --overwrite to rerun)"
else
    "$STAR_BIN" \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$H6_FASTQ" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$STAR_DIR/6h_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM MD \
        --slamQuantMode 1 \
        --slamSnpMaskIn "$MASK_BED" \
        --autoTrim variance \
        --trimSource "$H6_FASTQ" \
        --autoTrimDetectionReads 100000 \
        --slamQcReport "$QC_PREFIX" \
        > "$REPORT_DIR/star_6h_trim.log" 2>&1
fi

if [[ ! -f "$STAR_DIR/6h_SlamQuant.out" ]]; then
    echo "ERROR: STAR 6h output not created: $STAR_DIR/6h_SlamQuant.out"
    exit 1
fi
echo "✓ STAR 6h run complete"
echo

# Step 3: Parse trims from QC JSON
echo "[3/6] Parsing trims from QC report..."
if [[ ! -f "$QC_PREFIX.slam_qc.json" ]]; then
    echo "ERROR: QC JSON not found: $QC_PREFIX.slam_qc.json"
    exit 1
fi
TRIM5=$(python3 - "$QC_PREFIX.slam_qc.json" << 'PYTHON_EOF'
import json
import sys
try:
    with open(sys.argv[1]) as f:
        data = json.load(f)
    print(int(data.get("trim5p", 0)))
except Exception as e:
    print(f"Error reading trim: {e}", file=sys.stderr)
    print(0)
PYTHON_EOF
)

TRIM3=$(python3 - "$QC_PREFIX.slam_qc.json" << 'PYTHON_EOF'
import json
import sys
try:
    with open(sys.argv[1]) as f:
        data = json.load(f)
    print(int(data.get("trim3p", 0)))
except Exception as e:
    print(f"Error reading trim: {e}", file=sys.stderr)
    print(0)
PYTHON_EOF
)

# Parse error rate from QC JSON
SNP_ERR=$(python3 - "$QC_PREFIX.slam_qc.json" << 'PYTHON_EOF'
import json
import sys
try:
    with open(sys.argv[1]) as f:
        data = json.load(f)
    # Use snp_err_used if available, otherwise fallback to default
    err_used = data.get("snp_err_used", 0.001)
    err_est = data.get("snp_err_est", 0.0)
    fallback = data.get("snp_err_fallback_reason", "")
    print(f"{err_used:.6f}")
    if fallback:
        print(f"WARNING: Error rate fallback applied: {fallback}", file=sys.stderr)
        print(f"  Estimated: {err_est:.6f}, Used: {err_used:.6f}", file=sys.stderr)
except Exception as e:
    print(f"Error reading error rate: {e}", file=sys.stderr)
    print("0.001")  # Default fallback
PYTHON_EOF
)

echo "✓ Trim values detected:"
echo "  - trim5p: $TRIM5"
echo "  - trim3p: $TRIM3"
echo "✓ Error rate detected:"
echo "  - snp_err_used: $SNP_ERR"
echo

# Step 4: Run STAR-SLAM on 0h using detected trims
echo "[4/6] Running STAR-SLAM on 0h with detected trims..."
if [[ -f "$STAR_DIR/0h_SlamQuant.out" && -f "$STAR_DIR/0h_Aligned.sortedByCoord.out.bam" \
      && "$OVERWRITE" -eq 0 ]]; then
    echo "✓ STAR 0h outputs exist, skipping (use --overwrite to rerun)"
else
    "$STAR_BIN" \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$WT_FASTQ" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$STAR_DIR/0h_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM MD \
        --slamQuantMode 1 \
        --slamSnpMaskIn "$MASK_BED" \
        --autoTrim - \
        --slamCompatTrim5p "$TRIM5" \
        --slamCompatTrim3p "$TRIM3" \
        > "$REPORT_DIR/star_0h.log" 2>&1
fi

if [[ ! -f "$STAR_DIR/0h_SlamQuant.out" ]]; then
    echo "ERROR: STAR 0h output not created: $STAR_DIR/0h_SlamQuant.out"
    exit 1
fi
echo "✓ STAR 0h run complete"
echo

# Step 5: Run GEDI if available
if [[ "$RUN_GEDI" == "1" ]]; then
    echo "[5/6] Running GEDI/GRAND-SLAM on BAMs..."
    ensure_bam_index "$STAR_DIR/0h_Aligned.sortedByCoord.out.bam"
    ensure_bam_index "$STAR_DIR/6h_Aligned.sortedByCoord.out.bam"
    
    # 0h GEDI
    if [[ -f "$GEDI_DIR/0h.tsv.gz" && "$OVERWRITE" -eq 0 ]]; then
        echo "✓ GEDI 0h output exists, skipping (use --overwrite to rerun)"
    else
        echo "Running GEDI on 0h BAM (Sense strandness)..."
        "$GEDI_BIN" -e Slam \
            -reads "$STAR_DIR/0h_Aligned.sortedByCoord.out.bam" \
            -genomic "$GEDI_GENOME" \
            -prefix "$GEDI_DIR/0h" \
            -trim5p "$TRIM5" \
            -trim3p "$TRIM3" \
            -strandness Sense \
            -err "$SNP_ERR" \
            -D \
            > "$REPORT_DIR/gedi_0h.log" 2>&1
    fi
    
    # 6h GEDI
    if [[ -f "$GEDI_DIR/6h.tsv.gz" && "$OVERWRITE" -eq 0 ]]; then
        echo "✓ GEDI 6h output exists, skipping (use --overwrite to rerun)"
    else
        echo "Running GEDI on 6h BAM (Sense strandness)..."
        "$GEDI_BIN" -e Slam \
            -reads "$STAR_DIR/6h_Aligned.sortedByCoord.out.bam" \
            -genomic "$GEDI_GENOME" \
            -prefix "$GEDI_DIR/6h" \
            -trim5p "$TRIM5" \
            -trim3p "$TRIM3" \
            -strandness Sense \
            -err "$SNP_ERR" \
            -D \
            > "$REPORT_DIR/gedi_6h.log" 2>&1
    fi
    echo "✓ GEDI runs complete"
else
    echo "[5/6] Skipping GEDI (binary not available)"
fi
echo

# Step 6: Compare STAR vs GEDI
echo "[6/6] Comparing STAR vs GEDI correlations..."

COMPARE_SCRIPT="$SCRIPT_DIR/slam/compare_fixture.py"
if [[ ! -f "$COMPARE_SCRIPT" ]]; then
    echo "WARNING: compare_fixture.py not found at $COMPARE_SCRIPT"
    echo "Skipping comparison"
else
    if [[ "$RUN_GEDI" == "1" ]]; then
        # 0h comparison
        if [[ -f "$GEDI_DIR/0h.tsv.gz" ]]; then
            echo "Comparing 0h STAR vs GEDI..."
            if [[ -f "$REPORT_DIR/compare_0h.txt" && "$OVERWRITE" -eq 0 ]]; then
                echo "✓ compare_0h.txt exists, skipping (use --overwrite to recompute)"
            else
                python3 "$COMPARE_SCRIPT" \
                    --reference "$GEDI_DIR/0h.tsv.gz" \
                    --test "$STAR_DIR/0h_SlamQuant.out" \
                    --thresholds 20,50,100 \
                    | tee "$REPORT_DIR/compare_0h.txt"
            fi
        fi
        
        # 6h comparison
        if [[ -f "$GEDI_DIR/6h.tsv.gz" ]]; then
            echo "Comparing 6h STAR vs GEDI..."
            if [[ -f "$REPORT_DIR/compare_6h.txt" && "$OVERWRITE" -eq 0 ]]; then
                echo "✓ compare_6h.txt exists, skipping (use --overwrite to recompute)"
            else
                python3 "$COMPARE_SCRIPT" \
                    --reference "$GEDI_DIR/6h.tsv.gz" \
                    --test "$STAR_DIR/6h_SlamQuant.out" \
                    --thresholds 20,50,100 \
                    | tee "$REPORT_DIR/compare_6h.txt"
            fi
        fi
    else
        echo "GEDI not run, skipping comparison"
    fi
fi
echo

# Generate summary
echo "========================================================================"
echo "SUMMARY"
echo "========================================================================"
echo "Working directory: $WORK_DIR"
echo "SNP mask: $MASK_BED"
echo "Trim values (from 6h): trim5p=$TRIM5, trim3p=$TRIM3"
echo "Error rate (from 6h): snp_err_used=$SNP_ERR"
echo ""
echo "Output files:"
echo "  - STAR 0h: $STAR_DIR/0h_SlamQuant.out"
echo "  - STAR 6h: $STAR_DIR/6h_SlamQuant.out"
if [[ "$RUN_GEDI" == "1" ]]; then
    echo "  - GEDI 0h: $GEDI_DIR/0h.tsv.gz"
    echo "  - GEDI 6h: $GEDI_DIR/6h.tsv.gz"
fi
echo ""
echo "Comparison reports:"
echo "  - 0h: $REPORT_DIR/compare_0h.txt"
echo "  - 6h: $REPORT_DIR/compare_6h.txt"
echo ""
echo "QC report:"
echo "  - JSON: $QC_PREFIX.slam_qc.json"
echo "  - HTML: $QC_PREFIX.slam_qc.html"
echo ""
echo "Full log: $LOG_FILE"
echo "========================================================================"
echo "Run completed at $(date)"
echo "========================================================================"
