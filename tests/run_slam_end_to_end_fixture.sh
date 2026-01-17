#!/bin/bash
#
# Fixture smoke test for SLAM/GRAND-SLAM end-to-end flow (100K reads).
#
# Notes:
# - Uses the small fixture FASTQ for both 0h and 6h to exercise the pipeline.
# - Auto-trim detection is disabled (fixture is too small); trims default to 0.
# - GEDI is optional; if missing, comparison falls back to fixture reference.
#
# Usage:
#   bash tests/run_slam_end_to_end_fixture.sh [options]
#
# Options:
#   --work-dir DIR        Working directory (default: test/tmp_slam_fixture_e2e_$(date +%Y%m%d_%H%M%S))
#   --overwrite           Re-run steps even if outputs exist
#   --no-cleanup          Don't clean up intermediate files
#   --no-gedi             Force fixture reference comparison even if GEDI is available
#   --help                Show this help message

set -euo pipefail

WORK_DIR="${WORK_DIR:-}"
NO_CLEANUP=0
OVERWRITE=0
THREADS=4
TRIM5="${TRIM5:-0}"
TRIM3="${TRIM3:-0}"
FORCE_NO_GEDI=0

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
        --no-gedi)
            FORCE_NO_GEDI=1
            shift
            ;;
        --help)
            grep "^#" "$0" | tail -30
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [[ -z "$WORK_DIR" ]]; then
    WORK_DIR="$PROJECT_ROOT/test/tmp_slam_fixture_e2e_$(date +%Y%m%d_%H%M%S)"
fi

STAR_BIN="${STAR_BIN:-$PROJECT_ROOT/core/legacy/source/STAR}"
FIXTURE_ROOT="${FIXTURE_ROOT:-$PROJECT_ROOT/test/fixtures/slam}"
STAR_INDEX="${STAR_INDEX:-$FIXTURE_ROOT/ref/star_index}"
WT_FASTQ="${WT_FASTQ:-$FIXTURE_ROOT/raw/slam_100000_reads_SRR32576116.fastq.gz}"
H6_FASTQ="${H6_FASTQ:-$FIXTURE_ROOT/raw/slam_100000_reads_SRR32576116.fastq.gz}"
MASK_BED="${MASK_BED:-$FIXTURE_ROOT/ref/snps.bed}"
FIXTURE_REF="${FIXTURE_REF:-$FIXTURE_ROOT/expected/fixture_ref_human.tsv.gz}"

GEDI_BIN="${GEDI_BIN:-$PROJECT_ROOT/test/tmp_slam_fixture/gedi_install/bin/gedi}"
GEDI_GENOME="${GEDI_GENOME:-$HOME/.gedi/genomic/human_fixture.oml}"

if [[ ! -x "$STAR_BIN" ]]; then
    echo "ERROR: STAR binary not found: $STAR_BIN"
    exit 1
fi
if [[ ! -d "$STAR_INDEX" ]]; then
    echo "ERROR: STAR index not found: $STAR_INDEX"
    exit 1
fi
if [[ ! -f "$WT_FASTQ" ]]; then
    echo "ERROR: fixture FASTQ not found: $WT_FASTQ"
    exit 1
fi
if [[ ! -f "$MASK_BED" ]]; then
    echo "ERROR: fixture SNP BED not found: $MASK_BED"
    exit 1
fi
if [[ ! -f "$FIXTURE_REF" ]]; then
    echo "ERROR: fixture reference not found: $FIXTURE_REF"
    exit 1
fi

if [[ -x "$GEDI_BIN" && -f "$GEDI_GENOME" ]]; then
    RUN_GEDI=1
else
    RUN_GEDI=0
fi
if [[ "$FORCE_NO_GEDI" -eq 1 ]]; then
    RUN_GEDI=0
    GEDI_STATUS="disabled (--no-gedi)"
elif [[ "$RUN_GEDI" -eq 1 ]]; then
    GEDI_STATUS="enabled"
else
    GEDI_STATUS="missing binary/genome"
fi

mkdir -p "$WORK_DIR"/{star,gedi,qc,report}
STAR_DIR="$WORK_DIR/star"
GEDI_DIR="$WORK_DIR/gedi"
QC_DIR="$WORK_DIR/qc"
REPORT_DIR="$WORK_DIR/report"
QC_PREFIX="$QC_DIR/trim_fixture"
COMPARE_LABEL="fixture reference"
COMPARE_OUT="$REPORT_DIR/compare_fixture_ref.txt"

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
echo "SLAM Fixture Smoke Test"
echo "========================================================================"
echo "Date: $(date)"
echo "Working directory: $WORK_DIR"
echo "STAR binary: $STAR_BIN"
echo "STAR index: $STAR_INDEX"
echo "Fixture FASTQ (0h): $WT_FASTQ"
echo "Fixture FASTQ (6h): $H6_FASTQ"
echo "Fixture SNP BED: $MASK_BED"
echo "Trim values: trim5p=$TRIM5, trim3p=$TRIM3"
echo "GEDI: $GEDI_STATUS"
echo "========================================================================"
echo

echo "[1/4] Running STAR-SLAM on fixture (6h) with fixed trims..."
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
        --autoTrim - \
        --slamCompatTrim5p "$TRIM5" \
        --slamCompatTrim3p "$TRIM3" \
        --slamQcReport "$QC_PREFIX" \
        > "$REPORT_DIR/star_fixture_6h.log" 2>&1
fi
echo "✓ STAR 6h run complete"
echo

echo "[2/4] Running STAR-SLAM on fixture (0h) with fixed trims..."
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
        > "$REPORT_DIR/star_fixture_0h.log" 2>&1
fi
echo "✓ STAR 0h run complete"
echo

if [[ "$RUN_GEDI" == "1" ]]; then
    echo "[3/4] Running GEDI/GRAND-SLAM on fixture BAMs..."
    ensure_bam_index "$STAR_DIR/0h_Aligned.sortedByCoord.out.bam"
    ensure_bam_index "$STAR_DIR/6h_Aligned.sortedByCoord.out.bam"

    if [[ -f "$GEDI_DIR/0h.tsv.gz" && "$OVERWRITE" -eq 0 ]]; then
        echo "✓ GEDI 0h output exists, skipping (use --overwrite to rerun)"
    else
        "$GEDI_BIN" -e Slam \
            -reads "$STAR_DIR/0h_Aligned.sortedByCoord.out.bam" \
            -genomic "$GEDI_GENOME" \
            -prefix "$GEDI_DIR/0h" \
            -trim5p "$TRIM5" \
            -trim3p "$TRIM3" \
            -strandness Sense \
            -D \
            > "$REPORT_DIR/gedi_fixture_0h.log" 2>&1
    fi

    if [[ -f "$GEDI_DIR/6h.tsv.gz" && "$OVERWRITE" -eq 0 ]]; then
        echo "✓ GEDI 6h output exists, skipping (use --overwrite to rerun)"
    else
        "$GEDI_BIN" -e Slam \
            -reads "$STAR_DIR/6h_Aligned.sortedByCoord.out.bam" \
            -genomic "$GEDI_GENOME" \
            -prefix "$GEDI_DIR/6h" \
            -trim5p "$TRIM5" \
            -trim3p "$TRIM3" \
            -strandness Sense \
            -D \
            > "$REPORT_DIR/gedi_fixture_6h.log" 2>&1
    fi
    echo "✓ GEDI runs complete"
else
    echo "[3/4] Skipping GEDI ($GEDI_STATUS)"
fi
echo

echo "[4/4] Comparing STAR vs reference..."
COMPARE_SCRIPT="$SCRIPT_DIR/slam/compare_fixture.py"
if [[ ! -f "$COMPARE_SCRIPT" ]]; then
    echo "WARNING: compare_fixture.py not found at $COMPARE_SCRIPT"
    COMPARE_LABEL="not run (compare script missing)"
    COMPARE_OUT="N/A"
else
    COMPARE_REF="$FIXTURE_REF"
    if [[ "$RUN_GEDI" == "1" && -f "$GEDI_DIR/0h.tsv.gz" ]]; then
        COMPARE_REF="$GEDI_DIR/0h.tsv.gz"
        COMPARE_LABEL="GEDI 0h"
        COMPARE_OUT="$REPORT_DIR/compare_gedi_0h.txt"
    fi
    echo "Comparing STAR vs $COMPARE_LABEL..."
    python3 "$COMPARE_SCRIPT" \
        --reference "$COMPARE_REF" \
        --test "$STAR_DIR/0h_SlamQuant.out" \
        --thresholds 20,50,100 \
        | tee "$COMPARE_OUT"
fi

echo
echo "========================================================================"
echo "SUMMARY"
echo "========================================================================"
echo "Working directory: $WORK_DIR"
echo "Trim values: trim5p=$TRIM5, trim3p=$TRIM3"
echo "STAR 0h: $STAR_DIR/0h_SlamQuant.out"
echo "STAR 6h: $STAR_DIR/6h_SlamQuant.out"
if [[ "$RUN_GEDI" == "1" ]]; then
    echo "GEDI 0h: $GEDI_DIR/0h.tsv.gz"
    echo "GEDI 6h: $GEDI_DIR/6h.tsv.gz"
else
    echo "GEDI: $GEDI_STATUS"
fi
echo "Comparison: $COMPARE_LABEL -> $COMPARE_OUT"
echo "QC report JSON: $QC_PREFIX.slam_qc.json"
echo "QC report HTML: $QC_PREFIX.slam_qc.html"
echo "Log: $LOG_FILE"
echo "========================================================================"
