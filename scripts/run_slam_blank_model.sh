#!/usr/bin/env bash
set -euo pipefail

# Stage 1+2: SLAM-seq blank model and p_e estimation
# Uses Stage 0a artifacts as input

# Configuration (override via environment)
ARTIFACT_BASE=${ARTIFACT_BASE:-/mnt/pikachu/slam_blank_artifacts_20260201}
REF_FASTA=${REF_FASTA:-/storage/autoindex_110_44/cellranger_ref_cache/Homo_sapiens.GRCh38.dna.primary_assembly.fa}
TRIM5P=${TRIM5P:-11}
TRIM3P=${TRIM3P:-10}
MAX_READS=${MAX_READS:-0}  # 0 = process all reads

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BLANK_MODEL_SCRIPT="${SCRIPT_DIR}/../slam/tools/blank_model/build_blank_model.py"

echo "========================================================================"
echo "SLAM-seq Blank Model (Stage 1+2)"
echo "========================================================================"
echo "Date: $(date)"
echo "Artifact base: $ARTIFACT_BASE"
echo "Reference FASTA: $REF_FASTA"
echo "Trim: 5'=${TRIM5P}, 3'=${TRIM3P}"
echo "Max reads: $MAX_READS (0=all)"
echo "========================================================================"

# Validate inputs
if [[ ! -d "$ARTIFACT_BASE/star_internal" ]]; then
  echo "ERROR: Stage 0a artifacts not found: $ARTIFACT_BASE/star_internal"
  exit 1
fi

if [[ ! -f "$REF_FASTA" ]]; then
  echo "ERROR: Reference FASTA not found: $REF_FASTA"
  exit 1
fi

# BAM paths from Stage 0a
BLANK_BAM="$ARTIFACT_BASE/star_internal/0h_Aligned.sortedByCoord.out.bam"
QC_BLANK_BAM="$ARTIFACT_BASE/star_internal/blank_Aligned.sortedByCoord.out.bam"
LABELED_BAMS=(
  "$ARTIFACT_BASE/star_internal/6h_Aligned.sortedByCoord.out.bam"
  "$ARTIFACT_BASE/star_internal/24h_Aligned.sortedByCoord.out.bam"
)

# Validate BAMs exist
for bam in "$BLANK_BAM" "$QC_BLANK_BAM" "${LABELED_BAMS[@]}"; do
  if [[ ! -f "$bam" ]]; then
    echo "ERROR: BAM not found: $bam"
    exit 1
  fi
done

# Output directory
OUTPUT_DIR="$ARTIFACT_BASE/blank_model"
mkdir -p "$OUTPUT_DIR"

echo ""
echo "Input BAMs:"
echo "  Blank (0h): $BLANK_BAM"
echo "  QC Blank (no4su): $QC_BLANK_BAM"
echo "  Labeled: ${LABELED_BAMS[*]}"
echo "Output: $OUTPUT_DIR"
echo ""

# Check for pysam and numpy
python3 -c "import pysam; import numpy" 2>/dev/null || {
  echo "ERROR: Missing Python dependencies. Install with:"
  echo "  pip install pysam numpy"
  exit 1
}

# Run the blank model builder
echo "Running blank model builder..."
python3 "$BLANK_MODEL_SCRIPT" \
  --blank-bam "$BLANK_BAM" \
  --qc-blank-bam "$QC_BLANK_BAM" \
  --labeled-bams "${LABELED_BAMS[@]}" \
  --reference "$REF_FASTA" \
  --output-dir "$OUTPUT_DIR" \
  --trim5p "$TRIM5P" \
  --trim3p "$TRIM3P" \
  --max-reads "$MAX_READS" \
  --sample-names 6h 24h \
  --verbose

echo ""
echo "========================================================================"
echo "Blank model outputs:"
echo "========================================================================"
ls -la "$OUTPUT_DIR/"

echo ""
echo "========================================================================"
echo "Transition model summary:"
echo "========================================================================"
head -20 "$OUTPUT_DIR/blank_transition_model.tsv"

echo ""
echo "========================================================================"
echo "p_e estimates:"
echo "========================================================================"
cat "$OUTPUT_DIR/pe_estimates.tsv"

echo ""
echo "DONE. Outputs stored in: $OUTPUT_DIR"
