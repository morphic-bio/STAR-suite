#!/usr/bin/env bash
# Test script to validate CR-compat mode CRISPR feature calling integration
# This runs STAR with crMultiConfig and verifies crispr_analysis/ is generated

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Load external fixtures if available
EXTERNAL_ENV="${REPO_ROOT}/tests/external_fixtures_env.sh"
if [[ -f "${EXTERNAL_ENV}" ]]; then
  source "${EXTERNAL_ENV}"
fi

# Configuration
ROOT="${CR_MULTI_ROOT:-/storage/A375}"
FASTQ_ROOT="${CR_MULTI_FASTQ_ROOT:-${ROOT}/fastqs/1k_CRISPR_5p_gemx_fastqs}"
GEX_DIR="${CR_MULTI_GEX_DIR:-${FASTQ_ROOT}/gex}"
CRISPR_DIR="${CR_MULTI_CRISPR_DIR:-${FASTQ_ROOT}/crispr}"
FEATURE_REF="${CR_MULTI_FEATURE_REF:-${ROOT}/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv}"
WHITELIST="${CR_MULTI_WHITELIST:-${ROOT}/3M-5pgex-jan-2023.txt}"
GENOME_DIR="${CR_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
OUTPREFIX="${CR_COMPAT_TEST_OUTPREFIX:-/storage/A375/test_cr_compat_crispr_$(date +%Y%m%d_%H%M%S)/}"
THREADS="${A375_THREADS:-16}"
MIN_UMI="${CR_MIN_UMI:-10}"  # CR-compatible default

if [[ "${OUTPREFIX}" != */ ]]; then
  OUTPREFIX="${OUTPREFIX}/"
fi

echo "=========================================="
echo "CR-Compat CRISPR Feature Calling Test"
echo "=========================================="
echo "Output prefix: ${OUTPREFIX}"
echo "GEX FASTQ dir: ${GEX_DIR}"
echo "CRISPR FASTQ dir: ${CRISPR_DIR}"
echo "Feature ref: ${FEATURE_REF}"
echo "Genome dir: ${GENOME_DIR}"
echo "Min UMI threshold: ${MIN_UMI}"
echo ""

# Validate inputs
for f in "${WHITELIST}" "${FEATURE_REF}"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: Missing file: $f" >&2
    exit 1
  fi
done

for d in "${GENOME_DIR}" "${GEX_DIR}" "${CRISPR_DIR}"; do
  if [[ ! -d "$d" ]]; then
    echo "ERROR: Missing directory: $d" >&2
    exit 1
  fi
done

if [[ ! -x "${STAR_BIN}" ]]; then
  echo "ERROR: STAR binary not found at ${STAR_BIN}" >&2
  exit 1
fi

mkdir -p "${OUTPREFIX}"

# Get FASTQ files
GEX_R1_FILES=$(ls "${GEX_DIR}/"*R1*.fastq.gz | paste -sd, -)
GEX_R2_FILES=$(ls "${GEX_DIR}/"*R2*.fastq.gz | paste -sd, -)
CRISPR_R1_FILES=$(ls "${CRISPR_DIR}/"*R1*.fastq.gz | paste -sd, -)
CRISPR_R2_FILES=$(ls "${CRISPR_DIR}/"*R2*.fastq.gz | paste -sd, -)

# Create multi-config
MULTI_CONFIG="${OUTPREFIX}multi_config.csv"
cat > "${MULTI_CONFIG}" <<EOF
[libraries]
fastqs,sample,library_type,feature_types
${GEX_DIR},A375,Gene Expression,Gene Expression
${CRISPR_DIR},A375,CRISPR Guide Capture,CRISPR Guide Capture

[feature]
ref,${FEATURE_REF}
EOF

echo "Multi-config:"
cat "${MULTI_CONFIG}"
echo ""

# Run STAR with crMultiConfig
echo "Running STAR with CR-compat mode..."
"${STAR_BIN}" \
  --runThreadN "${THREADS}" \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${GEX_R2_FILES},${CRISPR_R2_FILES}" "${GEX_R1_FILES},${CRISPR_R1_FILES}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUTPREFIX}" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN gx gn sF ZG ZX sS sQ sM \
  --clipAdapterType CellRanger4 \
  --alignEndsType Local \
  --chimSegmentMin 1000000 \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 \
  --soloCBwhitelist "${WHITELIST}" \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloMultiMappers Unique \
  --soloCellFilter EmptyDrops_CR \
  --soloCbUbRequireTogether no \
  --soloStrand Unstranded \
  --soloAddTagsToUnsorted no \
  --soloFeatures GeneFull \
  --soloCrGexFeature genefull \
  --crMultiConfig "${MULTI_CONFIG}" \
  --crFeatureRef "${FEATURE_REF}" \
  --crWhitelist "${WHITELIST}" \
  --crMinUmi "${MIN_UMI}"

echo ""
echo "=========================================="
echo "Validating Outputs"
echo "=========================================="

# Check for expected outputs
EXPECTED_OUTPUTS=(
  "${OUTPREFIX}outs/raw_feature_bc_matrix/matrix.mtx.gz"
  "${OUTPREFIX}outs/filtered_feature_bc_matrix/matrix.mtx.gz"
  "${OUTPREFIX}outs/crispr_analysis/protospacer_calls_per_cell.csv"
  "${OUTPREFIX}outs/crispr_analysis/protospacer_calls_summary.csv"
  "${OUTPREFIX}outs/crispr_analysis/protospacer_umi_thresholds.csv"
)

PASS=true
for f in "${EXPECTED_OUTPUTS[@]}"; do
  if [[ -f "$f" ]]; then
    echo "  ✓ Found: $f"
  else
    echo "  ✗ Missing: $f"
    PASS=false
  fi
done

echo ""

if [[ "${PASS}" == "true" ]]; then
  echo "=========================================="
  echo "CRISPR Analysis Summary"
  echo "=========================================="
  cat "${OUTPREFIX}outs/crispr_analysis/protospacer_calls_summary.csv"
  echo ""
  echo "UMI Thresholds:"
  cat "${OUTPREFIX}outs/crispr_analysis/protospacer_umi_thresholds.csv"
  echo ""
  echo "First 10 calls:"
  head -11 "${OUTPREFIX}outs/crispr_analysis/protospacer_calls_per_cell.csv"
  echo ""
  echo "=========================================="
  echo "TEST PASSED: CR-compat CRISPR calling works!"
  echo "=========================================="
else
  echo "=========================================="
  echo "TEST FAILED: Missing expected outputs"
  echo "=========================================="
  exit 1
fi
