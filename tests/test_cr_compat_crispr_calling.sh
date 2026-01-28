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
DISK_POLL_SEC="${CR_COMPAT_DISK_POLL_SEC:-30}"

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

# Disk usage monitor (tracks peak usage during the run)
DISK_LOG="${OUTPREFIX}disk_usage.log"
(
  while true; do
    du -sb "${OUTPREFIX}" 2>/dev/null | awk '{print strftime("%F %T"), $1}' >> "${DISK_LOG}"
    sleep "${DISK_POLL_SEC}"
  done
) &
DISK_MON_PID=$!
cleanup_disk_mon() {
  if [[ -n "${DISK_MON_PID:-}" ]]; then
    kill "${DISK_MON_PID}" >/dev/null 2>&1 || true
  fi
}
trap cleanup_disk_mon EXIT

# Get FASTQ files (robust to empty globs)
mapfile -t GEX_R1_ARR < <(ls "${GEX_DIR}/"*R1*.fastq.gz 2>/dev/null || true)
mapfile -t GEX_R2_ARR < <(ls "${GEX_DIR}/"*R2*.fastq.gz 2>/dev/null || true)
mapfile -t CRISPR_R1_ARR < <(ls "${CRISPR_DIR}/"*R1*.fastq.gz 2>/dev/null || true)
mapfile -t CRISPR_R2_ARR < <(ls "${CRISPR_DIR}/"*R2*.fastq.gz 2>/dev/null || true)

join_by_comma() {
  local IFS=,
  echo "$*"
}

GEX_R1_FILES=$(join_by_comma "${GEX_R1_ARR[@]}")
GEX_R2_FILES=$(join_by_comma "${GEX_R2_ARR[@]}")
CRISPR_R1_FILES=$(join_by_comma "${CRISPR_R1_ARR[@]}")
CRISPR_R2_FILES=$(join_by_comma "${CRISPR_R2_ARR[@]}")

# Report FASTQ lists
echo "GEX R1 files:    ${#GEX_R1_ARR[@]}"
echo "GEX R2 files:    ${#GEX_R2_ARR[@]}"
echo "CRISPR R1 files: ${#CRISPR_R1_ARR[@]}"
echo "CRISPR R2 files: ${#CRISPR_R2_ARR[@]}"
echo "GEX R1 list:     ${GEX_R1_FILES:-<empty>}"
echo "GEX R2 list:     ${GEX_R2_FILES:-<empty>}"
echo "CRISPR R1 list:  ${CRISPR_R1_FILES:-<empty>}"
echo "CRISPR R2 list:  ${CRISPR_R2_FILES:-<empty>}"
echo ""

# Fail early if any FASTQ list is empty
if [[ -z "${GEX_R1_FILES}" || -z "${GEX_R2_FILES}" || -z "${CRISPR_R1_FILES}" || -z "${CRISPR_R2_FILES}" ]]; then
  echo "ERROR: FASTQ list is empty. Check FASTQ paths:" >&2
  echo "  GEX dir:    ${GEX_DIR}" >&2
  echo "  CRISPR dir: ${CRISPR_DIR}" >&2
  echo "  GEX R1:     ${GEX_R1_FILES:-<empty>}" >&2
  echo "  GEX R2:     ${GEX_R2_FILES:-<empty>}" >&2
  echo "  CRISPR R1:  ${CRISPR_R1_FILES:-<empty>}" >&2
  echo "  CRISPR R2:  ${CRISPR_R2_FILES:-<empty>}" >&2
  exit 1
fi

# Detect malformed lists (leading/trailing/double commas)
for lbl in GEX_R1_FILES GEX_R2_FILES CRISPR_R1_FILES CRISPR_R2_FILES; do
  val="${!lbl}"
  if [[ "${val}" == *,*","* || "${val}" == ,* || "${val}" == *, ]]; then
    echo "ERROR: malformed FASTQ list for ${lbl}: ${val}" >&2
    exit 1
  fi
done

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
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN gx gn sF sS sQ sM \
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
  --soloFeatures GeneFull \
  --soloCrGexFeature genefull \
  --crMultiConfig "${MULTI_CONFIG}" \
  --crFeatureRef "${FEATURE_REF}" \
  --crWhitelist "${WHITELIST}" \
  --crMinUmi "${MIN_UMI}"

cleanup_disk_mon

if [[ -f "${DISK_LOG}" ]]; then
  PEAK_BYTES=$(awk 'max<$2{max=$2} END{print max+0}' "${DISK_LOG}")
  if [[ "${PEAK_BYTES}" -gt 0 ]]; then
    echo "Peak disk usage (bytes): ${PEAK_BYTES}"
  fi
fi

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
