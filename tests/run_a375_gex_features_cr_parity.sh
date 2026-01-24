#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

EXTERNAL_ENV="${REPO_ROOT}/tests/external_fixtures_env.sh"
if [[ -f "${EXTERNAL_ENV}" ]]; then
  # shellcheck disable=SC1091
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

STAR_BIN="${STAR_BIN:-/mnt/pikachu/STAR-suite/core/legacy/source/STAR}"
OUTPREFIX="${A375_CR_PARITY_OUTPREFIX:-/storage/A375/star_gex_features_cr_parity_$(date +%Y%m%d_%H%M%S)/}"
THREADS="${A375_THREADS:-24}"
A375_SKIP_DOWNSAMPLE="${A375_SKIP_DOWNSAMPLE:-1}"

if [[ "${OUTPREFIX}" != */ ]]; then
  OUTPREFIX="${OUTPREFIX}/"
fi

# Validate inputs
if [[ ! -f "${WHITELIST}" ]]; then
  echo "Missing whitelist: ${WHITELIST}" >&2
  exit 1
fi

if [[ ! -d "${GENOME_DIR}" ]]; then
  echo "Missing genome index: ${GENOME_DIR}" >&2
  exit 1
fi

if [[ ! -x "${STAR_BIN}" ]]; then
  echo "STAR binary not found at ${STAR_BIN}. Build STAR first." >&2
  exit 1
fi

if [[ ! -f "${FEATURE_REF}" ]]; then
  echo "Missing feature reference: ${FEATURE_REF}" >&2
  exit 1
fi

# Setup FASTQ directories
FASTQ_DIR="${GEX_DIR}"
if [[ "${A375_SKIP_DOWNSAMPLE}" -eq 0 ]]; then
  FASTQ_DIR="${GEX_DIR}/downsampled"
fi

if [[ ! -d "${FASTQ_DIR}" ]] || ! ls "${FASTQ_DIR}"/*.fastq* >/dev/null 2>&1; then
  echo "Missing GEX FASTQ directory: ${FASTQ_DIR}" >&2
  exit 1
fi

mkdir -p "${OUTPREFIX}"
chmod 775 "${OUTPREFIX}"

GEX_R1_FILES=$(ls "${FASTQ_DIR}/"*R1*.fastq.gz | paste -sd, -)
GEX_R2_FILES=$(ls "${FASTQ_DIR}/"*R2*.fastq.gz | paste -sd, -)

# Get CRISPR FASTQ files
CRISPR_R1_FILES=$(ls "${CRISPR_DIR}/"*R1*.fastq.gz 2>/dev/null | paste -sd, - || echo "")
CRISPR_R2_FILES=$(ls "${CRISPR_DIR}/"*R2*.fastq.gz 2>/dev/null | paste -sd, - || echo "")

if [[ -z "${CRISPR_R1_FILES}" ]] || [[ -z "${CRISPR_R2_FILES}" ]]; then
  echo "Error: No CRISPR FASTQ files found in ${CRISPR_DIR}" >&2
  exit 1
fi

echo "=========================================="
echo "A375 GEX + Features CR-Parity Run"
echo "=========================================="
echo "Output prefix: ${OUTPREFIX}"
echo "GEX FASTQ dir: ${FASTQ_DIR}"
echo "CRISPR FASTQ dir: ${CRISPR_DIR}"
echo "Feature ref: ${FEATURE_REF}"
echo "Genome dir: ${GENOME_DIR}"
echo "Whitelist: ${WHITELIST}"
echo ""

# Create multi-config for GEX + Features
MULTI_CONFIG="${OUTPREFIX}multi_config.csv"
cat > "${MULTI_CONFIG}" <<EOF
[libraries]
fastqs,sample,library_type,feature_types
${FASTQ_DIR},A375,Gene Expression,Gene Expression
${CRISPR_DIR},A375,CRISPR Guide Capture,CRISPR Guide Capture

[feature]
ref,${FEATURE_REF}
EOF

echo "Multi-config created: ${MULTI_CONFIG}"
cat "${MULTI_CONFIG}"
echo ""

# Run STAR with both GEX and features
# GEX will use EmptyDrops_CR filtering
# Features will automatically use GEX filtered barcodes via --soloCrGexFeature gene
echo "Running STAR with GEX + Features..."
echo "  GEX: CR-parity mode with EmptyDrops_CR filtering"
echo "  Features: Will use GEX filtered barcodes (no EmptyDrops on features)"
echo ""

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
  --soloCbUbRequireTogether no \
  --soloCellFilter EmptyDrops_CR \
  --soloStrand Unstranded \
  --soloAddTagsToUnsorted no \
  --soloFeatures Gene \
  --soloCrGexFeature gene \
  --crMultiConfig "${MULTI_CONFIG}" \
  --crFeatureRef "${FEATURE_REF}" \
  --crWhitelist "${WHITELIST}"

echo ""
echo "=========================================="
echo "Output Summary"
echo "=========================================="

# GEX outputs
GEX_RAW_MEX="${OUTPREFIX}Solo.out/Gene/raw"
GEX_FILTERED_MEX="${OUTPREFIX}Solo.out/Gene/filtered"

echo "GEX outputs:"
echo "  Raw MEX:     ${GEX_RAW_MEX}"
echo "  Filtered MEX: ${GEX_FILTERED_MEX}"

# Features outputs (from crMultiConfig)
FEATURES_RAW_MEX="${OUTPREFIX}outs/raw_feature_bc_matrix"
FEATURES_FILTERED_MEX="${OUTPREFIX}outs/filtered_feature_bc_matrix"

echo ""
echo "Features outputs:"
echo "  Raw MEX:     ${FEATURES_RAW_MEX}"
echo "  Filtered MEX: ${FEATURES_FILTERED_MEX}"
echo "  (Filtered MEX uses GEX EmptyDrops barcodes, no EmptyDrops on features)"
echo ""

# Extract GEX filtered barcodes for reference
GEX_FILTERED_BARCODES="${OUTPREFIX}gex_filtered_barcodes.tsv"
if [[ -f "${GEX_FILTERED_MEX}/barcodes.tsv.gz" ]]; then
  zcat "${GEX_FILTERED_MEX}/barcodes.tsv.gz" > "${GEX_FILTERED_BARCODES}"
elif [[ -f "${GEX_FILTERED_MEX}/barcodes.tsv" ]]; then
  cp "${GEX_FILTERED_MEX}/barcodes.tsv" "${GEX_FILTERED_BARCODES}"
fi

if [[ -f "${GEX_FILTERED_BARCODES}" ]]; then
  echo "GEX filtered barcodes: ${GEX_FILTERED_BARCODES} ($(wc -l < "${GEX_FILTERED_BARCODES}") barcodes)"
fi
echo ""

# Step 5: Run comparison
echo "=========================================="
echo "Running Comparison"
echo "=========================================="
CR_MEX_DIR="${CR_MEX_DIR:-/storage/A375/outputs/unpacked/sample_filtered_feature_bc_matrix}"

if [[ ! -d "${CR_MEX_DIR}" ]]; then
  echo "Warning: Cell Ranger MEX directory not found: ${CR_MEX_DIR}" >&2
  echo "Skipping comparison. Set CR_MEX_DIR to run comparison."
else
  COMPARISON_OUT="${OUTPREFIX}comparison_report.txt"
  echo "Comparing STAR GEX filtered vs Cell Ranger filtered..."
  echo "Command: python3 ${REPO_ROOT}/tests/compare_a375_star_mex.py \\"
  echo "  ${CR_MEX_DIR} \\"
  echo "  ${GEX_FILTERED_MEX} \\"
  echo "  --feature-types \"Gene Expression\" \\"
  echo "  --min-counts 10 \\"
  echo "  --min-cells-pct 0.01"
  echo ""
  
  python3 "${REPO_ROOT}/tests/compare_a375_star_mex.py" \
    "${CR_MEX_DIR}" \
    "${GEX_FILTERED_MEX}" \
    --feature-types "Gene Expression" \
    --min-counts 10 \
    --min-cells-pct 0.01 \
    > "${COMPARISON_OUT}" 2>&1
  
  echo "Comparison report: ${COMPARISON_OUT}"
  echo ""
  echo "=== Correlation Results ==="
  cat "${COMPARISON_OUT}"
fi

echo ""
echo "Done!"
