#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

EXTERNAL_ENV="${REPO_ROOT}/tests/external_fixtures_env.sh"
if [[ -f "${EXTERNAL_ENV}" ]]; then
  # shellcheck disable=SC1091
  source "${EXTERNAL_ENV}"
fi

ROOT="${CR_MULTI_ROOT:-/storage/A375}"
FASTQ_ROOT="${CR_MULTI_FASTQ_ROOT:-${ROOT}/fastqs/1k_CRISPR_5p_gemx_fastqs}"
GEX_DIR="${CR_MULTI_GEX_DIR:-${FASTQ_ROOT}/gex}"
DOWN_SCRIPT="${CR_MULTI_DOWNSAMPLE_SCRIPT:-/mnt/pikachu/process_features/scripts/downsample_fastq_directory.sh}"
NREADS="${A375_NREADS:-10000}"
A375_SKIP_DOWNSAMPLE="${A375_SKIP_DOWNSAMPLE:-0}"

WHITELIST_GZ="${CR_MULTI_WHITELIST_GZ:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-5pgex-jan-2023.txt.gz}"
WHITELIST="${CR_MULTI_WHITELIST:-${ROOT}/3M-5pgex-jan-2023.txt}"
GENOME_DIR="${CR_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"

STAR_BIN="${STAR_BIN:-/mnt/pikachu/STAR-suite/core/legacy/source/STAR}"
OUTPREFIX="${A375_CRLIKE_OUTPREFIX:-/tmp/star_gex_cr_like/}"
THREADS="${A375_THREADS:-4}"
A375_WRITE_BAM="${A375_WRITE_BAM:-1}"
SOLO_MULTIMAPPERS="${A375_SOLO_MULTIMAPPERS:-EM}"
A375_REQUIRE_CBUB_TOGETHER="${A375_REQUIRE_CBUB_TOGETHER:-yes}"

if [[ "${OUTPREFIX}" != */ ]]; then
  OUTPREFIX="${OUTPREFIX}/"
fi

if [[ ! -x "${DOWN_SCRIPT}" ]]; then
  echo "Missing downsample script: ${DOWN_SCRIPT}" >&2
  exit 1
fi

if [[ ! -f "${WHITELIST}" ]]; then
  if [[ ! -f "${WHITELIST_GZ}" ]]; then
    echo "Missing whitelist: ${WHITELIST_GZ}" >&2
    exit 1
  fi
  zcat "${WHITELIST_GZ}" > "${WHITELIST}"
fi

if [[ ! -d "${GENOME_DIR}" ]]; then
  echo "Missing genome index: ${GENOME_DIR}" >&2
  exit 1
fi

if [[ ! -x "${STAR_BIN}" ]]; then
  echo "STAR binary not found at ${STAR_BIN}. Build STAR first." >&2
  exit 1
fi

FASTQ_DIR="${GEX_DIR}/downsampled"
if [[ "${A375_SKIP_DOWNSAMPLE}" -eq 1 ]]; then
  FASTQ_DIR="${GEX_DIR}"
else
  if [[ ! -d "${FASTQ_DIR}" ]] || ! ls "${FASTQ_DIR}"/*.fastq* >/dev/null 2>&1; then
    "${DOWN_SCRIPT}" "${NREADS}" "${GEX_DIR}"
  fi
fi

mkdir -p "${OUTPREFIX}"
chmod 775 "${OUTPREFIX}"

R1_FILES=$(ls "${FASTQ_DIR}/"*R1*.fastq.gz | paste -sd, -)
R2_FILES=$(ls "${FASTQ_DIR}/"*R2*.fastq.gz | paste -sd, -)

SAM_OUT_ARGS=()
if [[ "${A375_WRITE_BAM}" -eq 1 ]]; then
  SAM_OUT_ARGS+=(--outSAMtype BAM SortedByCoordinate)
else
  SAM_OUT_ARGS+=(--outSAMtype None)
fi

echo "Running A375 CR-like STARsolo (EM + EmptyDrops_CR)..."
"${STAR_BIN}" \
  --runThreadN "${THREADS}" \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${R2_FILES}" "${R1_FILES}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUTPREFIX}" \
  "${SAM_OUT_ARGS[@]}" \
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
  --soloMultiMappers "${SOLO_MULTIMAPPERS}" \
  --soloCbUbRequireTogether "${A375_REQUIRE_CBUB_TOGETHER}" \
  --soloCellFilter EmptyDrops_CR \
  --soloStrand Unstranded \
  --soloAddTagsToUnsorted no \
  --soloFeatures GeneFull

echo "Outputs:"
echo "  GeneFull raw:  ${OUTPREFIX}/Solo.out/GeneFull/raw"
echo "  GeneFull filt: ${OUTPREFIX}/Solo.out/GeneFull/filtered"
