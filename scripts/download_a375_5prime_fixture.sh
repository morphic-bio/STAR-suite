#!/usr/bin/env bash
set -euo pipefail
# Download the 10x Genomics A375 1K CRISPR 5' GEM-X public dataset.
#
# Source:  https://www.10xgenomics.com/datasets/1k-CRISPR-5p-gemx
# License: CC BY 4.0 (10x Genomics)
# Chemistry: 5' v3 GemX (TRU namespace, single-column whitelist 3M-5pgex-jan-2023)
# Libraries: GEX + CRISPR Guide Capture (11 guides, A375 melanoma cell line)
# Cells: ~1,000
# FASTQ tar: ~4.7 GB (GEX + CRISPR combined, uncompressed tar of gzipped FASTQs)
# CellRanger version used by 10x: 8.0.0
#
# Usage:
#   scripts/download_a375_5prime_fixture.sh [options]
#
# Options:
#   --outdir DIR       Output directory (default: /tmp/a375_5prime_fixture)
#   --downsample N     Create a downsampled tier with N reads per FASTQ (default: 100000)
#   --skip-downsample  Skip downsample tier creation
#   --force            Re-download even if files exist
#   -h, --help         Show help

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

DATASET_URL="https://www.10xgenomics.com/datasets/1k-CRISPR-5p-gemx"
FASTQS_TAR_URL="https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/1k_CRISPR_5p_gemx_Multiplex/1k_CRISPR_5p_gemx_Multiplex_fastqs.tar"
FEATURE_REF_URL="https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/1k_CRISPR_5p_gemx_Multiplex/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv"
FASTQS_TAR_BYTES=4961945600

OUTDIR="/tmp/a375_5prime_fixture"
DOWNSAMPLE_READS=100000
SKIP_DOWNSAMPLE=0
FORCE=0

usage() { sed -n '2,/^$/s/^# \?//p' "$0"; }

log() { printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*"; }
die() { echo "ERROR: $*" >&2; exit 2; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)          OUTDIR="$2";          shift 2 ;;
    --downsample)      DOWNSAMPLE_READS="$2"; shift 2 ;;
    --skip-downsample) SKIP_DOWNSAMPLE=1;    shift ;;
    --force)           FORCE=1;              shift ;;
    -h|--help)         usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

need_cmd curl
need_cmd tar

mkdir -p "${OUTDIR}"
OUTDIR="$(cd "${OUTDIR}" && pwd)"

FASTQS_TAR="${OUTDIR}/1k_CRISPR_5p_gemx_Multiplex_fastqs.tar"
FEATURE_REF="${OUTDIR}/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv"
FASTQ_ROOT="${OUTDIR}/fastqs"
GEX_DIR="${FASTQ_ROOT}/gex"
CRISPR_DIR="${FASTQ_ROOT}/crispr"
MANIFEST="${OUTDIR}/FIXTURE_MANIFEST.txt"

# ── Download FASTQ tar ──────────────────────────────────────────────
if [[ "${FORCE}" -eq 0 && -d "${GEX_DIR}" && -d "${CRISPR_DIR}" ]]; then
  log "FASTQs already extracted: ${FASTQ_ROOT}"
else
  if [[ "${FORCE}" -eq 0 && -f "${FASTQS_TAR}" ]]; then
    log "FASTQ tar already present: ${FASTQS_TAR}"
  else
    log "Downloading FASTQ tar (~4.7 GB)"
    curl -fSL --progress-bar -o "${FASTQS_TAR}" "${FASTQS_TAR_URL}"
  fi

  log "Extracting FASTQs"
  mkdir -p "${FASTQ_ROOT}"
  tar -xf "${FASTQS_TAR}" -C "${FASTQ_ROOT}" --strip-components=1
  log "Extracted to ${FASTQ_ROOT}"
fi

# ── Download feature reference ──────────────────────────────────────
if [[ "${FORCE}" -eq 0 && -f "${FEATURE_REF}" ]]; then
  log "Feature reference already present: ${FEATURE_REF}"
else
  log "Downloading feature reference CSV"
  curl -fSL -o "${FEATURE_REF}" "${FEATURE_REF_URL}"
fi

# ── Count reads ─────────────────────────────────────────────────────
GEX_R1="$(ls "${GEX_DIR}"/*_R1_*.fastq.gz 2>/dev/null | head -1)"
CRISPR_R1="$(ls "${CRISPR_DIR}"/*_R1_*.fastq.gz 2>/dev/null | head -1)"
if [[ -n "${GEX_R1}" ]]; then
  GEX_READS=$(zcat "${GEX_R1}" | awk 'END{print NR/4}')
else
  GEX_READS="unknown"
fi
if [[ -n "${CRISPR_R1}" ]]; then
  CRISPR_READS=$(zcat "${CRISPR_R1}" | awk 'END{print NR/4}')
else
  CRISPR_READS="unknown"
fi
log "GEX reads (lane 1 R1): ${GEX_READS}"
log "CRISPR reads (lane 1 R1): ${CRISPR_READS}"

# ── Downsample ──────────────────────────────────────────────────────
if [[ "${SKIP_DOWNSAMPLE}" -eq 0 && "${DOWNSAMPLE_READS}" -gt 0 ]]; then
  TIER_DIR="${OUTDIR}/downsampled_${DOWNSAMPLE_READS}"
  if [[ "${FORCE}" -eq 0 && -d "${TIER_DIR}/gex" && -d "${TIER_DIR}/crispr" ]]; then
    log "Downsample tier already exists: ${TIER_DIR}"
  else
    log "Creating ${DOWNSAMPLE_READS}-read downsample tier"
    MAX_LINES=$(( DOWNSAMPLE_READS * 4 ))
    for lib_name in gex crispr; do
      lib_src="${FASTQ_ROOT}/${lib_name}"
      lib_dst="${TIER_DIR}/${lib_name}"
      mkdir -p "${lib_dst}"
      for fq in "${lib_src}"/*.fastq.gz; do
        fname="$(basename "${fq}")"
        ( set +o pipefail; zcat "${fq}" | head -n "${MAX_LINES}" | gzip > "${lib_dst}/${fname}" )
      done
    done
    log "Downsample tier: ${TIER_DIR}"
  fi
fi

# ── Manifest ────────────────────────────────────────────────────────
{
  echo "A375 5' GEM-X CRISPR Public Fixture"
  echo "Generated (UTC): $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  echo "Script: ${REPO_ROOT}/scripts/download_a375_5prime_fixture.sh"
  echo "Dataset: ${DATASET_URL}"
  echo "FASTQ tar URL: ${FASTQS_TAR_URL}"
  echo "Feature ref URL: ${FEATURE_REF_URL}"
  echo "Chemistry: 5' v3 GemX (TRU)"
  echo "Libraries: GEX + CRISPR Guide Capture"
  echo "Whitelist: 3M-5pgex-jan-2023.txt"
  echo "GEX dir: ${GEX_DIR}"
  echo "CRISPR dir: ${CRISPR_DIR}"
  echo "Feature ref: ${FEATURE_REF}"
  echo "GEX reads (L001 R1): ${GEX_READS}"
  echo "CRISPR reads (L001 R1): ${CRISPR_READS}"
  [[ "${SKIP_DOWNSAMPLE}" -eq 0 ]] && echo "Downsample tier: ${TIER_DIR:-none} (${DOWNSAMPLE_READS} reads)"
} > "${MANIFEST}"

log "Wrote manifest: ${MANIFEST}"
log "Done. Output: ${OUTDIR}"
