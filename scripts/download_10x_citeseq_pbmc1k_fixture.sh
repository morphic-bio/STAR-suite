#!/usr/bin/env bash
set -euo pipefail
# Download and stage the public 10x PBMC 1k TotalSeq-B CITE-seq fixture.
#
# Source:
#   https://www.10xgenomics.com/datasets/1-k-human-pbm-cs-stained-with-a-panel-of-total-seq-b-antibodies-single-indexed-3-1-standard-4-0-0
#
# The staged layout is intended for process_features-only benchmarking:
#   fastqs/gex/        symlinks to Gene Expression FASTQs
#   fastqs/adt/        symlinks to Antibody Capture FASTQs
#   cellranger/filtered_feature_bc_matrix/
#   cellranger/raw_feature_bc_matrix/
#   SC3_v3_NextGem_SI_PBMC_CSP_1K_feature_ref.csv

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

DATASET_URL="https://www.10xgenomics.com/datasets/1-k-human-pbm-cs-stained-with-a-panel-of-total-seq-b-antibodies-single-indexed-3-1-standard-4-0-0"
SAMPLE="SC3_v3_NextGem_SI_PBMC_CSP_1K"
BASE_URL="https://cf.10xgenomics.com/samples/cell-exp/4.0.0/${SAMPLE}"
FASTQS_TAR_URL="${BASE_URL}/${SAMPLE}_fastqs.tar"
FEATURE_REF_URL="${BASE_URL}/${SAMPLE}_feature_ref.csv"
FILTERED_MEX_URL="${BASE_URL}/${SAMPLE}_filtered_feature_bc_matrix.tar.gz"
RAW_MEX_URL="${BASE_URL}/${SAMPLE}_raw_feature_bc_matrix.tar.gz"
METRICS_URL="${BASE_URL}/${SAMPLE}_metrics_summary.csv"
FASTQS_TAR_BYTES=6258360320

OUTDIR="${CITESEQ_PBMC1K_FIXTURE:-/mnt/pikachu/citeseq_pbmc1k_totalseqb_v31}"
FORCE=0
SKIP_FASTQS=0
SKIP_RAW_MEX=0

usage() {
  sed -n '2,/^$/s/^# \?//p' "$0"
  cat <<USAGE

Usage:
  scripts/download_10x_citeseq_pbmc1k_fixture.sh [options]

Options:
  --outdir DIR       Output root. Default: ${OUTDIR}
  --skip-fastqs      Download only metadata and Cell Ranger MEX outputs
  --skip-raw-mex     Skip the official raw_feature_bc_matrix tarball
  --force            Re-download and re-stage assets
  -h, --help         Show help
USAGE
}

log() { printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*"; }
die() { echo "ERROR: $*" >&2; exit 2; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)       OUTDIR="$2"; shift 2 ;;
    --skip-fastqs)  SKIP_FASTQS=1; shift ;;
    --skip-raw-mex) SKIP_RAW_MEX=1; shift ;;
    --force)        FORCE=1; shift ;;
    -h|--help)      usage; exit 0 ;;
    *)              die "Unknown argument: $1" ;;
  esac
done

need_cmd curl
need_cmd tar

mkdir -p "${OUTDIR}"
OUTDIR="$(cd "${OUTDIR}" && pwd)"

FASTQS_TAR="${OUTDIR}/${SAMPLE}_fastqs.tar"
FEATURE_REF="${OUTDIR}/${SAMPLE}_feature_ref.csv"
FILTERED_TAR="${OUTDIR}/${SAMPLE}_filtered_feature_bc_matrix.tar.gz"
RAW_TAR="${OUTDIR}/${SAMPLE}_raw_feature_bc_matrix.tar.gz"
METRICS_CSV="${OUTDIR}/${SAMPLE}_metrics_summary.csv"
RAW_FASTQ_ROOT="${OUTDIR}/fastqs_raw"
FASTQ_LINK_ROOT="${OUTDIR}/fastqs"
CELLRANGER_ROOT="${OUTDIR}/cellranger"
MANIFEST="${OUTDIR}/FIXTURE_MANIFEST.txt"

download_file() {
  local url="$1"
  local dest="$2"
  local label="$3"
  if [[ "${FORCE}" -eq 0 && -s "${dest}" ]]; then
    log "${label} already present: ${dest}"
    return
  fi
  log "Downloading ${label}: ${url}"
  curl -fSL --progress-bar -o "${dest}" "${url}"
}

extract_tar_once() {
  local tar_path="$1"
  local dest_dir="$2"
  local marker="$3"
  local label="$4"
  if [[ "${FORCE}" -eq 0 && -f "${marker}" ]]; then
    log "${label} already extracted: ${dest_dir}"
    return
  fi
  rm -rf "${dest_dir}"
  mkdir -p "${dest_dir}"
  log "Extracting ${label}"
  case "${tar_path}" in
    *.tar.gz|*.tgz) tar -xzf "${tar_path}" -C "${dest_dir}" ;;
    *.tar)          tar -xf "${tar_path}" -C "${dest_dir}" ;;
    *)              die "Unsupported tar suffix: ${tar_path}" ;;
  esac
  touch "${marker}"
}

stage_fastq_links() {
  rm -rf "${FASTQ_LINK_ROOT}"
  mkdir -p "${FASTQ_LINK_ROOT}/gex" "${FASTQ_LINK_ROOT}/adt"

  local gex_n=0
  local adt_n=0
  local ambiguous="${OUTDIR}/ambiguous_fastqs.txt"
  : > "${ambiguous}"

  while IFS= read -r fq; do
    local base lower dst
    base="$(basename "${fq}")"
    lower="${base,,}"
    if [[ "${lower}" == *antibody* || "${lower}" == *adt* || "${lower}" == *"cell_surface"* || "${lower}" == *"cell-surface"* ]]; then
      dst="${FASTQ_LINK_ROOT}/adt/${base}"
      ln -sfn "${fq}" "${dst}"
      adt_n=$((adt_n + 1))
    elif [[ "${lower}" == *gex* || "${lower}" == *"gene_expression"* || "${lower}" == *"gene-expression"* ]]; then
      dst="${FASTQ_LINK_ROOT}/gex/${base}"
      ln -sfn "${fq}" "${dst}"
      gex_n=$((gex_n + 1))
    else
      printf '%s\n' "${fq}" >> "${ambiguous}"
    fi
  done < <(find "${RAW_FASTQ_ROOT}" -type f -name '*.fastq.gz' | sort)

  if [[ "${adt_n}" -eq 0 ]]; then
    log "FASTQ files found:"
    find "${RAW_FASTQ_ROOT}" -type f -name '*.fastq.gz' | sort | sed 's/^/  /'
    die "Could not identify Antibody Capture FASTQs from file names"
  fi
  if [[ "${gex_n}" -eq 0 ]]; then
    log "NOTE: no GEX FASTQs were classified; process_features ADT benchmark can still run"
  fi

  log "Staged FASTQ symlinks: gex=${gex_n} adt=${adt_n}"
}

find_mex_dir() {
  local root="$1"
  find "${root}" -type f \( -name 'matrix.mtx.gz' -o -name 'matrix.mtx' \) -printf '%h\n' | head -n1
}

download_file "${FEATURE_REF_URL}" "${FEATURE_REF}" "feature reference"
download_file "${FILTERED_MEX_URL}" "${FILTERED_TAR}" "Cell Ranger filtered MEX"
download_file "${METRICS_URL}" "${METRICS_CSV}" "metrics summary"
if [[ "${SKIP_RAW_MEX}" -eq 0 ]]; then
  download_file "${RAW_MEX_URL}" "${RAW_TAR}" "Cell Ranger raw MEX"
fi

extract_tar_once "${FILTERED_TAR}" "${CELLRANGER_ROOT}/filtered_extract" \
  "${CELLRANGER_ROOT}/.filtered_extracted" "Cell Ranger filtered MEX"
FILTERED_MEX="$(find_mex_dir "${CELLRANGER_ROOT}/filtered_extract")"
[[ -n "${FILTERED_MEX}" ]] || die "Could not locate extracted filtered MEX"
rm -rf "${CELLRANGER_ROOT}/filtered_feature_bc_matrix"
ln -sfn "${FILTERED_MEX}" "${CELLRANGER_ROOT}/filtered_feature_bc_matrix"

if [[ "${SKIP_RAW_MEX}" -eq 0 ]]; then
  extract_tar_once "${RAW_TAR}" "${CELLRANGER_ROOT}/raw_extract" \
    "${CELLRANGER_ROOT}/.raw_extracted" "Cell Ranger raw MEX"
  RAW_MEX="$(find_mex_dir "${CELLRANGER_ROOT}/raw_extract")"
  [[ -n "${RAW_MEX}" ]] || die "Could not locate extracted raw MEX"
  rm -rf "${CELLRANGER_ROOT}/raw_feature_bc_matrix"
  ln -sfn "${RAW_MEX}" "${CELLRANGER_ROOT}/raw_feature_bc_matrix"
else
  RAW_MEX=""
fi

if [[ "${SKIP_FASTQS}" -eq 0 ]]; then
  download_file "${FASTQS_TAR_URL}" "${FASTQS_TAR}" "FASTQ tar (~$((FASTQS_TAR_BYTES / 1000000000)) GB)"
  extract_tar_once "${FASTQS_TAR}" "${RAW_FASTQ_ROOT}" "${OUTDIR}/.fastqs_extracted" "FASTQs"
  stage_fastq_links
else
  log "Skipping FASTQ download/extraction"
fi

{
  echo "10x PBMC 1k TotalSeq-B CITE-seq fixture"
  echo "Generated (UTC): $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  echo "Script: ${REPO_ROOT}/scripts/download_10x_citeseq_pbmc1k_fixture.sh"
  echo "Dataset: ${DATASET_URL}"
  echo "FASTQ tar URL: ${FASTQS_TAR_URL}"
  echo "Feature ref URL: ${FEATURE_REF_URL}"
  echo "Filtered MEX URL: ${FILTERED_MEX_URL}"
  echo "Raw MEX URL: ${RAW_MEX_URL}"
  echo "Feature ref: ${FEATURE_REF}"
  echo "ADT FASTQ dir: ${FASTQ_LINK_ROOT}/adt"
  echo "GEX FASTQ dir: ${FASTQ_LINK_ROOT}/gex"
  echo "Cell Ranger filtered MEX: ${CELLRANGER_ROOT}/filtered_feature_bc_matrix"
  [[ -n "${RAW_MEX:-}" ]] && echo "Cell Ranger raw MEX: ${CELLRANGER_ROOT}/raw_feature_bc_matrix"
  echo "Metrics summary: ${METRICS_CSV}"
} > "${MANIFEST}"

log "Wrote manifest: ${MANIFEST}"
log "Feature ref: ${FEATURE_REF}"
log "Filtered MEX: ${CELLRANGER_ROOT}/filtered_feature_bc_matrix"
if [[ "${SKIP_FASTQS}" -eq 0 ]]; then
  log "ADT FASTQs: ${FASTQ_LINK_ROOT}/adt"
fi
