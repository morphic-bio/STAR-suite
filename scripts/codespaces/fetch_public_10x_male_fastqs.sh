#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"
# shellcheck source=public_demo_sources.sh
source "${SCRIPT_DIR}/public_demo_sources.sh"

OUTDIR="${OUTDIR:-${DATA_ROOT}/public_10x_male_bcell}"
TAR_PATH="${TAR_PATH:-${CACHE_ROOT}/public_10x_male_bcell_fastqs.tar}"

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Download and extract a single public male 10x GEX FASTQ pair used to derive the
chr22+chrY Codespaces fixture.

Options:
  --outdir DIR     Extraction root. Default: ${DATA_ROOT}/public_10x_male_bcell
  --tar-path FILE  Cache path for the downloaded FASTQ tar. Default: ${CACHE_ROOT}/public_10x_male_bcell_fastqs.tar
  -h, --help       Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --tar-path) TAR_PATH="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

need_cmd curl
need_cmd tar

OUTDIR="$(mkdir -p "${OUTDIR}" && cd "${OUTDIR}" && pwd)"
RAW_DIR="${OUTDIR}/raw_fastqs"
mkdir -p "${RAW_DIR}" "$(dirname -- "${TAR_PATH}")"

if [[ ! -f "${TAR_PATH}" ]]; then
  log "Downloading ${CODESPACES_PUBLIC_10X_MALE_FASTQS_TAR_URL}"
  curl -L --fail -o "${TAR_PATH}" "${CODESPACES_PUBLIC_10X_MALE_FASTQS_TAR_URL}"
fi

if compgen -G "${RAW_DIR}/*_R1_001.fastq.gz" > /dev/null; then
  log "Using existing extracted FASTQ pair under ${RAW_DIR}"
  exit 0
fi

mapfile -t R1S < <(tar -tf "${TAR_PATH}" | grep '_R1_001.fastq.gz$' | sort)
[[ ${#R1S[@]} -gt 0 ]] || die "No R1 FASTQs found in ${TAR_PATH}"
FIRST_R1="${R1S[0]}"
FIRST_R2="${FIRST_R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
log "Extracting ${FIRST_R1} and ${FIRST_R2}"
tar -xf "${TAR_PATH}" -C "${RAW_DIR}" "${FIRST_R1}" "${FIRST_R2}"

# Flatten nested archive paths into RAW_DIR.
while IFS= read -r -d '' f; do
  base="$(basename -- "$f")"
  if [[ "$f" != "${RAW_DIR}/${base}" ]]; then
    mv "$f" "${RAW_DIR}/${base}"
  fi
done < <(find "${RAW_DIR}" -type f \( -name '*_R1_001.fastq.gz' -o -name '*_R2_001.fastq.gz' \) -print0)
find "${RAW_DIR}" -mindepth 1 -type d -empty -delete

log "Prepared ${RAW_DIR}"
