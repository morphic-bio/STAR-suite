#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

TINY_REPO_URL="${TINY_REPO_URL:-https://github.com/minoda-lab/universc.git}"
TINY_REPO_COMMIT="${TINY_REPO_COMMIT:-7cbd039613b45c64f4b6d8219906aafda28dd5f9}"
OUTDIR="${OUTDIR:-${DATA_ROOT}/public_10x_tiny}"
FORCE=0

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Fetch the public tiny 10x test surface used by universc.
This requires git-lfs because the FASTQ/reference assets are stored in LFS.

Options:
  --outdir DIR   Output directory. Default: ${DATA_ROOT}/public_10x_tiny
  --force        Re-stage the output directory
  -h, --help     Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

need_cmd git
need_cmd git-lfs

CACHE_REPO="${CACHE_ROOT}/universc-tiny"
OUTDIR="$(mkdir -p "${OUTDIR}" && cd "${OUTDIR}" && pwd)"

if [[ ! -d "${CACHE_REPO}/.git" ]]; then
  git clone "${TINY_REPO_URL}" "${CACHE_REPO}"
fi

(
  cd "${CACHE_REPO}"
  git fetch --all --tags
  git checkout "${TINY_REPO_COMMIT}"
  git lfs install --local
  git lfs pull --include="test/shared/cellranger-tiny-fastq/3.0.0/**,test/cellranger_reference/cellranger-tiny-ref/3.0.0/**"
)

FASTQ_SRC="${CACHE_REPO}/test/shared/cellranger-tiny-fastq/3.0.0"
REF_SRC="${CACHE_REPO}/test/cellranger_reference/cellranger-tiny-ref/3.0.0"
[[ -f "${FASTQ_SRC}/tinygex_S1_L001_R1_001.fastq.gz" ]] || die "Expected tiny FASTQ missing after git-lfs pull"
[[ -f "${REF_SRC}/fasta/genome.fa" ]] || die "Expected tiny reference FASTA missing after git-lfs pull"
if head -n 1 "${REF_SRC}/fasta/genome.fa" | grep -q 'git-lfs'; then
  die "git-lfs assets were not hydrated; install git-lfs and rerun"
fi

if [[ "${FORCE}" == "1" ]]; then
  perl -e 'use File::Path qw(remove_tree); for (@ARGV) { remove_tree($_) if -e $_ }' \
    "${OUTDIR}/cellranger-tiny-fastq" "${OUTDIR}/cellranger-tiny-ref"
fi

mkdir -p "${OUTDIR}"
ln -sfn "${FASTQ_SRC}" "${OUTDIR}/cellranger-tiny-fastq"
ln -sfn "${REF_SRC}" "${OUTDIR}/cellranger-tiny-ref"

cat > "${OUTDIR}/MANIFEST.txt" <<MANIFEST
Public tiny 10x fixture
Generated (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)
Source repo: ${TINY_REPO_URL}
Pinned commit: ${TINY_REPO_COMMIT}
FASTQ root: ${OUTDIR}/cellranger-tiny-fastq
Reference root: ${OUTDIR}/cellranger-tiny-ref
MANIFEST

log "FASTQ root: ${OUTDIR}/cellranger-tiny-fastq"
log "Reference root: ${OUTDIR}/cellranger-tiny-ref"
