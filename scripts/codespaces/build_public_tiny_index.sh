#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

TINY_ROOT="${TINY_ROOT:-${DATA_ROOT}/public_10x_tiny/cellranger-tiny-ref}"
OUTDIR="${OUTDIR:-${INDEX_ROOT}/public_10x_tiny_star}"
THREADS="${THREADS:-4}"
FORCE=0

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Build a tiny STAR index from the public 10x tiny reference.

Options:
  --tiny-root DIR  Tiny reference root. Default: ${TINY_ROOT}
  --outdir DIR     Output STAR genomeDir. Default: ${OUTDIR}
  --threads N      STAR genomeGenerate threads. Default: ${THREADS}
  --force          Rebuild even if the index already exists
  -h, --help       Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tiny-root) TINY_ROOT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

ensure_star_built
[[ -d "${TINY_ROOT}" ]] || die "Missing tiny reference root: ${TINY_ROOT}"
FASTA="${TINY_ROOT}/fasta/genome.fa"
GTF="$(resolve_tiny_gtf "${TINY_ROOT}")"
[[ -f "${FASTA}" ]] || die "Missing FASTA: ${FASTA}"

OUTDIR="$(mkdir -p "${OUTDIR}" && cd "${OUTDIR}" && pwd)"
if [[ -f "${OUTDIR}/Genome" && "${FORCE}" == "0" ]]; then
  log "Index already exists: ${OUTDIR}"
  exit 0
fi

TMP_GTF="${OUTDIR}/genes.gtf"
if [[ "${GTF}" == *.gz ]]; then
  gzip -cd "${GTF}" > "${TMP_GTF}"
else
  cp "${GTF}" "${TMP_GTF}"
fi

"${STAR_BIN}" \
  --runMode genomeGenerate \
  --runThreadN "${THREADS}" \
  --genomeDir "${OUTDIR}" \
  --genomeFastaFiles "${FASTA}" \
  --sjdbGTFfile "${TMP_GTF}" \
  --sjdbOverhang 90 \
  --genomeSAindexNbases 10

log "Built tiny STAR index: ${OUTDIR}"
