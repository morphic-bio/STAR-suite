#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

REF_ROOT="${REF_ROOT:-${DATA_ROOT}/public_human_chr22y_ref}"
OUTDIR="${OUTDIR:-${INDEX_ROOT}/public_human_chr22y_star}"
THREADS="${THREADS:-4}"
SA_INDEX_N_BASES="${SA_INDEX_N_BASES:-11}"

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Build a STAR genomeDir for the public chr22+chrY demo reference.

Options:
  --ref-root DIR            Reference root. Default: ${DATA_ROOT}/public_human_chr22y_ref
  --outdir DIR              STAR genomeDir output. Default: ${INDEX_ROOT}/public_human_chr22y_star
  --threads N               STAR threads. Default: ${THREADS}
  --sa-index-n-bases N      genomeSAindexNbases. Default: ${SA_INDEX_N_BASES}
  -h, --help                Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref-root) REF_ROOT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --sa-index-n-bases) SA_INDEX_N_BASES="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

ensure_star_built
"${SCRIPT_DIR}/fetch_public_chr22y_reference.sh" --outdir "${REF_ROOT}"

[[ -f "${REF_ROOT}/fasta/genome.fa" ]] || die "Missing genome.fa under ${REF_ROOT}"
[[ -f "${REF_ROOT}/genes/genes.gtf" ]] || die "Missing genes.gtf under ${REF_ROOT}"

OUTDIR="$(mkdir -p "${OUTDIR}" && cd "${OUTDIR}" && pwd)"
if [[ -f "${OUTDIR}/Genome" && -f "${OUTDIR}/SA" ]]; then
  log "Using existing STAR index ${OUTDIR}"
  exit 0
fi

"${STAR_BIN}" \
  --runMode genomeGenerate \
  --genomeDir "${OUTDIR}" \
  --runThreadN "${THREADS}" \
  --genomeFastaFiles "${REF_ROOT}/fasta/genome.fa" \
  --sjdbGTFfile "${REF_ROOT}/genes/genes.gtf" \
  --genomeSAindexNbases "${SA_INDEX_N_BASES}"

log "Built ${OUTDIR}"
