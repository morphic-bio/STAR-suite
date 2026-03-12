#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

ACCESSION="${ACCESSION:-SRR32576116}"
OUTDIR="${OUTDIR:-${DATA_ROOT}/public_slam_${ACCESSION}}"
MAX_READS="${MAX_READS:-100000}"
THREADS="${THREADS:-4}"
FORCE=0

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Download a small public SLAM-seq FASTQ subset.

Options:
  --accession ACC  SRA accession. Default: ${ACCESSION}
  --outdir DIR     Output directory. Default: ${OUTDIR}
  --max-reads N    Number of reads to keep. Default: ${MAX_READS}
  --threads N      pigz threads. Default: ${THREADS}
  --force          Re-download and overwrite
  -h, --help       Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --accession) ACCESSION="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --max-reads) MAX_READS="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

need_cmd fastq-dump
need_cmd pigz

OUTDIR="$(mkdir -p "${OUTDIR}" && cd "${OUTDIR}" && pwd)"
FASTQ_GZ="${OUTDIR}/${ACCESSION}.fastq.gz"
FASTQ_RAW="${OUTDIR}/${ACCESSION}.fastq"
PASS_FASTQ_RAW="${OUTDIR}/${ACCESSION}_pass.fastq"

if [[ -f "${FASTQ_GZ}" && "${FORCE}" == "0" ]]; then
  log "Fixture already present: ${FASTQ_GZ}"
else
  perl -e 'unlink @ARGV' "${FASTQ_GZ}" "${FASTQ_RAW}" "${PASS_FASTQ_RAW}"
  fastq-dump \
    --skip-technical \
    --readids \
    --read-filter pass \
    --dumpbase \
    --clip \
    -X "${MAX_READS}" \
    -O "${OUTDIR}" \
    "${ACCESSION}"
  if [[ -f "${PASS_FASTQ_RAW}" ]]; then
    mv "${PASS_FASTQ_RAW}" "${FASTQ_RAW}"
  fi
  [[ -f "${FASTQ_RAW}" ]] || die "Missing FASTQ after fastq-dump: ${FASTQ_RAW}"
  pigz -f -p "${THREADS}" "${FASTQ_RAW}"
fi

READ_COUNT="$(gzip -cd "${FASTQ_GZ}" | awk 'END {print NR/4}')"
cat > "${OUTDIR}/MANIFEST.txt" <<MANIFEST
Public SLAM fixture
Generated (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)
SRA accession: ${ACCESSION}
FASTQ: ${FASTQ_GZ}
Reads: ${READ_COUNT}
MANIFEST

log "FASTQ: ${FASTQ_GZ}"
log "Reads: ${READ_COUNT}"
