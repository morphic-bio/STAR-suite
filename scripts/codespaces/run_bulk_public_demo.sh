#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"
# shellcheck source=public_demo_sources.sh
source "${SCRIPT_DIR}/public_demo_sources.sh"

GENOME_DIR="${GENOME_DIR:-}"
OUTDIR="${OUTDIR:-${RUN_ROOT}/bulk_public_demo}"
THREADS="${THREADS:-4}"
RUN=0
DRY_RUN=0

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Fetch a small public bulk RNA-seq PE fixture and emit a STAR-core bulk command.
Pass --run to execute the command after the fixture is downloaded.

Options:
  --genome-dir DIR  STAR genomeDir for execution
  --outdir DIR      Output directory. Default: ${RUN_ROOT}/bulk_public_demo
  --threads N       STAR threads. Default: ${THREADS}
  --run             Execute STAR after writing RUN_COMMAND.sh
  --dry-run         Download fixture and write command without executing
  -h, --help        Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome-dir) GENOME_DIR="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --run) RUN=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

need_cmd pigz
need_cmd fastq-dump
ensure_star_built

FIXTURE_DIR="${DATA_ROOT}/bulk_${CODESPACES_PUBLIC_BULK_SRR}"
"${REPO_ROOT}/scripts/download_public_bulk_fixture.sh" \
  --accession "${CODESPACES_PUBLIC_BULK_SRR}" \
  --geo-series "${CODESPACES_PUBLIC_BULK_GSE}" \
  --geo-sample "${CODESPACES_PUBLIC_BULK_GSM}" \
  --outdir "${FIXTURE_DIR}"

R1="${FIXTURE_DIR}/${CODESPACES_PUBLIC_BULK_SRR}_1.fastq.gz"
R2="${FIXTURE_DIR}/${CODESPACES_PUBLIC_BULK_SRR}_2.fastq.gz"
OUTDIR="$(mkdir -p "${OUTDIR}" && cd "${OUTDIR}" && pwd)"

if [[ -z "${GENOME_DIR}" ]]; then
  GENOME_DIR="/path/to/star_bulk_index"
fi

CMD=(
  "${STAR_BIN}"
  --runMode alignReads
  --runThreadN "${THREADS}"
  --genomeDir "${GENOME_DIR}"
  --readFilesIn "${R2}" "${R1}"
  --readFilesCommand zcat
  --trimCutadapt Yes
  --outSAMtype BAM SortedByCoordinate
  --quantMode GeneCounts
  --outFileNamePrefix "${OUTDIR}/"
)

write_command_script "${OUTDIR}/RUN_COMMAND.sh" "${CMD[@]}"
cat > "${OUTDIR}/RUN_MANIFEST.txt" <<MANIFEST
Module: STAR-core bulk
Fixture: ${CODESPACES_PUBLIC_BULK_SRR}
Output: ${OUTDIR}
GenomeDir: ${GENOME_DIR}
R1: ${R1}
R2: ${R2}
MANIFEST

if [[ "${RUN}" == "1" && "${DRY_RUN}" == "0" ]]; then
  [[ -d "${GENOME_DIR}" ]] || die "Missing --genome-dir for --run: ${GENOME_DIR}"
  "${CMD[@]}"
else
  log "Wrote ${OUTDIR}/RUN_COMMAND.sh"
fi
