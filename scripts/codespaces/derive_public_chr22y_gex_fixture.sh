#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

OUTDIR="${OUTDIR:-${DATA_ROOT}/public_chr22y_gex_fixture}"
THREADS="${THREADS:-4}"
READ_LIMIT="${READ_LIMIT:-200000}"
INDEX_DIR="${INDEX_DIR:-${INDEX_ROOT}/public_human_chr22y_star}"
DRY_RUN=0

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Create a small public GEX fixture by aligning a bounded number of reads from a
public male 10x run to the chr22+chrY mini-reference and retaining only the
raw read pairs whose cDNA mate maps to chr22 or chrY.

Options:
  --outdir DIR      Output directory. Default: ${DATA_ROOT}/public_chr22y_gex_fixture
  --threads N       STAR threads. Default: ${THREADS}
  --read-limit N    Number of read pairs to inspect from the public input. Default: ${READ_LIMIT}
  --index-dir DIR   chr22+chrY STAR index. Default: ${INDEX_ROOT}/public_human_chr22y_star
  --dry-run         Prepare inputs/commands without executing STAR/filtering
  -h, --help        Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --read-limit) READ_LIMIT="$2"; shift 2 ;;
    --index-dir) INDEX_DIR="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

need_cmd samtools
ensure_star_built
"${SCRIPT_DIR}/fetch_public_10x_male_fastqs.sh"
"${SCRIPT_DIR}/build_public_chr22y_index.sh" --outdir "${INDEX_DIR}" --threads "${THREADS}"

RAW_DIR="${DATA_ROOT}/public_10x_male_bcell/raw_fastqs"
R1="$(find "${RAW_DIR}" -maxdepth 1 -type f -name '*_R1_001.fastq.gz' | sort | head -n1)"
[[ -n "${R1}" ]] || die "No extracted R1 FASTQ under ${RAW_DIR}"
R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
[[ -f "${R2}" ]] || die "Missing R2 mate for ${R1}"

OUTDIR="$(mkdir -p "${OUTDIR}" && cd "${OUTDIR}" && pwd)"
STAGE_DIR="${OUTDIR}/stage"
ALIGN_DIR="${OUTDIR}/align"
FIXTURE_DIR="${OUTDIR}/gex/public_chr22y_demo"
mkdir -p "${STAGE_DIR}" "${ALIGN_DIR}" "${FIXTURE_DIR}"

STAGE_R1="${STAGE_DIR}/public_stage_R1.fastq.gz"
STAGE_R2="${STAGE_DIR}/public_stage_R2.fastq.gz"
FIXTURE_R1="${FIXTURE_DIR}/public_chr22y_demo_S1_L001_R1_001.fastq.gz"
FIXTURE_R2="${FIXTURE_DIR}/public_chr22y_demo_S1_L001_R2_001.fastq.gz"
MAPPED_NAMES="${OUTDIR}/mapped_chr22y.names.txt"
MANIFEST="${OUTDIR}/MANIFEST.txt"

SAMPLE_COPY_CMD=(python3 "${SCRIPT_DIR}/copy_first_fastq_pairs.py" --r1 "${R1}" --r2 "${R2}" --out-r1 "${STAGE_R1}" --out-r2 "${STAGE_R2}" --read-limit "${READ_LIMIT}")
ALIGN_CMD=("${STAR_BIN}" --runThreadN "${THREADS}" --genomeDir "${INDEX_DIR}" --readFilesIn "${STAGE_R2}" --readFilesCommand zcat --outFileNamePrefix "${ALIGN_DIR}/" --outSAMtype BAM Unsorted --outSJtype None --clipAdapterType CellRanger4 --soloType None)
FILTER_CMD=(python3 "${SCRIPT_DIR}/filter_fastq_by_names.py" --names "${MAPPED_NAMES}" --r1 "${STAGE_R1}" --r2 "${STAGE_R2}" --out-r1 "${FIXTURE_R1}" --out-r2 "${FIXTURE_R2}")

write_command_script "${OUTDIR}/RUN_SAMPLE_COPY.sh" "${SAMPLE_COPY_CMD[@]}"
write_command_script "${OUTDIR}/RUN_ALIGN.sh" "${ALIGN_CMD[@]}"
write_command_script "${OUTDIR}/RUN_FILTER.sh" "${FILTER_CMD[@]}"

if [[ "${DRY_RUN}" == "1" ]]; then
  log "Prepared ${OUTDIR}"
  exit 0
fi

SAMPLED_COUNT="$(${SAMPLE_COPY_CMD[@]})"
if [[ ! -f "${ALIGN_DIR}/Aligned.out.bam" ]]; then
  "${ALIGN_CMD[@]}"
fi
samtools view -F 4 "${ALIGN_DIR}/Aligned.out.bam" | cut -f1 | LC_ALL=C sort -u > "${MAPPED_NAMES}"
KEPT_COUNT="$(${FILTER_CMD[@]})"
MAPPED_COUNT="$(wc -l < "${MAPPED_NAMES}")"

{
  echo "Public chr22+chrY GEX fixture"
  echo "Source R1: ${R1}"
  echo "Source R2: ${R2}"
  echo "Sampled read pairs: ${SAMPLED_COUNT}"
  echo "Mapped names: ${MAPPED_COUNT}"
  echo "Kept read pairs: ${KEPT_COUNT}"
  echo "Fixture dir: ${FIXTURE_DIR}"
  echo "Index dir: ${INDEX_DIR}"
} > "${MANIFEST}"

log "Prepared ${FIXTURE_DIR} (${KEPT_COUNT} read pairs)"
