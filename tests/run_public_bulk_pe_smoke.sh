#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

source "${SCRIPT_DIR}/external_fixtures_env.sh"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
SALMON_BIN="${SALMON_BIN:-salmon}"
GENOME_DIR="${PUBLIC_BULK_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
TRANSCRIPTOME="${PUBLIC_BULK_TRANSCRIPTOME:-/storage/autoindex_110_44/bulk_index/transcriptome.fa}"
R1="${PUBLIC_BULK_R1:-${PUBLIC_BULK_FIXTURE_ROOT}/SRR4422207_1.fastq.gz}"
R2="${PUBLIC_BULK_R2:-${PUBLIC_BULK_FIXTURE_ROOT}/SRR4422207_2.fastq.gz}"
THREADS="${PUBLIC_BULK_THREADS:-8}"
OUTDIR="${PUBLIC_BULK_OUTDIR:-${SCRIPT_DIR}/public_bulk_pe_smoke_output_$(date +%Y%m%d_%H%M%S)}"
STAR_MATE_ORDER="${PUBLIC_BULK_STAR_MATE_ORDER:-R2R1}"
EXPECTED_FORMAT="${PUBLIC_BULK_EXPECTED_FORMAT:-}"
MIN_SPEARMAN_ALL="${PUBLIC_BULK_MIN_SPEARMAN_ALL:-0.95}"
MIN_SPEARMAN_EXPRESSED="${PUBLIC_BULK_MIN_SPEARMAN_EXPRESSED:-0.99}"
MIN_PEARSON_EXPRESSED="${PUBLIC_BULK_MIN_PEARSON_EXPRESSED:-0.99}"
MAX_DROPPED_INCOMPAT="${PUBLIC_BULK_MAX_DROPPED_INCOMPAT:--1}"

Y_VALIDATOR="${SCRIPT_DIR}/validate_bulk_yremove_output.py"
VB_VALIDATOR="${SCRIPT_DIR}/transcriptvb/validate_pe_autodetect_output.py"

for cmd in "${STAR_BIN}" "${SALMON_BIN}" samtools python3; do
    if ! command -v "${cmd}" >/dev/null 2>&1 && [[ ! -x "${cmd}" ]]; then
        echo "ERROR: missing required command ${cmd}" >&2
        exit 2
    fi
done

[[ -f "${R1}" ]] || { echo "ERROR: missing fixture R1 ${R1}" >&2; exit 2; }
[[ -f "${R2}" ]] || { echo "ERROR: missing fixture R2 ${R2}" >&2; exit 2; }
[[ -d "${GENOME_DIR}" ]] || { echo "ERROR: missing genome dir ${GENOME_DIR}" >&2; exit 2; }
[[ -f "${TRANSCRIPTOME}" ]] || { echo "ERROR: missing transcriptome ${TRANSCRIPTOME}" >&2; exit 2; }

case "${STAR_MATE_ORDER}" in
    R1R2)
        STAR_READ_1="${R1}"
        STAR_READ_2="${R2}"
        ;;
    R2R1)
        STAR_READ_1="${R2}"
        STAR_READ_2="${R1}"
        ;;
    *)
        echo "ERROR: PUBLIC_BULK_STAR_MATE_ORDER must be R1R2 or R2R1" >&2
        exit 2
        ;;
esac

mkdir -p "${OUTDIR}"

echo "=== Public Bulk PE Smoke ==="
echo "STAR: ${STAR_BIN}"
echo "Salmon: ${SALMON_BIN}"
echo "Genome: ${GENOME_DIR}"
echo "Transcriptome: ${TRANSCRIPTOME}"
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "STAR readFilesIn order: ${STAR_MATE_ORDER}"
echo "Output: ${OUTDIR}"
echo

"${STAR_BIN}" \
    --runMode alignReads \
    --runThreadN "${THREADS}" \
    --genomeDir "${GENOME_DIR}" \
    --readFilesIn "${STAR_READ_1}" "${STAR_READ_2}" \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortMethod samtools \
    --limitIObufferSize 50000000 50000000 \
    --outFileNamePrefix "${OUTDIR}/" \
    --outTmpDir "${OUTDIR}/tmp" \
    --trimCutadapt Yes \
    --emitNoYBAM yes \
    --emitYNoYFastq yes \
    --emitYNoYFastqCompression gz \
    --quantMode TranscriptomeSAM TranscriptVB \
    --quantVBgcBias 1 \
    --quantVBLibType A

python3 "${Y_VALIDATOR}" --outdir "${OUTDIR}" --require-y-reads

"${SALMON_BIN}" quant \
    -t "${TRANSCRIPTOME}" \
    -l A \
    -a "${OUTDIR}/Aligned.toTranscriptome.out.bam" \
    -o "${OUTDIR}/salmon"

VALIDATE_ARGS=(
    --outdir "${OUTDIR}"
    --max-dropped-incompat "${MAX_DROPPED_INCOMPAT}"
    --salmon-quant "${OUTDIR}/salmon/quant.sf"
    --min-spearman-all "${MIN_SPEARMAN_ALL}"
    --min-spearman-expressed "${MIN_SPEARMAN_EXPRESSED}"
    --min-pearson-expressed "${MIN_PEARSON_EXPRESSED}"
)

if [[ -n "${EXPECTED_FORMAT}" ]]; then
    VALIDATE_ARGS+=(--expected-format "${EXPECTED_FORMAT}")
fi

python3 "${VB_VALIDATOR}" "${VALIDATE_ARGS[@]}"

echo
echo "PASS: public bulk PE smoke succeeded"
