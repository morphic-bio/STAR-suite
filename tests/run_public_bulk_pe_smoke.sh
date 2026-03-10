#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

source "${SCRIPT_DIR}/external_fixtures_env.sh"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
SALMON_BIN="${SALMON_BIN:-salmon}"
TRIMVALIDATE_BIN="${TRIMVALIDATE_BIN:-${REPO_ROOT}/core/features/vbem/tools/trimvalidate/trimvalidate}"
GENOME_DIR="${PUBLIC_BULK_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
TRANSCRIPTOME="${PUBLIC_BULK_TRANSCRIPTOME:-/storage/autoindex_110_44/bulk_index/transcriptome.fa}"
R1="${PUBLIC_BULK_R1:-${PUBLIC_BULK_FIXTURE_ROOT}/SRR4422207_1.fastq.gz}"
R2="${PUBLIC_BULK_R2:-${PUBLIC_BULK_FIXTURE_ROOT}/SRR4422207_2.fastq.gz}"
THREADS="${PUBLIC_BULK_THREADS:-8}"
OUTDIR="${PUBLIC_BULK_OUTDIR:-${SCRIPT_DIR}/public_bulk_pe_smoke_output_$(date +%Y%m%d_%H%M%S)}"
STAR_MATE_ORDER="${PUBLIC_BULK_STAR_MATE_ORDER:-R2R1}"
SALMON_LIBTYPE="${PUBLIC_BULK_SALMON_LIBTYPE:-IU}"
EXPECTED_FORMAT="${PUBLIC_BULK_EXPECTED_FORMAT:-}"
MIN_SPEARMAN_ALL="${PUBLIC_BULK_MIN_SPEARMAN_ALL:-0.93}"
MIN_SPEARMAN_EXPRESSED="${PUBLIC_BULK_MIN_SPEARMAN_EXPRESSED:-0.985}"
MIN_PEARSON_EXPRESSED="${PUBLIC_BULK_MIN_PEARSON_EXPRESSED:-0.99}"
MAX_DROPPED_INCOMPAT="${PUBLIC_BULK_MAX_DROPPED_INCOMPAT:--1}"
TRIM_QUALITY="${PUBLIC_BULK_TRIM_QUALITY:-20}"
TRIM_MIN_LENGTH="${PUBLIC_BULK_TRIM_MIN_LENGTH:-20}"
ADAPTER_R1="${PUBLIC_BULK_ADAPTER_R1:-AGATCGGAAGAGCACACGTCTGAACTCCAGTCA}"
ADAPTER_R2="${PUBLIC_BULK_ADAPTER_R2:-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT}"

Y_VALIDATOR="${SCRIPT_DIR}/validate_bulk_yremove_output.py"
VB_VALIDATOR="${SCRIPT_DIR}/transcriptvb/validate_pe_autodetect_output.py"
REFERENCE_SPLITTER="${SCRIPT_DIR}/reference_split_fastq_by_bam.py"

build_if_missing() {
    local target="$1"
    shift
    if [[ -x "${target}" ]]; then
        return
    fi
    "$@"
}

build_if_missing "${TRIMVALIDATE_BIN}" make -C "${REPO_ROOT}/core/features/vbem/tools/trimvalidate" -j8

for cmd in "${STAR_BIN}" "${SALMON_BIN}" "${TRIMVALIDATE_BIN}" samtools python3; do
    if ! command -v "${cmd}" >/dev/null 2>&1 && [[ ! -x "${cmd}" ]]; then
        echo "ERROR: missing required command ${cmd}" >&2
        exit 2
    fi
done

[[ -f "${R1}" ]] || { echo "ERROR: missing fixture R1 ${R1}" >&2; exit 2; }
[[ -f "${R2}" ]] || { echo "ERROR: missing fixture R2 ${R2}" >&2; exit 2; }
[[ -d "${GENOME_DIR}" ]] || { echo "ERROR: missing genome dir ${GENOME_DIR}" >&2; exit 2; }
[[ -f "${TRANSCRIPTOME}" ]] || { echo "ERROR: missing transcriptome ${TRANSCRIPTOME}" >&2; exit 2; }
[[ -f "${REFERENCE_SPLITTER}" ]] || { echo "ERROR: missing reference splitter ${REFERENCE_SPLITTER}" >&2; exit 2; }

case "${STAR_MATE_ORDER}" in
    R1R2|R2R1) ;;
    *)
        echo "ERROR: PUBLIC_BULK_STAR_MATE_ORDER must be R1R2 or R2R1" >&2
        exit 2
        ;;
esac

mkdir -p "${OUTDIR}"
REF_SPLIT_DIR="${OUTDIR}/reference_split"
RAW_R1="${REF_SPLIT_DIR}/raw_R1.fastq"
RAW_R2="${REF_SPLIT_DIR}/raw_R2.fastq"
TRIMMED_R1="${REF_SPLIT_DIR}/trimmed_R1.fastq"
TRIMMED_R2="${REF_SPLIT_DIR}/trimmed_R2.fastq"
REF_Y_R1="${REF_SPLIT_DIR}/trimmed_R1_Y.fastq"
REF_Y_R2="${REF_SPLIT_DIR}/trimmed_R2_Y.fastq"
REF_NOY_R1="${REF_SPLIT_DIR}/trimmed_R1_noY.fastq"
REF_NOY_R2="${REF_SPLIT_DIR}/trimmed_R2_noY.fastq"
mkdir -p "${REF_SPLIT_DIR}"

echo "=== Public Bulk PE Smoke ==="
echo "STAR: ${STAR_BIN}"
echo "Salmon: ${SALMON_BIN}"
echo "trimvalidate: ${TRIMVALIDATE_BIN}"
echo "reference splitter: ${REFERENCE_SPLITTER}"
echo "Genome: ${GENOME_DIR}"
echo "Transcriptome: ${TRANSCRIPTOME}"
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "STAR readFilesIn order: ${STAR_MATE_ORDER}"
echo "Salmon libType: ${SALMON_LIBTYPE}"
echo "Output: ${OUTDIR}"
echo

gzip -dc "${R1}" > "${RAW_R1}"
gzip -dc "${R2}" > "${RAW_R2}"

"${TRIMVALIDATE_BIN}" \
    -1 "${RAW_R1}" \
    -2 "${RAW_R2}" \
    -o1 "${TRIMMED_R1}" \
    -o2 "${TRIMMED_R2}" \
    --quality "${TRIM_QUALITY}" \
    --length "${TRIM_MIN_LENGTH}" \
    --adapter-r1 "${ADAPTER_R1}" \
    --adapter-r2 "${ADAPTER_R2}"

case "${STAR_MATE_ORDER}" in
    R1R2)
        STAR_READ_1="${TRIMMED_R1}"
        STAR_READ_2="${TRIMMED_R2}"
        ;;
    R2R1)
        STAR_READ_1="${TRIMMED_R2}"
        STAR_READ_2="${TRIMMED_R1}"
        ;;
esac

"${STAR_BIN}" \
    --runMode alignReads \
    --runThreadN "${THREADS}" \
    --genomeDir "${GENOME_DIR}" \
    --readFilesIn "${STAR_READ_1}" "${STAR_READ_2}" \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortMethod samtools \
    --limitIObufferSize 50000000 50000000 \
    --outFileNamePrefix "${OUTDIR}/" \
    --outTmpDir "${OUTDIR}/tmp" \
    --emitNoYBAM yes \
    --emitYNoYFastq yes \
    --emitYNoYFastqCompression gz \
    --quantMode TranscriptomeSAM TranscriptVB \
    --quantVBgcBias 1 \
    --quantVBLibType A

python3 "${REFERENCE_SPLITTER}" \
    --ybam "${OUTDIR}/Aligned.sortedByCoord.out_Y.bam" \
    --outdir "${REF_SPLIT_DIR}" \
    "${TRIMMED_R1}" "${TRIMMED_R2}"

python3 "${Y_VALIDATOR}" \
    --outdir "${OUTDIR}" \
    --require-y-reads \
    --reference-y-r1 "${REF_Y_R1}" \
    --reference-y-r2 "${REF_Y_R2}"

"${SALMON_BIN}" quant \
    -t "${TRANSCRIPTOME}" \
    -l "${SALMON_LIBTYPE}" \
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
