#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

ACCESSION="SRR4422207"
GEO_SERIES="GSE88509"
GEO_SAMPLE="GSM2344101"
SAMPLE_LABEL="male_esophagus_muscularis_mucosa"
OUTDIR="/tmp/public_bulk_fixture_SRR4422207"
MAX_SPOTS=500000
THREADS=8
FORCE=0
KEEP_UNCOMPRESSED=0

usage() {
    cat <<'EOF'
Usage:
  scripts/download_public_bulk_fixture.sh [options]

Downloads a deterministic public human male bulk RNA-seq PE fixture from GEO/SRA.
The default fixture is:
  GEO series:   GSE88509
  GEO sample:   GSM2344101
  SRA run:      SRR4422207
  Description:  Homo sapiens esophagus muscularis mucosa, male, stranded PE101

Options:
  --accession ACC            SRA run accession. Default: SRR4422207
  --geo-series ACC           GEO series accession. Default: GSE88509
  --geo-sample ACC           GEO sample accession. Default: GSM2344101
  --sample-label LABEL       Human-readable label for metadata.
  --outdir DIR               Output directory. Default: /tmp/public_bulk_fixture_SRR4422207
  --max-spots N              Number of spots to fetch with fastq-dump -X. Default: 500000
  --threads N                Compression threads. Default: 8
  --force                    Re-download and overwrite existing FASTQs.
  --keep-uncompressed        Keep intermediate .fastq files after gzip.
  -h, --help                 Show help.

Outputs:
  <outdir>/<accession>_1.fastq.gz
  <outdir>/<accession>_2.fastq.gz
  <outdir>/FIXTURE_MANIFEST.txt
EOF
}

log() {
    printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*"
}

die() {
    echo "ERROR: $*" >&2
    exit 2
}

need_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --accession) ACCESSION="$2"; shift 2 ;;
        --geo-series) GEO_SERIES="$2"; shift 2 ;;
        --geo-sample) GEO_SAMPLE="$2"; shift 2 ;;
        --sample-label) SAMPLE_LABEL="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --max-spots) MAX_SPOTS="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --force) FORCE=1; shift ;;
        --keep-uncompressed) KEEP_UNCOMPRESSED=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

need_cmd fastq-dump
need_cmd pigz

OUTDIR="$(mkdir -p "${OUTDIR}" && cd "${OUTDIR}" && pwd)"
MANIFEST="${OUTDIR}/FIXTURE_MANIFEST.txt"
R1_GZ="${OUTDIR}/${ACCESSION}_1.fastq.gz"
R2_GZ="${OUTDIR}/${ACCESSION}_2.fastq.gz"
R1_FASTQ="${OUTDIR}/${ACCESSION}_1.fastq"
R2_FASTQ="${OUTDIR}/${ACCESSION}_2.fastq"
R1_PASS_FASTQ="${OUTDIR}/${ACCESSION}_pass_1.fastq"
R2_PASS_FASTQ="${OUTDIR}/${ACCESSION}_pass_2.fastq"

if [[ "${FORCE}" -eq 0 && -f "${R1_GZ}" && -f "${R2_GZ}" ]]; then
    log "Fixture already present: ${R1_GZ} and ${R2_GZ}"
else
    if [[ "${FORCE}" -eq 0 && -f "${R1_PASS_FASTQ}" && -f "${R2_PASS_FASTQ}" ]]; then
        log "Resuming from partial fastq-dump output"
    elif [[ "${FORCE}" -eq 0 && -f "${R1_FASTQ}" && -f "${R2_FASTQ}" ]]; then
        log "Resuming from existing uncompressed FASTQs"
    else
        rm -f "${R1_FASTQ}" "${R2_FASTQ}" "${R1_PASS_FASTQ}" "${R2_PASS_FASTQ}" "${R1_GZ}" "${R2_GZ}"
        log "Downloading ${ACCESSION} with fastq-dump (max spots: ${MAX_SPOTS})"
        fastq-dump \
            --split-files \
            --skip-technical \
            --readids \
            --read-filter pass \
            --dumpbase \
            --clip \
            -X "${MAX_SPOTS}" \
            -O "${OUTDIR}" \
            "${ACCESSION}"
    fi

    if [[ -f "${R1_PASS_FASTQ}" && -f "${R2_PASS_FASTQ}" ]]; then
        mv "${R1_PASS_FASTQ}" "${R1_FASTQ}"
        mv "${R2_PASS_FASTQ}" "${R2_FASTQ}"
    fi

    [[ -f "${R1_FASTQ}" ]] || die "Missing output FASTQ: ${R1_FASTQ}"
    [[ -f "${R2_FASTQ}" ]] || die "Missing output FASTQ: ${R2_FASTQ}"

    log "Compressing FASTQs with pigz"
    pigz -f -p "${THREADS}" "${R1_FASTQ}" "${R2_FASTQ}"
fi

R1_COUNT="$(gzip -cd "${R1_GZ}" | awk 'END {print NR/4}')"
R2_COUNT="$(gzip -cd "${R2_GZ}" | awk 'END {print NR/4}')"
[[ "${R1_COUNT}" = "${R2_COUNT}" ]] || die "Paired FASTQ counts differ: R1=${R1_COUNT}, R2=${R2_COUNT}"

{
    echo "Public Bulk PE Fixture"
    echo "Generated (UTC): $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    echo "Script: ${REPO_ROOT}/scripts/download_public_bulk_fixture.sh"
    echo "SRA run: ${ACCESSION}"
    echo "GEO series: ${GEO_SERIES}"
    echo "GEO sample: ${GEO_SAMPLE}"
    echo "Sample label: ${SAMPLE_LABEL}"
    echo "Download method: fastq-dump --split-files -X ${MAX_SPOTS}"
    echo "R1: ${R1_GZ}"
    echo "R2: ${R2_GZ}"
    echo "Pairs: ${R1_COUNT}"
} > "${MANIFEST}"

if [[ "${KEEP_UNCOMPRESSED}" -eq 0 ]]; then
    rm -f "${R1_FASTQ}" "${R2_FASTQ}" "${R1_PASS_FASTQ}" "${R2_PASS_FASTQ}"
fi

log "Wrote fixture manifest: ${MANIFEST}"
log "R1: ${R1_GZ}"
log "R2: ${R2_GZ}"
log "Pairs: ${R1_COUNT}"
