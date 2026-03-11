#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REMOVE_Y_READS_BIN="${REMOVE_Y_READS_BIN:-${ROOT_DIR}/core/features/yremove_fastq/tools/remove_y_reads/remove_y_reads}"
YBAM_ROOT="${YBAM_ROOT:-/mnt/pikachu/slam_y}"
FASTQ_DIR="${FASTQ_DIR:-/mnt/pikachu/NW-5-21/SLAM-Seq}"
OUTDIR="${OUTDIR:-/mnt/pikachu/NW-5-21/SLAM-Seq_ysplit_$(date -u +%Y%m%d_%H%M%S)}"
JOBS="${JOBS:-4}"
GZIP_LEVEL="${GZIP_LEVEL:-6}"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --ybam-root DIR         Root containing per-sample *_Y.bam files
  --fastq-dir DIR         Directory containing full NW SLAM FASTQs
  --outdir DIR            Output root on /mnt/pikachu
  --jobs N                Parallel sample jobs (default: ${JOBS})
  --gzip-level N          Gzip compression level for remove_y_reads (default: ${GZIP_LEVEL})
  --remove-y-reads-bin P  Path to remove_y_reads binary
  -h, --help              Show help
EOF
}

log() {
    printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"
}

die() {
    log "ERROR: $*"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ybam-root) YBAM_ROOT="$2"; shift 2 ;;
        --fastq-dir) FASTQ_DIR="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --jobs) JOBS="$2"; shift 2 ;;
        --gzip-level) GZIP_LEVEL="$2"; shift 2 ;;
        --remove-y-reads-bin) REMOVE_Y_READS_BIN="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

[[ -d "${YBAM_ROOT}" ]] || die "Y BAM root not found: ${YBAM_ROOT}"
[[ -d "${FASTQ_DIR}" ]] || die "FASTQ directory not found: ${FASTQ_DIR}"
[[ -x "${REMOVE_Y_READS_BIN}" ]] || die "remove_y_reads binary not executable: ${REMOVE_Y_READS_BIN}"
command -v samtools >/dev/null 2>&1 || die "samtools not found on PATH"
[[ "${JOBS}" =~ ^[0-9]+$ ]] || die "--jobs must be an integer"
[[ "${GZIP_LEVEL}" =~ ^[0-9]+$ ]] || die "--gzip-level must be an integer"

mkdir -p "${OUTDIR}"/{logs,manifests}

MANIFEST="${OUTDIR}/manifests/samples.tsv"
MISSING="${OUTDIR}/manifests/missing_fastqs.tsv"
FAILED="${OUTDIR}/manifests/failed_samples.tsv"
: > "${MANIFEST}"
: > "${MISSING}"
: > "${FAILED}"

printf 'sample\tybam\tfastq\toutdir\ty_names\n' > "${MANIFEST}"
printf 'sample\tybam\treason\n' > "${MISSING}"
printf 'sample\tybam\tfastq\toutdir\tlog\n' > "${FAILED}"

map_fastq_for_sample() {
    local sample="$1"
    find "${FASTQ_DIR}" -maxdepth 1 -type f -name "${sample}_S*_R1_001.fastq.gz" | sort | head -n1
}

process_sample() {
    local ybam="$1"
    local sample
    local fastq
    local sample_out
    local names_out
    local log_out

    sample="$(basename "${ybam}")"
    sample="${sample%_Y.bam}"
    fastq="$(map_fastq_for_sample "${sample}")"

    if [[ -z "${fastq}" ]]; then
        printf '%s\t%s\tno matching FASTQ in %s\n' "${sample}" "${ybam}" "${FASTQ_DIR}" >> "${MISSING}"
        return 0
    fi

    sample_out="${OUTDIR}/${sample}"
    names_out="${sample_out}/${sample}_Y.names.txt"
    log_out="${OUTDIR}/logs/${sample}.log"
    mkdir -p "${sample_out}"

    samtools view "${ybam}" \
        | cut -f1 \
        | sed -E 's/[[:space:]].*$//' \
        | sed -E 's#(/1|/2|\\.1|\\.2)$##' \
        | sort -u > "${names_out}"

    {
        echo "sample=${sample}"
        echo "ybam=${ybam}"
        echo "fastq=${fastq}"
        echo "outdir=${sample_out}"
        echo "y_names=${names_out}"
        "${REMOVE_Y_READS_BIN}" \
            -y "${ybam}" \
            -o "${sample_out}" \
            --threads 1 \
            --gzip-level "${GZIP_LEVEL}" \
            "${fastq}"
    } > "${log_out}" 2>&1 || {
        printf '%s\t%s\t%s\t%s\t%s\n' "${sample}" "${ybam}" "${fastq}" "${sample_out}" "${log_out}" >> "${FAILED}"
        return 1
    }

    printf '%s\t%s\t%s\t%s\t%s\n' "${sample}" "${ybam}" "${fastq}" "${sample_out}" "${names_out}" >> "${MANIFEST}"
}

export ROOT_DIR REMOVE_Y_READS_BIN YBAM_ROOT FASTQ_DIR OUTDIR GZIP_LEVEL MANIFEST MISSING FAILED
export -f map_fastq_for_sample process_sample

log "Output root: ${OUTDIR}"
log "Scanning Y BAMs under ${YBAM_ROOT}"

find "${YBAM_ROOT}" -mindepth 2 -maxdepth 2 -type f -name '*_Y.bam' -print0 \
    | sort -z \
    | xargs -0 -n1 -P "${JOBS}" bash -lc 'process_sample "$@"' _

log "Completed. Sample manifest: ${MANIFEST}"
log "Missing FASTQs: ${MISSING}"
log "Failed samples: ${FAILED}"
