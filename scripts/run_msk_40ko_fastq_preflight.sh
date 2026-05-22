#!/usr/bin/env bash
# MSK 40KO FASTQ preflight bundle.
#
# Read-only with respect to FASTQs. It checks that each manifest row lands in
# the expected 10x whitelist family/namespace before any STAR production run is
# launched.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

MANIFEST="${MANIFEST:-${REPO_ROOT}/docs/MSK_40KO_FASTQ_MANIFEST.tsv}"
OUTDIR="${OUTDIR:-${REPO_ROOT}/plans/artifacts/msk_40ko_fastq_preflight_$(date +%Y%m%d_%H%M%S)}"
SAMPLE_READS="${SAMPLE_READS:-200000}"
MAX_FASTQS_PER_ROW="${MAX_FASTQS_PER_ROW:-4}"

FEB2018_TRU="${FEB2018_TRU:-/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt}"
FEB2018_NXT="${FEB2018_NXT:-/storage/scRNAseq_output/whitelists/3M-february-2018_NXT.txt}"
GEMX2023_TRU="${GEMX2023_TRU:-/storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_TRU.txt}"
GEMX2023_NXT="${GEMX2023_NXT:-/storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_NXT.txt}"

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --manifest TSV          FASTQ manifest (default: ${MANIFEST})
  --outdir DIR            Output directory (default: ${OUTDIR})
  --sample-reads N        Reads to sample per manifest row (default: ${SAMPLE_READS})
  --max-fastqs-per-row N  R1 FASTQs sampled per manifest row (default: ${MAX_FASTQS_PER_ROW})
  -h, --help              Show help

Environment overrides:
  FEB2018_TRU, FEB2018_NXT, GEMX2023_TRU, GEMX2023_NXT
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --sample-reads) SAMPLE_READS="$2"; shift 2 ;;
    --max-fastqs-per-row) MAX_FASTQS_PER_ROW="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for path in "${MANIFEST}" "${FEB2018_TRU}" "${FEB2018_NXT}" "${GEMX2023_TRU}" "${GEMX2023_NXT}"; do
  [[ -f "${path}" ]] || {
    echo "ERROR: missing required input: ${path}" >&2
    exit 1
  }
done

mkdir -p "${OUTDIR}"

python3 "${SCRIPT_DIR}/preflight_whitelist_family.py" \
  --manifest "${MANIFEST}" \
  --whitelist "feb2018:TRU:${FEB2018_TRU}" \
  --whitelist "feb2018:NXT:${FEB2018_NXT}" \
  --whitelist "may2023_gemx:TRU:${GEMX2023_TRU}" \
  --whitelist "may2023_gemx:NXT:${GEMX2023_NXT}" \
  --sample-reads "${SAMPLE_READS}" \
  --max-fastqs-per-manifest-row "${MAX_FASTQS_PER_ROW}" \
  --outdir "${OUTDIR}"

echo
echo "MSK 40KO FASTQ preflight outputs:"
echo "  ${OUTDIR}/whitelist_family_summary.tsv"
echo "  ${OUTDIR}/whitelist_family_rates.tsv"
echo "  ${OUTDIR}/whitelist_family_report.json"
