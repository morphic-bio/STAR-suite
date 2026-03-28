#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

OUTDIR="${OUTDIR:-${RUN_ROOT}/bulk_public_demo}"

usage() {
  cat <<USAGE
Usage: $(basename "$0") [--outdir DIR]

Verify the expected outputs from scripts/codespaces/run_bulk_public_demo.sh.
This checks that the command script, sorted BAM, and gene-count table exist.
If samtools is available, it also runs a lightweight BAM integrity check.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

OUTDIR="$(cd "$(dirname "${OUTDIR}")" && pwd)/$(basename "${OUTDIR}")"
CMD_SCRIPT="${OUTDIR}/RUN_COMMAND.sh"
BAM="${OUTDIR}/Aligned.sortedByCoord.out.bam"
GENE_COUNTS="${OUTDIR}/ReadsPerGene.out.tab"

[[ -x "${CMD_SCRIPT}" ]] || die "Missing executable command script: ${CMD_SCRIPT}"

if [[ ! -s "${BAM}" || ! -s "${GENE_COUNTS}" ]]; then
  cat >&2 <<EOF
ERROR: Expected bulk outputs are not present under ${OUTDIR}.

If this directory came from --dry-run, execute the saved command first:
  bash ${CMD_SCRIPT}
EOF
  [[ -s "${BAM}" ]] || echo "Missing BAM output: ${BAM}" >&2
  [[ -s "${GENE_COUNTS}" ]] || echo "Missing gene-count output: ${GENE_COUNTS}" >&2
  exit 2
fi

if command -v samtools >/dev/null 2>&1; then
  samtools quickcheck "${BAM}"
  BAM_READS="$(samtools view -c "${BAM}")"
else
  BAM_READS="samtools-not-installed"
fi

GENE_ROWS="$(wc -l < "${GENE_COUNTS}")"

cat <<EOF
Bulk public demo outputs look good.

Verified files:
  ${CMD_SCRIPT}
  ${BAM}
  ${GENE_COUNTS}

Quick summary:
  BAM reads: ${BAM_READS}
  ReadsPerGene rows: ${GENE_ROWS}

Useful inspection commands:
  bash ${CMD_SCRIPT}
  ls -lh "${BAM}" "${GENE_COUNTS}"
  head -n 8 "${GENE_COUNTS}"
EOF
