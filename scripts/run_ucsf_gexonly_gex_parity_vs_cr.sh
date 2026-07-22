#!/usr/bin/env bash
# GEX-only STAR Solo.out/GeneFull vs Cell Ranger 9 Gene Expression MEX.
# Use when STAR has no outs/crispr_analysis (run_ucsf_gexonly_no_bam_benchmark.sh).
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

STAR_RUN="${STAR_RUN:-}"
CR_RUN="${CR_RUN:-}"
OUT_REPORT="${OUT_REPORT:-}"

BARCODE_TRANSLATION="${BARCODE_TRANSLATION:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt.gz}"
TRANSLATION_DIRECTION="${TRANSLATION_DIRECTION:-left-to-right}"
GENE_CORR_MIN_COUNTS="${GENE_CORR_MIN_COUNTS:-20}"
GENE_CORR_MIN_CELLS_PCT="${GENE_CORR_MIN_CELLS_PCT:-0.01}"

usage() {
  cat <<EOF
Usage: STAR_RUN=/path/to/star_out $(basename "$0")

Optional env:
  CR_RUN=...              (required: matching UCSF EBs2_2 CR9 run)
  OUT_REPORT=path.txt     append stdout to this file (tee)
  BARCODE_TRANSLATION=... NXT->TRU map (default: Cell Ranger 9 path)
  TRANSLATION_DIRECTION=left-to-right|right-to-left
  GENE_CORR_MIN_COUNTS GENE_CORR_MIN_CELLS_PCT  (same as report_additional_parity_metrics.py)

Runs scripts/report_additional_parity_metrics.py with --skip-feature-call-parity.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ -z "${STAR_RUN}" ]]; then
  usage >&2
  exit 2
fi

if [[ -z "${CR_RUN}" ]]; then
  cat >&2 <<EOF
ERROR: CR_RUN is required.

The old default CR path used by this helper was a stale UCSF multi-library run
with only 233,828,574 reads and 7,325 filtered cells, which is not the
matching full EBs2_2 benchmark surface for current GEX-only parity.

Pass the actual matching Cell Ranger 9 run explicitly:
  CR_RUN=/path/to/cellranger_run STAR_RUN=/path/to/star_run $(basename "$0")
EOF
  exit 2
fi

for d in \
  "${STAR_RUN}/Solo.out/GeneFull/raw" \
  "${STAR_RUN}/Solo.out/GeneFull/filtered" \
  "${CR_RUN}/outs/raw_feature_bc_matrix" \
  "${CR_RUN}/outs/filtered_feature_bc_matrix"; do
  [[ -d "${d}" ]] || {
    echo "ERROR: missing directory: ${d}" >&2
    exit 2
  }
done

if [[ ! -f "${BARCODE_TRANSLATION}" ]]; then
  echo "ERROR: barcode translation file missing: ${BARCODE_TRANSLATION}" >&2
  exit 2
fi

CMD=(
  python3 "${REPO_ROOT}/scripts/report_additional_parity_metrics.py"
  --cr-run "${CR_RUN}"
  --star-run "${STAR_RUN}"
  --skip-feature-call-parity
  --barcode-translation "${BARCODE_TRANSLATION}"
  --translation-direction "${TRANSLATION_DIRECTION}"
  --translate both
  --gene-corr-min-counts "${GENE_CORR_MIN_COUNTS}"
  --gene-corr-min-cells-pct "${GENE_CORR_MIN_CELLS_PCT}"
)

if [[ -n "${OUT_REPORT}" ]]; then
  mkdir -p "$(dirname "${OUT_REPORT}")"
  "${CMD[@]}" | tee -a "${OUT_REPORT}"
else
  "${CMD[@]}"
fi
