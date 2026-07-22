#!/usr/bin/env bash
# Verify split-read discovery accepts SRR*_1.fastq.gz style names.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

FIXTURE_DIR="${CATATAC_GUIDE_FIXTURE:-/mnt/pikachu/catatac_gse288996/guide_redump/fixture}"
SPLIT_BIN="${REPO_ROOT}/core/features/process_features/tests/test_catatac_split_read"
ATAC_WHITELIST="${CATATAC_ATAC_WHITELIST:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt}"
FEATURE_REF="${CATATAC_GUIDE_FEATURE_REF:-${REPO_ROOT}/core/features/process_features/feature_lists/catatac_crispri_guide_capture.csv}"
OUT_DIR="${CATATAC_SRR_NAMING_OUT:-${REPO_ROOT}/tests/catatac_srr_naming_output}"

if [[ ! -x "${SPLIT_BIN}" ]]; then
  make -C "${REPO_ROOT}/core/features/process_features" tests/test_catatac_split_read
fi

TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/catatac_srr_naming.XXXXXX")"
trap 'rm -rf "${TMPDIR}"' EXIT

ln -s "${FIXTURE_DIR}/guide_R1.fastq.gz" "${TMPDIR}/SRR32265756_1.fastq.gz"
ln -s "${FIXTURE_DIR}/guide_R2.fastq.gz" "${TMPDIR}/SRR32265756_2.fastq.gz"
ln -s "${FIXTURE_DIR}/guide_R3.fastq.gz" "${TMPDIR}/SRR32265756_3.fastq.gz"

rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

export CATATAC_GUIDE_MAX_READS="${CATATAC_GUIDE_MAX_READS:-5000}"
"${SPLIT_BIN}" "${TMPDIR}" "${OUT_DIR}" "${ATAC_WHITELIST}" "${FEATURE_REF}"

test -s "${OUT_DIR}/sample/matrix.mtx"
test -s "${OUT_DIR}/split_read_metrics.tsv"
python3 - <<PY
import sys
from pathlib import Path
metrics = dict(
    line.split("\t")
    for line in Path("${OUT_DIR}/split_read_metrics.tsv").read_text().splitlines()
    if "\t" in line and not line.startswith("metric")
)
synth = int(metrics.get("barcode_synth_ok", "0"))
if synth <= 0:
    sys.exit("split-read produced no synthesized barcodes for SRR naming")
PY

echo "CAT-ATAC SRR fastq naming smoke OK: ${OUT_DIR}"
