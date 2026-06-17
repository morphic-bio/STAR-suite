#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONTRACT_DIR="${REPO_ROOT}/core/features/libchromap_contract"
CHROMAP_DIR="${CHROMAP_DIR:-/mnt/pikachu/Chromap-suite-signac-profile}"
PEAK_MEX="${PEAK_MEX:-${CONTRACT_DIR}/star_multiome_atac_peak_mex}"
OUT="${OUT:-$(mktemp -d /tmp/star_peak_mex_signac_profile.XXXXXX)}"

mkdir -p "${OUT}/mex" "${OUT}/tmp" "${OUT}/keep"

if [[ ! -x "${PEAK_MEX}" ]]; then
  make -C "${CONTRACT_DIR}" CHROMAP_DIR="${CHROMAP_DIR}" star_multiome_atac_peak_mex
fi

python3 - "${OUT}/sidecar.bin" <<'PY'
import struct
import sys

path = sys.argv[1]
records = []
for i in range(500):
    start = 1000 + i
    end = start + 50
    records.append((0, start, end, 1, 0))
for i in range(40):
    start = 10000 + i * 20
    end = start + 50
    records.append((0, start, end, 1, 0))

with open(path, "wb") as fh:
    fh.write(struct.pack("<4sIIIIIQ", b"AEV1", 1, 24, 16, 1, 0, len(records)))
    for record in records:
        fh.write(struct.pack("<iiiIQ", *record))

with open(path + ".chroms.tsv", "w", encoding="ascii") as fh:
    fh.write("0\tchr1\n")
PY

"${PEAK_MEX}" \
  --sidecar "${OUT}/sidecar.bin" \
  --call-peaks-from-sidecar \
  --peak-call-mode macs-bed \
  --macs-profile signac-atac \
  --peaks "${OUT}/peaks.narrowPeak" \
  --summits-out "${OUT}/summits.bed" \
  --out-dir "${OUT}/mex" \
  --metrics-tsv "${OUT}/metrics.tsv" \
  --temp-dir "${OUT}/tmp" \
  --keep-intermediates-dir "${OUT}/keep" \
  --force \
  >"${OUT}/peak_mex.stdout" \
  2>"${OUT}/peak_mex.stderr"

test -s "${OUT}/peaks.narrowPeak"
test -s "${OUT}/summits.bed"
test -s "${OUT}/mex/matrix.mtx.gz"
test -s "${OUT}/mex/barcodes.tsv.gz"
test -s "${OUT}/mex/features.tsv.gz"
test -s "${OUT}/metrics.tsv"
test -s "${OUT}/keep/sidecar_signac-atac_macs_bed.tsv"

grep -q "MACS BED profile=signac-atac" "${OUT}/peak_mex.stderr"

echo "PASS: STAR sidecar peak-MEX Signac MACS-BED profile smoke (${OUT})"
