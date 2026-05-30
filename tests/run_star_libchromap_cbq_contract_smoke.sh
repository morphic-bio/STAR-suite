#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

CHROMAP_BIN="${CHROMAP_BIN:-/mnt/pikachu/Chromap-suite/chromap}"
RUNNER="${RUNNER:-${ROOT_DIR}/core/features/libchromap_contract/star_libchromap_contract_runner}"
ENCODER="${ENCODER:-${ROOT_DIR}/core/legacy/source/cbq_ordered_encoder}"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_libchromap_cbq_contract_smoke}"

DATA_DIR="${OUT_ROOT}/data"
RUN_DIR="${OUT_ROOT}/run"
CBQ_DIR="${OUT_ROOT}/cbq"

rm -rf "${OUT_ROOT}"
mkdir -p "${DATA_DIR}" "${RUN_DIR}" "${CBQ_DIR}"

if [[ ! -x "${RUNNER}" ]]; then
  make -C "${ROOT_DIR}" star-libchromap-contract
fi
if [[ ! -x "${ENCODER}" ]]; then
  make -C "${ROOT_DIR}/core/legacy/source" cbq-ordered-encoder
fi
if [[ ! -x "${CHROMAP_BIN}" ]]; then
  echo "ERROR: chromap binary not found at ${CHROMAP_BIN}" >&2
  exit 1
fi

python3 - "${DATA_DIR}" <<'PY'
import random
import sys
from pathlib import Path

out = Path(sys.argv[1])
rng = random.Random(8675309)
alphabet = "ACGT"
genome = "".join(rng.choice(alphabet) for _ in range(12000))

def rc(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

barcodes = [
    "ACGTACGTACGTACGT",
    "TGCATGCATGCATGCA",
    "GGGGAAAACCCCTTTT",
    "AAAACCCCGGGGTTTT",
]
fragments = [
    ("read001", 400, 560, barcodes[0], "I"),
    ("read002", 2100, 2280, barcodes[1], "H"),
    ("read003", 5300, 5485, barcodes[2], "J"),
    ("read004", 8600, 8790, barcodes[3], "I"),
]

(out / "genome.fa").write_text(">chrSynthetic\n" + genome + "\n", encoding="ascii")
(out / "whitelist.txt").write_text("\n".join(barcodes) + "\n", encoding="ascii")

with (out / "reads_R1.fastq").open("wt", encoding="ascii") as r1, \
     (out / "reads_R2.fastq").open("wt", encoding="ascii") as r2, \
     (out / "barcodes.fastq").open("wt", encoding="ascii") as bc:
    for name, start, end, barcode, qchar in fragments:
        seq1 = genome[start:start + 90]
        seq2 = rc(genome[end - 90:end])
        r1.write(f"@{name}/1 1:N:0:ACGT\n{seq1}\n+\n{qchar * len(seq1)}\n")
        r2.write(f"@{name}/2 2:N:0:ACGT\n{seq2}\n+\n{qchar * len(seq2)}\n")
        bc.write(f"@{name} 3:N:0:ACGT\n{barcode}\n+\n{qchar * len(barcode)}\n")
PY

"${CHROMAP_BIN}" --build-index \
  -r "${DATA_DIR}/genome.fa" \
  -o "${DATA_DIR}/genome.idx" \
  -k 11 -w 5 \
  > "${RUN_DIR}/index.stdout" 2> "${RUN_DIR}/index.stderr"

"${ENCODER}" \
  --readFilesIn "${DATA_DIR}/reads_R1.fastq" "${DATA_DIR}/reads_R2.fastq" \
  --outFile "${CBQ_DIR}/reads_pair.cbq" \
  > "${RUN_DIR}/encode_reads.stdout" 2> "${RUN_DIR}/encode_reads.stderr"

"${ENCODER}" \
  --readFilesIn "${DATA_DIR}/barcodes.fastq" \
  --outFile "${CBQ_DIR}/barcodes.cbq" \
  > "${RUN_DIR}/encode_barcodes.stdout" 2> "${RUN_DIR}/encode_barcodes.stderr"

common_args=(
  --ref "${DATA_DIR}/genome.fa"
  --index "${DATA_DIR}/genome.idx"
  --barcode-whitelist "${DATA_DIR}/whitelist.txt"
  --read-format "bc:0:-1"
  --output-format BED
  --threads 1
)

"${RUNNER}" "${common_args[@]}" \
  --read1 "${DATA_DIR}/reads_R1.fastq" \
  --read2 "${DATA_DIR}/reads_R2.fastq" \
  --barcode "${DATA_DIR}/barcodes.fastq" \
  --output "${RUN_DIR}/fastq.fragments.bed" \
  > "${RUN_DIR}/fastq.stdout" 2> "${RUN_DIR}/fastq.stderr"

"${RUNNER}" "${common_args[@]}" \
  --input-format cbq \
  --read-pair-cbq "${CBQ_DIR}/reads_pair.cbq" \
  --barcode-cbq "${CBQ_DIR}/barcodes.cbq" \
  --output "${RUN_DIR}/cbq.fragments.bed" \
  > "${RUN_DIR}/cbq.stdout" 2> "${RUN_DIR}/cbq.stderr"

if [[ ! -s "${RUN_DIR}/fastq.fragments.bed" ]]; then
  echo "ERROR: FASTQ contract run produced no fragments" >&2
  exit 1
fi
if [[ ! -s "${RUN_DIR}/cbq.fragments.bed" ]]; then
  echo "ERROR: CBQ contract run produced no fragments" >&2
  exit 1
fi

LC_ALL=C sort "${RUN_DIR}/fastq.fragments.bed" > "${RUN_DIR}/fastq.fragments.sorted.bed"
LC_ALL=C sort "${RUN_DIR}/cbq.fragments.bed" > "${RUN_DIR}/cbq.fragments.sorted.bed"
if ! cmp -s "${RUN_DIR}/fastq.fragments.sorted.bed" "${RUN_DIR}/cbq.fragments.sorted.bed"; then
  echo "ERROR: FASTQ and CBQ contract fragments differ" >&2
  diff -u "${RUN_DIR}/fastq.fragments.sorted.bed" \
    "${RUN_DIR}/cbq.fragments.sorted.bed" | head -80 >&2 || true
  exit 1
fi

echo "PASS: STAR libchromap CBQ contract smoke output at ${OUT_ROOT}"
