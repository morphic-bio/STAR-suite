#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BENCH_DIR="${REPO_ROOT}/benchmarks/nfcore_rnaseq_compare"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
SALMON_BIN="${SALMON_BIN:-salmon}"
THREADS="${TRANSCRIPTVB_CHR22_THREADS:-32}"
SALMON_THREADS="${TRANSCRIPTVB_CHR22_SALMON_THREADS:-1}"
BASE_RUN="${TRANSCRIPTVB_CHR22_BASE_RUN:-${BENCH_DIR}/runs/chr22_20260626_003522}"
FIXTURE_DIR="${TRANSCRIPTVB_CHR22_FIXTURE_DIR:-${BENCH_DIR}/fixtures/chr22}"
OUTDIR="${TRANSCRIPTVB_CHR22_OUTDIR:-${SCRIPT_DIR}/transcriptvb_chr22_parity_output_$(date +%Y%m%d_%H%M%S)}"

MIN_NUMREADS_PEARSON="${TRANSCRIPTVB_CHR22_MIN_NUMREADS_PEARSON:-0.99998}"
MIN_TPM_PEARSON="${TRANSCRIPTVB_CHR22_MIN_TPM_PEARSON:-0.99995}"
MAX_TOTAL_DELTA="${TRANSCRIPTVB_CHR22_MAX_TOTAL_DELTA:-0.05}"
MAX_HALF_L1="${TRANSCRIPTVB_CHR22_MAX_HALF_L1:-20}"
MAX_GT1_TX="${TRANSCRIPTVB_CHR22_MAX_GT1_TX:-10}"
MAX_MAX_TX_DELTA="${TRANSCRIPTVB_CHR22_MAX_TX_DELTA:-5}"
EXPECTED_FORMAT="${TRANSCRIPTVB_CHR22_EXPECTED_FORMAT:-ISR}"

R1="${FIXTURE_DIR}/SRR4422207_chr22_1.fastq.gz"
R2="${FIXTURE_DIR}/SRR4422207_chr22_2.fastq.gz"
STAR_INDEX="${BASE_RUN}/reference/star_index"
TRANSCRIPTOME="${BASE_RUN}/reference/transcriptome.fa"

for cmd in "${STAR_BIN}" "${SALMON_BIN}" python3; do
    if ! command -v "${cmd}" >/dev/null 2>&1 && [[ ! -x "${cmd}" ]]; then
        echo "ERROR: missing required command ${cmd}" >&2
        exit 2
    fi
done

[[ -f "${R1}" ]] || { echo "ERROR: missing chr22 fixture R1: ${R1}" >&2; exit 2; }
[[ -f "${R2}" ]] || { echo "ERROR: missing chr22 fixture R2: ${R2}" >&2; exit 2; }
[[ -d "${STAR_INDEX}" ]] || { echo "ERROR: missing chr22 STAR index: ${STAR_INDEX}" >&2; exit 2; }
[[ -f "${TRANSCRIPTOME}" ]] || { echo "ERROR: missing chr22 transcriptome: ${TRANSCRIPTOME}" >&2; exit 2; }

STAR_OUT="${OUTDIR}/star_auto_A_noerr"
SALMON_OUT="${OUTDIR}/salmon_fixed_ISR_noerr_p1"
mkdir -p "${STAR_OUT}" "${SALMON_OUT}"

echo "=== TranscriptVB chr22 parity smoke ==="
echo "STAR: ${STAR_BIN}"
echo "Salmon: ${SALMON_BIN}"
echo "STAR threads: ${THREADS}"
echo "Salmon threads: ${SALMON_THREADS}"
echo "Fixture: ${FIXTURE_DIR}"
echo "Reference run: ${BASE_RUN}"
echo "Output: ${OUTDIR}"
echo

"${STAR_BIN}" \
    --runMode alignReads \
    --runThreadN "${THREADS}" \
    --genomeDir "${STAR_INDEX}" \
    --readFilesIn "${R1}" "${R2}" \
    --trimCutadapt Yes \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortMethod samtools \
    --quantMode TranscriptomeSAM TranscriptVB \
    --quantVBLibType A \
    --quantVBErrorModel off \
    --outFileNamePrefix "${STAR_OUT}/" \
    --outTmpDir "${STAR_OUT}/_tmp" \
    > "${STAR_OUT}/run.stdout.log" \
    2> "${STAR_OUT}/run.stderr.log"

if ! grep -q "Detected library format: ${EXPECTED_FORMAT}" "${STAR_OUT}/Log.out"; then
    echo "ERROR: expected detected library format ${EXPECTED_FORMAT}" >&2
    grep "Detected library format" "${STAR_OUT}/Log.out" >&2 || true
    exit 1
fi

"${SALMON_BIN}" quant \
    -t "${TRANSCRIPTOME}" \
    -l ISR \
    -p "${SALMON_THREADS}" \
    -a "${STAR_OUT}/Aligned.toTranscriptome.out.bam" \
    --noErrorModel \
    -o "${SALMON_OUT}" \
    > "${SALMON_OUT}/salmon.stdout.log" \
    2> "${SALMON_OUT}/salmon.stderr.log"

python3 - <<PY
import csv
import math
from pathlib import Path

import numpy as np

star_path = Path("${STAR_OUT}") / "quant.sf"
salmon_path = Path("${SALMON_OUT}") / "quant.sf"
summary_path = Path("${OUTDIR}") / "summary.tsv"

def load(path):
    out = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            out[row["Name"]] = (float(row["NumReads"]), float(row["TPM"]))
    return out

def pearson(a, b):
    if len(a) < 2 or float(np.std(a)) == 0.0 or float(np.std(b)) == 0.0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])

def spearman(a, b):
    if len(a) < 2:
        return float("nan")
    return pearson(np.argsort(np.argsort(a)).astype(float),
                   np.argsort(np.argsort(b)).astype(float))

star = load(star_path)
salmon = load(salmon_path)
keys = sorted(set(star) | set(salmon))
star_num = np.array([star.get(k, (0.0, 0.0))[0] for k in keys])
salmon_num = np.array([salmon.get(k, (0.0, 0.0))[0] for k in keys])
star_tpm = np.array([star.get(k, (0.0, 0.0))[1] for k in keys])
salmon_tpm = np.array([salmon.get(k, (0.0, 0.0))[1] for k in keys])
delta = star_num - salmon_num
abs_delta = np.abs(delta)
total_delta = float(delta.sum())
abs_delta_sum = float(abs_delta.sum())
half_l1 = 0.5 * (abs_delta_sum + abs(total_delta))

metrics = {
    "numreads_pearson": pearson(star_num, salmon_num),
    "numreads_spearman": spearman(star_num, salmon_num),
    "tpm_pearson": pearson(star_tpm, salmon_tpm),
    "tpm_spearman": spearman(star_tpm, salmon_tpm),
    "star_total": float(star_num.sum()),
    "salmon_total": float(salmon_num.sum()),
    "total_delta": total_delta,
    "abs_delta_sum": abs_delta_sum,
    "half_l1_moved": half_l1,
    "tx_absdiff_gt_0p1": int((abs_delta > 0.1).sum()),
    "tx_absdiff_gt_1": int((abs_delta > 1.0).sum()),
    "max_tx_delta": float(abs_delta.max()) if len(abs_delta) else 0.0,
}

thresholds = {
    "numreads_pearson": float("${MIN_NUMREADS_PEARSON}"),
    "tpm_pearson": float("${MIN_TPM_PEARSON}"),
    "max_abs_total_delta": float("${MAX_TOTAL_DELTA}"),
    "max_half_l1": float("${MAX_HALF_L1}"),
    "max_gt1_tx": int("${MAX_GT1_TX}"),
    "max_tx_delta": float("${MAX_MAX_TX_DELTA}"),
}

failures = []
if metrics["numreads_pearson"] < thresholds["numreads_pearson"]:
    failures.append(f"NumReads Pearson {metrics['numreads_pearson']:.12f} < {thresholds['numreads_pearson']}")
if metrics["tpm_pearson"] < thresholds["tpm_pearson"]:
    failures.append(f"TPM Pearson {metrics['tpm_pearson']:.12f} < {thresholds['tpm_pearson']}")
if abs(metrics["total_delta"]) > thresholds["max_abs_total_delta"]:
    failures.append(f"abs(total delta) {abs(metrics['total_delta']):.6f} > {thresholds['max_abs_total_delta']}")
if metrics["half_l1_moved"] > thresholds["max_half_l1"]:
    failures.append(f"half-L1 moved {metrics['half_l1_moved']:.6f} > {thresholds['max_half_l1']}")
if metrics["tx_absdiff_gt_1"] > thresholds["max_gt1_tx"]:
    failures.append(f"transcripts >1 read {metrics['tx_absdiff_gt_1']} > {thresholds['max_gt1_tx']}")
if metrics["max_tx_delta"] > thresholds["max_tx_delta"]:
    failures.append(f"max transcript delta {metrics['max_tx_delta']:.6f} > {thresholds['max_tx_delta']}")

with summary_path.open("w") as out:
    out.write("metric\tvalue\n")
    for key, value in metrics.items():
        out.write(f"{key}\t{value}\n")

for key, value in metrics.items():
    if isinstance(value, float):
        print(f"{key}: {value:.12f}")
    else:
        print(f"{key}: {value}")

if failures:
    print("\\nFAIL:")
    for failure in failures:
        print(f"  - {failure}")
    raise SystemExit(1)
print("\\nPASS: TranscriptVB chr22 parity smoke")
PY

echo
echo "Summary: ${OUTDIR}/summary.tsv"
