#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"

STAR_SUITE_ROOT="${STAR_SUITE_ROOT:-${ROOT_DIR}}"
MORPHIC_RECIPES_ROOT="${MORPHIC_RECIPES_ROOT:-/mnt/pikachu/morphic-recipes}"
RECIPE="${MORPHIC_RECIPES_ROOT}/scripts/run_jax_scrnaseq02_ocm_composite_smoke.sh"

OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_cbq_ocm_composite_smoke_${STAMP}}"
READ_PAIRS="${READ_PAIRS:-1000}"
THREADS="${THREADS:-8}"
RUN_FASTQ_PARITY="${RUN_FASTQ_PARITY:-1}"

RAW_DIR="${RAW_DIR:-/mnt/pikachu/JAX_scRNAseq02/raw}"
CONFIG="${CONFIG:-/mnt/pikachu/JAX_scRNAseq02/cellranger-logs/config.csv}"
GENOME_DIR="${GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
SOLO_CB_WHITELIST="${SOLO_CB_WHITELIST:-/storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_TRU.txt}"
STAR_BIN="${STAR_BIN:-${STAR_SUITE_ROOT}/core/legacy/source/STAR}"
CBQ_ORDERED_ENCODER_BIN="${CBQ_ORDERED_ENCODER_BIN:-${STAR_SUITE_ROOT}/core/legacy/source/cbq_ordered_encoder}"

skip() {
    echo "SKIP: $*"
    exit 0
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

[[ -f "${RECIPE}" ]] || skip "missing Morphic OCM recipe: ${RECIPE}"
[[ -d "${RAW_DIR}" ]] || skip "missing JAX OCM raw dir: ${RAW_DIR}"
[[ -f "${CONFIG}" ]] || skip "missing JAX OCM config: ${CONFIG}"
[[ -d "${GENOME_DIR}" ]] || skip "missing STAR genomeDir: ${GENOME_DIR}"
[[ -f "${SOLO_CB_WHITELIST}" ]] || skip "missing STAR whitelist: ${SOLO_CB_WHITELIST}"
[[ -x "${STAR_BIN}" ]] || skip "missing STAR binary: ${STAR_BIN}"
[[ -x "${CBQ_ORDERED_ENCODER_BIN}" ]] || skip "missing cbq_ordered_encoder: ${CBQ_ORDERED_ENCODER_BIN}"
[[ "${READ_PAIRS}" =~ ^[0-9]+$ && "${READ_PAIRS}" -gt 0 ]] || die "READ_PAIRS must be positive"
[[ "${THREADS}" =~ ^[0-9]+$ && "${THREADS}" -gt 0 ]] || die "THREADS must be positive"

rm -rf "${OUT_ROOT}"
mkdir -p "${OUT_ROOT}"

CBQ_ROOT="${OUT_ROOT}/cbq"
FASTQ_ROOT="${OUT_ROOT}/fastq"

echo "Running OCM CBQ smoke: ${CBQ_ROOT}"
STAR_SUITE_ROOT="${STAR_SUITE_ROOT}" \
RAW_DIR="${RAW_DIR}" \
CONFIG="${CONFIG}" \
GENOME_DIR="${GENOME_DIR}" \
SOLO_CB_WHITELIST="${SOLO_CB_WHITELIST}" \
STAR_BIN="${STAR_BIN}" \
CBQ_ORDERED_ENCODER_BIN="${CBQ_ORDERED_ENCODER_BIN}" \
bash "${RECIPE}" \
  --read-pairs "${READ_PAIRS}" \
  --out-root "${CBQ_ROOT}" \
  --prepare \
  --run-star \
  --force \
  --threads "${THREADS}" \
  --star-input-format cbq \
  --star-yremove no \
  --star-out-samtype None

test -f "${CBQ_ROOT}/STAR_COMPLETED.txt"
test -s "${CBQ_ROOT}/stage/star_composite_cbq/cbq_manifest.tsv"
grep -F -- '--readFilesType Binseq PE' "${CBQ_ROOT}/RUN_STAR_COMPOSITE.sh" >/dev/null
if grep -F -- '--readFilesCommand' "${CBQ_ROOT}/RUN_STAR_COMPOSITE.sh" >/dev/null; then
    die "CBQ STAR script unexpectedly contains --readFilesCommand"
fi
if grep -F -- '--emitYNoYFastq' "${CBQ_ROOT}/RUN_STAR_COMPOSITE.sh" >/dev/null; then
    die "CBQ STAR script unexpectedly contains --emitYNoYFastq"
fi
grep -F $'Number of input reads |\t'"${READ_PAIRS}" "${CBQ_ROOT}/star_composite/run/Log.final.out" >/dev/null \
    || die "CBQ Log.final.out does not report ${READ_PAIRS} input reads"

if [[ "${RUN_FASTQ_PARITY}" == "1" ]]; then
    echo "Running OCM FASTQ parity smoke: ${FASTQ_ROOT}"
    OCM_PREP_REUSE_ROOT="${CBQ_ROOT}" \
    STAR_SUITE_ROOT="${STAR_SUITE_ROOT}" \
    RAW_DIR="${RAW_DIR}" \
    CONFIG="${CONFIG}" \
    GENOME_DIR="${GENOME_DIR}" \
    SOLO_CB_WHITELIST="${SOLO_CB_WHITELIST}" \
    STAR_BIN="${STAR_BIN}" \
    bash "${RECIPE}" \
      --read-pairs "${READ_PAIRS}" \
      --out-root "${FASTQ_ROOT}" \
      --prepare \
      --run-star \
      --force \
      --threads "${THREADS}" \
      --star-input-format fastq \
      --star-yremove no \
      --star-out-samtype None

    test -f "${FASTQ_ROOT}/STAR_COMPLETED.txt"
    grep -F $'Number of input reads |\t'"${READ_PAIRS}" "${FASTQ_ROOT}/star_composite/run/Log.final.out" >/dev/null \
        || die "FASTQ Log.final.out does not report ${READ_PAIRS} input reads"

    diff -qr "${FASTQ_ROOT}/star_composite/run/Solo.out" "${CBQ_ROOT}/star_composite/run/Solo.out"
    diff -qr "${FASTQ_ROOT}/star_composite/outs" "${CBQ_ROOT}/star_composite/outs"
fi

set +e
GATE_ROOT="${OUT_ROOT}/cbq_yremove_gate"
STAR_SUITE_ROOT="${STAR_SUITE_ROOT}" \
RAW_DIR="${RAW_DIR}" \
CONFIG="${CONFIG}" \
GENOME_DIR="${GENOME_DIR}" \
SOLO_CB_WHITELIST="${SOLO_CB_WHITELIST}" \
STAR_BIN="${STAR_BIN}" \
CBQ_ORDERED_ENCODER_BIN="${CBQ_ORDERED_ENCODER_BIN}" \
bash "${RECIPE}" \
  --read-pairs 10 \
  --out-root "${GATE_ROOT}" \
  --prepare \
  --force \
  --threads 2 \
  --star-input-format cbq \
  --star-yremove yes \
  --star-out-samtype None >"${OUT_ROOT}/cbq_yremove_gate.stdout" 2>"${OUT_ROOT}/cbq_yremove_gate.stderr"
gate_status=$?
set -e
if [[ "${gate_status}" -eq 0 ]]; then
    die "CBQ/Y-removal gate unexpectedly passed"
fi
grep -F 'STAR_YREMOVE=yes is FASTQ-only' "${OUT_ROOT}/cbq_yremove_gate.stderr" >/dev/null \
    || die "CBQ/Y-removal gate did not emit the expected error"

echo "PASS: OCM CBQ smoke completed at ${OUT_ROOT}"
