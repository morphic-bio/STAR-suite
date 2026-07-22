#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
THREADS="${THREADS:-2}"
READ_LIMIT="${READ_LIMIT:-2000}"
WORKDIR="${WORKDIR:-/tmp/flex_tiny_public_smoke_$(date +%Y%m%d_%H%M%S)_$$}"
DATA_DIR="${WORKDIR}/data"
TINY_ROOT="${DATA_DIR}/public_10x_tiny"
ASSET_DIR="${WORKDIR}/assets"
FILTERED_REF_WORK="${WORKDIR}/refwork"
INDEX_DIR="${WORKDIR}/star_index"
OUT_BASE="${WORKDIR}/run_root"
RUN_ID="run"
DEMO_ROOT="${WORKDIR}/codespaces_demo"
CACHE_ROOT="${DEMO_ROOT}/cache"
INDEX_ROOT="${DEMO_ROOT}/indices"
RUN_ROOT_DIR="${DEMO_ROOT}/runs"

die() {
  echo "ERROR: $*" >&2
  exit 1
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "missing required command: $1"
}

need_cmd git
need_cmd git-lfs
need_cmd python3

[[ -x "${STAR_BIN}" ]] || die "STAR binary not found: ${STAR_BIN}"

mkdir -p "${WORKDIR}" "${DATA_DIR}" "${CACHE_ROOT}" "${INDEX_ROOT}" "${RUN_ROOT_DIR}"

DEMO_ROOT="${DEMO_ROOT}" \
CACHE_ROOT="${CACHE_ROOT}" \
DATA_ROOT="${DATA_DIR}" \
INDEX_ROOT="${INDEX_ROOT}" \
RUN_ROOT="${RUN_ROOT_DIR}" \
  "${REPO_ROOT}/scripts/codespaces/fetch_public_10x_tiny.sh" --outdir "${TINY_ROOT}"

TINY_FASTQ_DIR="${TINY_ROOT}/cellranger-tiny-fastq"
TINY_REF_ROOT="${TINY_ROOT}/cellranger-tiny-ref"
[[ -d "${TINY_FASTQ_DIR}" ]] || die "missing tiny FASTQ dir: ${TINY_FASTQ_DIR}"
[[ -d "${TINY_REF_ROOT}" ]] || die "missing tiny ref dir: ${TINY_REF_ROOT}"

python3 "${REPO_ROOT}/scripts/codespaces/generate_tiny_flex_demo.py" \
  --tiny-fastq-dir "${TINY_FASTQ_DIR}" \
  --tiny-ref-root "${TINY_REF_ROOT}" \
  --outdir "${ASSET_DIR}" \
  --read-limit "${READ_LIMIT}"

WHITELIST="${ASSET_DIR}/whitelist.txt"
CONFIG="${ASSET_DIR}/config.csv"
PROBE_CATALOG="${ASSET_DIR}/sample_probe_catalog.tsv"
[[ -f "${WHITELIST}" ]] || die "missing whitelist: ${WHITELIST}"
[[ -f "${CONFIG}" ]] || die "missing config: ${CONFIG}"
[[ -f "${PROBE_CATALOG}" ]] || die "missing sample probe catalog: ${PROBE_CATALOG}"

GTF_PATH="${TINY_REF_ROOT}/genes/genes.gtf"
if [[ ! -f "${GTF_PATH}" ]]; then
  GTF_PATH="${TINY_REF_ROOT}/genes/genes.gtf.gz"
fi
[[ -f "${GTF_PATH}" ]] || die "missing tiny GTF under ${TINY_REF_ROOT}/genes"

"${REPO_ROOT}/flex/scripts/build_filtered_reference.sh" \
  --probe-set "${ASSET_DIR}/probe_set.csv" \
  --base-fasta "${TINY_REF_ROOT}/fasta/genome.fa" \
  --base-gtf "${GTF_PATH}" \
  --work-dir "${FILTERED_REF_WORK}" \
  --skip-filter \
  --quiet

"${REPO_ROOT}/flex/scripts/make_filtered_star_index.sh" \
  --filtered-reference "${FILTERED_REF_WORK}/filtered_reference" \
  --output-dir "${INDEX_DIR}" \
  --threads "${THREADS}" \
  --sa-index-n-bases 11 \
  --star-bin "${STAR_BIN}"

"${REPO_ROOT}/scripts/run_flex_cr_config.sh" \
  --cr-config "${CONFIG}" \
  --genome-dir "${INDEX_DIR}" \
  --cb-whitelist "${WHITELIST}" \
  --solo-cb-start 1 \
  --solo-cb-len 16 \
  --solo-umi-start 17 \
  --solo-umi-len 10 \
  --sample-probe-catalog "${PROBE_CATALOG}" \
  --sample-probe-offset 68 \
  --out-base "${OUT_BASE}" \
  --run-id "${RUN_ID}" \
  --threads "${THREADS}"

RUN_ROOT="${OUT_BASE}/${RUN_ID}"
[[ -f "${RUN_ROOT}/RUN_MANIFEST.txt" ]] || die "missing RUN_MANIFEST.txt"
[[ -f "${RUN_ROOT}/Aligned.out.bam" ]] || die "missing BAM output"
[[ -f "${RUN_ROOT}/Log.final.out" ]] || die "missing Log.final.out"
[[ -f "${RUN_ROOT}/Solo.out/Barcodes.stats" ]] || die "missing Barcodes.stats"

grep -F "sample_probe_offset=68" "${RUN_ROOT}/RUN_MANIFEST.txt" >/dev/null || die "manifest missing sample_probe_offset"
grep -F "sample_probe_catalog=${PROBE_CATALOG}" "${RUN_ROOT}/RUN_MANIFEST.txt" >/dev/null || die "manifest missing sample_probe_catalog"
grep -F "Enabled Flex pipeline with production defaults" "${RUN_ROOT}/Log.out" >/dev/null || die "log missing Flex enablement"
grep -F "SampleDetector initialized successfully" "${RUN_ROOT}/Log.out" >/dev/null || die "log missing sample detector init"

echo "PASS: public tiny Flex smoke"
echo "Workdir: ${WORKDIR}"
