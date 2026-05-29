#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
ENCODER_BIN="${CBQ_ORDERED_ENCODER_BIN:-${REPO_ROOT}/core/legacy/source/cbq_ordered_encoder}"
THREADS="${THREADS:-2}"
READ_LIMIT="${READ_LIMIT:-500}"
DEFAULT_WORKDIR="/tmp/star_suite_cbq_flex_tiny_public_smoke_$(date +%Y%m%d_%H%M%S)_$$"
WORKDIR="${WORKDIR:-${OUT_ROOT:-${DEFAULT_WORKDIR}}}"
DATA_DIR="${WORKDIR}/data"
TINY_ROOT="${DATA_DIR}/public_10x_tiny"
ASSET_DIR="${WORKDIR}/assets"
FILTERED_REF_WORK="${WORKDIR}/refwork"
INDEX_DIR="${WORKDIR}/star_index"
OUT_BASE="${WORKDIR}/run_root"
CBQ_DIR="${WORKDIR}/cbq"
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
[[ -x "${ENCODER_BIN}" ]] || die "cbq_ordered_encoder not found: ${ENCODER_BIN}"

mkdir -p "${WORKDIR}" "${DATA_DIR}" "${CACHE_ROOT}" "${INDEX_ROOT}" "${RUN_ROOT_DIR}" "${CBQ_DIR}"

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
R1_FASTQ="${ASSET_DIR}/gex/tiny_flex/tinyflex_S1_L001_R1_001.fastq.gz"
R2_FASTQ="${ASSET_DIR}/gex/tiny_flex/tinyflex_S1_L001_R2_001.fastq.gz"
CBQ_FILE="${CBQ_DIR}/tinyflex_R2_R1.cbq"

[[ -f "${WHITELIST}" ]] || die "missing whitelist: ${WHITELIST}"
[[ -f "${CONFIG}" ]] || die "missing config: ${CONFIG}"
[[ -f "${PROBE_CATALOG}" ]] || die "missing sample probe catalog: ${PROBE_CATALOG}"
[[ -f "${R1_FASTQ}" ]] || die "missing R1 FASTQ: ${R1_FASTQ}"
[[ -f "${R2_FASTQ}" ]] || die "missing R2 FASTQ: ${R2_FASTQ}"

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

"${ENCODER_BIN}" \
  --readFilesIn "${R2_FASTQ}" "${R1_FASTQ}" \
  --outFile "${CBQ_FILE}"

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
  --run-id fastq \
  --threads "${THREADS}"

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
  --input-format cbq \
  --cbq-file "${CBQ_FILE}" \
  --out-base "${OUT_BASE}" \
  --run-id cbq \
  --threads "${THREADS}"

FASTQ_RUN="${OUT_BASE}/fastq"
CBQ_RUN="${OUT_BASE}/cbq"

for run_dir in "${FASTQ_RUN}" "${CBQ_RUN}"; do
  [[ -f "${run_dir}/RUN_MANIFEST.txt" ]] || die "missing RUN_MANIFEST.txt in ${run_dir}"
  [[ -f "${run_dir}/Aligned.out.bam" ]] || die "missing BAM output in ${run_dir}"
  [[ -f "${run_dir}/Log.final.out" ]] || die "missing Log.final.out in ${run_dir}"
  [[ -f "${run_dir}/Solo.out/Barcodes.stats" ]] || die "missing Barcodes.stats in ${run_dir}"
done

grep -F "input_format=cbq" "${CBQ_RUN}/RUN_MANIFEST.txt" >/dev/null || die "CBQ manifest missing input_format"
grep -F "cbq_file=${CBQ_FILE}" "${CBQ_RUN}/RUN_MANIFEST.txt" >/dev/null || die "CBQ manifest missing cbq_file"
grep -E "readFilesType[[:space:]]+Binseq[[:space:]]+PE" \
  "${CBQ_RUN}/Log.out" >/dev/null || die "CBQ log missing Binseq PE input type"
grep -F -- "--flex yes" "${CBQ_RUN}/Log.out" >/dev/null || die "CBQ log missing Flex mode"
grep -F "Flex pipeline: not active (CBQ/Binseq input uses the standard STAR CBQ adapter path)" \
  "${CBQ_RUN}/Log.out" >/dev/null || die "CBQ log missing FlexPipeline CBQ guard"

cmp -s "${FASTQ_RUN}/Solo.out/Barcodes.stats" "${CBQ_RUN}/Solo.out/Barcodes.stats" || {
  diff -u "${FASTQ_RUN}/Solo.out/Barcodes.stats" "${CBQ_RUN}/Solo.out/Barcodes.stats" >&2 || true
  die "FASTQ-vs-CBQ Barcodes.stats mismatch"
}

if command -v samtools >/dev/null 2>&1; then
  samtools view "${FASTQ_RUN}/Aligned.out.bam" > "${WORKDIR}/fastq.Aligned.body.sam"
  samtools view "${CBQ_RUN}/Aligned.out.bam" > "${WORKDIR}/cbq.Aligned.body.sam"
  cmp -s "${WORKDIR}/fastq.Aligned.body.sam" "${WORKDIR}/cbq.Aligned.body.sam" || {
    diff -u "${WORKDIR}/fastq.Aligned.body.sam" "${WORKDIR}/cbq.Aligned.body.sam" | head -200 >&2 || true
    die "FASTQ-vs-CBQ BAM body mismatch"
  }
else
  echo "WARN: samtools not found; skipped BAM body parity check" >&2
fi

echo "PASS: FLEX CBQ tiny public smoke"
echo "Workdir: ${WORKDIR}"
