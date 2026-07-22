#!/usr/bin/env bash
# STAR pf-multi smoke for CAT-ATAC trimodal guide arm (split-read + atac2gex map).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
SPLIT_BIN="${SPLIT_BIN:-${REPO_ROOT}/core/features/process_features/tests/test_catatac_split_read}"
OUT_DIR="${CATATAC_PF_MULTI_SMOKE_OUT:-${REPO_ROOT}/tests/catatac_guide_arm_pf_multi_smoke_output}"
FIXTURE_DIR="${CATATAC_GUIDE_FIXTURE:-/mnt/pikachu/catatac_gse288996/guide_redump/fixture}"
FEATURE_REF="${CATATAC_GUIDE_FEATURE_REF:-${REPO_ROOT}/core/features/process_features/feature_lists/catatac_crispri_guide_capture.csv}"
ATAC_WHITELIST="${CATATAC_ATAC_WHITELIST:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt}"
ATAC2GEX="${CATATAC_ATAC2GEX_MAP:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/atac2gex.tsv}"
GEX_WHITELIST="${CATATAC_GEX_WHITELIST:-/mnt/pikachu/GEX_whitelist/737K-arc-v1.txt}"
GENOME_DIR="${CATATAC_PF_MULTI_GENOME:-/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/star}"
LIBRARY_ID="${CATATAC_PF_MULTI_LIBRARY_ID:-catatac_guide_fixture}"
MAX_READS="${CATATAC_GUIDE_MAX_READS:-1000}"
THREADS="${CATATAC_PF_MULTI_THREADS:-8}"

for path in "${FIXTURE_DIR}" "${FEATURE_REF}" "${ATAC_WHITELIST}" "${ATAC2GEX}" "${GEX_WHITELIST}" "${GENOME_DIR}"; do
  if [[ ! -e "${path}" ]]; then
    echo "ERROR: missing required path: ${path}" >&2
    exit 1
  fi
done

for fq in guide_R1.fastq.gz guide_R2.fastq.gz guide_R3.fastq.gz; do
  if [[ ! -f "${FIXTURE_DIR}/${fq}" ]]; then
    echo "ERROR: missing fixture FASTQ: ${FIXTURE_DIR}/${fq}" >&2
    exit 1
  fi
done

if [[ ! -x "${STAR_BIN}" ]]; then
  make -C "${REPO_ROOT}/core/legacy/source" -j8 STAR
fi
if [[ ! -x "${SPLIT_BIN}" ]]; then
  make -C "${REPO_ROOT}/core/features/process_features" tests/test_catatac_split_read
fi

rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

NATIVE_OUT="${OUT_DIR}/native_baseline"
STAR_OUT="${OUT_DIR}/star"
PF_CONFIG="${OUT_DIR}/pf_multi_config.csv"
CS1="CAAGTTGATAACGGACTAGCC"
CS2="CAAGTTGTAAACGGACTAGCC"

cat > "${PF_CONFIG}" <<EOF
[libraries]
fastqs,sample,library_type,feature_types,star_library_id,star_feature_ref,star_whitelist,star_layout,star_barcode_read,star_barcode_format,star_umi_read,star_umi_start,star_umi_length,star_feature_read,star_capture_read,star_capture_sequences,star_capture_max_hamming,star_barcode_output_map,star_feature_search_mode
${FIXTURE_DIR},catatac_fixture,CRISPR Guide Capture,CRISPR Guide Capture,${LIBRARY_ID},${FEATURE_REF},${ATAC_WHITELIST},catatac_guide,R2,bc:8:23:-,R1,0,12,R3,R1,${CS1}|${CS2},0,${ATAC2GEX},free
EOF

echo "=== CAT-ATAC guide-arm pf-multi smoke ==="
echo "  STAR: ${STAR_BIN}"
echo "  Fixture: ${FIXTURE_DIR}"
echo "  Max reads: ${MAX_READS}"
echo "  Output: ${OUT_DIR}"
echo ""

export CATATAC_GUIDE_MAX_READS="${MAX_READS}"

echo "Running native split-read baseline..."
"${SPLIT_BIN}" "${FIXTURE_DIR}" "${NATIVE_OUT}" "${ATAC_WHITELIST}" "${FEATURE_REF}" \
  >"${OUT_DIR}/native_split_read.log" 2>&1

echo "Running STAR pf-multi..."
"${STAR_BIN}" \
  --runMode alignReads \
  --runThreadN "${THREADS}" \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${FIXTURE_DIR}/guide_R2.fastq.gz" "${FIXTURE_DIR}/guide_R1.fastq.gz" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${STAR_OUT}/" \
  --outSAMtype None \
  --soloType None \
  --readMapNumber "${MAX_READS}" \
  --pfMultiConfig "${PF_CONFIG}" \
  --crWhitelist "${ATAC_WHITELIST}" \
  --crAssignConsumerThreads -1 \
  --crAssignSearchThreads 1 \
  2>&1 | tee "${OUT_DIR}/star_stdout.log"

ASSIGN_BASE="${STAR_OUT}/cr_assign/CRISPR_Guide_Capture/${LIBRARY_ID}"
if [[ ! -d "${ASSIGN_BASE}" ]]; then
  echo "ERROR: missing pf-multi assign output: ${ASSIGN_BASE}" >&2
  exit 1
fi

python3 "${SCRIPT_DIR}/catatac_guide_pf_multi_verify.py" \
  --native-dir "${NATIVE_OUT}/sample" \
  --pf-multi-assign-dir "${ASSIGN_BASE}" \
  --log-out "${STAR_OUT}/Log.out" \
  --native-log "${OUT_DIR}/native_split_read.log" \
  --output-map "${ATAC2GEX}" \
  --gex-whitelist "${GEX_WHITELIST}" \
  --expected-features 54 \
  --report "${OUT_DIR}/pf_multi_verify.json"

echo "CAT-ATAC guide-arm pf-multi smoke completed: ${OUT_DIR}"
