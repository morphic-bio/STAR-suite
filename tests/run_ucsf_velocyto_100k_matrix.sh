#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

SAMPLE="${UCSF_VELOCYTO_100K_SAMPLE:-EBs2_2}"
SOURCE_SAMPLE_ROOT="${UCSF_VELOCYTO_100K_SOURCE_SAMPLE_ROOT:-${SCRIPT_DIR}/ucsf_corrected_production_100k_output_20260329_074313/samples/${SAMPLE}}"
STAGED_ROOT="${UCSF_VELOCYTO_100K_STAGED_ROOT:-${SOURCE_SAMPLE_ROOT}/staged_input}"
PF_CONFIG="${UCSF_VELOCYTO_100K_PF_CONFIG:-${SOURCE_SAMPLE_ROOT}/pf_multi_config.csv}"
STAR_BIN="${UCSF_VELOCYTO_100K_STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR.release}"
GENOME_DIR="${UCSF_VELOCYTO_100K_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
SOLO_CB_WHITELIST="${UCSF_VELOCYTO_100K_SOLO_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt}"
CR_WHITELIST="${UCSF_VELOCYTO_100K_CR_WHITELIST:-${SOLO_CB_WHITELIST}}"
FEATURE_REF="${UCSF_VELOCYTO_100K_FEATURE_REF:-/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv}"
THREADS="${UCSF_VELOCYTO_100K_THREADS:-24}"
OUT_ROOT="${UCSF_VELOCYTO_100K_OUT_ROOT:-${SCRIPT_DIR}/ucsf_velocyto_100k_matrix_output_$(date +%Y%m%d_%H%M%S)}"

die() {
  echo "ERROR: $*" >&2
  exit 1
}

join_by_comma() {
  local IFS=,
  echo "$*"
}

mtx_nnz() {
  local mtx="$1"
  if [[ -f "${mtx}" ]]; then
    awk 'NR==3{print $3}' "${mtx}"
  else
    echo "NA"
  fi
}

summary_metric() {
  local summary="$1"
  local key="$2"
  if [[ -f "${summary}" ]]; then
    awk -F, -v k="${key}" '$1==k{print $2}' "${summary}"
  else
    echo "NA"
  fi
}

find_fastqs() {
  local root="$1"
  local pattern="$2"
  find "${root}" -type f -name "${pattern}" | sort
}

[[ -x "${STAR_BIN}" ]] || die "Missing STAR binary: ${STAR_BIN}"
[[ -d "${STAGED_ROOT}/GEX" ]] || die "Missing staged GEX root: ${STAGED_ROOT}/GEX"
[[ -d "${STAGED_ROOT}/guides" ]] || die "Missing staged guide root: ${STAGED_ROOT}/guides"
[[ -f "${PF_CONFIG}" ]] || die "Missing pfMultiConfig: ${PF_CONFIG}"

mkdir -p "${OUT_ROOT}"
RESULTS_TSV="${OUT_ROOT}/results.tsv"

mapfile -t GEX_R2 < <(find_fastqs "${STAGED_ROOT}/GEX" '*_R2_001.fastq.gz')
mapfile -t GEX_R1 < <(find_fastqs "${STAGED_ROOT}/GEX" '*_R1_001.fastq.gz')
mapfile -t GUIDE_R2 < <(find_fastqs "${STAGED_ROOT}/guides" '*_R2_001.fastq.gz')
mapfile -t GUIDE_R1 < <(find_fastqs "${STAGED_ROOT}/guides" '*_R1_001.fastq.gz')

(( ${#GEX_R2[@]} > 0 )) || die "No staged GEX R2 FASTQs found"
(( ${#GEX_R1[@]} > 0 )) || die "No staged GEX R1 FASTQs found"
(( ${#GUIDE_R2[@]} > 0 )) || die "No staged guide R2 FASTQs found"
(( ${#GUIDE_R1[@]} > 0 )) || die "No staged guide R1 FASTQs found"

READS_R2="$(join_by_comma "${GEX_R2[@]}" "${GUIDE_R2[@]}")"
READS_R1="$(join_by_comma "${GEX_R1[@]}" "${GUIDE_R1[@]}")"

cat > "${RESULTS_TSV}" <<'EOF'
case_id	features	out_sam	y_mode	exit_code	velocyto_unique	velocyto_unique_multi	spliced_nnz	unspliced_nnz	ambiguous_nnz	run_dir
EOF

run_case() {
  local case_id="$1"
  local features="$2"
  local out_sam="$3"
  local y_mode="$4"

  local case_root="${OUT_ROOT}/${case_id}"
  local run_dir="${case_root}/run"
  local tmp_dir="${case_root}/tmp"
  local log_file="${case_root}/run.log"
  local cmd_file="${case_root}/RUN_COMMAND.sh"
  mkdir -p "${case_root}"
  rm -rf "${run_dir}" "${tmp_dir}"

  local -a cmd=(
    "${STAR_BIN}"
    --runThreadN "${THREADS}"
    --genomeDir "${GENOME_DIR}"
    --readFilesIn "${READS_R2}" "${READS_R1}"
    --readFilesCommand zcat
    --outFileNamePrefix "${run_dir}/"
    --outTmpDir "${tmp_dir}"
    --clipAdapterType CellRanger4
    --clip3pPolyG yes
    --alignEndsType Local
    --chimSegmentMin 1000000
    --soloType CB_UMI_Simple
    --soloCBstart 1
    --soloCBlen 16
    --soloUMIstart 17
    --soloUMIlen 12
    --soloBarcodeReadLength 0
    --soloCBwhitelist "${SOLO_CB_WHITELIST}"
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
    --soloUMIfiltering MultiGeneUMI_CR
    --soloUMIdedup 1MM_CR
    --soloMultiMappers Unique
    --soloCellFilter EmptyDrops_CR
    --soloCbUbRequireTogether no
    --soloStrand Forward
    --soloFeatures ${features}
    --soloCrGexFeature genefull
    --soloCrMultimapRescue yes
    --pfMultiConfig "${PF_CONFIG}"
    --crChemistry auto
    --crOutputChemistry TRU
    --crWhitelist "${CR_WHITELIST}"
    --crMinUmi 3
    --crAssignMaxHamming 1
    --crAssignFeatureOffset -1
    --crAssignLimitSearch -1
    --crAssignMinCounts 0
    --crAssignMaxBarcodeMismatches 5
    --crAssignFeatureN 0
    --crAssignBarcodeN 1
    --crAssignConsumerThreads -1
    --crAssignSearchThreads 1
    --crAssignSkipQcOutputs 1
    --dynamicThreadInterface 1
    --dynamicThreadConstMapPermits "${THREADS}"
    --dynamicThreadTelemetry 1
    --crFeatureRef "${FEATURE_REF}"
  )

  case "${out_sam}" in
    bam)
      cmd+=(--outSAMtype BAM Unsorted)
      ;;
    none)
      cmd+=(--outSAMtype None)
      ;;
    *)
      die "Unknown out_sam mode: ${out_sam}"
      ;;
  esac

  case "${y_mode}" in
    y)
      cmd+=(--emitNoYBAM yes --emitYNoYFastq yes --emitYNoYFastqCompression gz)
      ;;
    noy)
      cmd+=(--emitNoYBAM no --emitYNoYFastq no)
      ;;
    *)
      die "Unknown y_mode: ${y_mode}"
      ;;
  esac

  {
    echo '#!/usr/bin/env bash'
    echo 'set -euo pipefail'
    echo 'unset STAR_SOLO_NONFLEX_HASH_BRIDGE 2>/dev/null || true'
    printf '%q ' "${cmd[@]}"
    printf '\n'
  } > "${cmd_file}"
  chmod +x "${cmd_file}"

  local exit_code=0
  (
    cd "${REPO_ROOT}"
    unset STAR_SOLO_NONFLEX_HASH_BRIDGE 2>/dev/null || true
    "${cmd[@]}"
  ) > "${log_file}" 2>&1 || exit_code=$?

  local summary="${run_dir}/Solo.out/Velocyto/Summary.csv"
  local raw_root="${run_dir}/Solo.out/Velocyto/raw"
  local unique unique_multi spliced_nnz unspliced_nnz ambiguous_nnz
  unique="$(summary_metric "${summary}" 'Reads Mapped to Velocyto: Unique Velocyto')"
  unique_multi="$(summary_metric "${summary}" 'Reads Mapped to Velocyto: Unique+Multiple Velocyto')"
  spliced_nnz="$(mtx_nnz "${raw_root}/spliced.mtx")"
  unspliced_nnz="$(mtx_nnz "${raw_root}/unspliced.mtx")"
  ambiguous_nnz="$(mtx_nnz "${raw_root}/ambiguous.mtx")"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${case_id}" "${features}" "${out_sam}" "${y_mode}" "${exit_code}" \
    "${unique}" "${unique_multi}" "${spliced_nnz}" "${unspliced_nnz}" "${ambiguous_nnz}" "${run_dir}" \
    >> "${RESULTS_TSV}"
}

echo "=== UCSF 100K Velocyto matrix ==="
echo "Sample: ${SAMPLE}"
echo "Staged root: ${STAGED_ROOT}"
echo "Output root: ${OUT_ROOT}"
echo "Threads: ${THREADS}"

run_case "nogene_bam_y" "GeneFull Velocyto" "bam" "y"
run_case "gene_bam_y" "Gene GeneFull Velocyto" "bam" "y"
run_case "nogene_bam_noy" "GeneFull Velocyto" "bam" "noy"
run_case "gene_bam_noy" "Gene GeneFull Velocyto" "bam" "noy"
run_case "nogene_none_noy" "GeneFull Velocyto" "none" "noy"
run_case "gene_none_noy" "Gene GeneFull Velocyto" "none" "noy"

cat "${RESULTS_TSV}"
