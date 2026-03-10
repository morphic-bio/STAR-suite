#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

EXTERNAL_ENV="${SCRIPT_DIR}/external_fixtures_env.sh"
if [[ -f "${EXTERNAL_ENV}" ]]; then
  # shellcheck disable=SC1091
  source "${EXTERNAL_ENV}"
fi

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
SPLIT_HELPER="${SCRIPT_DIR}/split_fastq_consecutive_pe.py"
THREADS="${BULK_MULTISAMPLE_THREADS:-1}"
PAIRS_PER_SAMPLE="${BULK_MULTISAMPLE_PAIRS_PER_SAMPLE:-5000}"
GENOME_DIR="${BULK_MULTISAMPLE_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
INPUT_R1="${BULK_MULTISAMPLE_R1:-/storage/PE/downsampled/21033-09-01-13-01_S1_L007_R1_001.fastq.gz}"
INPUT_R2="${BULK_MULTISAMPLE_R2:-/storage/PE/downsampled/21033-09-01-13-01_S1_L007_R2_001.fastq.gz}"
OUTDIR="${BULK_MULTISAMPLE_OUTDIR:-${SCRIPT_DIR}/bulk_pe_multisample_equivalence_output_$(date +%Y%m%d_%H%M%S)}"
SAMPLE_A="${BULK_MULTISAMPLE_SAMPLE_A:-bulk_multi_A_S1_L001}"
SAMPLE_B="${BULK_MULTISAMPLE_SAMPLE_B:-bulk_multi_B_S2_L001}"
BATCH_SAMPLE_A="${BULK_MULTISAMPLE_BATCH_SAMPLE_A:-$(printf '%s\n' "${SAMPLE_A}" | sed -E 's/_S[0-9]+(_L[0-9]{3})?$//')}"
BATCH_SAMPLE_B="${BULK_MULTISAMPLE_BATCH_SAMPLE_B:-$(printf '%s\n' "${SAMPLE_B}" | sed -E 's/_S[0-9]+(_L[0-9]{3})?$//')}"
EXPECTED_FORMAT="${BULK_MULTISAMPLE_EXPECTED_FORMAT:-IU}"
QUALITY_CUTOFF="${BULK_MULTISAMPLE_QUALITY:-20}"
MIN_LENGTH="${BULK_MULTISAMPLE_MIN_LENGTH:-20}"
ADAPTER_R1="${BULK_MULTISAMPLE_ADAPTER_R1:-AGATCGGAAGAGCACACGTCTGAACTCCAGTCA}"
ADAPTER_R2="${BULK_MULTISAMPLE_ADAPTER_R2:-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT}"

log() {
  printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*"
}

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

require_file() {
  [[ -f "$1" ]] || fail "Missing file: $1"
}

find_one() {
  local dir="$1"
  local pattern="$2"
  mapfile -t matches < <(find "$dir" -maxdepth 1 -type f -name "$pattern" | sort)
  [[ "${#matches[@]}" -eq 1 ]] || fail "Expected exactly one match for ${pattern} under ${dir}, found ${#matches[@]}"
  printf '%s\n' "${matches[0]}"
}

compare_plain_file() {
  local label="$1"
  local left="$2"
  local right="$3"
  cmp -s "$left" "$right" || fail "Mismatch for ${label}: ${left} vs ${right}"
  printf '%s\tPASS\tplain\n' "${label}" >> "${OUTDIR}/comparison.tsv"
}

compare_optional_plain_file() {
  local label="$1"
  local left="$2"
  local right="$3"
  local left_exists=0
  local right_exists=0
  [[ -f "$left" ]] && left_exists=1
  [[ -f "$right" ]] && right_exists=1
  if [[ "${left_exists}" -ne "${right_exists}" ]]; then
    fail "Presence mismatch for ${label}: ${left} vs ${right}"
  fi
  if [[ "${left_exists}" -eq 1 ]]; then
    compare_plain_file "${label}" "${left}" "${right}"
  else
    printf '%s\tSKIP\toptional-missing\n' "${label}" >> "${OUTDIR}/comparison.tsv"
  fi
}

compare_gzip_contents() {
  local label="$1"
  local left="$2"
  local right="$3"
  cmp -s <(gzip -dc "$left") <(gzip -dc "$right") || fail "Mismatch for ${label}: ${left} vs ${right}"
  printf '%s\tPASS\tgzip-decoded\n' "${label}" >> "${OUTDIR}/comparison.tsv"
}

compare_bam_canonical() {
  local label="$1"
  local left="$2"
  local right="$3"
  diff -u <(samtools view -h "$left" | awk '!/^@PG\t/ && !/^@CO\t/') \
          <(samtools view -h "$right" | awk '!/^@PG\t/ && !/^@CO\t/') >/dev/null \
    || fail "Mismatch for ${label}: ${left} vs ${right}"
  printf '%s\tPASS\tbam-canonical\n' "${label}" >> "${OUTDIR}/comparison.tsv"
}

assert_y_reads_present() {
  local sample="$1"
  local fastq="$2"
  local count
  count="$(gzip -dc "${fastq}" | awk 'END{print NR/4}')"
  [[ "${count}" -gt 0 ]] || fail "${sample} emitted zero Y reads; increase BULK_MULTISAMPLE_PAIRS_PER_SAMPLE"
  printf '%s:Y-count\tPASS\t%s\n' "${sample}" "${count}" >> "${OUTDIR}/comparison.tsv"
}

compare_sample() {
  local batch_sample="$1"
  local seq_sample="$2"
  local batch_counts="${OUTDIR}/batch/counts/${batch_sample}"
  local batch_align="${OUTDIR}/batch/alignments/${batch_sample}"
  local batch_y="${OUTDIR}/batch/y_separated/${batch_sample}"
  local seq_dir="${OUTDIR}/sequential/${seq_sample}"
  local seq_y="${seq_dir}/y_separated"
  local batch_y_r1 batch_y_r2 batch_noy_r1 batch_noy_r2 seq_y_r1 seq_y_r2 seq_noy_r1 seq_noy_r2
  local batch_log seq_log

  require_file "${batch_counts}/${batch_sample}.quant.sf"
  require_file "${batch_counts}/${batch_sample}.quant.genes.sf"
  require_file "${seq_dir}/quant.sf"
  require_file "${seq_dir}/quant.genes.sf"

  compare_plain_file "${seq_sample}:quant.sf" "${batch_counts}/${batch_sample}.quant.sf" "${seq_dir}/quant.sf"
  compare_plain_file "${seq_sample}:quant.genes.sf" "${batch_counts}/${batch_sample}.quant.genes.sf" "${seq_dir}/quant.genes.sf"
  compare_optional_plain_file "${seq_sample}:quant.genes.tximport.sf" \
    "${batch_counts}/${batch_sample}.quant.genes.tximport.sf" \
    "${seq_dir}/quant.genes.tximport.sf"

  compare_bam_canonical "${seq_sample}:Aligned.sortedByCoord.out.bam" \
    "${batch_align}/${batch_sample}_Aligned.sortedByCoord.out.bam" \
    "${seq_dir}/Aligned.sortedByCoord.out.bam"
  compare_bam_canonical "${seq_sample}:Aligned.toTranscriptome.out.bam" \
    "${batch_align}/${batch_sample}_Aligned.toTranscriptome.out.bam" \
    "${seq_dir}/Aligned.toTranscriptome.out.bam"
  compare_bam_canonical "${seq_sample}:Y.bam" \
    "${batch_y}/${batch_sample}_Y.bam" \
    "${seq_dir}/Aligned.sortedByCoord.out_Y.bam"
  compare_bam_canonical "${seq_sample}:noY.bam" \
    "${batch_y}/${batch_sample}_noY.bam" \
    "${seq_dir}/Aligned.sortedByCoord.out_noY.bam"

  batch_y_r1="$(find_one "${batch_y}" '*_Y_R1_001.fastq.gz')"
  batch_y_r2="$(find_one "${batch_y}" '*_Y_R2_001.fastq.gz')"
  batch_noy_r1="$(find_one "${batch_y}" '*_noY_R1_001.fastq.gz')"
  batch_noy_r2="$(find_one "${batch_y}" '*_noY_R2_001.fastq.gz')"
  seq_y_r1="$(find_one "${seq_y}" '*_Y_R1_001.fastq.gz')"
  seq_y_r2="$(find_one "${seq_y}" '*_Y_R2_001.fastq.gz')"
  seq_noy_r1="$(find_one "${seq_y}" '*_noY_R1_001.fastq.gz')"
  seq_noy_r2="$(find_one "${seq_y}" '*_noY_R2_001.fastq.gz')"

  assert_y_reads_present "${seq_sample}" "${batch_y_r1}"

  compare_gzip_contents "${seq_sample}:Y_R1" "${batch_y_r1}" "${seq_y_r1}"
  compare_gzip_contents "${seq_sample}:Y_R2" "${batch_y_r2}" "${seq_y_r2}"
  compare_gzip_contents "${seq_sample}:noY_R1" "${batch_noy_r1}" "${seq_noy_r1}"
  compare_gzip_contents "${seq_sample}:noY_R2" "${batch_noy_r2}" "${seq_noy_r2}"

  compare_plain_file "${seq_sample}:SJ.out.tab" \
    "${batch_align}/${batch_sample}_SJ.out.tab" \
    "${seq_dir}/SJ.out.tab"

  batch_log="$(find_one "${batch_align}" '*_Log.out')"
  seq_log="${seq_dir}/Log.out"
  require_file "${batch_log}"
  require_file "${seq_log}"
  grep -q "Detected library format: ${EXPECTED_FORMAT}" "${batch_log}" \
    || fail "${seq_sample} batch log did not detect ${EXPECTED_FORMAT}"
  grep -q "Detected library format: ${EXPECTED_FORMAT}" "${seq_log}" \
    || fail "${seq_sample} sequential log did not detect ${EXPECTED_FORMAT}"
  printf '%s:autodetect\tPASS\t%s\n' "${seq_sample}" "${EXPECTED_FORMAT}" >> "${OUTDIR}/comparison.tsv"

}

for cmd in python3 samtools gzip; do
  command -v "${cmd}" >/dev/null 2>&1 || fail "Missing required command: ${cmd}"
done
if [[ ! -x "${STAR_BIN}" ]]; then
  fail "Missing required STAR binary: ${STAR_BIN}"
fi
require_file "${INPUT_R1}"
require_file "${INPUT_R2}"
require_file "${SPLIT_HELPER}"
[[ -d "${GENOME_DIR}" ]] || fail "Missing genome dir: ${GENOME_DIR}"

rm -rf "${OUTDIR}"
mkdir -p "${OUTDIR}/inputs" "${OUTDIR}/batch" "${OUTDIR}/sequential"
printf 'label\tstatus\tmode\n' > "${OUTDIR}/comparison.tsv"

log "Creating consecutive pseudo-samples under ${OUTDIR}/inputs"
python3 "${SPLIT_HELPER}" \
  --r1 "${INPUT_R1}" \
  --r2 "${INPUT_R2}" \
  --outdir "${OUTDIR}/inputs" \
  --sample-a "${SAMPLE_A}" \
  --sample-b "${SAMPLE_B}" \
  --pairs-per-sample "${PAIRS_PER_SAMPLE}"

A_R1="${OUTDIR}/inputs/${SAMPLE_A}_R1_001.fastq.gz"
A_R2="${OUTDIR}/inputs/${SAMPLE_A}_R2_001.fastq.gz"
B_R1="${OUTDIR}/inputs/${SAMPLE_B}_R1_001.fastq.gz"
B_R2="${OUTDIR}/inputs/${SAMPLE_B}_R2_001.fastq.gz"

log "Running multisample batch STAR"
"${STAR_BIN}" \
  --runMode alignReads \
  --runThreadN "${THREADS}" \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${A_R1},${B_R1}" "${A_R2},${B_R2}" \
  --readFilesCommand zcat \
  --trimCutadapt Yes \
  --trimCutadaptQuality "${QUALITY_CUTOFF}" \
  --trimCutadaptMinLength "${MIN_LENGTH}" \
  --trimCutadaptAdapter "${ADAPTER_R1}" "${ADAPTER_R2}" \
  --outSAMtype BAM SortedByCoordinate \
  --outBAMsortMethod samtools \
  --emitNoYBAM yes \
  --keepBAM yes \
  --emitYNoYFastq yes \
  --emitYNoYFastqCompression gz \
  --quantMode TranscriptomeSAM TranscriptVB \
  --quantVBgcBias 1 \
  --quantVBLibType A \
  --batchMode 1 \
  --outFileNamePrefixAuto 1 \
  --outFileNamePrefix "${OUTDIR}/batch/" \
  --outTmpDir "${OUTDIR}/tmp_batch"

for sample in "${SAMPLE_A}" "${SAMPLE_B}"; do
  sample_r1="${OUTDIR}/inputs/${sample}_R1_001.fastq.gz"
  sample_r2="${OUTDIR}/inputs/${sample}_R2_001.fastq.gz"
  sample_out="${OUTDIR}/sequential/${sample}"
  mkdir -p "${sample_out}"
  log "Running sequential STAR for ${sample}"
  "${STAR_BIN}" \
    --runMode alignReads \
    --runThreadN "${THREADS}" \
    --genomeDir "${GENOME_DIR}" \
    --readFilesIn "${sample_r1}" "${sample_r2}" \
    --readFilesCommand zcat \
    --trimCutadapt Yes \
    --trimCutadaptQuality "${QUALITY_CUTOFF}" \
    --trimCutadaptMinLength "${MIN_LENGTH}" \
    --trimCutadaptAdapter "${ADAPTER_R1}" "${ADAPTER_R2}" \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortMethod samtools \
    --emitNoYBAM yes \
    --keepBAM yes \
    --emitYNoYFastq yes \
    --emitYNoYFastqCompression gz \
    --quantMode TranscriptomeSAM TranscriptVB \
    --quantVBgcBias 1 \
    --quantVBLibType A \
    --outFileNamePrefix "${sample_out}/" \
    --outTmpDir "${OUTDIR}/tmp_${sample}"
done

compare_sample "${BATCH_SAMPLE_A}" "${SAMPLE_A}"
compare_sample "${BATCH_SAMPLE_B}" "${SAMPLE_B}"

cat > "${OUTDIR}/SUMMARY.txt" <<EOF
Bulk PE multisample equivalence: PASS
Input R1: ${INPUT_R1}
Input R2: ${INPUT_R2}
Pairs per sample: ${PAIRS_PER_SAMPLE}
Threads: ${THREADS}
Sample A: ${SAMPLE_A}
Sample B: ${SAMPLE_B}
Expected TranscriptVB format: ${EXPECTED_FORMAT}
Comparison table: ${OUTDIR}/comparison.tsv
EOF

log "PASS: multisample batch and sequential bulk PE outputs matched"
