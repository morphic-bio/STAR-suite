#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
CR_INPUT_HELPER="${CR_INPUT_HELPER:-${SCRIPT_DIR}/ucsf_parity/render_star_inputs_from_cr_config.py}"
TRIM_QC_HELPER="${TRIM_QC_HELPER:-${SCRIPT_DIR}/lib/star_trim_qc.sh}"

# shellcheck disable=SC1090
source "${TRIM_QC_HELPER}"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR.release}"
PREPARE_VELOCYTO_MEX="${PREPARE_VELOCYTO_MEX:-${REPO_ROOT}/scripts/prepare_velocyto_mex.py}"
ALLOW_LEGACY_PREPARE_VELOCYTO="${ALLOW_LEGACY_PREPARE_VELOCYTO:-0}"
DOWNSTREAM_WRAPPER="${DOWNSTREAM_WRAPPER:-${REPO_ROOT}/scripts/run_scrna_downstream_gene_full_velocyto.sh}"
DATASET_ROOT="${DATASET_ROOT:-/mnt/pikachu/ucsf-perturb-seq-corrected}"
SAMPLE_MAP="${SAMPLE_MAP:-${DATASET_ROOT}/sample_fastq_guide_map.csv}"
FEATURE_REF="${FEATURE_REF:-/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv}"
GENOME_DIR="${GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
SOLO_CB_WHITELIST="${SOLO_CB_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-february-2018_TRU.txt}"
CR_WHITELIST="${CR_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt}"
OUT_ROOT="${OUT_ROOT:-/mnt/pikachu/ucsf-perturb-yremove_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-24}"
SAMPLES_CSV=""
ALL_SAMPLES=0
CR_CONFIG=""
CR_SAMPLE_ID=""
DOWNSAMPLE_READS=0
DOWNSAMPLE_SEED="${DOWNSAMPLE_SEED:-1}"
GLOBUS_SRC_ENDPOINT="${GLOBUS_SRC_ENDPOINT:-}"
GLOBUS_DST_ENDPOINT="${GLOBUS_DST_ENDPOINT:-}"
GLOBUS_DST_ROOT="${GLOBUS_DST_ROOT:-}"
GLOBUS_POLL_SECONDS="${GLOBUS_POLL_SECONDS:-30}"
STAR_TRIM_QC_ENABLE="${STAR_TRIM_QC_ENABLE:-0}"
STAR_TRIM_QC_MAX_READS="${STAR_TRIM_QC_MAX_READS:-250000}"
RUN_DOWNSTREAM=0
RUN_CELLBENDER=0
ADAPTIVE_FILTER=1
DOWNSTREAM_OUTPUT_NAME=""
CELLBENDER_CPU_CORES=""
CELLBENDER_GPU=0
DRY_RUN=0

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --samples CSV            Comma-separated sample IDs from sample_fastq_guide_map.csv
  --all-samples            Run every sample in the sample map
  --cr-config FILE         Use a Cell Ranger config.csv as the exact UCSF input source for one sample
  --cr-sample-id ID        Override the logical STAR sample ID when using --cr-config
  --downsample-reads N     Downsample each library (GEX and guides) to N read pairs before STAR
  --threads N              STAR threads per sample run (default: ${THREADS})
  --out-root DIR           Output root (default: ${OUT_ROOT})
  --dataset-root DIR       UCSF dataset root (default: ${DATASET_ROOT})
  --sample-map FILE        Sample map CSV (default: ${SAMPLE_MAP})
  --feature-ref FILE       CRISPR feature ref CSV (default: ${FEATURE_REF})
  --genome-dir DIR         STAR genomeDir (default: ${GENOME_DIR})
  --solo-cb-whitelist P    Solo CB whitelist (default: ${SOLO_CB_WHITELIST})
  --cr-whitelist P         CR whitelist (default: ${CR_WHITELIST})
  --star-bin PATH          STAR binary (default: ${STAR_BIN})
  --globus-src-endpoint ID Source Globus collection ID for pikachu outputs
  --globus-dst-endpoint ID Destination Globus collection ID
  --globus-dst-root PATH   Destination root path (e.g. /UCSF-perturb/run_name)
  --globus-poll-seconds N  Poll interval while waiting for transfer cleanup
  --run-downstream         Build local h5ad outputs after each sample run
  --adaptive-filter        Use adaptive n_genes and MT percentage downstream QC (default: on)
  --run-cellbender         Run CellBender in downstream processing
  --downstream-output-name NAME
                           Downstream output directory name under each sample root
  --cellbender-cpu-cores N CPU cores passed to CellBender (default: threads)
  --cellbender-gpu         Run CellBender with GPU support
  --trim-qc                Emit STAR read-level FastQC-like HTML/JSON reports
  --trim-qc-max-reads N    Limit reads sampled by trim-QC reporting
  --dry-run                Write manifests/commands but do not run STAR
  -h, --help               Show help
EOF
}

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"
}

die() {
  log "ERROR: $*" >&2
  exit 1
}

require_file() {
  local path="$1"
  [[ -f "${path}" ]] || die "Missing file: ${path}"
}

require_dir() {
  local path="$1"
  [[ -d "${path}" ]] || die "Missing directory: ${path}"
}

require_exec() {
  local path="$1"
  [[ -x "${path}" ]] || die "Missing executable: ${path}"
}

join_by_comma() {
  local IFS=,
  echo "$*"
}

realpath_if_possible() {
  local path="$1"
  if command -v realpath >/dev/null 2>&1; then
    realpath "${path}"
  else
    readlink -f "${path}"
  fi
}

free_bytes() {
  local path="$1"
  df -PB1 "${path}" | awk 'NR==2 {print $4}'
}

sum_bytes() {
  local total=0
  local path
  for path in "$@"; do
    total=$((total + $(stat -c %s "${path}")))
  done
  echo "${total}"
}

count_reads_gz() {
  local path="$1"
  pigz -dc "${path}" | awk 'END {print NR/4}'
}

find_fastq_files() {
  local dir="$1"
  local pattern="$2"
  find -L "${dir}" -maxdepth 1 -type f -name "${pattern}" | sort
}

validate_pairs() {
  local dir="$1"
  local label="$2"
  mapfile -t _r1_files < <(find_fastq_files "${dir}" '*_R1_001.fastq.gz')
  [[ "${#_r1_files[@]}" -gt 0 ]] || die "No R1 FASTQs found for ${label} under ${dir}"
  local r1
  local expected_r2
  for r1 in "${_r1_files[@]}"; do
    expected_r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    [[ -f "${expected_r2}" ]] || die "Missing R2 pair for ${r1}"
  done
}

lookup_sample_row() {
  local sample="$1"
  if [[ -f "${SAMPLE_MAP}" ]]; then
  python3 - <<'PY' "${SAMPLE_MAP}" "${sample}"
import csv, sys
sample_map, sample = sys.argv[1], sys.argv[2]
with open(sample_map, newline='') as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        if row["sample"] == sample:
            print("\t".join([
                row["sample"],
                row["gex_fastq_dir"],
                row["guide_fastq_dir"],
                row.get("feature_ref", ""),
            ]))
            break
    else:
        sys.exit(1)
PY
    return
  fi

  local sample_root="${DATASET_ROOT}/${sample}"
  local gex_dir="${sample_root}/GEX"
  local guide_dir="${sample_root}/guides"
  local sample_pf="${sample_root}/pf_multi_config.csv"
  [[ -d "${gex_dir}" && -d "${guide_dir}" ]] || return 1
  printf '%s\t%s\t%s\t%s\t%s\n' "${sample}" "${gex_dir}" "${guide_dir}" "${FEATURE_REF}" "${sample_pf}"
}

render_from_cr_config() {
  local config="$1"
  local sample_override="$2"
  local pf_out="$3"
  local args=(
    python3 "${CR_INPUT_HELPER}"
    --config "${config}"
    --pf-multi-out "${pf_out}"
  )
  if [[ -n "${sample_override}" ]]; then
    args+=(--sample-id "${sample_override}")
  fi
  args+=(--emit-env)
  "${args[@]}"
}

enumerate_samples() {
  if [[ -n "${CR_CONFIG}" ]]; then
    if [[ -n "${CR_SAMPLE_ID}" ]]; then
      printf '%s\n' "${CR_SAMPLE_ID}"
      return
    fi
    python3 "${CR_INPUT_HELPER}" --config "${CR_CONFIG}" --pf-multi-out /dev/null --emit-env \
      | awk -F= '/^SAMPLE_ID=/{sub(/^SAMPLE_ID='\''?/,""); gsub(/^'\''|'\''$/,""); print $0}'
    return
  fi
  if [[ "${ALL_SAMPLES}" == "1" ]]; then
    if [[ -f "${SAMPLE_MAP}" ]]; then
    python3 - <<'PY' "${SAMPLE_MAP}"
import csv, sys
with open(sys.argv[1], newline='') as handle:
    for row in csv.DictReader(handle):
        print(row["sample"])
PY
    else
      find "${DATASET_ROOT}" -mindepth 1 -maxdepth 1 -type d | while read -r sample_dir; do
        [[ -d "${sample_dir}/GEX" && -d "${sample_dir}/guides" ]] || continue
        basename "${sample_dir}"
      done | sort
    fi
    return
  fi
  [[ -n "${SAMPLES_CSV}" ]] || die "Specify --samples or --all-samples"
  tr ',' '\n' <<< "${SAMPLES_CSV}" | sed '/^$/d'
}

compute_lane_targets() {
  local target="$1"
  shift
  python3 - <<'PY' "${target}" "$@"
import sys
target = int(sys.argv[1])
counts = [int(x) for x in sys.argv[2:]]
total = sum(counts)
if total == 0:
    print("\n".join("0" for _ in counts))
    raise SystemExit
raw = [target * c / total for c in counts]
floors = [int(x) for x in raw]
remainder = target - sum(floors)
order = sorted(range(len(counts)), key=lambda i: (raw[i] - floors[i]), reverse=True)
for i in order[:remainder]:
    floors[i] += 1
for value in floors:
    print(value)
PY
}

write_globus_helper() {
  local helper="$1"
  local batch="$2"
  cat > "${helper}" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

SRC_ENDPOINT="${1:?Usage: $0 <src-endpoint> <dst-endpoint> <dst-root>}"
DST_ENDPOINT="${2:?Usage: $0 <src-endpoint> <dst-endpoint> <dst-root>}"
DST_ROOT="${3:?Usage: $0 <src-endpoint> <dst-endpoint> <dst-root>}"
BATCH_FILE="${4:-}"

if [[ -z "${BATCH_FILE}" ]]; then
  BATCH_FILE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/globus_batch.tsv"
fi

command -v globus >/dev/null 2>&1 || {
  echo "globus CLI not found on PATH" >&2
  exit 1
}
[[ -f "${BATCH_FILE}" ]] || {
  echo "Missing batch file: ${BATCH_FILE}" >&2
  exit 1
}

TMP_BATCH="$(mktemp)"
trap 'rm -f "${TMP_BATCH}"' EXIT
awk -F '\t' -v root="${DST_ROOT}" '{print $1 "\t" root "/" $2}' "${BATCH_FILE}" > "${TMP_BATCH}"
globus transfer "${SRC_ENDPOINT}" "${DST_ENDPOINT}" --batch "${TMP_BATCH}"
EOF
  chmod +x "${helper}"
}

globus_enabled() {
  [[ -n "${GLOBUS_SRC_ENDPOINT}" && -n "${GLOBUS_DST_ENDPOINT}" && -n "${GLOBUS_DST_ROOT}" ]]
}

globus_task_status() {
  local task_id="$1"
  globus task show "${task_id}" | awk -F': *' '/^Status:/ {print $2}'
}

cleanup_sample_large_outputs() {
  local sample="$1"
  local sample_root="$2"
  local run_dir="$3"
  local cleanup_file="${sample_root}/TRANSFER_CLEANED.ok"

  rm -f \
    "${run_dir}/Aligned.out_Y.bam" \
    "${run_dir}/Aligned.out_noY.bam"
  rm -rf "${run_dir}/y_separated"
  rm -rf "${sample_root}/staged_input"
  rm -rf "${sample_root}/tmp"

  {
    printf 'cleaned_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'sample=%s\n' "${sample}"
  } > "${cleanup_file}"
  log "Cleaned transferred BAM/YFASTQ outputs for ${sample}"
}

submit_sample_transfer() {
  local sample="$1"
  local sample_root="$2"
  local transfer_helper="$3"
  local submitted_file="${sample_root}/TRANSFER_SUBMITTED.ok"
  local task_file="${sample_root}/TRANSFER_TASK.txt"

  globus_enabled || return 0
  [[ -x "${transfer_helper}" ]] || die "Missing transfer helper for ${sample}: ${transfer_helper}"

  if [[ -f "${submitted_file}" ]]; then
    log "Transfer already submitted for ${sample}"
    return 0
  fi

  local task_output task_id
  task_output="$("${transfer_helper}" "${GLOBUS_SRC_ENDPOINT}" "${GLOBUS_DST_ENDPOINT}" "${GLOBUS_DST_ROOT}")"
  printf '%s\n' "${task_output}" > "${task_file}"
  task_id="$(awk '/Task ID:/ {print $3}' "${task_file}" | tail -n1)"
  [[ -n "${task_id}" ]] || die "Could not parse Globus task ID for ${sample}"
  {
    printf 'submitted_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'task_id=%s\n' "${task_id}"
    printf 'dst_root=%s\n' "${GLOBUS_DST_ROOT}"
  } > "${submitted_file}"
  log "Submitted Globus transfer for ${sample}: ${task_id}"
}

reap_completed_transfers() {
  globus_enabled || return 0

  local sample_dir sample sample_root run_dir submitted_file cleanup_file task_file task_id status
  while IFS= read -r sample_dir; do
    [[ -n "${sample_dir}" ]] || continue
    sample="$(basename "${sample_dir}")"
    sample_root="${sample_dir}"
    run_dir="${sample_root}/run"
    submitted_file="${sample_root}/TRANSFER_SUBMITTED.ok"
    cleanup_file="${sample_root}/TRANSFER_CLEANED.ok"
    task_file="${sample_root}/TRANSFER_TASK.txt"

    [[ -f "${submitted_file}" ]] || continue
    [[ -f "${cleanup_file}" ]] && continue
    [[ -f "${task_file}" ]] || die "Missing TRANSFER_TASK.txt for ${sample}"
    task_id="$(awk '/Task ID:/ {print $3}' "${task_file}" | tail -n1)"
    [[ -n "${task_id}" ]] || die "Could not parse task ID for ${sample}"

    status="$(globus_task_status "${task_id}")"
    case "${status}" in
      SUCCEEDED)
        cleanup_sample_large_outputs "${sample}" "${sample_root}" "${run_dir}"
        ;;
      ACTIVE|INACTIVE)
        ;;
      *)
        die "Globus transfer for ${sample} entered unexpected status ${status} (task ${task_id})"
        ;;
    esac
  done < <(find "${OUT_ROOT}/samples" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
}

pending_transfer_count() {
  find "${OUT_ROOT}/samples" -mindepth 1 -maxdepth 1 -type d \
    -exec test -f '{}/TRANSFER_SUBMITTED.ok' ';' \
    ! -exec test -f '{}/TRANSFER_CLEANED.ok' ';' \
    -print 2>/dev/null | wc -l
}

ensure_space_for_sample() {
  local sample="$1"
  local required_bytes="$2"
  local avail_bytes

  while true; do
    reap_completed_transfers
    avail_bytes="$(free_bytes "${OUT_ROOT}")"
    if (( avail_bytes >= required_bytes )); then
      printf '%s\n' "${avail_bytes}"
      return 0
    fi

    if ! globus_enabled; then
      die "Insufficient free space for ${sample}: need ${required_bytes} bytes (4x source FASTQ bytes), have ${avail_bytes}"
    fi

    if [[ "$(pending_transfer_count)" == "0" ]]; then
      die "Insufficient free space for ${sample}: need ${required_bytes} bytes (4x source FASTQ bytes), have ${avail_bytes}, and no pending transfers remain to free space"
    fi

    log "Waiting for Globus transfers to complete before starting ${sample}: need ${required_bytes} bytes, have ${avail_bytes}"
    sleep "${GLOBUS_POLL_SECONDS}"
  done
}

wait_for_all_transfers_and_cleanup() {
  globus_enabled || return 0
  while [[ "$(pending_transfer_count)" != "0" ]]; do
    reap_completed_transfers
    if [[ "$(pending_transfer_count)" == "0" ]]; then
      break
    fi
    log "Waiting for remaining Globus transfers to complete"
    sleep "${GLOBUS_POLL_SECONDS}"
  done
  reap_completed_transfers
}

finalize_sample_outputs() {
  local sample="$1"
  local sample_root="$2"
  local run_dir="$3"
  local transfer_batch="$4"
  local transfer_helper="$5"
  local done_file="$6"

  [[ -f "${run_dir}/Aligned.out_Y.bam" ]] || die "Missing Y BAM for ${sample}"
  [[ -f "${run_dir}/Aligned.out_noY.bam" ]] || die "Missing noY BAM for ${sample}"
  [[ -d "${run_dir}/y_separated" ]] || die "Missing y_separated FASTQs for ${sample}"

  : > "${transfer_batch}"
  {
    printf '%s\t%s\n' "${run_dir}/Aligned.out_Y.bam" "${sample}/Aligned.out_Y.bam"
    printf '%s\t%s\n' "${run_dir}/Aligned.out_noY.bam" "${sample}/Aligned.out_noY.bam"
  } >> "${transfer_batch}"
  while IFS= read -r fastq; do
    printf '%s\t%s\n' "${fastq}" "${sample}/y_separated/$(basename "${fastq}")" >> "${transfer_batch}"
  done < <(find "${run_dir}/y_separated" -maxdepth 1 -type f -name '*.fastq.gz' | sort)
  while IFS= read -r qc_file; do
    [[ -f "${qc_file}" ]] || continue
    printf '%s\t%s\n' "${qc_file}" "${sample}/$(basename "${qc_file}")" >> "${transfer_batch}"
  done < <(star_trim_qc_list_outputs "${run_dir}")
  write_globus_helper "${transfer_helper}" "${transfer_batch}"

  printf 'completed_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "${done_file}"
  log "Completed sample ${sample}: ${sample_root}"
}

package_velocyto_outputs() {
  local sample="$1"
  local run_dir="$2"
  local raw_dir="${run_dir}/outs/raw_velocyto_feature_bc_matrix"
  local filtered_dir="${run_dir}/outs/filtered_velocyto_feature_bc_matrix"
  local manifest="${run_dir}/outs/velocyto_feature_bc_matrix_manifest.json"

  if [[ -d "${raw_dir}" && -d "${filtered_dir}" && -f "${manifest}" ]]; then
    log "Native Velocyto packaged outputs present for ${sample}"
    return
  fi

  if [[ "${ALLOW_LEGACY_PREPARE_VELOCYTO}" == "1" ]]; then
    [[ -f "${PREPARE_VELOCYTO_MEX}" ]] || die "Missing helper: ${PREPARE_VELOCYTO_MEX}"
    log "WARNING: using legacy prepare_velocyto_mex.py fallback for ${sample}"
    python3 "${PREPARE_VELOCYTO_MEX}" --run-dir "${run_dir}" >/dev/null
  fi

  [[ -d "${raw_dir}" ]] || die "Missing packaged raw velocyto MEX for ${sample}"
  [[ -d "${filtered_dir}" ]] || die "Missing packaged filtered velocyto MEX for ${sample}"
  [[ -f "${manifest}" ]] || die "Missing velocyto manifest for ${sample}"
}

run_downstream_for_sample() {
  local sample="$1"
  local sample_root="$2"
  local run_dir="$3"

  [[ "${RUN_DOWNSTREAM}" == "1" ]] || return 0
  [[ -x "${DOWNSTREAM_WRAPPER}" ]] || die "Missing downstream wrapper: ${DOWNSTREAM_WRAPPER}"

  local output_name="${DOWNSTREAM_OUTPUT_NAME}"
  if [[ -z "${output_name}" ]]; then
    if [[ "${RUN_CELLBENDER}" == "1" ]]; then
      output_name="downstream_genefull_velocyto_cellbender"
    else
      output_name="downstream_genefull_velocyto"
    fi
  fi

  local output_dir="${sample_root}/${output_name}"
  if [[ -f "${output_dir}/summary.txt" ]]; then
    log "Downstream outputs already present for ${sample}: ${output_dir}"
    return 0
  fi

  local cellbender_cores="${CELLBENDER_CPU_CORES:-${THREADS}}"
  local args=(
    "${DOWNSTREAM_WRAPPER}"
    --run-dir "${run_dir}"
    --output-dir "${output_dir}"
  )
  if [[ "${ADAPTIVE_FILTER}" == "1" ]]; then
    args+=(--adaptive-filter)
  fi
  if [[ "${RUN_CELLBENDER}" == "1" ]]; then
    args+=(--run-cellbender --cellbender-cpu-cores "${cellbender_cores}")
    if [[ "${CELLBENDER_GPU}" == "1" ]]; then
      args+=(--cellbender-gpu)
    fi
  fi

  "${args[@]}"
}

run_downstream_batch() {
  [[ "${RUN_DOWNSTREAM}" == "1" ]] || return 0

  local sample sample_root run_dir
  while IFS= read -r sample; do
    [[ -n "${sample}" ]] || continue
    sample_root="${OUT_ROOT}/samples/${sample}"
    run_dir="${sample_root}/run"
    [[ -d "${run_dir}" ]] || die "Missing run directory for downstream sample ${sample}: ${run_dir}"
    log "Starting downstream ${sample}"
    run_downstream_for_sample "${sample}" "${sample_root}" "${run_dir}"
  done < <(enumerate_samples)
}

downsample_library_dir() {
  local src_dir="$1"
  local dst_dir="$2"
  local target_reads="$3"
  local seed_base="$4"
  local label="$5"

  mkdir -p "${dst_dir}"
  mapfile -t r1_files < <(find_fastq_files "${src_dir}" '*_R1_001.fastq.gz')
  [[ "${#r1_files[@]}" -gt 0 ]] || die "No R1 FASTQs found under ${src_dir}"

  local counts_file="${dst_dir}/.lane_counts.tsv"
  local targets_file="${dst_dir}/.lane_targets.tsv"
  local done_file="${dst_dir}/.downsample_complete"
  if [[ -f "${done_file}" ]]; then
    log "Downsample already complete for ${label}: ${dst_dir}"
    return
  fi

  : > "${counts_file}"
  local total_reads=0
  local lane_idx=0
  local r1
  local r2
  local count
  for r1 in "${r1_files[@]}"; do
    r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    count="$(count_reads_gz "${r1}")"
    printf '%s\t%s\t%s\n' "$(basename "${r1}")" "$(basename "${r2}")" "${count}" >> "${counts_file}"
    total_reads=$((total_reads + count))
    lane_idx=$((lane_idx + 1))
  done

  if (( total_reads <= target_reads )); then
    log "${label}: total reads ${total_reads} <= target ${target_reads}; hardlinking originals"
    for r1 in "${r1_files[@]}"; do
      r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
      ln -f "$(realpath_if_possible "${r1}")" "${dst_dir}/$(basename "${r1}")"
      ln -f "$(realpath_if_possible "${r2}")" "${dst_dir}/$(basename "${r2}")"
    done
    printf 'source_dir=%s\ntotal_reads=%s\ntarget_reads=%s\nmode=hardlink\n' "${src_dir}" "${total_reads}" "${target_reads}" > "${done_file}"
    return
  fi

  mapfile -t counts < <(awk -F '\t' '{print $3}' "${counts_file}")
  mapfile -t targets < <(compute_lane_targets "${target_reads}" "${counts[@]}")
  : > "${targets_file}"
  for idx in "${!targets[@]}"; do
    printf '%s\t%s\n' "$(basename "${r1_files[$idx]}")" "${targets[$idx]}" >> "${targets_file}"
  done

  for idx in "${!r1_files[@]}"; do
    r1="${r1_files[$idx]}"
    r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    local target_lane="${targets[$idx]}"
    local lane_seed=$((seed_base + idx))
    if (( target_lane <= 0 )); then
      : | gzip -c > "${dst_dir}/$(basename "${r1}")"
      : | gzip -c > "${dst_dir}/$(basename "${r2}")"
      continue
    fi
    seqtk sample -s"${lane_seed}" "${r1}" "${target_lane}" | pigz -p 8 > "${dst_dir}/$(basename "${r1}")"
    seqtk sample -s"${lane_seed}" "${r2}" "${target_lane}" | pigz -p 8 > "${dst_dir}/$(basename "${r2}")"
  done

  printf 'source_dir=%s\ntotal_reads=%s\ntarget_reads=%s\nmode=seqtk-count\n' "${src_dir}" "${total_reads}" "${target_reads}" > "${done_file}"
}

run_sample() {
  local sample="$1"
  local gex_dir guide_dir feature_ref map_sample sample_pf_config=""
  local cr_gene_expression_reference=""
  local cr_gene_expression_chemistry=""
  local -a gex_dirs=()
  local -a feature_dirs=()
  local -a all_input_dirs=()

  if [[ -n "${CR_CONFIG}" ]]; then
    local cr_env_file="${OUT_ROOT}/samples/${sample}/cr_config_inputs.env"
    local pre_pf="${OUT_ROOT}/samples/${sample}/pf_multi_config.prestage.csv"
    mkdir -p "$(dirname "${cr_env_file}")"
    render_from_cr_config "${CR_CONFIG}" "${CR_SAMPLE_ID:-${sample}}" "${pre_pf}" > "${cr_env_file}"
    # shellcheck disable=SC1090
    source "${cr_env_file}"
    map_sample="${SAMPLE_ID}"
    IFS=',' read -r -a gex_dirs <<< "${GEX_FASTQ_DIRS:-}"
    IFS=',' read -r -a feature_dirs <<< "${FEATURE_FASTQ_DIRS:-}"
    [[ "${#gex_dirs[@]}" -gt 0 ]] || die "CR config rendered no GEX fastq dirs for ${sample}"
    [[ "${#feature_dirs[@]}" -gt 0 ]] || die "CR config rendered no feature fastq dirs for ${sample}"
    gex_dir="$(python3 - <<'PY' "${GEX_R1}"
import sys
print(sys.argv[1].split(',')[0].rsplit('/', 1)[0])
PY
)"
    guide_dir="${feature_dirs[0]}"
    feature_ref="${CR_FEATURE_REF}"
    cr_gene_expression_reference="${CR_GENE_EXPRESSION_REFERENCE:-}"
    cr_gene_expression_chemistry="${CR_GENE_EXPRESSION_CHEMISTRY:-}"
  else
    local row
    row="$(lookup_sample_row "${sample}")" || die "Sample not found in map: ${sample}"
    IFS=$'\t' read -r map_sample gex_dir guide_dir feature_ref sample_pf_config <<< "${row}"
    [[ -n "${feature_ref}" ]] || feature_ref="${FEATURE_REF}"
    gex_dirs=("${gex_dir}")
    feature_dirs=("${guide_dir}")
  fi

  all_input_dirs=("${gex_dirs[@]}" "${feature_dirs[@]}")
  local dir_index=0
  for gex_dir in "${gex_dirs[@]}"; do
    require_dir "${gex_dir}"
    validate_pairs "${gex_dir}" "${sample} GEX library ${dir_index}"
    dir_index=$((dir_index + 1))
  done
  dir_index=0
  for guide_dir in "${feature_dirs[@]}"; do
    require_dir "${guide_dir}"
    validate_pairs "${guide_dir}" "${sample} feature library ${dir_index}"
    dir_index=$((dir_index + 1))
  done
  if [[ -n "${feature_ref}" ]]; then
    require_file "${feature_ref}"
  fi

  local input_dir
  mapfile -t src_fastqs < <(
    for input_dir in "${all_input_dirs[@]}"; do
      find_fastq_files "${input_dir}" '*.fastq.gz'
    done
  )
  [[ "${#src_fastqs[@]}" -gt 0 ]] || die "No FASTQs found for ${sample}"

  local source_bytes
  source_bytes="$(sum_bytes "${src_fastqs[@]}")"
  local required_bytes=$((source_bytes * 4))
  local avail_bytes
  avail_bytes="$(ensure_space_for_sample "${sample}" "${required_bytes}")"

  local sample_root="${OUT_ROOT}/samples/${sample}"
  local stage_root="${sample_root}/staged_input"
  local stage_gex="${stage_root}/GEX/${sample}"
  local stage_guides="${stage_root}/guides/${sample}"
  local run_dir="${sample_root}/run"
  local pf_config="${sample_root}/pf_multi_config.csv"
  local sample_manifest="${sample_root}/RUN_MANIFEST.txt"
  local command_script="${sample_root}/RUN_COMMAND.sh"
  local transfer_batch="${sample_root}/globus_batch.tsv"
  local transfer_helper="${sample_root}/RUN_GLOBUS_TRANSFER.sh"
  local done_file="${sample_root}/RUN_COMPLETE.ok"
  mkdir -p "${sample_root}" "${run_dir}" "$(dirname "${stage_gex}")" "$(dirname "${stage_guides}")"

  if [[ -f "${done_file}" ]]; then
    log "Sample already complete, skipping: ${sample}"
    return
  fi

  if (( DOWNSAMPLE_READS > 0 )); then
    log "Downsampling ${sample} libraries to ${DOWNSAMPLE_READS} read pairs each"
    if [[ -n "${CR_CONFIG}" ]]; then
      local -a rewrite_args=()
      local idx=0
      for gex_dir in "${gex_dirs[@]}"; do
        local stage_dir="${stage_root}/GEX/$(printf '%02d' "${idx}")_$(basename "${gex_dir}")"
        downsample_library_dir "${gex_dir}" "${stage_dir}" "${DOWNSAMPLE_READS}" "$((DOWNSAMPLE_SEED + idx))" "${sample} GEX ${idx}"
        rewrite_args+=(--rewrite-fastq-dir "${gex_dir}=${stage_dir}")
        idx=$((idx + 1))
      done
      idx=0
      for guide_dir in "${feature_dirs[@]}"; do
        local stage_dir="${stage_root}/FEATURE/$(printf '%02d' "${idx}")_$(basename "${guide_dir}")"
        downsample_library_dir "${guide_dir}" "${stage_dir}" "${DOWNSAMPLE_READS}" "$((DOWNSAMPLE_SEED + 1000 + idx))" "${sample} feature ${idx}"
        rewrite_args+=(--rewrite-fastq-dir "${guide_dir}=${stage_dir}")
        idx=$((idx + 1))
      done
      python3 "${CR_INPUT_HELPER}" \
        --config "${CR_CONFIG}" \
        --sample-id "${CR_SAMPLE_ID:-${sample}}" \
        --pf-multi-out "${pf_config}" \
        "${rewrite_args[@]}" \
        --emit-env > "${cr_env_file}"
      # shellcheck disable=SC1090
      source "${cr_env_file}"
      stage_gex="${stage_root}/GEX"
      stage_guides="${stage_root}/FEATURE"
    else
      downsample_library_dir "${gex_dir}" "${stage_gex}" "${DOWNSAMPLE_READS}" "${DOWNSAMPLE_SEED}" "${sample} GEX"
      downsample_library_dir "${guide_dir}" "${stage_guides}" "${DOWNSAMPLE_READS}" "$((DOWNSAMPLE_SEED + 1000))" "${sample} guides"
    fi
  else
    if [[ -n "${CR_CONFIG}" ]]; then
      python3 "${CR_INPUT_HELPER}" \
        --config "${CR_CONFIG}" \
        --sample-id "${CR_SAMPLE_ID:-${sample}}" \
        --pf-multi-out "${pf_config}" \
        --emit-env > "${cr_env_file}"
      # shellcheck disable=SC1090
      source "${cr_env_file}"
      stage_gex="$(IFS=,; echo "${GEX_FASTQ_DIRS}")"
      stage_guides="$(IFS=,; echo "${FEATURE_FASTQ_DIRS}")"
    else
      stage_gex="${gex_dir}"
      stage_guides="${guide_dir}"
    fi
  fi

  if [[ -z "${CR_CONFIG}" ]]; then
    if [[ -n "${sample_pf_config}" && -f "${sample_pf_config}" && "${DOWNSAMPLE_READS}" -eq 0 ]]; then
      cp -f "${sample_pf_config}" "${pf_config}"
    else
    cat > "${pf_config}" <<EOF
[libraries]
fastqs,sample,library_type,feature_types
${stage_gex},${sample},Gene Expression,Gene Expression
${stage_guides},${sample},CRISPR Guide Capture,CRISPR Guide Capture

[feature]
ref,${feature_ref}
EOF
    fi
  fi

  local reads_r2 reads_r1
  if [[ -n "${CR_CONFIG}" ]]; then
    reads_r2="${GEX_R2}"
    reads_r1="${GEX_R1}"
    if [[ -n "${FEATURE_R2:-}" ]]; then
      reads_r2="${reads_r2},${FEATURE_R2}"
      reads_r1="${reads_r1},${FEATURE_R1}"
    fi
  else
    mapfile -t gex_r2 < <(find_fastq_files "${stage_gex}" '*_R2_001.fastq.gz')
    mapfile -t gex_r1 < <(find_fastq_files "${stage_gex}" '*_R1_001.fastq.gz')
    mapfile -t guide_r2 < <(find_fastq_files "${stage_guides}" '*_R2_001.fastq.gz')
    mapfile -t guide_r1 < <(find_fastq_files "${stage_guides}" '*_R1_001.fastq.gz')
    local all_r2=("${gex_r2[@]}" "${guide_r2[@]}")
    local all_r1=("${gex_r1[@]}" "${guide_r1[@]}")
    reads_r2="$(join_by_comma "${all_r2[@]}")"
    reads_r1="$(join_by_comma "${all_r1[@]}")"
  fi

  local cmd=(
    "${STAR_BIN}"
    --runThreadN "${THREADS}"
    --genomeDir "${GENOME_DIR}"
    --readFilesIn "${reads_r2}" "${reads_r1}"
    --readFilesCommand zcat
    --outFileNamePrefix "${run_dir}/"
    --outTmpDir "${sample_root}/tmp"
    --outSAMtype BAM Unsorted
    --emitNoYBAM yes
    --emitYNoYFastq yes
    --emitYNoYFastqCompression gz
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
    --soloFeatures GeneFull Velocyto
    --soloCrGexFeature genefull
    --soloCrMultimapRescue yes
    --pfMultiConfig "${pf_config}"
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
  )
  if [[ -n "${feature_ref}" ]]; then
    cmd+=(--crFeatureRef "${feature_ref}")
  fi
  star_trim_qc_append_args cmd "${run_dir}"

  {
    printf 'sample=%s\n' "${sample}"
    printf 'date_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'star_bin=%s\n' "${STAR_BIN}"
    printf 'genome_dir=%s\n' "${GENOME_DIR}"
    printf 'feature_ref=%s\n' "${feature_ref}"
    printf 'solo_cb_whitelist=%s\n' "${SOLO_CB_WHITELIST}"
    printf 'cr_whitelist=%s\n' "${CR_WHITELIST}"
    printf 'threads=%s\n' "${THREADS}"
    printf 'star_solo_nonflex_hash_bridge=%s\n' "disabled_for_bam_velocyto"
    printf 'downsample_reads=%s\n' "${DOWNSAMPLE_READS}"
    printf 'source_fastq_bytes=%s\n' "${source_bytes}"
    printf 'required_free_bytes=%s\n' "${required_bytes}"
    printf 'free_bytes_at_start=%s\n' "${avail_bytes}"
    printf 'gex_dir=%s\n' "${gex_dir}"
    printf 'guide_dir=%s\n' "${guide_dir}"
    printf 'gex_dirs=%s\n' "$(join_by_comma "${gex_dirs[@]}")"
    printf 'feature_dirs=%s\n' "$(join_by_comma "${feature_dirs[@]}")"
    printf 'staged_gex_dir=%s\n' "${stage_gex}"
    printf 'staged_guide_dir=%s\n' "${stage_guides}"
    printf 'run_dir=%s\n' "${run_dir}"
    if [[ -n "${CR_CONFIG}" ]]; then
      printf 'cr_config=%s\n' "${CR_CONFIG}"
      printf 'cr_sample_id=%s\n' "${CR_SAMPLE_ID:-${sample}}"
      printf 'cr_gene_expression_reference=%s\n' "${cr_gene_expression_reference}"
      printf 'cr_gene_expression_chemistry=%s\n' "${cr_gene_expression_chemistry}"
    fi
  } > "${sample_manifest}"
  star_trim_qc_write_manifest "${sample_manifest}" "${run_dir}"

  {
    echo '#!/usr/bin/env bash'
    echo 'set -euo pipefail'
    echo 'unset STAR_SOLO_NONFLEX_HASH_BRIDGE 2>/dev/null || true'
    printf '%q ' "${cmd[@]}"
    printf '\n'
  } > "${command_script}"
  chmod +x "${command_script}"

  if (( DRY_RUN == 1 )); then
    log "DRY_RUN: prepared ${sample} at ${sample_root}"
    return
  fi

  if [[ -f "${run_dir}/Aligned.out_Y.bam" && -f "${run_dir}/Aligned.out_noY.bam" && -d "${run_dir}/y_separated" ]]; then
    log "Existing STAR outputs found for ${sample}; regenerating transfer helpers only"
    package_velocyto_outputs "${sample}" "${run_dir}"
    finalize_sample_outputs "${sample}" "${sample_root}" "${run_dir}" "${transfer_batch}" "${transfer_helper}" "${done_file}"
    submit_sample_transfer "${sample}" "${sample_root}" "${transfer_helper}"
    return
  fi

  (
    unset STAR_SOLO_NONFLEX_HASH_BRIDGE 2>/dev/null || true
    "${cmd[@]}"
  )
  package_velocyto_outputs "${sample}" "${run_dir}"
  finalize_sample_outputs "${sample}" "${sample_root}" "${run_dir}" "${transfer_batch}" "${transfer_helper}" "${done_file}"
  submit_sample_transfer "${sample}" "${sample_root}" "${transfer_helper}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples) SAMPLES_CSV="$2"; shift 2 ;;
    --all-samples) ALL_SAMPLES=1; shift ;;
    --cr-config) CR_CONFIG="$2"; shift 2 ;;
    --cr-sample-id) CR_SAMPLE_ID="$2"; shift 2 ;;
    --downsample-reads) DOWNSAMPLE_READS="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --out-root) OUT_ROOT="$2"; shift 2 ;;
    --dataset-root) DATASET_ROOT="$2"; shift 2 ;;
    --sample-map) SAMPLE_MAP="$2"; shift 2 ;;
    --feature-ref) FEATURE_REF="$2"; shift 2 ;;
    --genome-dir) GENOME_DIR="$2"; shift 2 ;;
    --solo-cb-whitelist) SOLO_CB_WHITELIST="$2"; shift 2 ;;
    --cr-whitelist) CR_WHITELIST="$2"; shift 2 ;;
    --star-bin) STAR_BIN="$2"; shift 2 ;;
    --globus-src-endpoint) GLOBUS_SRC_ENDPOINT="$2"; shift 2 ;;
    --globus-dst-endpoint) GLOBUS_DST_ENDPOINT="$2"; shift 2 ;;
    --globus-dst-root) GLOBUS_DST_ROOT="$2"; shift 2 ;;
    --globus-poll-seconds) GLOBUS_POLL_SECONDS="$2"; shift 2 ;;
    --run-downstream) RUN_DOWNSTREAM=1; shift ;;
    --adaptive-filter) ADAPTIVE_FILTER=1; shift ;;
    --run-cellbender) RUN_CELLBENDER=1; RUN_DOWNSTREAM=1; shift ;;
    --downstream-output-name) DOWNSTREAM_OUTPUT_NAME="$2"; shift 2 ;;
    --cellbender-cpu-cores) CELLBENDER_CPU_CORES="$2"; shift 2 ;;
    --cellbender-gpu) CELLBENDER_GPU=1; RUN_CELLBENDER=1; RUN_DOWNSTREAM=1; shift ;;
    --trim-qc) STAR_TRIM_QC_ENABLE=1; shift ;;
    --trim-qc-max-reads) STAR_TRIM_QC_MAX_READS="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

require_exec "${STAR_BIN}"
require_dir "${DATASET_ROOT}"
if [[ -n "${CR_CONFIG}" ]]; then
  require_file "${CR_CONFIG}"
  require_exec "${CR_INPUT_HELPER}"
else
  if [[ -f "${SAMPLE_MAP}" ]]; then
    require_file "${SAMPLE_MAP}"
  fi
  require_file "${FEATURE_REF}"
fi
require_dir "${GENOME_DIR}"
require_file "${SOLO_CB_WHITELIST}"
require_file "${CR_WHITELIST}"
if [[ "${ALLOW_LEGACY_PREPARE_VELOCYTO}" == "1" ]]; then
  require_file "${PREPARE_VELOCYTO_MEX}"
fi
command -v seqtk >/dev/null 2>&1 || die "seqtk not found on PATH"
command -v pigz >/dev/null 2>&1 || die "pigz not found on PATH"
if globus_enabled; then
  command -v globus >/dev/null 2>&1 || die "globus CLI not found on PATH"
fi
if [[ "${RUN_DOWNSTREAM}" == "1" ]]; then
  [[ -x "${DOWNSTREAM_WRAPPER}" ]] || die "Downstream wrapper not executable: ${DOWNSTREAM_WRAPPER}"
fi

mkdir -p "${OUT_ROOT}"
if globus_enabled && (( DRY_RUN == 0 )); then
  globus mkdir "${GLOBUS_DST_ENDPOINT}:${GLOBUS_DST_ROOT}" >/dev/null 2>&1 || true
fi

TOP_MANIFEST="${OUT_ROOT}/RUNS.tsv"
if [[ ! -f "${TOP_MANIFEST}" ]]; then
  printf 'sample\tstatus\tsample_root\n' > "${TOP_MANIFEST}"
fi

log "Output root: ${OUT_ROOT}"
log "Dataset root: ${DATASET_ROOT}"
log "STAR binary: ${STAR_BIN}"
log "Threads: ${THREADS}"
log "Downsample reads per library: ${DOWNSAMPLE_READS}"
if [[ "${RUN_DOWNSTREAM}" == "1" ]]; then
  log "Downstream wrapper: ${DOWNSTREAM_WRAPPER}"
  if [[ "${RUN_CELLBENDER}" == "1" ]]; then
    log "CellBender CPU cores: ${CELLBENDER_CPU_CORES:-${THREADS}}"
  fi
fi
if globus_enabled; then
  log "Globus source endpoint: ${GLOBUS_SRC_ENDPOINT}"
  log "Globus destination endpoint: ${GLOBUS_DST_ENDPOINT}"
  log "Globus destination root: ${GLOBUS_DST_ROOT}"
fi
if [[ -n "${CR_CONFIG}" ]]; then
  log "Cell Ranger config: ${CR_CONFIG}"
fi

while IFS= read -r sample; do
  [[ -n "${sample}" ]] || continue
  log "Starting sample ${sample}"
  run_sample "${sample}"
  printf '%s\t%s\t%s\n' "${sample}" "done" "${OUT_ROOT}/samples/${sample}" >> "${TOP_MANIFEST}"
done < <(enumerate_samples)

run_downstream_batch
wait_for_all_transfers_and_cleanup
log "All requested samples complete"
