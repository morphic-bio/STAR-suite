#!/usr/bin/env bash
set -euo pipefail

# Production runner for the NW/SLAM-Seq paired-end panel.
# Runs one STAR invocation per sample so Y/noY BAM/FASTQ outputs can be uploaded
# after each sample without waiting for the whole panel.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"

FASTQ_DIR="${FASTQ_DIR:-/mnt/pikachu/SLAM-Seq}"
OUT_BASE="${OUT_BASE:-/mnt/pikachu/SLAM-Seq-PE-results/prod_run_${STAMP}}"
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
STAR_INDEX="${STAR_INDEX:-/storage/autoindex_110_44/bulk_index}"
SNP_MASK="${SNP_MASK:-/mnt/pikachu/slam_blank_artifacts_20260201/mask/snps_from_vcf.bed.gz}"
THREADS="${THREADS:-16}"
OUT_SAM_ORDER="${OUT_SAM_ORDER:-SortedByCoordinate}"
SLAM_STRANDNESS="${SLAM_STRANDNESS:-Sense}"
SLAM_CB_FORMAT="${SLAM_CB_FORMAT:-star}"
SLAM_MIN_CALLABLE_LENGTH="${SLAM_MIN_CALLABLE_LENGTH:-30}"

# Fixed trims from the noSU 100K smoke on 2026-05-11.
SE_TRIM5P="${SE_TRIM5P:-8}"
SE_TRIM3P="${SE_TRIM3P:-12}"
PE_TRIM5P_M1="${PE_TRIM5P_M1:-8}"
PE_TRIM3P_M1="${PE_TRIM3P_M1:-13}"
PE_TRIM5P_M2="${PE_TRIM5P_M2:-19}"
PE_TRIM3P_M2="${PE_TRIM3P_M2:-14}"

GLOBUS_SRC_ENDPOINT="${GLOBUS_SRC_ENDPOINT:-}"
GLOBUS_DST_ENDPOINT="${GLOBUS_DST_ENDPOINT:-61fb8b9a-9b52-456e-928c-30c0fb0140bf}"
GLOBUS_DST_ROOT="${GLOBUS_DST_ROOT:-SLAM-seq-PE-results}"
GLOBUS_POLL_SECONDS="${GLOBUS_POLL_SECONDS:-30}"
SUBMIT_GLOBUS=1
WAIT_FOR_GLOBUS=0
CLEANUP_AFTER_GLOBUS=1

DRY_RUN=0
SE_ONLY=0
RESUME=0
LIMIT=0
PANEL_REGEX=""
declare -a SAMPLE_FILTERS=()

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Default action: run paired-end STAR-SLAM per sample from ${FASTQ_DIR}.

Selection:
  --sample NAME              Restrict to one sample. May be repeated.
  --samples CSV              Restrict to comma-separated sample names.
  --panel REGEX              Restrict to samples whose names match REGEX.
  --pilot                    Shortcut for ARID1A-no4su_S50 and ARID1A-6h-1_S43.
  --limit N                  Run at most N selected samples.
  --se-only                  Force R1-only mode even when R2 exists.

Paths:
  --fastq-dir DIR            Input FASTQ directory (default: ${FASTQ_DIR})
  --outdir DIR               Output root (default: ${OUT_BASE})
  --star-bin PATH            STAR binary (default: ${STAR_BIN})
  --genome-dir DIR           STAR genomeDir (default: ${STAR_INDEX})
  --snp-mask PATH            SLAM SNP mask BED.GZ (default: ${SNP_MASK})

Runtime:
  --threads N                STAR threads per sample (default: ${THREADS})
  --out-sam-order MODE       BAM order for --outSAMtype BAM MODE (default: ${OUT_SAM_ORDER})
  --resume                   Skip samples with an existing successful Log.final.out.
  --dry-run                  Write manifests and commands, but do not run STAR or Globus.

SLAM:
  --slam-cb-format FMT       Count-binomial output format: star|ezbakr (default: ${SLAM_CB_FORMAT})
  --slam-min-callable N      Minimum callable positions per read/pair; 0 disables (default: ${SLAM_MIN_CALLABLE_LENGTH})

Globus:
  --globus-src-endpoint ID   Source Globus collection ID for local outputs.
  --globus-dst-endpoint ID   Destination Globus collection ID
                             (default: ${GLOBUS_DST_ENDPOINT})
  --globus-dst-root PATH     Destination subdirectory/root
                             (default: ${GLOBUS_DST_ROOT})
  --no-globus                Do not submit Globus transfers; still write batch manifests.
  --wait-for-globus          Wait for each submitted transfer before starting next sample.
  --no-cleanup-after-globus  Keep local Y/noY BAM and FASTQ files after Globus succeeds.
  --globus-poll-seconds N    Poll interval when waiting (default: ${GLOBUS_POLL_SECONDS})

Examples:
  # Command/manifest dry-run for the two-sample pilot.
  $(basename "$0") --pilot --dry-run

  # Full ARID1A paired-end run with background Globus uploads.
  $(basename "$0") --panel '^ARID1A' --globus-src-endpoint <pikachu-src-uuid>
EOF
}

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"
}

die() {
  log "ERROR: $*" >&2
  exit 1
}

trap 'die "Command failed at line ${LINENO}: ${BASH_COMMAND}"' ERR

require_file() {
  [[ -f "$1" ]] || die "Missing file: $1"
}

require_dir() {
  [[ -d "$1" ]] || die "Missing directory: $1"
}

require_exec() {
  [[ -x "$1" ]] || die "Missing executable: $1"
}

quote_cmd() {
  printf '%q ' "$@"
}

sample_name_from_r1() {
  local fq="$1"
  local base
  base="$(basename "$fq")"
  base="${base%.fastq.gz}"
  base="${base%.fq.gz}"
  base="${base%.fastq}"
  base="${base%.fq}"
  base="${base%_R1_001}"
  base="${base%_R1}"
  base="${base%.R1}"
  echo "$base"
}

r2_from_r1() {
  local r1="$1"
  local direct=""
  local sample_id=""
  local dir
  local base
  local -a candidates=()
  if [[ "$r1" == *"_R1_"* ]]; then
    direct="${r1/_R1_/_R2_}"
  elif [[ "$r1" == *".R1."* ]]; then
    direct="${r1/.R1./.R2.}"
  elif [[ "$r1" == *"_R1."* ]]; then
    direct="${r1/_R1./_R2.}"
  elif [[ "$r1" == *".R1.fastq.gz" ]]; then
    direct="${r1%.R1.fastq.gz}.R2.fastq.gz"
  else
    die "Cannot infer R2 filename from R1: ${r1}"
  fi
  if [[ -f "$direct" ]]; then
    echo "$direct"
    return 0
  fi

  # Some delivered R2 files have a mismatched gene label but the same sample
  # index, e.g. *_S35_R1_001.fastq.gz paired with a different *_S35_R2_001.
  dir="$(dirname "$r1")"
  base="$(basename "$r1")"
  sample_id="$(sed -nE 's/^.*_(S[0-9]+)_R1_001\.fastq\.gz$/\1/p' <<< "$base")"
  if [[ -n "$sample_id" ]]; then
    mapfile -t candidates < <(find "$dir" -maxdepth 1 -type f -name "*_${sample_id}_R2_001.fastq.gz" | sort)
    if [[ "${#candidates[@]}" -eq 1 ]]; then
      echo "${candidates[0]}"
      return 0
    fi
  fi
  echo "$direct"
}

csv_to_samples() {
  tr ',' '\n' <<< "$1" | sed '/^$/d'
}

sample_filter_matches() {
  local sample="$1"
  local found=0
  local filter
  if ((${#SAMPLE_FILTERS[@]} > 0)); then
    for filter in "${SAMPLE_FILTERS[@]}"; do
      if [[ "$sample" == "$filter" ]]; then
        found=1
        break
      fi
    done
    ((found == 1)) || return 1
  fi
  if [[ -n "$PANEL_REGEX" && ! "$sample" =~ $PANEL_REGEX ]]; then
    return 1
  fi
  return 0
}

globus_enabled() {
  [[ "$SUBMIT_GLOBUS" == "1" && -n "$GLOBUS_SRC_ENDPOINT" && -n "$GLOBUS_DST_ENDPOINT" ]]
}

detect_globus_src_endpoint() {
  if [[ "$SUBMIT_GLOBUS" == "1" && -z "$GLOBUS_SRC_ENDPOINT" ]] && command -v globus >/dev/null 2>&1; then
    GLOBUS_SRC_ENDPOINT="$(globus endpoint local-id 2>/dev/null || true)"
  fi
}

globus_task_status() {
  local output
  output="$(globus task show "$1" 2>&1)" || return $?
  awk -F': *' '
    /^Status:/ {
      print $2
      found = 1
      exit
    }
    END {
      if (!found) {
        exit 1
      }
    }
  ' <<< "$output"
}

wait_for_globus_task() {
  local task_id="$1"
  local status=""
  [[ -n "$task_id" ]] || return 0
  while true; do
    status="$(globus_task_status "$task_id" || true)"
    case "$status" in
      SUCCEEDED)
        log "Globus transfer completed: ${task_id}"
        return 0
        ;;
      ACTIVE|INACTIVE|QUEUED|PENDING|"")
        sleep "$GLOBUS_POLL_SECONDS"
        ;;
      *)
        die "Globus transfer ${task_id} entered unexpected status: ${status}"
        ;;
    esac
  done
}

cleanup_sample_large_outputs() {
  local sample="$1"
  local sample_root="$2"
  local run_dir="${sample_root}/run"
  local cleanup_file="${sample_root}/TRANSFER_CLEANED.ok"

  [[ "$CLEANUP_AFTER_GLOBUS" == "1" ]] || return 0
  [[ -d "$run_dir" ]] || return 0
  rm -f "${run_dir}"/*_Y.bam "${run_dir}"/*_noY.bam
  rm -rf "${run_dir}/y_separated" "${run_dir}/tmp"
  {
    printf 'cleaned_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'sample=%s\n' "$sample"
  } > "$cleanup_file"
  log "Cleaned transferred BAM/FASTQ outputs for ${sample}"
}

reap_completed_transfers() {
  globus_enabled || return 0
  [[ "$CLEANUP_AFTER_GLOBUS" == "1" ]] || return 0
  [[ -d "${OUT_BASE}/samples" ]] || return 0

  local sample_root sample submitted_file cleaned_file task_file task_id status
  while IFS= read -r sample_root; do
    sample="$(basename "$sample_root")"
    submitted_file="${sample_root}/TRANSFER_SUBMITTED.ok"
    cleaned_file="${sample_root}/TRANSFER_CLEANED.ok"
    task_file="${sample_root}/TRANSFER_TASK.txt"
    [[ -f "$submitted_file" ]] || continue
    [[ ! -f "$cleaned_file" ]] || continue
    [[ -f "$task_file" ]] || die "Missing transfer task file for ${sample}: ${task_file}"
    task_id="$(awk -F': *' '/Task ID/ {print $2; exit}' "$task_file")"
    [[ -n "$task_id" ]] || die "Could not parse Globus task ID for ${sample}"
    status="$(globus_task_status "$task_id" || true)"
    case "$status" in
      SUCCEEDED)
        cleanup_sample_large_outputs "$sample" "$sample_root"
        ;;
      ACTIVE|INACTIVE|QUEUED|PENDING|"")
        if [[ -z "$status" ]]; then
          log "Globus status unavailable for ${sample} (${task_id}); will retry later"
        fi
        ;;
      *)
        die "Globus transfer for ${sample} entered unexpected status ${status} (task ${task_id})"
        ;;
    esac
  done < <(find "${OUT_BASE}/samples" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
}

write_top_level_manifests() {
  mkdir -p "${OUT_BASE}/manifests" "${OUT_BASE}/logs"
  cat > "${OUT_BASE}/manifests/run_config.env" <<EOF
FASTQ_DIR=${FASTQ_DIR}
OUT_BASE=${OUT_BASE}
STAR_BIN=${STAR_BIN}
STAR_INDEX=${STAR_INDEX}
SNP_MASK=${SNP_MASK}
THREADS=${THREADS}
OUT_SAM_ORDER=${OUT_SAM_ORDER}
SLAM_STRANDNESS=${SLAM_STRANDNESS}
SLAM_CB_FORMAT=${SLAM_CB_FORMAT}
SLAM_MIN_CALLABLE_LENGTH=${SLAM_MIN_CALLABLE_LENGTH}
SE_ONLY=${SE_ONLY}
SE_TRIM5P=${SE_TRIM5P}
SE_TRIM3P=${SE_TRIM3P}
PE_TRIM5P_M1=${PE_TRIM5P_M1}
PE_TRIM3P_M1=${PE_TRIM3P_M1}
PE_TRIM5P_M2=${PE_TRIM5P_M2}
PE_TRIM3P_M2=${PE_TRIM3P_M2}
GLOBUS_SRC_ENDPOINT=${GLOBUS_SRC_ENDPOINT}
GLOBUS_DST_ENDPOINT=${GLOBUS_DST_ENDPOINT}
GLOBUS_DST_ROOT=${GLOBUS_DST_ROOT}
CLEANUP_AFTER_GLOBUS=${CLEANUP_AFTER_GLOBUS}
EOF
  printf 'sample\tmode\tr1\tr2\tsample_root\trun_dir\tcounts_dir\tqc_dir\n' \
    > "${OUT_BASE}/manifests/samples.tsv"
  printf 'sample\tr1\tr2\treason\n' > "${OUT_BASE}/manifests/pairing_warnings.tsv"
  printf 'sample\tcommand\n' > "${OUT_BASE}/manifests/commands.tsv"
  : > "${OUT_BASE}/manifests/r1_files.txt"
  : > "${OUT_BASE}/manifests/r2_files.txt"
  printf 'sample\tbatch_file\ttask_id\tstatus\n' > "${OUT_BASE}/manifests/globus_tasks.tsv"
}

append_sample_manifest() {
  local sample="$1"
  local mode="$2"
  local r1="$3"
  local r2="$4"
  local sample_root="$5"
  local run_dir="$6"
  local counts_dir="$7"
  local qc_dir="$8"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$mode" "$r1" "$r2" "$sample_root" "$run_dir" "$counts_dir" "$qc_dir" \
    >> "${OUT_BASE}/manifests/samples.tsv"
  printf '%s\n' "$r1" >> "${OUT_BASE}/manifests/r1_files.txt"
  if [[ -n "$r2" ]]; then
    printf '%s\n' "$r2" >> "${OUT_BASE}/manifests/r2_files.txt"
  fi
}

build_star_command() {
  local sample="$1"
  local r1="$2"
  local r2="$3"
  local mode="$4"
  local run_dir="$5"
  local counts_dir="$6"
  local qc_dir="$7"
  local -n out_cmd="$8"

  out_cmd=(
    "$STAR_BIN"
    --runThreadN "$THREADS"
    --genomeDir "$STAR_INDEX"
    --readFilesIn "$r1"
  )
  if [[ "$mode" == "PE" ]]; then
    out_cmd+=("$r2")
  fi
  out_cmd+=(
    --readFilesCommand zcat
    --outFileNamePrefix "${run_dir}/star_"
    --outTmpDir "${run_dir}/tmp"
    --outSAMtype BAM "$OUT_SAM_ORDER"
    --outSAMattributes NH HI AS nM MD
    --outSAMprimaryFlag OneBestScore
    --emitNoYBAM yes
    --keepBAM no
    --emitYNoYFastq yes
    --emitYNoYFastqCompression gz
    --quantMode TranscriptVB
    --quantTranscriptomeSAMoutput BanSingleEnd
    --quantVBgcBias 1
    --quantVBgenes 1
    --quantVBgenesMode Tximport
    --quantVBLibType A
    --slamQuantMode 1
    --slamGrandSlamOut 1
    --slamCbOut 1
    --slamCbFormat "$SLAM_CB_FORMAT"
    --slamCbOutFile "${counts_dir}/${sample}.SlamQuant.cB.tsv"
    --slamStrandness "$SLAM_STRANDNESS"
    --slamOutFile "${counts_dir}/${sample}.SlamQuant.out"
    --slamQcReport "${qc_dir}/${sample}"
    --slamMinCallableLength "$SLAM_MIN_CALLABLE_LENGTH"
  )
  if [[ -n "$SNP_MASK" && "$SNP_MASK" != "-" ]]; then
    out_cmd+=(--slamSnpMaskIn "$SNP_MASK")
  fi
  if [[ "$mode" == "PE" ]]; then
    out_cmd+=(
      --slamCompatTrim5pMate1 "$PE_TRIM5P_M1"
      --slamCompatTrim3pMate1 "$PE_TRIM3P_M1"
      --slamCompatTrim5pMate2 "$PE_TRIM5P_M2"
      --slamCompatTrim3pMate2 "$PE_TRIM3P_M2"
    )
  else
    out_cmd+=(
      --slamCompatTrim5p "$SE_TRIM5P"
      --slamCompatTrim3p "$SE_TRIM3P"
    )
  fi
}

write_sample_command() {
  local sample="$1"
  local sample_root="$2"
  shift 2
  {
    printf '#!/usr/bin/env bash\nset -euo pipefail\n\n'
    quote_cmd "$@"
    printf '\n'
  } > "${sample_root}/RUN_COMMAND.sh"
  chmod +x "${sample_root}/RUN_COMMAND.sh"
  printf '%s\t' "$sample" >> "${OUT_BASE}/manifests/commands.tsv"
  quote_cmd "$@" >> "${OUT_BASE}/manifests/commands.tsv"
  printf '\n' >> "${OUT_BASE}/manifests/commands.tsv"
}

write_globus_batch() {
  local sample="$1"
  local sample_root="$2"
  local batch="$3"
  local dest_run_root
  local file
  local rel
  dest_run_root="${GLOBUS_DST_ROOT%/}/$(basename "$OUT_BASE")/samples/${sample}"
  : > "$batch"
  while IFS= read -r -d '' file; do
    rel="${file#${sample_root}/}"
    case "$rel" in
      tmp/*|run/tmp/*|*_STARtmp/*) continue ;;
    esac
    printf '%s\t%s/%s\n' "$file" "$dest_run_root" "$rel" >> "$batch"
  done < <(find "$sample_root" -type f -print0 | sort -z)
}

write_globus_helper() {
  local helper="$1"
  local batch="$2"
  cat > "$helper" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

SRC_ENDPOINT="${GLOBUS_SRC_ENDPOINT:?Set GLOBUS_SRC_ENDPOINT}"
DST_ENDPOINT="${GLOBUS_DST_ENDPOINT:?Set GLOBUS_DST_ENDPOINT}"
BATCH_FILE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/globus_batch.tsv"
LABEL="${GLOBUS_LABEL:-STAR-suite SLAM PE sample upload}"

globus transfer "${SRC_ENDPOINT}" "${DST_ENDPOINT}" --batch "${BATCH_FILE}" --label "${LABEL}"
EOF
  chmod +x "$helper"
}

submit_globus_batch() {
  local sample="$1"
  local batch="$2"
  local sample_root="$3"
  local output
  local task_id
  if ! globus_enabled; then
    printf '%s\t%s\t%s\t%s\n' "$sample" "$batch" "" "not-submitted" \
      >> "${OUT_BASE}/manifests/globus_tasks.tsv"
    return 0
  fi
  command -v globus >/dev/null 2>&1 || die "globus CLI not found on PATH"
  output="$(globus transfer "$GLOBUS_SRC_ENDPOINT" "$GLOBUS_DST_ENDPOINT" \
    --batch "$batch" --label "STAR-suite SLAM PE ${sample}" 2>&1)"
  printf '%s\n' "$output" > "${OUT_BASE}/logs/${sample}.globus_submit.log"
  printf '%s\n' "$output" > "${sample_root}/TRANSFER_TASK.txt"
  task_id="$(awk -F': *' '/Task ID/ {print $2; exit}' "${OUT_BASE}/logs/${sample}.globus_submit.log" || true)"
  [[ -n "$task_id" ]] || die "Could not parse Globus task ID for ${sample}; see ${OUT_BASE}/logs/${sample}.globus_submit.log"
  {
    printf 'submitted_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'task_id=%s\n' "$task_id"
    printf 'batch=%s\n' "$batch"
  } > "${sample_root}/TRANSFER_SUBMITTED.ok"
  printf '%s\t%s\t%s\t%s\n' "$sample" "$batch" "$task_id" "submitted" \
    >> "${OUT_BASE}/manifests/globus_tasks.tsv"
  log "Submitted Globus transfer for ${sample}: ${task_id:-unknown-task-id}"
  if [[ "$WAIT_FOR_GLOBUS" == "1" ]]; then
    wait_for_globus_task "$task_id"
    cleanup_sample_large_outputs "$sample" "$sample_root"
  fi
}

sample_done() {
  local run_dir="$1"
  local sample_root
  sample_root="$(cd -- "${run_dir}/.." && pwd)"
  [[ -s "${run_dir}/star_Log.final.out" ]] &&
    rg -q 'finished successfully' "${sample_root}/logs/star.stdout.log" 2>/dev/null
}

run_sample() {
  local sample="$1"
  local r1="$2"
  local r2="$3"
  local mode="$4"
  local sample_root="${OUT_BASE}/samples/${sample}"
  local run_dir="${sample_root}/run"
  local counts_dir="${sample_root}/counts"
  local qc_dir="${sample_root}/qc"
  local logs_dir="${sample_root}/logs"
  local manifest_dir="${sample_root}/manifests"
  local -a cmd
  local batch
  local helper

  mkdir -p "$run_dir" "$counts_dir" "$qc_dir" "$logs_dir" "$manifest_dir"
  reap_completed_transfers
  append_sample_manifest "$sample" "$mode" "$r1" "$r2" "$sample_root" "$run_dir" "$counts_dir" "$qc_dir"
  build_star_command "$sample" "$r1" "$r2" "$mode" "$run_dir" "$counts_dir" "$qc_dir" cmd
  write_sample_command "$sample" "$sample_root" "${cmd[@]}"

  if [[ "$DRY_RUN" == "1" ]]; then
    log "DRY-RUN ${sample} (${mode})"
    quote_cmd "${cmd[@]}"
    printf '\n'
    return 0
  fi

  if [[ "$RESUME" == "1" ]] && sample_done "$run_dir"; then
    log "Skipping ${sample}; existing successful STAR output found"
  else
    log "Running ${sample} (${mode})"
    /usr/bin/time -v -o "${logs_dir}/star.time.log" \
      "${cmd[@]}" \
      > "${logs_dir}/star.stdout.log" \
      2> "${logs_dir}/star.stderr.log"
  fi

  batch="${manifest_dir}/globus_batch.tsv"
  helper="${manifest_dir}/submit_globus.sh"
  write_globus_batch "$sample" "$sample_root" "$batch"
  write_globus_helper "$helper" "$batch"
  submit_globus_batch "$sample" "$batch" "$sample_root"
  reap_completed_transfers
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE_FILTERS+=("$2"); shift 2 ;;
    --samples)
      while IFS= read -r sample; do SAMPLE_FILTERS+=("$sample"); done < <(csv_to_samples "$2")
      shift 2
      ;;
    --panel) PANEL_REGEX="$2"; shift 2 ;;
    --pilot)
      SAMPLE_FILTERS+=("ARID1A-no4su_S50" "ARID1A-6h-1_S43")
      shift
      ;;
    --limit) LIMIT="$2"; shift 2 ;;
    --se-only) SE_ONLY=1; shift ;;
    --fastq-dir) FASTQ_DIR="$2"; shift 2 ;;
    --outdir|--out-base) OUT_BASE="$2"; shift 2 ;;
    --star-bin) STAR_BIN="$2"; shift 2 ;;
    --genome-dir) STAR_INDEX="$2"; shift 2 ;;
    --snp-mask) SNP_MASK="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --out-sam-order) OUT_SAM_ORDER="$2"; shift 2 ;;
    --resume) RESUME=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    --slam-cb-format) SLAM_CB_FORMAT="$2"; shift 2 ;;
    --slam-min-callable) SLAM_MIN_CALLABLE_LENGTH="$2"; shift 2 ;;
    --globus-src-endpoint) GLOBUS_SRC_ENDPOINT="$2"; shift 2 ;;
    --globus-dst-endpoint) GLOBUS_DST_ENDPOINT="$2"; shift 2 ;;
    --globus-dst-root) GLOBUS_DST_ROOT="$2"; shift 2 ;;
    --globus-poll-seconds) GLOBUS_POLL_SECONDS="$2"; shift 2 ;;
    --no-globus) SUBMIT_GLOBUS=0; shift ;;
    --wait-for-globus) WAIT_FOR_GLOBUS=1; shift ;;
    --no-cleanup-after-globus) CLEANUP_AFTER_GLOBUS=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -gt 0 ]] || die "--threads must be a positive integer"
[[ "$LIMIT" =~ ^[0-9]+$ ]] || die "--limit must be a non-negative integer"
[[ "$SLAM_MIN_CALLABLE_LENGTH" =~ ^[0-9]+$ ]] || die "--slam-min-callable must be a non-negative integer"
[[ "$SLAM_CB_FORMAT" == "star" || "$SLAM_CB_FORMAT" == "ezbakr" ]] || die "--slam-cb-format must be star or ezbakr"
[[ "$OUT_SAM_ORDER" == "Unsorted" || "$OUT_SAM_ORDER" == "SortedByCoordinate" ]] || die "--out-sam-order must be Unsorted or SortedByCoordinate"

require_dir "$FASTQ_DIR"
require_exec "$STAR_BIN"
require_dir "$STAR_INDEX"
if [[ -n "$SNP_MASK" && "$SNP_MASK" != "-" ]]; then
  require_file "$SNP_MASK"
fi

detect_globus_src_endpoint
if [[ "$SUBMIT_GLOBUS" == "1" && -z "$GLOBUS_SRC_ENDPOINT" ]]; then
  log "Globus destination is configured (${GLOBUS_DST_ENDPOINT}:${GLOBUS_DST_ROOT}), but no source endpoint was provided; transfers will not be submitted."
fi

write_top_level_manifests

log "STAR-SLAM PE production runner"
log "FASTQ dir: ${FASTQ_DIR}"
log "Output root: ${OUT_BASE}"
log "Globus destination: ${GLOBUS_DST_ENDPOINT}:${GLOBUS_DST_ROOT}"

mapfile -t R1_FILES < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name '*_R1_*.fastq.gz' | sort)
((${#R1_FILES[@]} > 0)) || die "No R1 FASTQs found under ${FASTQ_DIR}"

selected=0
for r1 in "${R1_FILES[@]}"; do
  sample="$(sample_name_from_r1 "$r1")"
  sample_filter_matches "$sample" || continue
  r2=""
  mode="SE"
  if [[ "$SE_ONLY" == "0" ]]; then
    r2="$(r2_from_r1 "$r1")"
    [[ -f "$r2" ]] || die "Missing R2 pair for ${r1}; expected ${r2}"
    r2_sample="$(basename "$r2")"
    r2_sample="${r2_sample%_R2_001.fastq.gz}"
    if [[ "$r2_sample" != "$sample" ]]; then
      printf '%s\t%s\t%s\t%s\n' "$sample" "$r1" "$r2" "paired_by_sample_index_name_mismatch" \
        >> "${OUT_BASE}/manifests/pairing_warnings.tsv"
    fi
    mode="PE"
  fi
  run_sample "$sample" "$r1" "$r2" "$mode"
  selected=$((selected + 1))
  if [[ "$LIMIT" -gt 0 && "$selected" -ge "$LIMIT" ]]; then
    break
  fi
done

((selected > 0)) || die "No samples matched selection"

log "Done. Selected samples: ${selected}"
if [[ "$DRY_RUN" == "0" ]]; then
  reap_completed_transfers
fi
log "Sample manifest: ${OUT_BASE}/manifests/samples.tsv"
log "Command manifest: ${OUT_BASE}/manifests/commands.tsv"
log "Globus task manifest: ${OUT_BASE}/manifests/globus_tasks.tsv"
