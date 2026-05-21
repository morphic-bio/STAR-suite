#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

RAW_DIR="${RAW_DIR:-/mnt/pikachu/JAX_scRNAseq02/raw}"
OUT_PARENT="${OUT_PARENT:-/mnt/pikachu/JAX_scRNAseq02_processed}"
BATCH_ROOT="${BATCH_ROOT:-${OUT_PARENT}/ocm_prod_batch_$(date -u +%Y%m%dT%H%M%SZ)}"
FIRST_RUN_ROOT="${FIRST_RUN_ROOT:-}"
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
GENOME_DIR="${GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
SOLO_CB_WHITELIST="${SOLO_CB_WHITELIST:-/storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_TRU.txt}"
THREADS="${THREADS:-16}"

REMOTE_HOST="${REMOTE_HOST:-10.159.4.53}"
REMOTE_ROOT="${REMOTE_ROOT:-/home/lhhung/jax_scrnaseq02_remote_downstream}"
DOWNSTREAM_NAME="${DOWNSTREAM_NAME:-downstream_genefull_velocyto_cellbender_remote}"

SRC_EP="${SRC_EP:-07446cad-33b8-11f0-8c0c-0afffb017b7d}"
DST_EP="${DST_EP:-61fb8b9a-9b52-456e-928c-30c0fb0140bf}"
DST_PREFIX="${DST_PREFIX:-/JAX_scRNAseq02_processed/large_files}"

SKIP_FIRST_POST="${SKIP_FIRST_POST:-0}"
SKIP_REMAINING_STAR="${SKIP_REMAINING_STAR:-0}"
KEEP_REMOTE="${KEEP_REMOTE:-0}"

mkdir -p "${BATCH_ROOT}/logs" "${BATCH_ROOT}/configs"
GLOBAL_STATUS="${BATCH_ROOT}/logs/batch_status.tsv"
POST_LOCK="${BATCH_ROOT}/logs/post_phase.lock"
: > "${GLOBAL_STATUS}"

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"
}

mark() {
  printf '%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$1" "$2" >> "${GLOBAL_STATUS}"
}

die() {
  log "ERROR: $*" >&2
  mark "FAILED" "$*"
  exit 1
}

require_file() {
  [[ -f "$1" ]] || die "missing file: $1"
}

require_dir() {
  [[ -d "$1" ]] || die "missing directory: $1"
}

LIBRARIES=(
  "25E32-L3"
  "25E32-L4"
  "25E34-L3"
  "25E34-L4"
  "25E35-L3"
  "25E35-L4"
)

declare -A STEMS=(
  ["25E32-L3"]="25E32-L3_GT25-03394_ACCTCGAGCT-ATCGAACACA_S44"
  ["25E32-L4"]="25E32-L4_GT25-03395_CGAAGTATAC-CTCCAAGTTC_S45"
  ["25E34-L3"]="25E34-L3_GT25-03396_GCACTGAGAA-TTCACGCATA_S40"
  ["25E34-L4"]="25E34-L4_GT25-03397_GCTACAAAGC-AGGGCACGTG_S41"
  ["25E35-L3"]="25E35-L3_GT25-03398_CGCTGAAATC-GCAGACACCT_S42"
  ["25E35-L4"]="25E35-L4_GT25-03399_GAGCAAGGGC-CCAAGTCAAT_S43"
)

declare -A SAMPLE_ROWS=(
  ["25E32-L3"]="GCM1-Day-4,OB1,iPSCs;GRHL1-Day-4,OB2,iPSCs;OVOL1-Day-4,OB3,iPSCs;WT-PrS-20pct-Day-4,OB4,iPSCs"
  ["25E32-L4"]="EPAS1-Day-4,OB1|OB2,iPSCs;WT-PrS-3pct-Day-4,OB3|OB4,iPSCs"
  ["25E34-L3"]="GCM1-Day-5,OB1,iPSCs;GRHL1-Day-5,OB2,iPSCs;OVOL1-Day-5,OB3,iPSCs;WT-PrS-20pct-Day-5,OB4,iPSCs"
  ["25E34-L4"]="ISL1-Day-5,OB1,iPSCs;EPAS1-Day-5,OB2,iPSCs;WT-PrS-3pct-Day-5,OB3,iPSCs;WT-ExM-Day-5,OB4,iPSCs"
  ["25E35-L3"]="GCM1-Day-6,OB1,iPSCs;GRHL1-Day-6,OB2,iPSCs;OVOL1-Day-6,OB3,iPSCs;WT-PrS-20pct-Day-6,OB4,iPSCs"
  ["25E35-L4"]="ISL1-Day-6,OB1,iPSCs;EPAS1-Day-6,OB2,iPSCs;WT-PrS-3pct-Day-6,OB3,iPSCs;WT-ExM-Day-6,OB4,iPSCs"
)

samples_for_library() {
  local lib="$1"
  tr ';' '\n' <<< "${SAMPLE_ROWS[${lib}]}" | cut -d, -f1
}

config_for_library() {
  printf '%s/configs/%s/config.csv' "${BATCH_ROOT}" "$1"
}

run_root_for_library() {
  printf '%s/%s' "${BATCH_ROOT}" "$1"
}

write_ocm_config() {
  local lib="$1"
  local cfg
  cfg="$(config_for_library "${lib}")"
  mkdir -p "$(dirname "${cfg}")"
  {
    printf '[gene-expression]\n'
    printf 'reference,/sc/service/pipelines/references/10x-rna/refdata-gex-GRCh38-2024-A\n'
    printf 'create-bam,true\n'
    printf 'include-introns,true\n\n'
    printf '[libraries]\n'
    printf 'fastq_id,fastqs,feature_types\n'
    printf '%s,%s,Gene Expression\n\n' "${STEMS[${lib}]%_S*}" "${RAW_DIR}"
    printf '[samples]\n'
    printf 'sample_id,ocm_barcode_ids,description\n'
    tr ';' '\n' <<< "${SAMPLE_ROWS[${lib}]}"
  } > "${cfg}"
}

source_fastq() {
  local lib="$1"
  local lane="$2"
  local read="$3"
  printf '%s/%s_%s_%s_001.fastq.gz' "${RAW_DIR}" "${STEMS[${lib}]}" "${lane}" "${read}"
}

stage_fastqs() {
  local lib="$1"
  local run_root="$2"
  local stage_dir="${run_root}/stage_fastqs"
  local manifest="${run_root}/downsample_manifest.tsv"
  mkdir -p "${stage_dir}"
  printf 'sample_id\tlane\tread_pairs\tr1\tr2\n' > "${manifest}"
  local lane src_r1 src_r2 dst_r1 dst_r2
  for lane in L007 L008; do
    src_r1="$(source_fastq "${lib}" "${lane}" R1)"
    src_r2="$(source_fastq "${lib}" "${lane}" R2)"
    require_file "${src_r1}"
    require_file "${src_r2}"
    dst_r1="${stage_dir}/$(basename "${src_r1}")"
    dst_r2="${stage_dir}/$(basename "${src_r2}")"
    ln -sfn "${src_r1}" "${dst_r1}"
    ln -sfn "${src_r2}" "${dst_r2}"
    printf '%s\t%s\tfull\t%s\t%s\n' "${lib}" "${lane}" "${dst_r1}" "${dst_r2}" >> "${manifest}"
  done
}

join_by_comma() {
  local IFS=,
  printf '%s' "$*"
}

render_star_script() {
  local lib="$1"
  local run_root="$2"
  local cfg="$3"
  local sample_dir="${run_root}/samples/${lib}"
  local run_dir="${sample_dir}/run"
  local tmp_dir="${sample_dir}/tmp"
  local log_dir="${run_root}/logs"
  local stage_dir="${run_root}/stage_fastqs"
  local run_script="${run_root}/RUN_STAR.sh"
  mkdir -p "${run_dir}" "${log_dir}"

  local r1_files=()
  local r2_files=()
  local lane
  for lane in L007 L008; do
    r1_files+=("${stage_dir}/$(basename "$(source_fastq "${lib}" "${lane}" R1)")")
    r2_files+=("${stage_dir}/$(basename "$(source_fastq "${lib}" "${lane}" R2)")")
  done

  local r1_csv r2_csv
  r1_csv="$(join_by_comma "${r1_files[@]}")"
  r2_csv="$(join_by_comma "${r2_files[@]}")"

  local cmd=(
    "${STAR_BIN}"
    --runThreadN "${THREADS}"
    --dynamicThreadInterface 1
    --genomeDir "${GENOME_DIR}"
    --readFilesIn "${r2_csv}" "${r1_csv}"
    --readFilesCommand zcat
    --outFileNamePrefix "${run_dir}/"
    --outTmpDir "${tmp_dir}"
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
    --soloInlineCBCorrection yes
    --soloUMIfiltering MultiGeneUMI_CR
    --soloUMIdedup 1MM_CR
    --soloMultiMappers Unique
    --soloCellFilter EmptyDrops_CR
    --soloCbUbRequireTogether no
    --soloStrand Forward
    --soloFeatures GeneFull Velocyto
    --soloCrGexFeature genefull
    --soloCrMultimapRescue yes
    --soloInlineHashMode no
    --ocmMultiEnable yes
    --ocmMultiConfig "${cfg}"
    --ocmMultiBarcodeMode flex
    --ocmMultiOutputCompat cellranger
  )

  {
    printf '#!/usr/bin/env bash\n'
    printf 'set -euo pipefail\n\n'
    printf 'cd %q\n' "${REPO_ROOT}"
    printf ': "${STAR_VELOCYTO_LOW_MEM:=1}"\n'
    printf ': "${STAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS:=8192}"\n'
    printf ': "${STAR_VELOCYTO_UMI_RESERVE_CAP:=32}"\n'
    printf ': "${STAR_SOLO_BINARY_SPOOL:=1}"\n'
    printf ': "${MALLOC_ARENA_MAX:=2}"\n'
    printf ': "${MALLOC_TRIM_THRESHOLD_:=131072}"\n'
    printf 'export STAR_VELOCYTO_LOW_MEM STAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS STAR_VELOCYTO_UMI_RESERVE_CAP STAR_SOLO_BINARY_SPOOL MALLOC_ARENA_MAX MALLOC_TRIM_THRESHOLD_\n'
    printf 'mkdir -p %q %q\n' "${run_dir}" "${log_dir}"
    printf 'printf "started_utc=%%s\\n" "$(date -u +%%Y-%%m-%%dT%%H:%%M:%%SZ)" > %q\n' "${run_root}/STAR_STARTED.txt"
    printf 'printf "STAR_VELOCYTO_LOW_MEM=%%s\\nSTAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS=%%s\\nSTAR_VELOCYTO_UMI_RESERVE_CAP=%%s\\nSTAR_SOLO_BINARY_SPOOL=%%s\\nMALLOC_ARENA_MAX=%%s\\nMALLOC_TRIM_THRESHOLD_=%%s\\n" "$STAR_VELOCYTO_LOW_MEM" "$STAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS" "$STAR_VELOCYTO_UMI_RESERVE_CAP" "$STAR_SOLO_BINARY_SPOOL" "$MALLOC_ARENA_MAX" "$MALLOC_TRIM_THRESHOLD_" > %q\n' "${log_dir}/star.env.txt"
    printf 'cmd=('
    local arg
    for arg in "${cmd[@]}"; do
      printf ' %q' "${arg}"
    done
    printf ' )\n'
    printf 'printf "STAR command:" | tee %q\n' "${log_dir}/star.command.txt"
    printf 'printf " %%q" "${cmd[@]}" | tee -a %q\n' "${log_dir}/star.command.txt"
    printf 'printf "\\n" | tee -a %q\n' "${log_dir}/star.command.txt"
    printf '"${cmd[@]}" 2>&1 | tee %q\n' "${log_dir}/star.log"
    printf 'rm -rf %q\n' "${tmp_dir}"
    printf 'printf "completed_utc=%%s\\n" "$(date -u +%%Y-%%m-%%dT%%H:%%M:%%SZ)" > %q\n' "${run_root}/STAR_COMPLETED.txt"
  } > "${run_script}"
  chmod +x "${run_script}"
}

star_finished() {
  local lib="$1"
  local run_root="$2"
  local run_dir="${run_root}/samples/${lib}/run"
  [[ -f "${run_root}/STAR_COMPLETED.txt" ]] && rg -q "ALL DONE!|finished successfully" "${run_dir}/Log.out"
}

run_star_library() {
  local lib="$1"
  local run_root="$2"
  local cfg="$3"
  local log_dir="${run_root}/logs"
  mkdir -p "${log_dir}"
  if star_finished "${lib}" "${run_root}"; then
    log "${lib}: STAR already complete; skipping mapping"
    mark "STAR_SKIP" "library=${lib};run_root=${run_root}"
    return 0
  fi

  log "${lib}: staging FASTQ symlinks and rendering STAR command"
  write_ocm_config "${lib}"
  stage_fastqs "${lib}" "${run_root}"
  render_star_script "${lib}" "${run_root}" "${cfg}"

  mark "STAR_START" "library=${lib};run_root=${run_root}"
  "${run_root}/RUN_STAR.sh"
  star_finished "${lib}" "${run_root}" || die "${lib}: STAR did not finish successfully"
  mark "STAR_PASS" "library=${lib};run_root=${run_root}"
}

init_qc_summary() {
  local run_root="$1"
  mkdir -p "${run_root}/logs"
  printf 'sample\tfinal_h5ad\tn_obs\tn_vars\tlayers\tmt_pct_threshold\tmt_pct_median\tmt_pct_mad\tmt_pct_n_mad\tmt_pct_flag\tcellbender_counts\n' > "${run_root}/logs/per_sample_h5ad_qc.tsv"
}

run_remote_downstream_for_sample() {
  local run_root="$1"
  local sample="$2"
  local sample_dir="${run_root}/samples/${sample}"
  local out_dir="${sample_dir}/${DOWNSTREAM_NAME}"
  local local_log="${run_root}/logs/${sample}.remote_downstream.log"
  local qc_summary="${run_root}/logs/per_sample_h5ad_qc.tsv"

  if [[ -f "${out_dir}/final_counts.h5ad" && -f "${out_dir}/adaptive_qc_threshold.json" ]]; then
    log "${sample}: downstream already complete; skipping"
    mark "REMOTE_DOWNSTREAM_SKIP" "sample=${sample};run_root=${run_root}"
    return 0
  fi

  require_dir "${sample_dir}/run/outs/filtered_feature_bc_matrix"
  require_dir "${sample_dir}/run/outs/raw_feature_bc_matrix"
  require_dir "${sample_dir}/run/outs/raw_velocyto_feature_bc_matrix"

  mark "REMOTE_DOWNSTREAM_START" "sample=${sample};run_root=${run_root}"
  local remote_args=()
  if [[ "${KEEP_REMOTE}" == "1" ]]; then
    remote_args+=(--keep-remote)
  fi
  "${REPO_ROOT}/scripts/run_remote_scrna_downstream_rsync.sh" \
    --sample-dir "${sample_dir}" \
    --remote-host "${REMOTE_HOST}" \
    --remote-root "${REMOTE_ROOT}" \
    --output-name "${DOWNSTREAM_NAME}" \
    --run-cellbender \
    --adaptive-filter \
    --n-mad 3 \
    --cellbender-gpu \
    --cellbender-cpu-cores 16 \
    --no-sync-images \
    --local-log "${local_log}" \
    "${remote_args[@]}"

  require_file "${out_dir}/final_counts.h5ad"
  require_file "${out_dir}/filtered_counts.h5ad"
  require_file "${out_dir}/unfiltered_counts.h5ad"
  require_file "${out_dir}/adaptive_qc_threshold.json"
  require_file "${out_dir}/summary.txt"

  python3 - <<'PY' "${sample}" "${out_dir}/final_counts.h5ad" "${out_dir}/adaptive_qc_threshold.json" "${out_dir}/cellbender/cellbender_counts.h5" "${qc_summary}"
import json
import sys
from pathlib import Path

import anndata as ad

sample, h5ad_s, qc_json_s, cellbender_s, summary_s = sys.argv[1:]
h5ad = Path(h5ad_s)
qc_json = Path(qc_json_s)
cellbender = Path(cellbender_s)
summary = Path(summary_s)

with qc_json.open("r", encoding="utf-8") as handle:
    qc = json.load(handle)
required_json = [
    "mt_pct_median",
    "mt_pct_mad",
    "mt_pct_n_mad",
    "mt_pct_floor",
    "mt_pct_raw_threshold",
    "mt_pct_threshold",
    "filter_cells_mt_adaptive",
    "singlet_filtered_cells_mt_adaptive",
]
missing_json = [key for key in required_json if key not in qc]
if missing_json:
    raise SystemExit(f"{sample}: adaptive QC JSON missing keys: {missing_json}")

adata = ad.read_h5ad(h5ad, backed="r")
try:
    obs_cols = set(adata.obs.columns)
    layers = list(adata.layers.keys())
    required_obs = [
        "mt_pct",
        "n_genes",
        "filter",
        "singlet_filtered",
        "filter_strict_mt5",
        "singlet_filtered_strict_mt5",
    ]
    missing_obs = [col for col in required_obs if col not in obs_cols]
    if missing_obs:
        raise SystemExit(f"{sample}: final h5ad missing obs columns: {missing_obs}")
    if cellbender.exists() and "denoised" not in layers:
        raise SystemExit(f"{sample}: CellBender output exists but denoised layer is missing")
    line = "\t".join(
        [
            sample,
            str(h5ad),
            str(adata.n_obs),
            str(adata.n_vars),
            ",".join(layers),
            str(qc["mt_pct_threshold"]),
            str(qc["mt_pct_median"]),
            str(qc["mt_pct_mad"]),
            str(qc["mt_pct_n_mad"]),
            str(qc.get("mt_pct_flag", False)),
            str(cellbender.exists()),
        ]
    )
finally:
    adata.file.close()

with summary.open("a", encoding="utf-8") as handle:
    handle.write(line + "\n")
PY

  mark "REMOTE_DOWNSTREAM_PASS" "sample=${sample};final=${out_dir}/final_counts.h5ad"
}

add_large_file() {
  local run_root="$1"
  local path="$2"
  local large_batch="$3"
  local large_inventory="$4"
  [[ -f "${path}" ]] || return 0
  local dst_root="${DST_PREFIX}/$(basename "${run_root}")"
  local rel="${path#${run_root}/}"
  local dest="${dst_root}/${rel}"
  local size
  size="$(stat -c '%s' "${path}")"
  printf '%s %s\n' "${path}" "${dest}" >> "${large_batch}"
  printf '%s\t%s\t%s\n' "${path}" "${size}" "${dest}" >> "${large_inventory}"
}

build_large_file_inventory() {
  local lib="$1"
  local run_root="$2"
  local run_dir="${run_root}/samples/${lib}/run"
  local large_batch="${run_root}/logs/globus_large_files.batch"
  local large_inventory="${run_root}/logs/globus_large_files_inventory.tsv"
  : > "${large_batch}"
  printf 'source_path\tsize_bytes\tdestination_path\n' > "${large_inventory}"

  while IFS= read -r -d '' file; do
    add_large_file "${run_root}" "${file}" "${large_batch}" "${large_inventory}"
  done < <(find "${run_root}" -type f -name '*.bam' -print0 | sort -z)

  if [[ -d "${run_dir}/y_separated" ]]; then
    while IFS= read -r -d '' file; do
      add_large_file "${run_root}" "${file}" "${large_batch}" "${large_inventory}"
    done < <(find "${run_dir}/y_separated" -type f -name '*.fastq.gz' -print0 | sort -z)
  fi

  [[ -s "${large_batch}" ]]
}

transfer_large_files_and_cleanup() {
  local lib="$1"
  local run_root="$2"
  local run_dir="${run_root}/samples/${lib}/run"
  local large_batch="${run_root}/logs/globus_large_files.batch"
  local large_inventory="${run_root}/logs/globus_large_files_inventory.tsv"
  local deleted_inventory="${run_root}/logs/deleted_generated_large_files.tsv"
  local dst_root="${DST_PREFIX}/$(basename "${run_root}")"

  if ! build_large_file_inventory "${lib}" "${run_root}"; then
    log "${lib}: no generated large files found for Globus transfer"
    mark "GLOBUS_LARGE_SKIP" "library=${lib};run_root=${run_root}"
    return 0
  fi

  mark "GLOBUS_LARGE_START" "library=${lib};batch=${large_batch}"
  globus mkdir "${DST_EP}:${dst_root}" >/dev/null 2>&1 || true
  local task_json task_id task_status
  task_json="$(globus transfer "${SRC_EP}" "${DST_EP}" \
    --batch "${large_batch}" \
    --sync-level checksum \
    --label "JAX scRNAseq02 OCM generated large files $(basename "${run_root}")" \
    --notify off \
    --format json)"
  printf '%s\n' "${task_json}" > "${run_root}/logs/globus_large_files_task.json"
  task_id="$(python3 - <<'PY' "${run_root}/logs/globus_large_files_task.json"
import json, sys
with open(sys.argv[1], "r", encoding="utf-8") as handle:
    print(json.load(handle)["task_id"])
PY
)"
  echo "${task_id}" > "${run_root}/logs/globus_large_files_task_id.txt"
  globus task wait --polling-interval 60 "${task_id}"
  globus task show "${task_id}" --format json > "${run_root}/logs/globus_large_files_task_final.json"
  task_status="$(python3 - <<'PY' "${run_root}/logs/globus_large_files_task_final.json"
import json, sys
with open(sys.argv[1], "r", encoding="utf-8") as handle:
    print(json.load(handle).get("status"))
PY
)"
  [[ "${task_status}" == "SUCCEEDED" ]] || die "${lib}: large-file Globus task status ${task_status}"

  printf 'deleted_utc\tsource_path\tsize_bytes\tdestination_path\n' > "${deleted_inventory}"
  tail -n +2 "${large_inventory}" | while IFS=$'\t' read -r source_path size_bytes destination_path; do
    if [[ -f "${source_path}" ]]; then
      rm -f "${source_path}"
      printf '%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${source_path}" "${size_bytes}" "${destination_path}" >> "${deleted_inventory}"
    fi
  done
  find "${run_dir}/y_separated" -type d -empty -delete 2>/dev/null || true
  mark "GLOBUS_LARGE_PASS" "library=${lib};task_id=${task_id};deleted_inventory=${deleted_inventory}"
}

post_phase() {
  local lib="$1"
  local run_root="$2"
  (
    flock 9
    log "${lib}: starting serialized post-STAR phase"
    mark "POST_START" "library=${lib};run_root=${run_root}"
    init_qc_summary "${run_root}"
    local samples=()
    mapfile -t samples < <(samples_for_library "${lib}")
    local sample
    for sample in "${samples[@]}"; do
      [[ -n "${sample}" ]] || continue
      run_remote_downstream_for_sample "${run_root}" "${sample}"
    done
    transfer_large_files_and_cleanup "${lib}" "${run_root}"
    mark "POST_PASS" "library=${lib};run_root=${run_root}"
    log "${lib}: post-STAR phase complete"
  ) 9>"${POST_LOCK}"
}

preflight() {
  require_dir "${RAW_DIR}"
  require_file "${STAR_BIN}"
  [[ -x "${STAR_BIN}" ]] || die "STAR is not executable: ${STAR_BIN}"
  require_dir "${GENOME_DIR}"
  require_file "${SOLO_CB_WHITELIST}"
  command -v globus >/dev/null 2>&1 || die "globus CLI not found"
  command -v ssh >/dev/null 2>&1 || die "ssh not found"
  command -v rsync >/dev/null 2>&1 || die "rsync not found"
  command -v docker >/dev/null 2>&1 || die "docker not found"
  command -v python3 >/dev/null 2>&1 || die "python3 not found"

  if ! ssh -n -o BatchMode=yes -o ConnectTimeout=10 -o StrictHostKeyChecking=accept-new "${REMOTE_HOST}" \
    "hostname; docker image inspect biodepot/scrna-matrices:latest >/dev/null; docker image inspect biodepot/cellbender:0.3.2 >/dev/null; nvidia-smi --query-gpu=index,name,memory.free,memory.total --format=csv,noheader" \
    > "${BATCH_ROOT}/logs/remote_gpu_preflight.txt" 2>&1; then
    die "remote GPU preflight failed; see ${BATCH_ROOT}/logs/remote_gpu_preflight.txt"
  fi
  if ! globus endpoint local-id > "${BATCH_ROOT}/logs/globus_local_endpoint.txt" 2>&1; then
    die "Globus local endpoint check failed; see ${BATCH_ROOT}/logs/globus_local_endpoint.txt"
  fi
  mark "PREFLIGHT_PASS" "batch_root=${BATCH_ROOT};remote_host=${REMOTE_HOST}"
}

main() {
  log "Starting JAX scRNAseq02 OCM production batch"
  mark "START" "batch_root=${BATCH_ROOT}"
  preflight

  local post_pids=()

  if [[ "${SKIP_FIRST_POST}" != "1" && -n "${FIRST_RUN_ROOT}" ]]; then
    require_dir "${FIRST_RUN_ROOT}"
    if ! star_finished "25E32-L3" "${FIRST_RUN_ROOT}"; then
      die "first run root is not a completed STAR run: ${FIRST_RUN_ROOT}"
    fi
    log "25E32-L3: launching post-STAR phase from existing completed run"
    post_phase "25E32-L3" "${FIRST_RUN_ROOT}" > "${FIRST_RUN_ROOT}/logs/post_phase.log" 2>&1 &
    post_pids+=("$!")
    mark "POST_QUEUED" "library=25E32-L3;pid=${post_pids[-1]};run_root=${FIRST_RUN_ROOT}"
  fi

  if [[ "${SKIP_REMAINING_STAR}" != "1" ]]; then
    local lib run_root cfg
    for lib in "${LIBRARIES[@]}"; do
      if [[ "${SKIP_FIRST_POST}" != "1" && -n "${FIRST_RUN_ROOT}" && "${lib}" == "25E32-L3" ]]; then
        log "Skipping fresh ${lib}; using FIRST_RUN_ROOT=${FIRST_RUN_ROOT}"
        continue
      fi
      run_root="$(run_root_for_library "${lib}")"
      cfg="$(config_for_library "${lib}")"
      mkdir -p "${run_root}/logs"
      write_ocm_config "${lib}"
      run_star_library "${lib}" "${run_root}" "${cfg}"
      log "${lib}: launching post-STAR phase"
      post_phase "${lib}" "${run_root}" > "${run_root}/logs/post_phase.log" 2>&1 &
      post_pids+=("$!")
      mark "POST_QUEUED" "library=${lib};pid=${post_pids[-1]};run_root=${run_root}"
    done
  fi

  local pid
  for pid in "${post_pids[@]}"; do
    wait "${pid}"
  done

  mark "COMPLETE" "batch_root=${BATCH_ROOT}"
  log "PASS: JAX scRNAseq02 OCM production batch complete"
}

main "$@"
