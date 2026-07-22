#!/usr/bin/env bash

# Full-sample FLEX no-align seed + inline-hash snapshot replay harness.
#
# FLEX_SNAPSHOT_COMPARE:
#   count (default) — replay sets STAR_SOLO_FLEX_SNAPSHOT_STOP_AFTER_COUNT=1 (no raw/filtered
#     outputs on replay). After both runs, assert snapshot is non-empty (STAR also refuses
#     0-entry seed unless STAR_SOLO_FLEX_HASH_SNAPSHOT_ALLOW_EMPTY=1) and seed/replay hash
#     entry counts in Log.out match.
#   full — replay runs full Solo output; diff seed vs replay per_sample (includes Gene/raw and
#     Gene/filtered under each BC*) + Solo.out (fails on any difference). Omit STOP_AFTER_COUNT on replay.
#
# STAR_SOLO_FLEX_HASH_SNAPSHOT_ALLOW_EMPTY=1 — allow writing/loading a 0-entry flex snapshot
#   (default is refuse; avoids silent false-positive harness passes).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${SCRIPT_DIR}/.."
STAR_BIN="${STAR_BIN:-${REPO_DIR}/core/legacy/source/STAR}"

MODE="${1:-both}"  # seed | replay | both

OUT_ROOT="${OUT_ROOT:-/storage/flex_snapshot_noalign_$(date -u +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-32}"
TMP_ROOT="${TMP_ROOT:-/storage/flex_snapshot_noalign_tmp}"
SNAPSHOT_PATH="${SNAPSHOT_PATH:-${OUT_ROOT}/flex_snapshot.bin}"
FLEX_SNAPSHOT_COMPARE="${FLEX_SNAPSHOT_COMPARE:-count}"

GENOME_DIR="${GENOME_DIR:-/storage/flex_filtered_reference_2024/star_index}"
CB_WHITELIST="${CB_WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}"
SAMPLE_WHITELIST="${SAMPLE_WHITELIST:-/mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv}"
PROBE_LIST="${PROBE_LIST:-/storage/flex_filtered_reference_2024/star_index/flex_probe_artifacts/probe_list.txt}"
SAMPLE_PROBES="${SAMPLE_PROBES:-/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}"
HASH_CACHE="${HASH_CACHE:-/storage/downsampled_100K/SC2300771/results/flex_h01_2024_20260320_081246/h01_cache.bin}"

READS_R2="${READS_R2:-$(find /storage/JAX_sequences -maxdepth 1 -type f -name 'SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L*_R2_001.fastq.gz' | sort | paste -sd, -)}"
READS_R1="${READS_R1:-$(find /storage/JAX_sequences -maxdepth 1 -type f -name 'SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L*_R1_001.fastq.gz' | sort | paste -sd, -)}"

if [[ ! -x "${STAR_BIN}" ]]; then
  echo "ERROR: STAR binary not found: ${STAR_BIN}" >&2
  exit 1
fi

for path in "${GENOME_DIR}" "${CB_WHITELIST}" "${SAMPLE_WHITELIST}" "${PROBE_LIST}" "${SAMPLE_PROBES}" "${HASH_CACHE}"; do
  [[ -e "${path}" ]] || { echo "ERROR: Missing required input: ${path}" >&2; exit 1; }
done

[[ -n "${READS_R1}" && -n "${READS_R2}" ]] || {
  echo "ERROR: Failed to resolve full-sample FASTQ lists under /storage/JAX_sequences" >&2
  exit 1
}

mkdir -p "${OUT_ROOT}"

run_mode() {
  local label="$1"
  local out_dir="${OUT_ROOT}/${label}"
  local tmp_dir="${TMP_ROOT}/${label}"

  rm -rf "${out_dir}" "${tmp_dir}"
  mkdir -p "${out_dir}"
  mkdir -p "$(dirname "${tmp_dir}")"

  local -a cmd=(
    "${STAR_BIN}"
    --runThreadN "${THREADS}"
    --outTmpDir "${tmp_dir}"
    --genomeDir "${GENOME_DIR}"
    --soloType CB_UMI_Simple
    --soloCBstart 1 --soloUMIstart 17 --soloCBlen 16 --soloUMIlen 12
    --soloBarcodeReadLength 0
    --soloCBwhitelist "${CB_WHITELIST}"
    --flex yes
    --soloFlexExpectedCellsPerTag 3000
    --soloSampleWhitelist "${SAMPLE_WHITELIST}"
    --soloProbeList "${PROBE_LIST}"
    --soloSampleProbes "${SAMPLE_PROBES}"
    --soloSampleProbeOffset 68
    --soloFlexAllowedTags "${SAMPLE_WHITELIST}"
    --soloFlexOutputPrefix "${out_dir}/per_sample"
    --limitIObufferSize 50000000 50000000
    --outSJtype None
    --outSAMtype None
    --soloFeatures Gene
    --soloCellFilter None
    --soloMultiMappers Rescue
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
    --soloUMIfiltering MultiGeneUMI_CR
    --soloUMIdedup 1MM_CR
    --soloStrand Unstranded
    --clipAdapterType CellRanger4
    --alignEndsType Local
    --chimSegmentMin 1000000
    --soloKeysCompat cr
    --outSAMattributes None
    --soloSampleSearchNearby no
    --soloHashScreenFile "${HASH_CACHE}"
    --flexPipeline yes
    --flexPipelineNTriage 0
    --flexPipelineNSolo 0
    --flexNoAlign 1
    --readFilesIn "${READS_R2}" "${READS_R1}"
    --outFileNamePrefix "${out_dir}/"
  )

  echo "=== ${label} ==="
  echo "Output: ${out_dir}"
  echo "Tmp:    ${tmp_dir}"
  echo "Snap:   ${SNAPSHOT_PATH}"
  echo "Compare mode (replay): ${FLEX_SNAPSHOT_COMPARE}"

  if [[ "${label}" == "seed" ]]; then
    STAR_SOLO_FLEX_HASH_SNAPSHOT_OUT="${SNAPSHOT_PATH}" \
    /usr/bin/time -v "${cmd[@]}" \
      >"${out_dir}/stdout.log" 2>"${out_dir}/stderr.log"
  else
    if [[ "${FLEX_SNAPSHOT_COMPARE}" == "full" ]]; then
      STAR_SOLO_FLEX_HASH_SNAPSHOT_IN="${SNAPSHOT_PATH}" \
      STAR_SOLO_FLEX_HASH_SNAPSHOT_REPLAY_SKIP_READS=1 \
      /usr/bin/time -v "${cmd[@]}" \
        >"${out_dir}/stdout.log" 2>"${out_dir}/stderr.log"
    else
      STAR_SOLO_FLEX_HASH_SNAPSHOT_IN="${SNAPSHOT_PATH}" \
      STAR_SOLO_FLEX_HASH_SNAPSHOT_REPLAY_SKIP_READS=1 \
      STAR_SOLO_FLEX_SNAPSHOT_STOP_AFTER_COUNT=1 \
      /usr/bin/time -v "${cmd[@]}" \
        >"${out_dir}/stdout.log" 2>"${out_dir}/stderr.log"
    fi
  fi
}

# FLEX_SNAPSHOT_COMPARE=count: replay stops before processRecords output trees; only compare
# snapshot + log-reported hash sizes (STAR refuses 0-entry flex snapshots unless ALLOW_EMPTY=1).
compare_count_surface() {
  local seed_log="${OUT_ROOT}/seed/Log.out"
  local replay_log="${OUT_ROOT}/replay/Log.out"
  local report="${OUT_ROOT}/compare.log"

  [[ -f "${SNAPSHOT_PATH}" ]] || {
    echo "ERROR: Snapshot missing: ${SNAPSHOT_PATH}" >&2
    return 1
  }
  local sz
  sz=$(stat -c%s "${SNAPSHOT_PATH}" 2>/dev/null || stat -f%z "${SNAPSHOT_PATH}")
  if [[ "${sz}" -lt 64 ]]; then
    echo "ERROR: Snapshot file implausibly small (${sz} bytes): ${SNAPSHOT_PATH}" >&2
    return 1
  fi

  local seed_n replay_n
  seed_n=$(grep 'STAR_SOLO_FLEX_HASH_SNAPSHOT_OUT wrote' "${seed_log}" | sed -n 's/.*wrote \([0-9][0-9]*\) flex hash entries.*/\1/p' | tail -1)
  replay_n=$(grep 'STAR_SOLO_FLEX_HASH_SNAPSHOT_IN replay: loaded' "${replay_log}" | sed -n 's/.*loaded \([0-9][0-9]*\) hash entries.*/\1/p' | tail -1)

  {
    echo "FLEX_SNAPSHOT_COMPARE=count (replay stopped after countCBgeneUMI; no output-tree parity here)"
    echo "snapshot_bytes=${sz}"
    echo "seed_flex_hash_entries_log=${seed_n:-}"
    echo "replay_flex_hash_entries_log=${replay_n:-}"
  } | tee "${report}"

  [[ -n "${seed_n}" ]] || { echo "ERROR: Could not parse seed flex snapshot write count from ${seed_log}" >&2; return 1; }
  [[ "${seed_n}" -gt 0 ]] || { echo "ERROR: Seed log reports 0 flex hash entries" >&2; return 1; }
  [[ -n "${replay_n}" ]] || { echo "ERROR: Could not parse replay flex snapshot load count from ${replay_log}" >&2; return 1; }
  [[ "${seed_n}" == "${replay_n}" ]] || {
    echo "ERROR: Flex hash entry count mismatch seed=${seed_n} replay=${replay_n}" >&2
    return 1
  }

  grep -q 'inline-hash snapshot replay: stopping after countCBgeneUMI' "${replay_log}" || {
    echo "ERROR: Expected replay early-stop log line missing in ${replay_log}" >&2
    return 1
  }

  echo "Count-surface check OK (seed/replay flex hash entries match: ${seed_n})"
}

# FLEX_SNAPSHOT_COMPARE=full: replay runs full output; diff trees (strict).
# Filtered per-sample MEX lives under per_sample/<BC>/Gene/filtered/ (soloFlexOutputPrefix + sample dirs).
compare_outputs_full() {
  local seed_dir="${OUT_ROOT}/seed"
  local replay_dir="${OUT_ROOT}/replay"
  local report="${OUT_ROOT}/compare.log"

  {
    echo "FLEX_SNAPSHOT_COMPARE=full — diff per_sample (raw+filtered) and Solo.out"
    echo "Comparing per_sample outputs"
    diff -qr "${seed_dir}/per_sample" "${replay_dir}/per_sample"
    echo
    echo "Comparing Solo.out"
    diff -qr "${seed_dir}/Solo.out" "${replay_dir}/Solo.out"
  } | tee "${report}"

  echo "Full output diff OK"
}

case "${MODE}" in
  seed)
    run_mode seed
    ;;
  replay)
    [[ -f "${SNAPSHOT_PATH}" ]] || { echo "ERROR: Snapshot missing for replay: ${SNAPSHOT_PATH}" >&2; exit 1; }
    run_mode replay
    ;;
  both)
    run_mode seed
    run_mode replay
    if [[ "${FLEX_SNAPSHOT_COMPARE}" == "full" ]]; then
      compare_outputs_full
    else
      compare_count_surface
    fi
    ;;
  *)
    echo "Usage: $(basename "$0") [seed|replay|both]" >&2
    exit 2
    ;;
esac

echo "Outputs: ${OUT_ROOT}"
