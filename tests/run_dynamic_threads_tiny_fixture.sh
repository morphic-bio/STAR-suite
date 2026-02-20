#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAR_BIN="${STAR_BIN:-${ROOT_DIR}/core/legacy/source/STAR}"
PARSER="${ROOT_DIR}/tests/dynamic_threads/mock_consumer_report.py"

TINY_REF_DIR="${TINY_REF_DIR:-/home/lhhung/cellranger-9.0.1/external/cellranger_tiny_ref}"
TINY_FASTQ_DIR="${TINY_FASTQ_DIR:-/home/lhhung/cellranger-9.0.1/external/cellranger_tiny_fastq}"

RUN_THREADS="${RUN_THREADS:-4}"
MAP_PERMITS="${MAP_PERMITS:-2}"
READ_MAP_NUMBER="${READ_MAP_NUMBER:-120000}"
OUT_SAMTYPE="${OUT_SAMTYPE:-None}"
VARIABLE_THREADS="${VARIABLE_THREADS:-0}"
VARIABLE_THREADS_RETUNE_EVERY_ACQUIRES="${VARIABLE_THREADS_RETUNE_EVERY_ACQUIRES:-0}"
VARIABLE_THREADS_PERMIT_SEQUENCE="${VARIABLE_THREADS_PERMIT_SEQUENCE:-}"
MIN_RETUNES="${MIN_RETUNES:-0}"
CHECK_LOG_FINAL_PARITY="${CHECK_LOG_FINAL_PARITY:-1}"
CHECK_BAM_PARITY="${CHECK_BAM_PARITY:-auto}"
REQUIRE_BAM_BODY_MD5_MATCH="${REQUIRE_BAM_BODY_MD5_MATCH:-1}"
ASSERT_RETUNE_TRACE_CYCLE="${ASSERT_RETUNE_TRACE_CYCLE:-1}"
KEEP_OUTPUTS="${KEEP_OUTPUTS:-1}"

OUT_BASE="${OUT_BASE:-/tmp/dynamic_threads_tiny_$(date +%Y%m%d_%H%M%S)_$$}"
OUT_OFF="${OUT_BASE}/off"
OUT_ON="${OUT_BASE}/on"
REPORT_JSON="${OUT_BASE}/dynamic_thread_report.json"
REPORT_TXT="${OUT_BASE}/dynamic_thread_report.txt"
LOG_FINAL_DIFF="${OUT_BASE}/log_final.diff"
BAM_PARITY_SUMMARY="${OUT_BASE}/bam_parity_summary.txt"

GENOME_DIR="${TINY_REF_DIR}/star"
FASTA_FILE="${TINY_REF_DIR}/fasta/genome.fa"
GTF_FILE_GZ="${TINY_REF_DIR}/genes/genes.gtf.gz"
GTF_FILE_UNZIPPED="${OUT_BASE}/genes.gtf"
R2_L001="${TINY_FASTQ_DIR}/tinygex_S1_L001_R2_001.fastq.gz"
R2_L002="${TINY_FASTQ_DIR}/tinygex_S1_L002_R2_001.fastq.gz"

extract_metric() {
    local metric_name="$1"
    local log_final="$2"
    awk -F'|' -v metric="${metric_name}" '
        $1 ~ metric {
            gsub(/^[ \t]+|[ \t]+$/, "", $2);
            gsub(/,/, "", $2);
            print $2;
            exit;
        }
    ' "${log_final}"
}

canonicalize_log_final_metrics() {
    local log_final="$1"
    awk -F'|' '
        NF >= 2 {
            key = $1
            value = $2
            gsub(/^[ \t]+|[ \t]+$/, "", key)
            gsub(/^[ \t]+|[ \t]+$/, "", value)
            gsub(/,/, "", value)
            if (key == "Started job on" || key == "Started mapping on" || key == "Finished on" || key == "Mapping speed, Million of reads per hour") {
                next
            }
            if (key != "") {
                print key "=" value
            }
        }
    ' "${log_final}" | LC_ALL=C sort
}

resolve_bam_path() {
    local out_dir="$1"
    if [[ "${OUT_SAMTYPE}" == "BAM SortedByCoordinate" ]]; then
        echo "${out_dir}/Aligned.sortedByCoord.out.bam"
    elif [[ "${OUT_SAMTYPE}" == "BAM Unsorted" ]]; then
        echo "${out_dir}/Aligned.out.bam"
    else
        echo ""
    fi
}

require_file() {
    local path="$1"
    if [[ ! -f "${path}" ]]; then
        echo "ERROR: missing required file: ${path}" >&2
        exit 1
    fi
}

run_star_mode() {
    local label="$1"
    local out_dir="$2"
    local dynamic_interface="$3"
    local dynamic_telemetry="$4"
    local dynamic_permits="$5"
    local variable_threads="$6"
    local retune_every_acquires="$7"
    local permit_sequence="$8"

    rm -rf "${out_dir}"
    mkdir -p "${out_dir}"

    echo "=== Running mode: ${label} ==="
    local extra_args=()
    local out_samtype_spaced="${OUT_SAMTYPE//,/ }"
    # shellcheck disable=SC2206
    local out_samtype_array=(${out_samtype_spaced})

    if [[ "${variable_threads}" == "1" && "${retune_every_acquires}" -gt 0 && -n "${permit_sequence}" ]]; then
        local permit_sequence_spaced="${permit_sequence//,/ }"
        # shellcheck disable=SC2206
        local permit_sequence_array=(${permit_sequence_spaced})
        extra_args+=(
            --variableThreadsRetuneEveryAcquires "${retune_every_acquires}"
            --variableThreadsPermitSequence
            "${permit_sequence_array[@]}"
        )
    fi

    # shellcheck disable=SC2086
    "${STAR_BIN}" ${STAR_EXTRA_ARGS:-} \
        --runMode alignReads \
        --runThreadN "${RUN_THREADS}" \
        --genomeDir "${GENOME_DIR}" \
        --readFilesIn "${R2_L001},${R2_L002}" \
        --readFilesCommand zcat \
        --readMapNumber "${READ_MAP_NUMBER}" \
        --outSAMtype "${out_samtype_array[@]}" \
        --outFileNamePrefix "${out_dir}/" \
        --outTmpDir "${out_dir}/_STARtmp" \
        --dynamicThreadInterface "${dynamic_interface}" \
        --dynamicThreadConstMapPermits "${dynamic_permits}" \
        --dynamicThreadTelemetry "${dynamic_telemetry}" \
        --variableThreads "${variable_threads}" \
        "${extra_args[@]}"
}

if [[ ! -x "${STAR_BIN}" ]]; then
    echo "ERROR: STAR binary not found at ${STAR_BIN}. Build with: make core" >&2
    exit 1
fi

require_file "${FASTA_FILE}"
require_file "${GTF_FILE_GZ}"
require_file "${R2_L001}"
require_file "${R2_L002}"
require_file "${PARSER}"

echo "=== Dynamic Thread Tiny Fixture Smoke ==="
echo "STAR_BIN=${STAR_BIN}"
echo "TINY_REF_DIR=${TINY_REF_DIR}"
echo "TINY_FASTQ_DIR=${TINY_FASTQ_DIR}"
echo "OUT_BASE=${OUT_BASE}"
echo "runThreadN=${RUN_THREADS}, mapPermits=${MAP_PERMITS}, readMapNumber=${READ_MAP_NUMBER}, outSAMtype='${OUT_SAMTYPE}', variableThreads=${VARIABLE_THREADS}, retuneEveryAcquires=${VARIABLE_THREADS_RETUNE_EVERY_ACQUIRES}, permitSequence='${VARIABLE_THREADS_PERMIT_SEQUENCE}', minRetunes=${MIN_RETUNES}, checkLogFinalParity=${CHECK_LOG_FINAL_PARITY}, checkBamParity=${CHECK_BAM_PARITY}"

mkdir -p "${OUT_BASE}"

GENOME_DIR="${OUT_BASE}/genome_index"
echo "=== Building tiny index for current STAR binary ==="
zcat "${GTF_FILE_GZ}" > "${GTF_FILE_UNZIPPED}"
"${STAR_BIN}" \
    --runMode genomeGenerate \
    --runThreadN "${RUN_THREADS}" \
    --outFileNamePrefix "${OUT_BASE}/genome_generate_" \
    --genomeDir "${GENOME_DIR}" \
    --genomeFastaFiles "${FASTA_FILE}" \
    --sjdbGTFfile "${GTF_FILE_UNZIPPED}" \
    --sjdbOverhang 91 \
    --genomeSAindexNbases 10 >/dev/null

run_star_mode "baseline_off" "${OUT_OFF}" 0 0 0 0 0 ""
run_star_mode "dynamic_on" "${OUT_ON}" 1 1 "${MAP_PERMITS}" "${VARIABLE_THREADS}" "${VARIABLE_THREADS_RETUNE_EVERY_ACQUIRES}" "${VARIABLE_THREADS_PERMIT_SEQUENCE}"

if [[ "${CHECK_BAM_PARITY}" == "auto" ]]; then
    if [[ "${OUT_SAMTYPE}" == BAM* ]]; then
        CHECK_BAM_PARITY=1
    else
        CHECK_BAM_PARITY=0
    fi
fi

python3 "${PARSER}" \
    --log "${OUT_ON}/Log.out" \
    --json-out "${REPORT_JSON}" \
    --summary-out "${REPORT_TXT}" \
    --require-telemetry \
    --min-acquires 1 \
    --min-work-units 1 \
    --min-retunes "${MIN_RETUNES}"

off_reads="$(extract_metric "Number of input reads" "${OUT_OFF}/Log.final.out")"
on_reads="$(extract_metric "Number of input reads" "${OUT_ON}/Log.final.out")"

if [[ -z "${off_reads}" || -z "${on_reads}" ]]; then
    echo "ERROR: failed to extract 'Number of input reads' from Log.final.out" >&2
    exit 1
fi

if [[ "${off_reads}" != "${on_reads}" ]]; then
    echo "ERROR: baseline and dynamic modes processed different read counts (${off_reads} vs ${on_reads})" >&2
    exit 1
fi

if [[ "${CHECK_LOG_FINAL_PARITY}" == "1" ]]; then
    OFF_FINAL_METRICS="${OUT_BASE}/off_log_final_metrics.txt"
    ON_FINAL_METRICS="${OUT_BASE}/on_log_final_metrics.txt"
    canonicalize_log_final_metrics "${OUT_OFF}/Log.final.out" > "${OFF_FINAL_METRICS}"
    canonicalize_log_final_metrics "${OUT_ON}/Log.final.out" > "${ON_FINAL_METRICS}"
    if ! diff -u "${OFF_FINAL_METRICS}" "${ON_FINAL_METRICS}" > "${LOG_FINAL_DIFF}"; then
        echo "ERROR: Log.final.out parity failed. See ${LOG_FINAL_DIFF}" >&2
        exit 1
    fi
fi

if [[ "${ASSERT_RETUNE_TRACE_CYCLE}" == "1" && "${VARIABLE_THREADS}" == "1" && "${VARIABLE_THREADS_RETUNE_EVERY_ACQUIRES}" -gt 0 && -n "${VARIABLE_THREADS_PERMIT_SEQUENCE}" ]]; then
    python3 - "${REPORT_JSON}" "${VARIABLE_THREADS_PERMIT_SEQUENCE}" "${RUN_THREADS}" <<'PY'
import json
import re
import sys
from pathlib import Path

report_path = Path(sys.argv[1])
sequence_raw = sys.argv[2]
run_threads = int(sys.argv[3])

report = json.loads(report_path.read_text(encoding="utf-8"))
telemetry = report.get("telemetry") or {}
trace = telemetry.get("retune_trace") or []
dropped = int(telemetry.get("retune_trace_dropped") or 0)
retunes = int(telemetry.get("retunes") or 0)

expected_raw = [int(x) for x in re.split(r"[\s,]+", sequence_raw.strip()) if x]
if not expected_raw:
    raise SystemExit("expected retune sequence is empty")
expected = [max(1, min(x, run_threads)) for x in expected_raw]

if not trace:
    raise SystemExit(f"retune trace is empty (retunes={retunes})")
if dropped != 0:
    raise SystemExit(f"retune trace dropped entries={dropped}, cannot assert full cycle")
if retunes != len(trace):
    raise SystemExit(f"retunes ({retunes}) != trace length ({len(trace)})")

for idx, target in enumerate(trace):
    exp = expected[idx % len(expected)]
    if target != exp:
        raise SystemExit(
            f"retune trace mismatch at index {idx}: got {target}, expected {exp}, trace={trace}, expected_cycle={expected}"
        )

print(f"retune_trace_cycle_ok: trace={trace}, expected_cycle={expected}")
PY
fi

if [[ "${CHECK_BAM_PARITY}" == "1" ]]; then
    if ! command -v samtools >/dev/null 2>&1; then
        echo "ERROR: CHECK_BAM_PARITY=1 requires samtools in PATH" >&2
        exit 1
    fi

    OFF_BAM="$(resolve_bam_path "${OUT_OFF}")"
    ON_BAM="$(resolve_bam_path "${OUT_ON}")"
    if [[ -z "${OFF_BAM}" || -z "${ON_BAM}" ]]; then
        echo "ERROR: unsupported OUT_SAMTYPE for BAM parity checks: ${OUT_SAMTYPE}" >&2
        exit 1
    fi
    require_file "${OFF_BAM}"
    require_file "${ON_BAM}"

    samtools quickcheck "${OFF_BAM}" "${ON_BAM}"
    OFF_BAM_TOTAL_READS="$(samtools view -c "${OFF_BAM}")"
    ON_BAM_TOTAL_READS="$(samtools view -c "${ON_BAM}")"
    OFF_BAM_PRIMARY_READS="$(samtools view -F 0x900 -c "${OFF_BAM}")"
    ON_BAM_PRIMARY_READS="$(samtools view -F 0x900 -c "${ON_BAM}")"
    OFF_BAM_BODY_MD5="$(samtools view "${OFF_BAM}" | md5sum | awk '{print $1}')"
    ON_BAM_BODY_MD5="$(samtools view "${ON_BAM}" | md5sum | awk '{print $1}')"

    {
        echo "OUT_SAMTYPE=${OUT_SAMTYPE}"
        echo "OFF_BAM=${OFF_BAM}"
        echo "ON_BAM=${ON_BAM}"
        echo "OFF_BAM_TOTAL_READS=${OFF_BAM_TOTAL_READS}"
        echo "ON_BAM_TOTAL_READS=${ON_BAM_TOTAL_READS}"
        echo "OFF_BAM_PRIMARY_READS=${OFF_BAM_PRIMARY_READS}"
        echo "ON_BAM_PRIMARY_READS=${ON_BAM_PRIMARY_READS}"
        echo "OFF_BAM_BODY_MD5=${OFF_BAM_BODY_MD5}"
        echo "ON_BAM_BODY_MD5=${ON_BAM_BODY_MD5}"
    } > "${BAM_PARITY_SUMMARY}"

    if [[ "${OFF_BAM_TOTAL_READS}" != "${ON_BAM_TOTAL_READS}" ]]; then
        echo "ERROR: BAM total read count mismatch (${OFF_BAM_TOTAL_READS} vs ${ON_BAM_TOTAL_READS})" >&2
        exit 1
    fi
    if [[ "${OFF_BAM_PRIMARY_READS}" != "${ON_BAM_PRIMARY_READS}" ]]; then
        echo "ERROR: BAM primary read count mismatch (${OFF_BAM_PRIMARY_READS} vs ${ON_BAM_PRIMARY_READS})" >&2
        exit 1
    fi
    if [[ "${REQUIRE_BAM_BODY_MD5_MATCH}" == "1" && "${OFF_BAM_BODY_MD5}" != "${ON_BAM_BODY_MD5}" ]]; then
        echo "ERROR: BAM body MD5 mismatch (${OFF_BAM_BODY_MD5} vs ${ON_BAM_BODY_MD5})" >&2
        exit 1
    fi
fi

echo "PASS: dynamic-thread tiny fixture smoke succeeded"
echo "Report JSON: ${REPORT_JSON}"
echo "Report TXT: ${REPORT_TXT}"
if [[ "${CHECK_LOG_FINAL_PARITY}" == "1" ]]; then
    echo "Log.final.out parity diff (expected empty): ${LOG_FINAL_DIFF}"
fi
if [[ "${CHECK_BAM_PARITY}" == "1" ]]; then
    echo "BAM parity summary: ${BAM_PARITY_SUMMARY}"
fi
echo "Run outputs: ${OUT_BASE}"

if [[ "${KEEP_OUTPUTS}" == "0" ]]; then
    rm -rf "${OUT_BASE}"
    echo "Cleanup complete (KEEP_OUTPUTS=0)."
fi
