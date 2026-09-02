#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-}"
case "${MODE}" in gzip|bgzf|cbq) ;; *) echo "Usage: $0 gzip|bgzf|cbq" >&2; exit 2 ;; esac

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
STAR_BIN="${STAR_BIN:-${ROOT_DIR}/core/legacy/source/STAR}"
THREADS="${THREADS:-32}"
FASTQ_DIR="${SCFFPE_320K_FASTQ_DIR:-/mnt/pikachu/tenx_320k_scFFPE/320k_scFFPE_16-plex_GEM-X_FLEX_fastqs}"
CBQ_DIR="${SCFFPE_320K_CBQ_DIR:-/mnt/pikachu/tenx_320k_scFFPE/320k_scFFPE_16-plex_GEM-X_FLEX_cbq}"
ARTIFACT_ROOT="${ARTIFACT_ROOT:-/home/lhhung/STAR-suite-bgzf-ingest-benchmark-20260901}"
LOG_DIR="${ROOT_DIR}/docs/benchmarks/bgzf_ingest_20260901"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
RUN_DIR="${ARTIFACT_ROOT}/scffpe320k_${MODE}_${STAMP}"
LOG_PREFIX="${LOG_DIR}/scffpe320k_${MODE}_${STAMP}"
SOLO_BUCKET_MODE="${SOLO_BUCKET_MODE:-auto}"
SOLO_BUCKET_MEM_GB="${SOLO_BUCKET_MEM_GB:-32}"
SOLO_BUCKET_COUNT="${SOLO_BUCKET_COUNT:-256}"
SOLO_BUCKET_SPILL_DIR="${SOLO_BUCKET_SPILL_DIR:-${RUN_DIR}/bucket_spill}"

case "${SOLO_BUCKET_MODE}" in
    off|ram|spill|auto) ;;
    *) echo "ERROR: SOLO_BUCKET_MODE must be off, ram, spill, or auto" >&2; exit 2 ;;
esac

GENOME_DIR="${GENOME_DIR:-/storage/flex_filtered_reference_2024/star_index}"
CB_WHITELIST="${CB_WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}"
SAMPLE_WHITELIST="${SAMPLE_WHITELIST:-/mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv}"
PROBE_LIST="${PROBE_LIST:-/storage/flex_filtered_reference_2024/star_index/flex_probe_artifacts/probe_list.txt}"
SAMPLE_PROBES="${SAMPLE_PROBES:-/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}"
HASH_CACHE="${HASH_CACHE:-/storage/downsampled_100K/SC2300771/results/flex_h01_2024_20260320_081246/h01_cache.bin}"

for path in "${STAR_BIN}" "${GENOME_DIR}" "${CB_WHITELIST}" "${SAMPLE_WHITELIST}" \
            "${PROBE_LIST}" "${SAMPLE_PROBES}" "${HASH_CACHE}"; do
    [[ -e "${path}" ]] || { echo "ERROR: required input is absent: ${path}" >&2; exit 1; }
done
[[ ! -e "${RUN_DIR}" ]] || { echo "ERROR: fresh run directory already exists: ${RUN_DIR}" >&2; exit 1; }

join_comma() { local IFS=,; echo "$*"; }
case "${MODE}" in
    gzip|bgzf)
        mapfile -t reads_r2 < <(find "${FASTQ_DIR}" -maxdepth 1 -type f \
            -name '16-plex_GEM-X_FLEX_S1_L*_R2_001.fastq.gz' | sort)
        mapfile -t reads_r1 < <(find "${FASTQ_DIR}" -maxdepth 1 -type f \
            -name '16-plex_GEM-X_FLEX_S1_L*_R1_001.fastq.gz' | sort)
        [[ "${#reads_r1[@]}" -eq 4 && "${#reads_r2[@]}" -eq 4 ]] || {
            echo "ERROR: expected 4 320K scFFPE FASTQs per mate" >&2; exit 1;
        }
        reads_r2_csv="$(join_comma "${reads_r2[@]}")"
        reads_r1_csv="$(join_comma "${reads_r1[@]}")"
        ;;
    cbq)
        mapfile -t cbq_files < <(find "${CBQ_DIR}" -maxdepth 1 -type f \
            -name 'lane_*.cbq' | sort)
        [[ "${#cbq_files[@]}" -eq 4 ]] || {
            echo "ERROR: expected 4 ordered CBQ lanes" >&2; exit 1;
        }
        cbq_csv="$(join_comma "${cbq_files[@]}")"
        ;;
esac

mkdir -p "${RUN_DIR}" "${LOG_DIR}"
{
    date -u '+date_utc=%Y-%m-%dT%H:%M:%SZ'
    printf 'mode=%s\nthreads=%s\nrun_dir=%s\ncommit=%s\n' \
        "${MODE}" "${THREADS}" "${RUN_DIR}" "$(git -C "${ROOT_DIR}" rev-parse HEAD)"
    printf 'solo_bucket_mode=%s\nsolo_bucket_mem_gb=%s\nsolo_bucket_count=%s\nsolo_bucket_spill_dir=%s\n' \
        "${SOLO_BUCKET_MODE}" "${SOLO_BUCKET_MEM_GB}" "${SOLO_BUCKET_COUNT}" \
        "${SOLO_BUCKET_SPILL_DIR}"
    uptime
    ps -eo pid,comm,%cpu,%mem,etime --sort=-%cpu | sed -n '1,20p'
} > "${LOG_PREFIX}.preflight.txt"

common=(
    "${STAR_BIN}"
    --runThreadN "${THREADS}"
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
    --soloFlexOutputPrefix "${RUN_DIR}/per_sample"
    --limitIObufferSize 50000000 50000000
    --outSJtype None
    --outSAMtype None
    --outSAMattributes None
    --soloFeatures Gene
    --soloCellFilter None
    --soloMultiMappers Rescue
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
    --soloUMIfiltering MultiGeneUMI_CR
    --soloUMIdedup 1MM_CR
    --soloStrand Unstranded
    --clipAdapterType CellRanger4
    --alignEndsType Local
    --chimSegmentMin 0
    --soloKeysCompat cr
    --soloSampleSearchNearby no
    --soloHashScreenFile "${HASH_CACHE}"
    --soloBucketMode "${SOLO_BUCKET_MODE}"
    --soloBucketMemGB "${SOLO_BUCKET_MEM_GB}"
    --soloBucketCount "${SOLO_BUCKET_COUNT}"
    --soloBucketSpillDir "${SOLO_BUCKET_SPILL_DIR}"
    --flexPipeline yes
    --flexPipelineNTriage 0
    --flexPipelineNSolo 0
    --flexNoAlign 1
    --dynamicThreadInterface 1
    --crAssignConsumerThreads -1
    --crAssignSearchThreads 1
    --outFileNamePrefix "${RUN_DIR}/"
)

case "${MODE}" in
    gzip) input=(--readFilesBgzfMode off --readFilesIn "${reads_r2_csv}" "${reads_r1_csv}") ;;
    bgzf) input=(--readFilesBgzfMode range --bgzfReaderThreads "${THREADS}" \
                 --bgzfCrcCheck 1 --readFilesIn "${reads_r2_csv}" "${reads_r1_csv}") ;;
    cbq) input=(--readFilesType Binseq PE --readFilesCbqRangeMode range \
                --readFilesIn "${cbq_csv}") ;;
esac
cmd=("${common[@]}" "${input[@]}")
printf '%q ' "${cmd[@]}" > "${LOG_PREFIX}.command.txt"
printf '\n' >> "${LOG_PREFIX}.command.txt"

# Keep paper timings comparable: reject hosts that cannot provide a cold cache
# instead of silently recording a warm-cache result.
sudo -n sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches'
/usr/bin/time -v -o "${LOG_PREFIX}.time.txt" "${cmd[@]}" \
    > "${LOG_PREFIX}.stdout.txt" 2> "${LOG_PREFIX}.stderr.txt"

if ! grep -Fq 'Flex count-only no-genome: active' "${RUN_DIR}/Log.out"; then
    echo "ERROR: cache-only/no-genome path did not activate; benchmark is invalid" >&2
    exit 1
fi

for log in Log.out Log.final.out Log.progress.out SJ.out.tab; do
    [[ -f "${RUN_DIR}/${log}" ]] && cp "${RUN_DIR}/${log}" "${LOG_PREFIX}.${log}"
done
printf 'status=success\nrun_dir=%s\n' "${RUN_DIR}" > "${LOG_PREFIX}.completion.txt"
echo "PASS: ${MODE} benchmark completed at ${RUN_DIR}"
