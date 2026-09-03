#!/usr/bin/env bash
# Help-path neutrality control for the JAX matrix.
#
# The fused pipeline now aligns hash misses inside the producer whenever alignQ
# is full. The staged pipeline (--flexPipelineNTriage > 0) is a different
# arrangement of the same work: separate reader, triage and solo threads, plus
# dedicated alignment consumers that keep the original blocking push. If both
# arrangements produce the same per-sample outputs on the same reads, then
# aligning inside a producer is output-neutral.
#
# gzip input only: BGZF range readers require the fully-fused pipeline, so the
# staged control cannot use them. Same staged inputs, references, thread count
# and bucket mode as the matrix.
set -uo pipefail

LOG_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${LOG_DIR}/../../.." && pwd)"
STAR_BIN="${STAR_BIN:-${ROOT_DIR}/core/legacy/source/STAR}"
THREADS="${THREADS:-32}"
N_TRIAGE="${N_TRIAGE:-4}"
N_SOLO="${N_SOLO:-4}"
BUCKET_MODE="${BUCKET_MODE:-ram}"
BUCKET_COUNT="${BUCKET_COUNT:-256}"

STAGE_ROOT="${JAX_STAGE_ROOT:-/home/lhhung/jax_stage_20260903}"
STAGE_FASTQ="${STAGE_ROOT}/fastq"
STAGE_REF="${STAGE_ROOT}/ref"
ARTIFACT_ROOT="${ARTIFACT_ROOT:-/home/lhhung/jax_matrix_20260903}"
export TMPDIR="${ARTIFACT_ROOT}/tmp"
GENOME_DIR="${STAGE_REF}/star_index"
PROBE_LIST="${GENOME_DIR}/flex_probe_artifacts/probe_list.txt"
CB_WHITELIST="${STAGE_REF}/737K-fixed-rna-profiling.txt"
HASH_CACHE="${STAGE_REF}/h01_cache.bin"
SAMPLE_WHITELIST="${STAGE_REF}/sample_whitelist_full_16.tsv"
SAMPLE_PROBES="${STAGE_REF}/probe-barcodes-fixed-rna-profiling-rna.txt"

STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
RUN_DIR="${ARTIFACT_ROOT}/staged_control_${STAMP}"
PREFIX="${LOG_DIR}/staged_control"

die() { echo "FATAL: $*" >&2; exit 1; }
for p in "${STAR_BIN}" "${GENOME_DIR}/Genome" "${PROBE_LIST}" "${CB_WHITELIST}" \
         "${HASH_CACHE}" "${SAMPLE_WHITELIST}" "${SAMPLE_PROBES}" "${STAGE_FASTQ}"; do
    [[ -e "${p}" ]] || die "required input is absent: ${p} (run the matrix first so the stage exists)"
done

mapfile -t reads_r2 < <(find "${STAGE_FASTQ}" -maxdepth 1 -type f -name '*_R2_001.fastq.gz' | sort)
mapfile -t reads_r1 < <(find "${STAGE_FASTQ}" -maxdepth 1 -type f -name '*_R1_001.fastq.gz' | sort)
[[ "${#reads_r1[@]}" -eq 8 && "${#reads_r2[@]}" -eq 8 ]] || die "expected 8 staged FASTQs per mate"
join_comma() { local IFS=,; echo "$*"; }

mkdir -p "${RUN_DIR}" "${TMPDIR}"
cmd=(
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
    --outSJtype None --outSAMtype None --outSAMattributes None
    --soloFeatures Gene --soloCellFilter None --soloMultiMappers Rescue
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
    --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR
    --soloStrand Unstranded --clipAdapterType CellRanger4
    --alignEndsType Local --chimSegmentMin 0
    --soloKeysCompat cr --soloSampleSearchNearby no
    --soloHashScreenFile "${HASH_CACHE}"
    --soloBucketMode "${BUCKET_MODE}" --soloBucketCount "${BUCKET_COUNT}"
    --flexPipeline yes
    --flexPipelineNTriage "${N_TRIAGE}"
    --flexPipelineNSolo "${N_SOLO}"
    --flexNoAlign 0
    --readFilesBgzfMode off
    --dynamicThreadInterface 1
    --crAssignConsumerThreads -1 --crAssignSearchThreads 1
    --outFileNamePrefix "${RUN_DIR}/"
    --readFilesIn "$(join_comma "${reads_r2[@]}")" "$(join_comma "${reads_r1[@]}")"
)
{ printf '%q ' "${cmd[@]}"; printf '\n'; } > "${PREFIX}.command.txt"

sudo -n sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches' || die "cache drop unavailable"
echo "[$(date -u +%FT%TZ)] staged control starting (triage=${N_TRIAGE}, solo=${N_SOLO}, dedicated align workers)"
/usr/bin/time -v -o "${PREFIX}.time.txt" "${cmd[@]}" \
    > "${PREFIX}.stdout.txt" 2> "${PREFIX}.stderr.txt"
rc=$?
for log in Log.out Log.final.out Log.progress.out; do
    [[ -f "${RUN_DIR}/${log}" ]] && cp "${RUN_DIR}/${log}" "${PREFIX}.${log}"
done

python3 - "${RUN_DIR}" "${PREFIX}.manifest.tsv" <<'PY'
import gzip, hashlib, sys
from pathlib import Path
root, out = Path(sys.argv[1]), Path(sys.argv[2])
names = {"matrix.mtx", "matrix.mtx.gz", "barcodes.tsv", "barcodes.tsv.gz",
         "features.tsv", "features.tsv.gz", "Summary.csv", "Barcodes.stats"}
rows = []
for path in sorted(root.rglob("*")):
    if path.is_file() and path.name in names:
        payload = gzip.open(path, "rb").read() if path.suffix == ".gz" else path.read_bytes()
        rows.append(f"{path.relative_to(root)}\t{hashlib.sha256(payload).hexdigest()}")
out.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")
PY

printf 'exit=%s\nrun_dir=%s\n' "${rc}" "${RUN_DIR}" > "${PREFIX}.completion.txt"
echo "[$(date -u +%FT%TZ)] staged control finished with exit ${rc}; run dir ${RUN_DIR}"
exit "${rc}"
