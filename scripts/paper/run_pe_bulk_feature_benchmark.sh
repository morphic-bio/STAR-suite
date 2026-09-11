#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
GATHER_SCRIPT="${SCRIPT_DIR}/gather_pe_bulk_external_tools.sh"

OUT_ROOT=""
TOOL_ENV=""
THREADS="16"
STAGES="downsampled,storage,pikachu"
DRY_RUN=0
SALMON_LIBTYPE="${SALMON_LIBTYPE:-IU}"

STAR_INDEX="${STAR_INDEX:-/storage/autoindex_110_44/bulk_index}"
TRANSCRIPTOME="${TRANSCRIPTOME:-}"
GTF="${GTF:-}"

DOWNSAMPLED_FASTQ_DIR="/storage/PE/downsampled"
STORAGE_FASTQ_DIR="/storage/PE"
PIKACHU_FASTQ_DIR="/mnt/pikachu/paper_pe_bulk_benchmark_inputs/PE"

QUALITY_CUTOFF="20"
MIN_LENGTH="20"
ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
YREMOVE=1
FORCE_COPY=0
SKIP_COPY=0
SKIP_INTEGRATED=0
SKIP_EXTERNAL=0
SKIP_COMPARE=0
PARITY_QC=0
TRIM_QC=1
TRIM_QC_BASENAME="${TRIM_QC_BASENAME:-read_qc}"
TRIM_QC_MAX_READS="${TRIM_QC_MAX_READS:-250000}"
EXTERNAL_TRIMMER="${EXTERNAL_TRIMMER:-trim_galore}"
TRIMGALORE_CORES="${TRIMGALORE_CORES:-}"
EXTERNAL_YREMOVE="${EXTERNAL_YREMOVE:-awk}"

usage() {
    cat <<'EOF'
Usage:
  scripts/paper/run_pe_bulk_feature_benchmark.sh --out-root DIR [options]

Required:
  --out-root DIR                 Root output directory for all benchmark stages.

Optional:
  --tool-env FILE                Existing tool_env.sh produced by gather_pe_bulk_external_tools.sh.
  --threads N                    Threads for STAR and Salmon (default: 16).
  --salmon-libtype STR           Salmon alignment-mode libtype for QC quant (default: IU).
  --stages LIST                  Comma-separated: downsampled,storage,pikachu (default: all).
  --star-index DIR               STAR index (default: /storage/autoindex_110_44/bulk_index).
  --transcriptome FILE           Transcriptome FASTA (default: <star-index>/transcriptome.fa).
  --gtf FILE                     GTF for tx2gene generation (default: <star-index>/cellranger_ref/genes.gtf).
  --downsampled-fastq-dir DIR    Downsampled PE directory (default: /storage/PE/downsampled).
  --storage-fastq-dir DIR        Full PE directory (default: /storage/PE).
  --pikachu-fastq-dir DIR        HDD-local PE directory (default: /mnt/pikachu/paper_pe_bulk_benchmark_inputs/PE).
  --quality Q                    Adapter trim quality cutoff for trimvalidate (default: 20).
  --min-length N                 Adapter trim minimum length for trimvalidate (default: 20).
  --adapter-r1 SEQ               R1 adapter for trimvalidate.
  --adapter-r2 SEQ               R2 adapter for trimvalidate.
  --external-trimmer MODE        External trimming arm: trim_galore or trimvalidate
                                  (default: trim_galore).
  --trimgalore-cores N           Cores passed to Trim Galore --cores
                                  (default: --threads value).
  --external-yremove MODE        External Y-removal arm: awk or remove_y_reads
                                  (default: awk).
  --yremove                      Enable Y-chromosome removal (default).
  --no-yremove                   Disable Y-chromosome removal on both arms.
  --skip-copy                    Do not copy /storage PE inputs into --pikachu-fastq-dir.
  --force-copy                   Re-copy into --pikachu-fastq-dir even if files already exist.
  --skip-integrated              Skip integrated STAR-suite arm.
  --skip-external                Skip external stepwise arm.
  --integrated-only              Run only the STAR-suite production arm
                                  (implies --skip-external --skip-compare).
  --skip-compare                 Skip quant/Y-split comparison stage.
  --parity-qc                    Also emit integrated TranscriptomeSAM and run
                                  Salmon QC/comparison. Excluded from the
                                  production benchmark by default.
  --trim-qc                      Emit FastQC-like trim-QC reports on both arms
                                  and include external trim-QC timing (default).
  --no-trim-qc                   Disable FastQC-like trim-QC reports on both arms.
  --trim-qc-max-reads N          Limit reads sampled by trim-QC reporting
                                  (default: 250000; 0 = no limit).
  --dry-run                      Print resolved commands only.
  -h, --help                     Show this help.

Modes:
  --yremove (default):
    Integrated arm:
      raw FASTQ -> STAR (trimCutadapt + trim-QC + emitNoYBAM
                + emitYNoYFastq + TranscriptVB)

    External arm:
      raw FASTQ -> Trim Galore + FastQC -> STAR (full BAM + TranscriptomeSAM)
                -> awk/samtools Y-removal from BAM-derived ynames
                -> gzip-compressed Y/no-Y FASTQs
                -> Salmon alignment-mode QC on no-Y transcriptome BAM

  --no-yremove:
    Integrated arm:
      raw FASTQ -> STAR (trimCutadapt + trim-QC + TranscriptVB)

    External arm:
      raw FASTQ -> Trim Galore + FastQC -> STAR (TranscriptomeSAM)
                -> Salmon alignment-mode QC on transcriptome BAM

Note:
  Salmon auto libtype detection (-l A) mis-detects the current PE benchmark sample as ISR.
  This script therefore defaults Salmon QC to explicit IU for reproducible comparison.
  STAR-suite internal TranscriptVB is the production quantifier. Integrated
  TranscriptomeSAM and integrated Salmon QC are disabled by default; use
  --parity-qc when parity artifacts are required.
EOF
}

log() {
    printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*"
}

die() {
    echo "ERROR: $*" >&2
    exit 2
}

render_cmd() {
    local out=""
    local arg
    for arg in "$@"; do
        out+=$(printf '%q ' "$arg")
    done
    printf '%s\n' "${out% }"
}

run_or_echo() {
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        log "DRY_RUN: $(render_cmd "$@")"
    else
        "$@"
    fi
}

run_timed_cmd() {
    local time_log="$1"
    local stdout_log="$2"
    shift 2
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        log "DRY_RUN: /usr/bin/time -v -o ${time_log} $(render_cmd "$@")"
    else
        /usr/bin/time -v -o "${time_log}" "$@" 2>&1 | tee "${stdout_log}"
    fi
}

require_file() {
    [[ -f "$1" ]] || die "Missing file: $1"
}

require_dir() {
    [[ -d "$1" ]] || die "Missing directory: $1"
}

extract_time_metric() {
    local log_file="$1"
    local pattern="$2"
    if [[ ! -f "${log_file}" ]]; then
        echo "NA"
        return
    fi
    grep -F "${pattern}" "${log_file}" | head -n1 \
        | sed -E 's/^.*\):[[:space:]]*//; t; s/^.*:[[:space:]]*//'
}

write_cmd_script() {
    local script_path="$1"
    shift
    {
        echo '#!/usr/bin/env bash'
        echo 'set -euo pipefail'
        render_cmd "$@"
    } > "${script_path}"
    chmod +x "${script_path}"
}

find_single_pair() {
    local dir="$1"
    mapfile -t r1_files < <(find "${dir}" -maxdepth 1 -type f -name '*_R1_001.fastq.gz' | sort)
    [[ "${#r1_files[@]}" -gt 0 ]] || die "No *_R1_001.fastq.gz files found in ${dir}"
    [[ "${#r1_files[@]}" -eq 1 ]] || die "Expected exactly one R1 FASTQ in ${dir}, found ${#r1_files[@]}"
    local r1="${r1_files[0]}"
    local r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    [[ -f "${r2}" ]] || die "Matching R2 FASTQ not found for ${r1}"
    printf '%s\n%s\n' "${r1}" "${r2}"
}

ensure_tx2gene() {
    local out_tsv="$1"
    if [[ -s "${out_tsv}" ]]; then
        return
    fi
    run_or_echo "${MAKE_GENE_MAP_SCRIPT}" --gtf "${GTF}" --out "${out_tsv}"
    if [[ "${DRY_RUN}" -eq 0 ]]; then
        require_file "${out_tsv}"
    fi
}

count_fastq_reads() {
    local fastq="$1"
    local magic=""
    if [[ ! -f "${fastq}" ]]; then
        echo "NA"
        return
    fi
    magic="$(head -c 2 "${fastq}" 2>/dev/null | od -An -tx1 | tr -d ' \n')"
    if [[ "${fastq}" == *.gz || "${magic}" == "1f8b" ]]; then
        gzip -dc "${fastq}" | awk 'END {print NR/4}'
    else
        awk 'END {print NR/4}' "${fastq}"
    fi
}

write_awk_y_remove_script() {
    local script_path="$1"
    cat > "${script_path}" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 9 ]]; then
    echo "Usage: $0 SAMTOOLS AWK GENOME_BAM TRANSCRIPTOME_BAM R1_FASTQ R2_FASTQ BAM_OUT_DIR FASTQ_OUT_DIR SAMPLE_PREFIX" >&2
    exit 2
fi

SAMTOOLS_BIN="$1"
AWK_BIN="$2"
GENOME_BAM="$3"
TRANSCRIPTOME_BAM="$4"
R1_FASTQ="$5"
R2_FASTQ="$6"
BAM_OUT_DIR="$7"
FASTQ_OUT_DIR="$8"
SAMPLE_PREFIX="$9"

[[ -x "${SAMTOOLS_BIN}" ]] || { echo "samtools not executable: ${SAMTOOLS_BIN}" >&2; exit 2; }
[[ -x "${AWK_BIN}" ]] || { echo "awk not executable: ${AWK_BIN}" >&2; exit 2; }
[[ -f "${GENOME_BAM}" ]] || { echo "genome BAM not found: ${GENOME_BAM}" >&2; exit 2; }
[[ -f "${TRANSCRIPTOME_BAM}" ]] || { echo "transcriptome BAM not found: ${TRANSCRIPTOME_BAM}" >&2; exit 2; }
[[ -f "${R1_FASTQ}" ]] || { echo "R1 FASTQ not found: ${R1_FASTQ}" >&2; exit 2; }
[[ -f "${R2_FASTQ}" ]] || { echo "R2 FASTQ not found: ${R2_FASTQ}" >&2; exit 2; }
mkdir -p "${BAM_OUT_DIR}" "${FASTQ_OUT_DIR}"

Y_NAMES="${BAM_OUT_DIR}/ynames.txt"
GENOME_Y_BAM="${BAM_OUT_DIR}/Aligned.sortedByCoord.out_Y.bam"
GENOME_NOY_BAM="${BAM_OUT_DIR}/Aligned.sortedByCoord.out_noY.bam"
TX_Y_BAM="${BAM_OUT_DIR}/Aligned.toTranscriptome.out_Y.bam"
TX_NOY_BAM="${BAM_OUT_DIR}/Aligned.toTranscriptome.out_noY.bam"
Y_R1="${FASTQ_OUT_DIR}/${SAMPLE_PREFIX}_R1_Y.fastq.gz"
Y_R2="${FASTQ_OUT_DIR}/${SAMPLE_PREFIX}_R2_Y.fastq.gz"
NOY_R1="${FASTQ_OUT_DIR}/${SAMPLE_PREFIX}_R1_noY.fastq.gz"
NOY_R2="${FASTQ_OUT_DIR}/${SAMPLE_PREFIX}_R2_noY.fastq.gz"

"${SAMTOOLS_BIN}" view "${GENOME_BAM}" \
  | "${AWK_BIN}" '$3 == "chrY" {
        name=$1
        sub(/\/[12]$/, "", name)
        unique[name]=1
    }
    END {
        for (name in unique) print name
    }' > "${Y_NAMES}"

"${SAMTOOLS_BIN}" view -h -b -N "${Y_NAMES}" \
    -o "${GENOME_Y_BAM}" \
    -U "${GENOME_NOY_BAM}" \
    "${GENOME_BAM}"

"${SAMTOOLS_BIN}" view -h -b -N "${Y_NAMES}" \
    -o "${TX_Y_BAM}" \
    -U "${TX_NOY_BAM}" \
    "${TRANSCRIPTOME_BAM}"

split_one_fastq() {
    local in_fastq="$1"
    local y_out="$2"
    local noy_out="$3"
    local tmp_dir="${FASTQ_OUT_DIR}/.awk_filter.$(basename "${in_fastq}").$$"
    local y_fifo="${tmp_dir}/y.fifo"
    local noy_fifo="${tmp_dir}/noy.fifo"
    mkdir -p "${tmp_dir}"
    mkfifo "${y_fifo}" "${noy_fifo}"

    gzip -c < "${y_fifo}" > "${y_out}" &
    local y_pid=$!
    gzip -c < "${noy_fifo}" > "${noy_out}" &
    local noy_pid=$!

    if [[ "${in_fastq}" == *.gz ]]; then
        gzip -dc -- "${in_fastq}"
    else
        cat -- "${in_fastq}"
    fi | "${AWK_BIN}" -v y_names="${Y_NAMES}" -v y_out="${y_fifo}" -v noy_out="${noy_fifo}" '
        BEGIN {
            while ((getline line < y_names) > 0) {
                if (line != "") y[line]=1
            }
            close(y_names)
        }
        NR % 4 == 1 {
            read=$0
            sub(/^@/, "", read)
            sub(/[[:space:]].*/, "", read)
            sub(/\/[12]$/, "", read)
            dest = (read in y) ? y_out : noy_out
        }
        { print >> dest }
        END {
            close(y_out)
            close(noy_out)
        }'

    wait "${y_pid}"
    wait "${noy_pid}"
    rm -rf "${tmp_dir}"
}

split_one_fastq "${R1_FASTQ}" "${Y_R1}" "${NOY_R1}"
split_one_fastq "${R2_FASTQ}" "${Y_R2}" "${NOY_R2}"

for fastq in "${Y_R1}" "${Y_R2}" "${NOY_R1}" "${NOY_R2}"; do
    reads="$(gzip -dc -- "${fastq}" | "${AWK_BIN}" 'END {print NR/4}')"
    printf '%s: reads=%s\n' "${fastq}" "${reads}"
done
EOF
    chmod +x "${script_path}"
}

append_quant_compare_row() {
    local stage_name="$1"
    local label="$2"
    local file_a="$3"
    local file_b="$4"
    local out_tsv="$5"
    python3 - "${stage_name}" "${label}" "${file_a}" "${file_b}" "${COMPARE_SALMON_STAR_SCRIPT}" "${out_tsv}" <<'PY'
import importlib.util
import math
import os
import sys
from pathlib import Path

stage_name = sys.argv[1]
label = sys.argv[2]
file_a = Path(sys.argv[3])
file_b = Path(sys.argv[4])
compare_script = Path(sys.argv[5])
out_tsv = Path(sys.argv[6])

columns = [
    "stage",
    "label",
    "file_a",
    "file_b",
    "status",
    "spearman_all",
    "pearson_all",
    "spearman_expressed",
    "pearson_expressed",
    "total_transcripts",
    "jointly_expressed",
    "a_only",
    "b_only",
    "a_total_reads",
    "b_total_reads",
    "abs_total_reads_diff",
    "mean_abs_diff",
    "max_abs_diff",
]

def fmt(value):
    if value is None:
        return "NA"
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return "NA"
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)

if not out_tsv.exists():
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_tsv.write_text("\t".join(columns) + "\n")

if not file_a.exists() or not file_b.exists():
    row = {
        "stage": stage_name,
        "label": label,
        "file_a": str(file_a),
        "file_b": str(file_b),
        "status": "MISSING",
    }
else:
    spec = importlib.util.spec_from_file_location("compare_salmon_star", compare_script)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    try:
        a_df = mod.load_quant_sf(str(file_a))
        b_df = mod.load_quant_sf(str(file_b))
        result = mod.compare_quantifications(a_df, b_df, verbose=False)
        row = {
            "stage": stage_name,
            "label": label,
            "file_a": str(file_a),
            "file_b": str(file_b),
            "status": "OK",
            "spearman_all": result.get("spearman_all"),
            "pearson_all": result.get("pearson_all"),
            "spearman_expressed": result.get("spearman_expressed"),
            "pearson_expressed": result.get("pearson_expressed"),
            "total_transcripts": result.get("total_transcripts"),
            "jointly_expressed": result.get("jointly_expressed"),
            "a_only": result.get("salmon_only"),
            "b_only": result.get("star_only"),
            "a_total_reads": result.get("salmon_total_reads"),
            "b_total_reads": result.get("star_total_reads"),
            "abs_total_reads_diff": abs(result.get("salmon_total_reads", 0.0) - result.get("star_total_reads", 0.0)),
            "mean_abs_diff": result.get("mean_diff"),
            "max_abs_diff": result.get("max_diff"),
        }
    except Exception as exc:
        row = {
            "stage": stage_name,
            "label": label,
            "file_a": str(file_a),
            "file_b": str(file_b),
            "status": f"ERROR:{exc}",
        }

with out_tsv.open("a") as handle:
    handle.write("\t".join(fmt(row.get(col)) for col in columns) + "\n")
PY
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out-root) OUT_ROOT="$2"; shift 2 ;;
        --tool-env) TOOL_ENV="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --salmon-libtype) SALMON_LIBTYPE="$2"; shift 2 ;;
        --stages) STAGES="$2"; shift 2 ;;
        --star-index) STAR_INDEX="$2"; shift 2 ;;
        --transcriptome) TRANSCRIPTOME="$2"; shift 2 ;;
        --gtf) GTF="$2"; shift 2 ;;
        --downsampled-fastq-dir) DOWNSAMPLED_FASTQ_DIR="$2"; shift 2 ;;
        --storage-fastq-dir) STORAGE_FASTQ_DIR="$2"; shift 2 ;;
        --pikachu-fastq-dir) PIKACHU_FASTQ_DIR="$2"; shift 2 ;;
        --quality) QUALITY_CUTOFF="$2"; shift 2 ;;
        --min-length) MIN_LENGTH="$2"; shift 2 ;;
        --adapter-r1) ADAPTER_R1="$2"; shift 2 ;;
        --adapter-r2) ADAPTER_R2="$2"; shift 2 ;;
        --external-trimmer) EXTERNAL_TRIMMER="$2"; shift 2 ;;
        --trimgalore-cores) TRIMGALORE_CORES="$2"; shift 2 ;;
        --external-yremove) EXTERNAL_YREMOVE="$2"; shift 2 ;;
        --yremove) YREMOVE=1; shift ;;
        --no-yremove) YREMOVE=0; shift ;;
        --skip-copy) SKIP_COPY=1; shift ;;
        --force-copy) FORCE_COPY=1; shift ;;
        --skip-integrated) SKIP_INTEGRATED=1; shift ;;
        --skip-external) SKIP_EXTERNAL=1; shift ;;
        --integrated-only) SKIP_EXTERNAL=1; SKIP_COMPARE=1; shift ;;
        --skip-compare) SKIP_COMPARE=1; shift ;;
        --parity-qc) PARITY_QC=1; shift ;;
        --trim-qc) TRIM_QC=1; shift ;;
        --no-trim-qc) TRIM_QC=0; shift ;;
        --trim-qc-max-reads) TRIM_QC_MAX_READS="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

[[ -n "${OUT_ROOT}" ]] || die "--out-root is required"
[[ "${THREADS}" =~ ^[0-9]+$ ]] || die "--threads must be an integer"
[[ "${THREADS}" -gt 0 ]] || die "--threads must be > 0"
[[ -n "${SALMON_LIBTYPE}" ]] || die "--salmon-libtype must be non-empty"
[[ "${QUALITY_CUTOFF}" =~ ^[0-9]+$ ]] || die "--quality must be an integer"
[[ "${MIN_LENGTH}" =~ ^[0-9]+$ ]] || die "--min-length must be an integer"
[[ "${TRIM_QC_MAX_READS}" =~ ^[0-9]+$ ]] || die "--trim-qc-max-reads must be an integer"
case "${EXTERNAL_TRIMMER}" in
    trim_galore|trimvalidate) ;;
    *) die "--external-trimmer must be trim_galore or trimvalidate" ;;
esac
case "${EXTERNAL_YREMOVE}" in
    awk|remove_y_reads) ;;
    *) die "--external-yremove must be awk or remove_y_reads" ;;
esac
if [[ -z "${TRIMGALORE_CORES}" ]]; then
    TRIMGALORE_CORES="${THREADS}"
fi
[[ "${TRIMGALORE_CORES}" =~ ^[0-9]+$ ]] || die "--trimgalore-cores must be an integer"
[[ "${TRIMGALORE_CORES}" -gt 0 ]] || die "--trimgalore-cores must be > 0"

if [[ -z "${TRANSCRIPTOME}" ]]; then
    TRANSCRIPTOME="${STAR_INDEX}/transcriptome.fa"
fi
if [[ -z "${GTF}" ]]; then
    GTF="${STAR_INDEX}/cellranger_ref/genes.gtf"
fi

require_dir "${STAR_INDEX}"
require_file "${TRANSCRIPTOME}"
require_file "${GTF}"
require_file "${GATHER_SCRIPT}"
command -v python3 >/dev/null 2>&1 || die "python3 is required"
command -v gzip >/dev/null 2>&1 || die "gzip is required"
command -v samtools >/dev/null 2>&1 || die "samtools is required"
command -v awk >/dev/null 2>&1 || die "awk is required"
require_file /usr/bin/time

OUT_ROOT="$(mkdir -p "${OUT_ROOT}" && cd "${OUT_ROOT}" && pwd)"

if [[ -z "${TOOL_ENV}" ]]; then
    TOOL_DIR="${OUT_ROOT}/tools"
    mkdir -p "${TOOL_DIR}"
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        "${GATHER_SCRIPT}" --outdir "${TOOL_DIR}" --no-build-missing
    else
        "${GATHER_SCRIPT}" --outdir "${TOOL_DIR}"
    fi
    TOOL_ENV="${TOOL_DIR}/tool_env.sh"
fi

require_file "${TOOL_ENV}"
# shellcheck source=/dev/null
source "${TOOL_ENV}"

[[ -x "${STAR_BIN}" ]] || die "Resolved STAR_BIN not executable: ${STAR_BIN}"
[[ -x "${SALMON_BIN}" ]] || die "Resolved SALMON_BIN not executable: ${SALMON_BIN}"
[[ -x "${TRIMVALIDATE_BIN}" ]] || die "Resolved TRIMVALIDATE_BIN not executable: ${TRIMVALIDATE_BIN}"
if [[ "${EXTERNAL_TRIMMER}" == "trim_galore" ]]; then
    [[ -x "${TRIM_GALORE_BIN:-}" ]] || die "Resolved TRIM_GALORE_BIN not executable: ${TRIM_GALORE_BIN:-unset}"
    [[ -x "${FASTQC_BIN:-}" ]] || die "Resolved FASTQC_BIN not executable: ${FASTQC_BIN:-unset}"
fi
if [[ "${EXTERNAL_YREMOVE}" == "awk" ]]; then
    AWK_BIN="${AWK_BIN:-$(command -v awk)}"
    [[ -x "${AWK_BIN}" ]] || die "Resolved AWK_BIN not executable: ${AWK_BIN}"
fi
if [[ "${TRIM_QC}" -eq 1 ]]; then
    [[ -x "${TRIM_QC_FASTQ_BIN:-}" ]] || die "Resolved TRIM_QC_FASTQ_BIN not executable: ${TRIM_QC_FASTQ_BIN:-unset}"
fi
[[ -x "${REMOVE_Y_READS_BIN}" ]] || die "Resolved REMOVE_Y_READS_BIN not executable: ${REMOVE_Y_READS_BIN}"
[[ -x "${MAKE_GENE_MAP_SCRIPT}" ]] || die "Resolved MAKE_GENE_MAP_SCRIPT not executable: ${MAKE_GENE_MAP_SCRIPT}"
require_file "${COMPARE_SALMON_STAR_SCRIPT}"

TX2GENE="${OUT_ROOT}/inputs/tx2gene.tsv"
mkdir -p "${OUT_ROOT}/inputs"
ensure_tx2gene "${TX2GENE}"

RUN_MANIFEST="${OUT_ROOT}/RUN_MANIFEST.txt"
{
    echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "repo_root=${REPO_ROOT}"
    echo "repo_commit=$(git -C "${REPO_ROOT}" rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "star_index=${STAR_INDEX}"
    echo "transcriptome=${TRANSCRIPTOME}"
    echo "gtf=${GTF}"
    echo "tx2gene=${TX2GENE}"
    echo "threads=${THREADS}"
    echo "salmon_libtype=${SALMON_LIBTYPE}"
    echo "stages=${STAGES}"
    echo "yremove=${YREMOVE}"
    echo "dry_run=${DRY_RUN}"
    echo "skip_copy=${SKIP_COPY}"
    echo "skip_integrated=${SKIP_INTEGRATED}"
    echo "skip_external=${SKIP_EXTERNAL}"
    echo "skip_compare=${SKIP_COMPARE}"
    echo "parity_qc=${PARITY_QC}"
    echo "trim_qc=${TRIM_QC}"
    echo "trim_qc_basename=${TRIM_QC_BASENAME}"
    echo "trim_qc_max_reads=${TRIM_QC_MAX_READS}"
    echo "external_trimmer=${EXTERNAL_TRIMMER}"
    echo "trimgalore_cores=${TRIMGALORE_CORES}"
    echo "external_yremove=${EXTERNAL_YREMOVE}"
} > "${RUN_MANIFEST}"

run_stage() {
    local stage_name="$1"
    local source_dir=""
    local stage_fastq_dir=""
    local copy_time_log=""
    local copy_stdout_log=""

    case "${stage_name}" in
        downsampled)
            source_dir="${DOWNSAMPLED_FASTQ_DIR}"
            stage_fastq_dir="${DOWNSAMPLED_FASTQ_DIR}"
            ;;
        storage)
            source_dir="${STORAGE_FASTQ_DIR}"
            stage_fastq_dir="${STORAGE_FASTQ_DIR}"
            ;;
        pikachu)
            source_dir="${STORAGE_FASTQ_DIR}"
            stage_fastq_dir="${PIKACHU_FASTQ_DIR}"
            ;;
        *)
            die "Unsupported stage: ${stage_name}"
            ;;
    esac

    require_dir "${source_dir}"

    local stage_dir="${OUT_ROOT}/${stage_name}"
    local inputs_dir="${stage_dir}/inputs"
    local integrated_dir="${stage_dir}/integrated"
    local external_dir="${stage_dir}/external"
    local comparison_dir="${stage_dir}/comparison"
    local stage_manifest="${stage_dir}/BENCHMARK_MANIFEST.txt"
    local stage_summary="${stage_dir}/BENCHMARK_SUMMARY.txt"
    local stage_compare_tsv="${comparison_dir}/comparison_metrics.tsv"
    mkdir -p "${inputs_dir}" "${integrated_dir}" "${external_dir}" "${comparison_dir}"

    if [[ "${stage_name}" == "pikachu" && "${SKIP_COPY}" -eq 0 ]]; then
        mkdir -p "${stage_fastq_dir}"
        mapfile -t src_pair < <(find_single_pair "${source_dir}")
        local src_r1="${src_pair[0]}"
        local src_r2="${src_pair[1]}"
        local dest_r1="${stage_fastq_dir}/$(basename "${src_r1}")"
        local dest_r2="${stage_fastq_dir}/$(basename "${src_r2}")"
        copy_time_log="${stage_dir}/copy.time_v.log"
        copy_stdout_log="${stage_dir}/copy.log"
        if [[ "${FORCE_COPY}" -eq 1 || ! -f "${dest_r1}" || ! -f "${dest_r2}" ]]; then
            log "Stage ${stage_name}: copying FASTQs to ${stage_fastq_dir}"
            run_timed_cmd "${copy_time_log}" "${copy_stdout_log}" cp -fp "${src_r1}" "${src_r2}" "${stage_fastq_dir}/"
        else
            log "Stage ${stage_name}: reusing existing copied FASTQs in ${stage_fastq_dir}"
        fi
    fi

    require_dir "${stage_fastq_dir}"
    mapfile -t pair < <(find_single_pair "${stage_fastq_dir}")
    local raw_r1="${pair[0]}"
    local raw_r2="${pair[1]}"
    local sample_name
    sample_name="$(basename "${raw_r1}" _R1_001.fastq.gz)"

    {
        echo "stage=${stage_name}"
        echo "sample=${sample_name}"
        echo "fastq_dir=${stage_fastq_dir}"
        echo "raw_r1=${raw_r1}"
        echo "raw_r2=${raw_r2}"
        echo "integrated_dir=${integrated_dir}"
        echo "external_dir=${external_dir}"
        echo "comparison_dir=${comparison_dir}"
        echo "salmon_libtype=${SALMON_LIBTYPE}"
        echo "yremove=${YREMOVE}"
        echo "parity_qc=${PARITY_QC}"
    } > "${stage_manifest}"

    local integrated_quant="${integrated_dir}/quant.sf"
    local integrated_gene_quant="${integrated_dir}/quant.genes.sf"
    local integrated_transcriptome_bam="${integrated_dir}/Aligned.toTranscriptome.out.bam"
    local integrated_salmon_dir="${integrated_dir}/salmon_qc"
    local integrated_trim_qc_prefix="${integrated_dir}/${TRIM_QC_BASENAME}"
    local external_trim_qc_prefix="${external_dir}/${TRIM_QC_BASENAME}"
    local external_trim_qc_r1_prefix="${external_trim_qc_prefix}_R1"
    local external_trim_qc_r2_prefix="${external_trim_qc_prefix}_R2"

    local external_raw_r1="${external_dir}/${sample_name}.raw_R1.fastq"
    local external_raw_r2="${external_dir}/${sample_name}.raw_R2.fastq"
    local external_trim_r1="${external_dir}/${sample_name}.trimmed_R1.fastq"
    local external_trim_r2="${external_dir}/${sample_name}.trimmed_R2.fastq"
    if [[ "${EXTERNAL_TRIMMER}" == "trim_galore" ]]; then
        external_trim_r1="${external_dir}/${sample_name}.trimmed_val_1.fq.gz"
        external_trim_r2="${external_dir}/${sample_name}.trimmed_val_2.fq.gz"
    fi
    local external_genome_bam="${external_dir}/Aligned.sortedByCoord.out.bam"
    local external_transcriptome_bam="${external_dir}/Aligned.toTranscriptome.out.bam"
    local external_y_bam="${external_dir}/Aligned.sortedByCoord.out_Y.bam"
    local external_noy_bam="${external_dir}/Aligned.sortedByCoord.out_noY.bam"
    local external_y_transcriptome_bam="${external_dir}/Aligned.toTranscriptome.out_Y.bam"
    local external_noy_transcriptome_bam="${external_dir}/Aligned.toTranscriptome.out_noY.bam"
    local external_salmon_bam="${external_transcriptome_bam}"
    local external_salmon_dir="${external_dir}/salmon_qc"

    if [[ "${SKIP_INTEGRATED}" -eq 0 ]]; then
        local integrated_quant_modes=(TranscriptVB)
        if [[ "${PARITY_QC}" -eq 1 ]]; then
            integrated_quant_modes=(TranscriptomeSAM TranscriptVB)
        fi
        local integrated_cmd=(
            "${STAR_BIN}"
            --runMode alignReads
            --runThreadN "${THREADS}"
            --genomeDir "${STAR_INDEX}"
            --readFilesIn "${raw_r1}" "${raw_r2}"
            --readFilesCommand zcat
            --trimCutadapt Yes
            --trimCutadaptQuality "${QUALITY_CUTOFF}"
            --trimCutadaptMinLength "${MIN_LENGTH}"
            --trimCutadaptAdapter "${ADAPTER_R1}" "${ADAPTER_R2}"
            --outSAMtype BAM SortedByCoordinate
            --outBAMsortMethod samtools
            --keepBAM yes
            --quantMode "${integrated_quant_modes[@]}"
            --transcriptomeFasta "${TRANSCRIPTOME}"
            --quantVBgcBias 1
            --outFileNamePrefix "${integrated_dir}/"
        )
        if [[ "${YREMOVE}" -eq 1 ]]; then
            integrated_cmd+=(
                --emitNoYBAM yes
                --emitYNoYFastq yes
                --emitYNoYFastqCompression gz
            )
        fi
        if [[ "${TRIM_QC}" -eq 1 ]]; then
            integrated_cmd+=(
                --trimQcReport "${integrated_trim_qc_prefix}"
                --trimQcMaxReads "${TRIM_QC_MAX_READS}"
            )
        fi
        write_cmd_script "${integrated_dir}/run_integrated_star.sh" "${integrated_cmd[@]}"
        log "Stage ${stage_name}: integrated STAR-suite arm (yremove=${YREMOVE})"
        run_timed_cmd "${integrated_dir}/star.time_v.log" "${integrated_dir}/star.log" "${integrated_cmd[@]}"
        if [[ "${DRY_RUN}" -eq 0 ]]; then
            require_file "${integrated_quant}"
            if [[ "${PARITY_QC}" -eq 1 ]]; then
                require_file "${integrated_transcriptome_bam}"
            fi
            if [[ "${YREMOVE}" -eq 1 ]]; then
                require_file "${integrated_dir}/Aligned.sortedByCoord.out_Y.bam"
            fi
            if [[ "${TRIM_QC}" -eq 1 ]]; then
                require_file "${integrated_trim_qc_prefix}.trim_qc.json"
                require_file "${integrated_trim_qc_prefix}.trim_qc.html"
            fi
        fi

        if [[ "${PARITY_QC}" -eq 1 ]]; then
            local integrated_salmon_cmd=(
                "${SALMON_BIN}" quant
                -t "${TRANSCRIPTOME}"
                -l "${SALMON_LIBTYPE}"
                -a "${integrated_transcriptome_bam}"
                -g "${TX2GENE}"
                --gcBias
                -p "${THREADS}"
                -o "${integrated_salmon_dir}"
            )
            write_cmd_script "${integrated_dir}/run_integrated_salmon_qc.sh" "${integrated_salmon_cmd[@]}"
            log "Stage ${stage_name}: Salmon QC on integrated transcriptome BAM"
            run_timed_cmd "${integrated_dir}/salmon.time_v.log" "${integrated_dir}/salmon.log" "${integrated_salmon_cmd[@]}"
            if [[ "${DRY_RUN}" -eq 0 ]]; then
                require_file "${integrated_salmon_dir}/quant.sf"
            fi
        else
            log "Stage ${stage_name}: integrated Salmon QC disabled (production benchmark mode)"
        fi
    else
        log "Stage ${stage_name}: skipping integrated arm"
    fi

    if [[ "${SKIP_EXTERNAL}" -eq 0 ]]; then
        if [[ "${EXTERNAL_TRIMMER}" == "trim_galore" ]]; then
            local trim_cmd=(
                "${TRIM_GALORE_BIN}"
                --paired
                --quality "${QUALITY_CUTOFF}"
                --length "${MIN_LENGTH}"
                --adapter "${ADAPTER_R1}"
                --adapter2 "${ADAPTER_R2}"
                --cores "${TRIMGALORE_CORES}"
                --basename "${sample_name}.trimmed"
                --output_dir "${external_dir}"
            )
            if [[ "${TRIM_QC}" -eq 1 ]]; then
                trim_cmd+=(--fastqc --fastqc_args "--threads ${THREADS}")
            fi
            trim_cmd+=("${raw_r1}" "${raw_r2}")
            write_cmd_script "${external_dir}/run_trimgalore.sh" "${trim_cmd[@]}"
            log "Stage ${stage_name}: external Trim Galore adapter-removal arm (cores=${TRIMGALORE_CORES}, fastqc=${TRIM_QC})"
            run_timed_cmd "${external_dir}/trimgalore.time_v.log" "${external_dir}/trimgalore.log" "${trim_cmd[@]}"
            if [[ "${DRY_RUN}" -eq 0 ]]; then
                require_file "${external_trim_r1}"
                require_file "${external_trim_r2}"
            fi
        else
            local external_decompress_cmd=(
                bash -lc
                "gzip -dc $(printf '%q' "${raw_r1}") > $(printf '%q' "${external_raw_r1}") && gzip -dc $(printf '%q' "${raw_r2}") > $(printf '%q' "${external_raw_r2}")"
            )
            write_cmd_script "${external_dir}/run_decompress_raw_fastq.sh" "${external_decompress_cmd[@]}"
            log "Stage ${stage_name}: external FASTQ decompression for trimvalidate"
            run_timed_cmd "${external_dir}/decompress.time_v.log" "${external_dir}/decompress.log" "${external_decompress_cmd[@]}"
            if [[ "${DRY_RUN}" -eq 0 ]]; then
                require_file "${external_raw_r1}"
                require_file "${external_raw_r2}"
            fi

            if [[ "${TRIM_QC}" -eq 1 ]]; then
                local external_trim_qc_cmd=(
                    bash -lc
                    "$(printf '%q' "${TRIM_QC_FASTQ_BIN}") --input $(printf '%q' "${external_raw_r1}") --report $(printf '%q' "${external_trim_qc_r1_prefix}") --stage external_fastq_r1 --mate-count 2 --mate-index 1 --max-reads $(printf '%q' "${TRIM_QC_MAX_READS}") && $(printf '%q' "${TRIM_QC_FASTQ_BIN}") --input $(printf '%q' "${external_raw_r2}") --report $(printf '%q' "${external_trim_qc_r2_prefix}") --stage external_fastq_r2 --mate-count 2 --mate-index 2 --max-reads $(printf '%q' "${TRIM_QC_MAX_READS}")"
                )
                write_cmd_script "${external_dir}/run_external_trim_qc.sh" "${external_trim_qc_cmd[@]}"
                log "Stage ${stage_name}: external FastQC-like trim-QC reports"
                run_timed_cmd "${external_dir}/trim_qc.time_v.log" "${external_dir}/trim_qc.log" "${external_trim_qc_cmd[@]}"
                if [[ "${DRY_RUN}" -eq 0 ]]; then
                    require_file "${external_trim_qc_r1_prefix}.trim_qc.json"
                    require_file "${external_trim_qc_r1_prefix}.trim_qc.html"
                    require_file "${external_trim_qc_r2_prefix}.trim_qc.json"
                    require_file "${external_trim_qc_r2_prefix}.trim_qc.html"
                fi
            else
                log "Stage ${stage_name}: skipping external trim-QC reports"
            fi

            local trim_cmd=(
                "${TRIMVALIDATE_BIN}"
                -1 "${external_raw_r1}"
                -2 "${external_raw_r2}"
                -o1 "${external_trim_r1}"
                -o2 "${external_trim_r2}"
                --quality "${QUALITY_CUTOFF}"
                --length "${MIN_LENGTH}"
                --adapter-r1 "${ADAPTER_R1}"
                --adapter-r2 "${ADAPTER_R2}"
            )
            write_cmd_script "${external_dir}/run_trimvalidate.sh" "${trim_cmd[@]}"
            log "Stage ${stage_name}: external trimvalidate adapter-removal arm"
            run_timed_cmd "${external_dir}/trimvalidate.time_v.log" "${external_dir}/trimvalidate.log" "${trim_cmd[@]}"
            if [[ "${DRY_RUN}" -eq 0 ]]; then
                require_file "${external_trim_r1}"
                require_file "${external_trim_r2}"
            fi
        fi

        local external_star_cmd=(
            "${STAR_BIN}"
            --runMode alignReads
            --runThreadN "${THREADS}"
            --genomeDir "${STAR_INDEX}"
            --readFilesIn "${external_trim_r1}" "${external_trim_r2}"
            --outSAMtype BAM SortedByCoordinate
            --outBAMsortMethod samtools
            --keepBAM yes
            --quantMode TranscriptomeSAM
            --outFileNamePrefix "${external_dir}/"
        )
        if [[ "${external_trim_r1}" == *.gz ]]; then
            external_star_cmd+=(--readFilesCommand zcat)
        fi
        if [[ "${YREMOVE}" -eq 1 && "${EXTERNAL_YREMOVE}" != "awk" ]]; then
            external_star_cmd+=(--emitNoYBAM yes)
        fi
        write_cmd_script "${external_dir}/run_external_star.sh" "${external_star_cmd[@]}"
        log "Stage ${stage_name}: external STAR alignment arm (yremove=${YREMOVE})"
        run_timed_cmd "${external_dir}/star.time_v.log" "${external_dir}/star.log" "${external_star_cmd[@]}"
        if [[ "${DRY_RUN}" -eq 0 ]]; then
            require_file "${external_transcriptome_bam}"
            if [[ "${YREMOVE}" -eq 1 ]]; then
                require_file "${external_genome_bam}"
                if [[ "${EXTERNAL_YREMOVE}" != "awk" ]]; then
                    require_file "${external_y_bam}"
                fi
            fi
        fi

        if [[ "${YREMOVE}" -eq 1 ]]; then
            mkdir -p "${external_dir}/y_fastq_split"
            if [[ "${EXTERNAL_YREMOVE}" == "awk" ]]; then
                local awk_y_split_script="${external_dir}/awk_y_remove.sh"
                local awk_y_prefix="${sample_name}.trimmed"
                write_awk_y_remove_script "${awk_y_split_script}"
                local awk_y_cmd=(
                    "${awk_y_split_script}"
                    "${SAMTOOLS_BIN:-$(command -v samtools)}"
                    "${AWK_BIN:-$(command -v awk)}"
                    "${external_genome_bam}"
                    "${external_transcriptome_bam}"
                    "${external_trim_r1}"
                    "${external_trim_r2}"
                    "${external_dir}"
                    "${external_dir}/y_fastq_split"
                    "${awk_y_prefix}"
                )
                write_cmd_script "${external_dir}/run_awk_y_remove.sh" "${awk_y_cmd[@]}"
                log "Stage ${stage_name}: external awk/samtools Y-removal for BAM and FASTQ artifacts"
                run_timed_cmd "${external_dir}/awk_y_split.time_v.log" "${external_dir}/awk_y_split.log" "${awk_y_cmd[@]}"
                external_salmon_bam="${external_noy_transcriptome_bam}"
                if [[ "${DRY_RUN}" -eq 0 ]]; then
                    require_file "${external_y_bam}"
                    require_file "${external_noy_bam}"
                    require_file "${external_y_transcriptome_bam}"
                    require_file "${external_noy_transcriptome_bam}"
                    require_file "${external_salmon_bam}"
                fi
            else
                local remove_y_cmd=(
                    "${REMOVE_Y_READS_BIN}"
                    -y "${external_y_bam}"
                    --threads "${THREADS}"
                    -o "${external_dir}/y_fastq_split"
                    "${external_trim_r1}"
                    "${external_trim_r2}"
                )
                write_cmd_script "${external_dir}/run_remove_y_reads.sh" "${remove_y_cmd[@]}"
                log "Stage ${stage_name}: external remove_y_reads FASTQ split"
                run_timed_cmd "${external_dir}/remove_y_reads.time_v.log" "${external_dir}/remove_y_reads.log" "${remove_y_cmd[@]}"
            fi
        else
            log "Stage ${stage_name}: skipping external Y/no-Y FASTQ split (yremove disabled)"
        fi

        local external_salmon_cmd=(
            "${SALMON_BIN}" quant
            -t "${TRANSCRIPTOME}"
            -l "${SALMON_LIBTYPE}"
            -a "${external_salmon_bam}"
            -g "${TX2GENE}"
            --gcBias
            -p "${THREADS}"
            -o "${external_salmon_dir}"
        )
        write_cmd_script "${external_dir}/run_external_salmon_qc.sh" "${external_salmon_cmd[@]}"
        log "Stage ${stage_name}: external Salmon QC"
        run_timed_cmd "${external_dir}/salmon.time_v.log" "${external_dir}/salmon.log" "${external_salmon_cmd[@]}"
        if [[ "${DRY_RUN}" -eq 0 ]]; then
            require_file "${external_salmon_dir}/quant.sf"
        fi
    else
        log "Stage ${stage_name}: skipping external arm"
    fi

    if [[ "${SKIP_COMPARE}" -eq 0 && "${PARITY_QC}" -eq 1 && "${DRY_RUN}" -eq 0 ]]; then
        log "Stage ${stage_name}: collecting comparison metrics"
        append_quant_compare_row "${stage_name}" "integrated_transcriptvb_vs_integrated_salmon" \
            "${integrated_salmon_dir}/quant.sf" "${integrated_quant}" "${stage_compare_tsv}"
        append_quant_compare_row "${stage_name}" "integrated_salmon_vs_external_salmon" \
            "${integrated_salmon_dir}/quant.sf" "${external_salmon_dir}/quant.sf" "${stage_compare_tsv}"
        append_quant_compare_row "${stage_name}" "integrated_transcriptvb_vs_external_salmon" \
            "${external_salmon_dir}/quant.sf" "${integrated_quant}" "${stage_compare_tsv}"

        if [[ -f "${integrated_gene_quant}" && -f "${integrated_salmon_dir}/quant.genes.sf" ]]; then
            append_quant_compare_row "${stage_name}" "integrated_gene_transcriptvb_vs_integrated_gene_salmon" \
                "${integrated_salmon_dir}/quant.genes.sf" "${integrated_gene_quant}" "${stage_compare_tsv}"
        fi
        if [[ -f "${integrated_gene_quant}" && -f "${external_salmon_dir}/quant.genes.sf" ]]; then
            append_quant_compare_row "${stage_name}" "integrated_gene_transcriptvb_vs_external_gene_salmon" \
                "${external_salmon_dir}/quant.genes.sf" "${integrated_gene_quant}" "${stage_compare_tsv}"
        fi
    elif [[ "${SKIP_COMPARE}" -eq 0 && "${PARITY_QC}" -eq 1 ]]; then
        log "Stage ${stage_name}: DRY_RUN active, skipping comparison execution"
    elif [[ "${SKIP_COMPARE}" -eq 0 ]]; then
        log "Stage ${stage_name}: skipping comparison stage (production benchmark mode; use --parity-qc for Salmon comparison)"
    else
        log "Stage ${stage_name}: skipping comparison stage"
    fi

    local integrated_y_r1=""
    local integrated_y_r2=""
    local integrated_noy_r1=""
    local integrated_noy_r2=""
    local external_y_r1=""
    local external_y_r2=""
    local external_noy_r1=""
    local external_noy_r2=""
    if [[ "${YREMOVE}" -eq 1 ]]; then
        integrated_y_r1="$(find "${integrated_dir}" -maxdepth 2 -type f -name '*_Y_R1*.fastq*' | head -n1 || true)"
        integrated_y_r2="$(find "${integrated_dir}" -maxdepth 2 -type f -name '*_Y_R2*.fastq*' | head -n1 || true)"
        integrated_noy_r1="$(find "${integrated_dir}" -maxdepth 2 -type f -name '*_noY_R1*.fastq*' | head -n1 || true)"
        integrated_noy_r2="$(find "${integrated_dir}" -maxdepth 2 -type f -name '*_noY_R2*.fastq*' | head -n1 || true)"
        external_y_r1="$(find "${external_dir}/y_fastq_split" -maxdepth 1 -type f -name '*R1*_Y.fastq*' | head -n1 || true)"
        external_y_r2="$(find "${external_dir}/y_fastq_split" -maxdepth 1 -type f -name '*R2*_Y.fastq*' | head -n1 || true)"
        external_noy_r1="$(find "${external_dir}/y_fastq_split" -maxdepth 1 -type f -name '*R1*_noY.fastq*' | head -n1 || true)"
        external_noy_r2="$(find "${external_dir}/y_fastq_split" -maxdepth 1 -type f -name '*R2*_noY.fastq*' | head -n1 || true)"
    fi

    local copy_wall="NA"
    local integrated_star_wall="NA"
    local integrated_salmon_wall="NA"
    local external_decompress_wall="NA"
    local external_trim_qc_wall="NA"
    local external_trim_wall="NA"
    local external_star_wall="NA"
    local external_ysplit_wall="NA"
    local external_salmon_wall="NA"
    copy_wall="$(extract_time_metric "${copy_time_log}" "Elapsed (wall clock) time")"
    integrated_star_wall="$(extract_time_metric "${integrated_dir}/star.time_v.log" "Elapsed (wall clock) time")"
    if [[ "${PARITY_QC}" -eq 1 ]]; then
        integrated_salmon_wall="$(extract_time_metric "${integrated_dir}/salmon.time_v.log" "Elapsed (wall clock) time")"
    else
        integrated_salmon_wall="(disabled)"
    fi
    if [[ "${EXTERNAL_TRIMMER}" == "trim_galore" ]]; then
        external_decompress_wall="(folded into Trim Galore)"
        external_trim_qc_wall="(FastQC folded into Trim Galore)"
        external_trim_wall="$(extract_time_metric "${external_dir}/trimgalore.time_v.log" "Elapsed (wall clock) time")"
    else
        external_decompress_wall="$(extract_time_metric "${external_dir}/decompress.time_v.log" "Elapsed (wall clock) time")"
        if [[ "${TRIM_QC}" -eq 1 ]]; then
        external_trim_qc_wall="$(extract_time_metric "${external_dir}/trim_qc.time_v.log" "Elapsed (wall clock) time")"
        else
        external_trim_qc_wall="(disabled)"
        fi
        external_trim_wall="$(extract_time_metric "${external_dir}/trimvalidate.time_v.log" "Elapsed (wall clock) time")"
    fi
    external_star_wall="$(extract_time_metric "${external_dir}/star.time_v.log" "Elapsed (wall clock) time")"
    if [[ "${YREMOVE}" -eq 1 ]]; then
        if [[ "${EXTERNAL_YREMOVE}" == "awk" ]]; then
            external_ysplit_wall="$(extract_time_metric "${external_dir}/awk_y_split.time_v.log" "Elapsed (wall clock) time")"
        else
            external_ysplit_wall="$(extract_time_metric "${external_dir}/remove_y_reads.time_v.log" "Elapsed (wall clock) time")"
        fi
    else
        external_ysplit_wall="(skipped)"
    fi
    external_salmon_wall="$(extract_time_metric "${external_dir}/salmon.time_v.log" "Elapsed (wall clock) time")"

    {
        echo "=============================================================="
        echo "PE Bulk Feature Benchmark Summary"
        echo "=============================================================="
        echo "Stage:                  ${stage_name}"
        echo "Y-removal:              ${YREMOVE}"
        echo "Parity QC:              ${PARITY_QC}"
        echo "Trim QC:                ${TRIM_QC}"
        echo "Trim QC max reads:      ${TRIM_QC_MAX_READS}"
        echo "External trimmer:       ${EXTERNAL_TRIMMER}"
        echo "Trim Galore cores:      ${TRIMGALORE_CORES}"
        echo "External Y-removal:     ${EXTERNAL_YREMOVE}"
        echo "Sample:                 ${sample_name}"
        echo "FASTQ dir:              ${stage_fastq_dir}"
        echo "STAR index:             ${STAR_INDEX}"
        echo "Transcriptome:          ${TRANSCRIPTOME}"
        echo "GTF:                    ${GTF}"
        echo "tx2gene:                ${TX2GENE}"
        echo "Salmon libtype:         ${SALMON_LIBTYPE}"
        echo
        echo "Timings:"
        echo "  Copy to pikachu:      ${copy_wall}"
        echo "  Integrated STAR:      ${integrated_star_wall} (production total)"
        echo "  Integrated Salmon QC: ${integrated_salmon_wall}"
        echo "  External decompress:  ${external_decompress_wall}"
        echo "  External trim QC:     ${external_trim_qc_wall}"
        echo "  External trimmer:     ${external_trim_wall}"
        echo "  External STAR:        ${external_star_wall}"
        echo "  External remove_y:    ${external_ysplit_wall}"
        echo "  External Salmon:      ${external_salmon_wall}"
        echo
        echo "Integrated outputs:"
        echo "  quant.sf:             ${integrated_quant}"
        echo "  quant.genes.sf:       ${integrated_gene_quant}"
        if [[ "${PARITY_QC}" -eq 1 ]]; then
            echo "  transcriptome BAM:    ${integrated_transcriptome_bam}"
            echo "  Salmon QC dir:        ${integrated_salmon_dir}"
        else
            echo "  transcriptome BAM:    NA (production mode)"
            echo "  Salmon QC dir:        NA (production mode)"
        fi
        echo "  Integrated Y R1:      ${integrated_y_r1:-NA}"
        echo "  Integrated Y R2:      ${integrated_y_r2:-NA}"
        echo "  Integrated noY R1:    ${integrated_noy_r1:-NA}"
        echo "  Integrated noY R2:    ${integrated_noy_r2:-NA}"
        if [[ "${TRIM_QC}" -eq 1 ]]; then
            echo "  Trim QC JSON:         ${integrated_trim_qc_prefix}.trim_qc.json"
            echo "  Trim QC HTML:         ${integrated_trim_qc_prefix}.trim_qc.html"
        else
            echo "  Trim QC JSON:         NA (disabled)"
            echo "  Trim QC HTML:         NA (disabled)"
        fi
        echo
        echo "External outputs:"
        if [[ "${SKIP_EXTERNAL}" -eq 1 ]]; then
            echo "  NA (skipped)"
        else
            if [[ "${EXTERNAL_TRIMMER}" == "trimvalidate" ]]; then
                echo "  raw R1 (plain):       ${external_raw_r1}"
                echo "  raw R2 (plain):       ${external_raw_r2}"
            else
                echo "  raw R1 source:        ${raw_r1}"
                echo "  raw R2 source:        ${raw_r2}"
            fi
            echo "  trimmed R1:           ${external_trim_r1}"
            echo "  trimmed R2:           ${external_trim_r2}"
            echo "  genome BAM:           ${external_genome_bam}"
            echo "  transcriptome BAM:    ${external_transcriptome_bam}"
            if [[ "${YREMOVE}" -eq 1 ]]; then
                echo "  noY genome BAM:       ${external_noy_bam}"
                echo "  noY transcriptome BAM:${external_noy_transcriptome_bam}"
                echo "  Salmon input BAM:     ${external_salmon_bam}"
            fi
            echo "  Salmon QC dir:        ${external_salmon_dir}"
            echo "  External Y R1:        ${external_y_r1:-NA}"
            echo "  External Y R2:        ${external_y_r2:-NA}"
            echo "  External noY R1:      ${external_noy_r1:-NA}"
            echo "  External noY R2:      ${external_noy_r2:-NA}"
            if [[ "${EXTERNAL_TRIMMER}" == "trim_galore" && "${TRIM_QC}" -eq 1 ]]; then
                echo "  FastQC outputs:       ${external_dir}/*fastqc.{html,zip}"
            elif [[ "${TRIM_QC}" -eq 1 ]]; then
                echo "  Trim QC R1 JSON:      ${external_trim_qc_r1_prefix}.trim_qc.json"
                echo "  Trim QC R1 HTML:      ${external_trim_qc_r1_prefix}.trim_qc.html"
                echo "  Trim QC R2 JSON:      ${external_trim_qc_r2_prefix}.trim_qc.json"
                echo "  Trim QC R2 HTML:      ${external_trim_qc_r2_prefix}.trim_qc.html"
            else
                echo "  Trim QC reports:      NA (disabled)"
            fi
        fi
        echo
        if [[ "${DRY_RUN}" -eq 0 && "${YREMOVE}" -eq 1 ]]; then
            echo "FASTQ counts (Y-split):"
            echo "  Integrated Y R1 reads:    $(count_fastq_reads "${integrated_y_r1:-/dev/null}")"
            echo "  Integrated noY R1 reads:  $(count_fastq_reads "${integrated_noy_r1:-/dev/null}")"
            echo "  External Y R1 reads:      $(count_fastq_reads "${external_y_r1:-/dev/null}")"
            echo "  External noY R1 reads:    $(count_fastq_reads "${external_noy_r1:-/dev/null}")"
            echo
        fi
        echo "Comparison TSV:         ${stage_compare_tsv}"
        echo "Manifest:               ${stage_manifest}"
        echo "Summary:                ${stage_summary}"
        echo "=============================================================="
    } | tee "${stage_summary}"
}

IFS=',' read -r -a stage_list <<< "${STAGES}"
for stage_name in "${stage_list[@]}"; do
    stage_name="${stage_name//[[:space:]]/}"
    [[ -n "${stage_name}" ]] || continue
    run_stage "${stage_name}"
done

log "Run manifest: ${RUN_MANIFEST}"
