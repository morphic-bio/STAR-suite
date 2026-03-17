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
  --yremove                      Enable Y-chromosome removal (default).
  --no-yremove                   Disable Y-chromosome removal on both arms.
  --skip-copy                    Do not copy /storage PE inputs into --pikachu-fastq-dir.
  --force-copy                   Re-copy into --pikachu-fastq-dir even if files already exist.
  --skip-integrated              Skip integrated STAR-suite arm.
  --skip-external                Skip external stepwise arm.
  --skip-compare                 Skip quant/Y-split comparison stage.
  --dry-run                      Print resolved commands only.
  -h, --help                     Show this help.

Modes:
  --yremove (default):
    Integrated arm:
      raw FASTQ -> STAR (trimCutadapt + emitNoYBAM + emitYNoYFastq + TranscriptomeSAM + TranscriptVB)
                -> Salmon alignment-mode QC on integrated transcriptome BAM

    External arm:
      raw FASTQ -> trimvalidate -> STAR (TranscriptomeSAM + emitNoYBAM)
                -> remove_y_reads on trimmed FASTQs
                -> Salmon alignment-mode QC on transcriptome BAM

  --no-yremove:
    Integrated arm:
      raw FASTQ -> STAR (trimCutadapt + TranscriptomeSAM + TranscriptVB)
                -> Salmon alignment-mode QC on integrated transcriptome BAM

    External arm:
      raw FASTQ -> trimvalidate -> STAR (TranscriptomeSAM)
                -> Salmon alignment-mode QC on transcriptome BAM

Note:
  Salmon auto libtype detection (-l A) mis-detects the current PE benchmark sample as ISR.
  This script therefore defaults Salmon QC to explicit IU for reproducible comparison.
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
    grep -F "${pattern}" "${log_file}" | head -n1 | sed -E 's/^.*:[[:space:]]*//'
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
        --yremove) YREMOVE=1; shift ;;
        --no-yremove) YREMOVE=0; shift ;;
        --skip-copy) SKIP_COPY=1; shift ;;
        --force-copy) FORCE_COPY=1; shift ;;
        --skip-integrated) SKIP_INTEGRATED=1; shift ;;
        --skip-external) SKIP_EXTERNAL=1; shift ;;
        --skip-compare) SKIP_COMPARE=1; shift ;;
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
    } > "${stage_manifest}"

    local integrated_quant="${integrated_dir}/quant.sf"
    local integrated_gene_quant="${integrated_dir}/quant.genes.sf"
    local integrated_transcriptome_bam="${integrated_dir}/Aligned.toTranscriptome.out.bam"
    local integrated_salmon_dir="${integrated_dir}/salmon_qc"

    local external_raw_r1="${external_dir}/${sample_name}.raw_R1.fastq"
    local external_raw_r2="${external_dir}/${sample_name}.raw_R2.fastq"
    local external_trim_r1="${external_dir}/${sample_name}.trimmed_R1.fastq"
    local external_trim_r2="${external_dir}/${sample_name}.trimmed_R2.fastq"
    local external_transcriptome_bam="${external_dir}/Aligned.toTranscriptome.out.bam"
    local external_y_bam="${external_dir}/Aligned.sortedByCoord.out_Y.bam"
    local external_salmon_dir="${external_dir}/salmon_qc"

    if [[ "${SKIP_INTEGRATED}" -eq 0 ]]; then
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
            --quantMode TranscriptomeSAM TranscriptVB
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
        write_cmd_script "${integrated_dir}/run_integrated_star.sh" "${integrated_cmd[@]}"
        log "Stage ${stage_name}: integrated STAR-suite arm (yremove=${YREMOVE})"
        run_timed_cmd "${integrated_dir}/star.time_v.log" "${integrated_dir}/star.log" "${integrated_cmd[@]}"
        if [[ "${DRY_RUN}" -eq 0 ]]; then
            require_file "${integrated_quant}"
            require_file "${integrated_transcriptome_bam}"
            if [[ "${YREMOVE}" -eq 1 ]]; then
                require_file "${integrated_dir}/Aligned.sortedByCoord.out_Y.bam"
            fi
        fi

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
        log "Stage ${stage_name}: skipping integrated arm"
    fi

    if [[ "${SKIP_EXTERNAL}" -eq 0 ]]; then
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
        if [[ "${YREMOVE}" -eq 1 ]]; then
            external_star_cmd+=(--emitNoYBAM yes)
        fi
        write_cmd_script "${external_dir}/run_external_star.sh" "${external_star_cmd[@]}"
        log "Stage ${stage_name}: external STAR alignment arm (yremove=${YREMOVE})"
        run_timed_cmd "${external_dir}/star.time_v.log" "${external_dir}/star.log" "${external_star_cmd[@]}"
        if [[ "${DRY_RUN}" -eq 0 ]]; then
            require_file "${external_transcriptome_bam}"
            if [[ "${YREMOVE}" -eq 1 ]]; then
                require_file "${external_y_bam}"
            fi
        fi

        if [[ "${YREMOVE}" -eq 1 ]]; then
            mkdir -p "${external_dir}/y_fastq_split"
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
        else
            log "Stage ${stage_name}: skipping external remove_y_reads (yremove disabled)"
        fi

        local external_salmon_cmd=(
            "${SALMON_BIN}" quant
            -t "${TRANSCRIPTOME}"
            -l "${SALMON_LIBTYPE}"
            -a "${external_transcriptome_bam}"
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

    if [[ "${SKIP_COMPARE}" -eq 0 && "${DRY_RUN}" -eq 0 ]]; then
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
    elif [[ "${SKIP_COMPARE}" -eq 0 ]]; then
        log "Stage ${stage_name}: DRY_RUN active, skipping comparison execution"
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
    local external_trim_wall="NA"
    local external_star_wall="NA"
    local external_ysplit_wall="NA"
    local external_salmon_wall="NA"
    copy_wall="$(extract_time_metric "${copy_time_log}" "Elapsed (wall clock) time")"
    integrated_star_wall="$(extract_time_metric "${integrated_dir}/star.time_v.log" "Elapsed (wall clock) time")"
    integrated_salmon_wall="$(extract_time_metric "${integrated_dir}/salmon.time_v.log" "Elapsed (wall clock) time")"
    external_decompress_wall="$(extract_time_metric "${external_dir}/decompress.time_v.log" "Elapsed (wall clock) time")"
    external_trim_wall="$(extract_time_metric "${external_dir}/trimvalidate.time_v.log" "Elapsed (wall clock) time")"
    external_star_wall="$(extract_time_metric "${external_dir}/star.time_v.log" "Elapsed (wall clock) time")"
    if [[ "${YREMOVE}" -eq 1 ]]; then
        external_ysplit_wall="$(extract_time_metric "${external_dir}/remove_y_reads.time_v.log" "Elapsed (wall clock) time")"
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
        echo "  Integrated STAR:      ${integrated_star_wall}"
        echo "  Integrated Salmon:    ${integrated_salmon_wall}"
        echo "  External decompress:  ${external_decompress_wall}"
        echo "  External trimvalidate:${external_trim_wall}"
        echo "  External STAR:        ${external_star_wall}"
        echo "  External remove_y:    ${external_ysplit_wall}"
        echo "  External Salmon:      ${external_salmon_wall}"
        echo
        echo "Integrated outputs:"
        echo "  quant.sf:             ${integrated_quant}"
        echo "  quant.genes.sf:       ${integrated_gene_quant}"
        echo "  transcriptome BAM:    ${integrated_transcriptome_bam}"
        echo "  Salmon QC dir:        ${integrated_salmon_dir}"
        echo "  Integrated Y R1:      ${integrated_y_r1:-NA}"
        echo "  Integrated Y R2:      ${integrated_y_r2:-NA}"
        echo "  Integrated noY R1:    ${integrated_noy_r1:-NA}"
        echo "  Integrated noY R2:    ${integrated_noy_r2:-NA}"
        echo
        echo "External outputs:"
        echo "  raw R1 (plain):       ${external_raw_r1}"
        echo "  raw R2 (plain):       ${external_raw_r2}"
        echo "  trimmed R1:           ${external_trim_r1}"
        echo "  trimmed R2:           ${external_trim_r2}"
        echo "  transcriptome BAM:    ${external_transcriptome_bam}"
        echo "  Salmon QC dir:        ${external_salmon_dir}"
        echo "  External Y R1:        ${external_y_r1:-NA}"
        echo "  External Y R2:        ${external_y_r2:-NA}"
        echo "  External noY R1:      ${external_noy_r1:-NA}"
        echo "  External noY R2:      ${external_noy_r2:-NA}"
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
