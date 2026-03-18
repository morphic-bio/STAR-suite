#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

OUTDIR=""
DRY_RUN=0
BUILD_MISSING=1

STAR_BIN="${REPO_ROOT}/core/legacy/source/STAR"
SALMON_BIN_DEFAULT="/mnt/pikachu/salmon/build/src/salmon"
SALMON_BIN=""
FASTQ_DUMP_BIN="${FASTQ_DUMP_BIN:-}"
FASTERQ_DUMP_BIN="${FASTERQ_DUMP_BIN:-}"
PREFETCH_BIN="${PREFETCH_BIN:-}"
PIGZ_BIN="${PIGZ_BIN:-}"
SEQTK_BIN="${SEQTK_BIN:-}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-}"
TRIMVALIDATE_BIN="${REPO_ROOT}/core/features/vbem/tools/trimvalidate/trimvalidate"
REMOVE_Y_READS_BIN="${REPO_ROOT}/core/features/yremove_fastq/tools/remove_y_reads/remove_y_reads"
TXIMPORT_BIN="${REPO_ROOT}/core/features/vbem/tools/tximport_compat/tximport_compat"
SAMPLE_FLD_BIN="${REPO_ROOT}/core/features/vbem/tools/sample_fld/sample_fld"
COMPUTE_EXPECTED_GC_BIN="${REPO_ROOT}/core/features/vbem/tools/compute_expected_gc/compute_expected_gc"
COMPUTE_GC_BIAS_BIN="${REPO_ROOT}/core/features/vbem/tools/compute_gc_bias/compute_gc_bias"
COMPARE_SCRIPT="${REPO_ROOT}/tests/transcriptvb/compare_salmon_star.py"
MAKE_GENE_MAP_SCRIPT="${REPO_ROOT}/tests/transcriptvb/make_gene_map_from_gtf.sh"

usage() {
    cat <<'EOF'
Usage:
  scripts/paper/gather_pe_bulk_external_tools.sh --outdir DIR [options]

Required:
  --outdir DIR                 Output directory for tool manifest/env files.

Optional:
  --star-bin PATH              STAR binary (default: repo STAR).
  --salmon-bin PATH            Salmon binary (default: /mnt/pikachu/salmon/build/src/salmon, fallback PATH).
  --fastq-dump-bin PATH        fastq-dump binary (fallback PATH).
  --fasterq-dump-bin PATH      fasterq-dump binary (fallback PATH).
  --prefetch-bin PATH          prefetch binary (fallback PATH).
  --pigz-bin PATH              pigz binary (fallback PATH).
  --seqtk-bin PATH             seqtk binary (fallback PATH).
  --samtools-bin PATH          samtools binary (fallback PATH).
  --trimvalidate-bin PATH      trimvalidate binary (default: repo path).
  --remove-y-reads-bin PATH    remove_y_reads binary (default: repo path).
  --tximport-bin PATH          tximport_compat binary (default: repo path).
  --sample-fld-bin PATH        sample_fld binary (default: repo path).
  --compute-expected-gc-bin PATH  compute_expected_gc binary (default: repo path).
  --compute-gc-bias-bin PATH   compute_gc_bias binary (default: repo path).
  --compare-script PATH        compare_salmon_star.py path.
  --make-gene-map-script PATH  make_gene_map_from_gtf.sh path.
  --no-build-missing           Do not build missing repo-local tools.
  --dry-run                    Print build actions without executing them.
  -h, --help                   Show this help.

Outputs:
  TOOL_MANIFEST.tsv
  tool_env.sh
  SUMMARY.txt
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

require_file() {
    [[ -f "$1" ]] || die "Missing file: $1"
}

sha_or_na() {
    local path="$1"
    if [[ ! -e "${path}" ]]; then
        echo "NA"
        return
    fi
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "${path}" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "${path}" | awk '{print $1}'
    else
        cksum "${path}" | awk '{print $1 ":" $2}'
    fi
}

version_or_na() {
    local path="$1"
    if [[ ! -x "${path}" ]]; then
        echo "NA"
        return
    fi
    case "$(basename "${path}")" in
        STAR|salmon)
            "${path}" --version 2>/dev/null | head -n1 || echo "NA"
            ;;
        *)
            echo "sha256=$(sha_or_na "${path}")"
            ;;
    esac
}

build_if_missing() {
    local target="$1"
    shift
    if [[ -x "${target}" ]]; then
        return
    fi
    if [[ "${BUILD_MISSING}" -ne 1 ]]; then
        die "Required binary missing and --no-build-missing set: ${target}"
    fi
    log "Building missing tool: ${target}"
    run_or_echo "$@"
}

resolve_salmon() {
    if [[ -n "${SALMON_BIN}" ]]; then
        [[ -x "${SALMON_BIN}" ]] || die "Salmon binary not executable: ${SALMON_BIN}"
        return
    fi
    if [[ -x "${SALMON_BIN_DEFAULT}" ]]; then
        SALMON_BIN="${SALMON_BIN_DEFAULT}"
        return
    fi
    if command -v salmon >/dev/null 2>&1; then
        SALMON_BIN="$(command -v salmon)"
        return
    fi
    die "Could not resolve Salmon binary. Checked ${SALMON_BIN_DEFAULT} and PATH."
}

resolve_from_path() {
    local var_name="$1"
    local current_value="$2"
    local command_name="$3"
    if [[ -n "${current_value}" ]]; then
        [[ -x "${current_value}" ]] || die "${command_name} binary not executable: ${current_value}"
        printf -v "${var_name}" '%s' "${current_value}"
        return
    fi
    if command -v "${command_name}" >/dev/null 2>&1; then
        printf -v "${var_name}" '%s' "$(command -v "${command_name}")"
        return
    fi
    die "Could not resolve ${command_name} from PATH"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir) OUTDIR="$2"; shift 2 ;;
        --star-bin) STAR_BIN="$2"; shift 2 ;;
        --salmon-bin) SALMON_BIN="$2"; shift 2 ;;
        --fastq-dump-bin) FASTQ_DUMP_BIN="$2"; shift 2 ;;
        --fasterq-dump-bin) FASTERQ_DUMP_BIN="$2"; shift 2 ;;
        --prefetch-bin) PREFETCH_BIN="$2"; shift 2 ;;
        --pigz-bin) PIGZ_BIN="$2"; shift 2 ;;
        --seqtk-bin) SEQTK_BIN="$2"; shift 2 ;;
        --samtools-bin) SAMTOOLS_BIN="$2"; shift 2 ;;
        --trimvalidate-bin) TRIMVALIDATE_BIN="$2"; shift 2 ;;
        --remove-y-reads-bin) REMOVE_Y_READS_BIN="$2"; shift 2 ;;
        --tximport-bin) TXIMPORT_BIN="$2"; shift 2 ;;
        --sample-fld-bin) SAMPLE_FLD_BIN="$2"; shift 2 ;;
        --compute-expected-gc-bin) COMPUTE_EXPECTED_GC_BIN="$2"; shift 2 ;;
        --compute-gc-bias-bin) COMPUTE_GC_BIAS_BIN="$2"; shift 2 ;;
        --compare-script) COMPARE_SCRIPT="$2"; shift 2 ;;
        --make-gene-map-script) MAKE_GENE_MAP_SCRIPT="$2"; shift 2 ;;
        --no-build-missing) BUILD_MISSING=0; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

[[ -n "${OUTDIR}" ]] || die "--outdir is required"

OUTDIR="$(mkdir -p "${OUTDIR}" && cd "${OUTDIR}" && pwd)"
MANIFEST="${OUTDIR}/TOOL_MANIFEST.tsv"
ENV_FILE="${OUTDIR}/tool_env.sh"
SUMMARY="${OUTDIR}/SUMMARY.txt"

build_if_missing "${STAR_BIN}" make -C "${REPO_ROOT}/core/legacy/source" -j8 STAR
build_if_missing "${TRIMVALIDATE_BIN}" make -C "${REPO_ROOT}/core/features/vbem/tools/trimvalidate" -j8
build_if_missing "${REMOVE_Y_READS_BIN}" make -C "${REPO_ROOT}/core/features/yremove_fastq/tools/remove_y_reads" -j8
build_if_missing "${TXIMPORT_BIN}" make -C "${REPO_ROOT}/core/features/vbem/tools/tximport_compat" -j8
build_if_missing "${SAMPLE_FLD_BIN}" make -C "${REPO_ROOT}/core/features/vbem/tools/sample_fld" -j8
build_if_missing "${COMPUTE_EXPECTED_GC_BIN}" make -C "${REPO_ROOT}/core/features/vbem/tools/compute_expected_gc" -j8
build_if_missing "${COMPUTE_GC_BIAS_BIN}" make -C "${REPO_ROOT}/core/features/vbem/tools/compute_gc_bias" -j8

resolve_salmon
resolve_from_path FASTQ_DUMP_BIN "${FASTQ_DUMP_BIN}" fastq-dump
resolve_from_path FASTERQ_DUMP_BIN "${FASTERQ_DUMP_BIN}" fasterq-dump
resolve_from_path PREFETCH_BIN "${PREFETCH_BIN}" prefetch
resolve_from_path PIGZ_BIN "${PIGZ_BIN}" pigz
resolve_from_path SEQTK_BIN "${SEQTK_BIN}" seqtk
resolve_from_path SAMTOOLS_BIN "${SAMTOOLS_BIN}" samtools

require_file "${COMPARE_SCRIPT}"
require_file "${MAKE_GENE_MAP_SCRIPT}"
[[ -x "${STAR_BIN}" ]] || die "STAR binary not executable: ${STAR_BIN}"
[[ -x "${SALMON_BIN}" ]] || die "Salmon binary not executable: ${SALMON_BIN}"
[[ -x "${FASTQ_DUMP_BIN}" ]] || die "fastq-dump binary not executable: ${FASTQ_DUMP_BIN}"
[[ -x "${FASTERQ_DUMP_BIN}" ]] || die "fasterq-dump binary not executable: ${FASTERQ_DUMP_BIN}"
[[ -x "${PREFETCH_BIN}" ]] || die "prefetch binary not executable: ${PREFETCH_BIN}"
[[ -x "${PIGZ_BIN}" ]] || die "pigz binary not executable: ${PIGZ_BIN}"
[[ -x "${SEQTK_BIN}" ]] || die "seqtk binary not executable: ${SEQTK_BIN}"
[[ -x "${SAMTOOLS_BIN}" ]] || die "samtools binary not executable: ${SAMTOOLS_BIN}"
[[ -x "${TRIMVALIDATE_BIN}" ]] || die "trimvalidate binary not executable: ${TRIMVALIDATE_BIN}"
[[ -x "${REMOVE_Y_READS_BIN}" ]] || die "remove_y_reads binary not executable: ${REMOVE_Y_READS_BIN}"
[[ -x "${TXIMPORT_BIN}" ]] || die "tximport_compat binary not executable: ${TXIMPORT_BIN}"
[[ -x "${SAMPLE_FLD_BIN}" ]] || die "sample_fld binary not executable: ${SAMPLE_FLD_BIN}"
[[ -x "${COMPUTE_EXPECTED_GC_BIN}" ]] || die "compute_expected_gc binary not executable: ${COMPUTE_EXPECTED_GC_BIN}"
[[ -x "${COMPUTE_GC_BIAS_BIN}" ]] || die "compute_gc_bias binary not executable: ${COMPUTE_GC_BIAS_BIN}"
[[ -x "${MAKE_GENE_MAP_SCRIPT}" ]] || die "make_gene_map script not executable: ${MAKE_GENE_MAP_SCRIPT}"

{
    printf 'tool\tpath\torigin\tkind\tidentity\n'
    printf 'STAR\t%s\trepo\tbinary\t%s\n' "${STAR_BIN}" "$(version_or_na "${STAR_BIN}")"
    printf 'salmon\t%s\texternal\tbinary\t%s\n' "${SALMON_BIN}" "$(version_or_na "${SALMON_BIN}")"
    printf 'fastq-dump\t%s\texternal\tbinary\t%s\n' "${FASTQ_DUMP_BIN}" "$(version_or_na "${FASTQ_DUMP_BIN}")"
    printf 'fasterq-dump\t%s\texternal\tbinary\t%s\n' "${FASTERQ_DUMP_BIN}" "$(version_or_na "${FASTERQ_DUMP_BIN}")"
    printf 'prefetch\t%s\texternal\tbinary\t%s\n' "${PREFETCH_BIN}" "$(version_or_na "${PREFETCH_BIN}")"
    printf 'pigz\t%s\texternal\tbinary\t%s\n' "${PIGZ_BIN}" "$(version_or_na "${PIGZ_BIN}")"
    printf 'seqtk\t%s\texternal\tbinary\t%s\n' "${SEQTK_BIN}" "$(version_or_na "${SEQTK_BIN}")"
    printf 'samtools\t%s\texternal\tbinary\t%s\n' "${SAMTOOLS_BIN}" "$(version_or_na "${SAMTOOLS_BIN}")"
    printf 'trimvalidate\t%s\trepo\tbinary\t%s\n' "${TRIMVALIDATE_BIN}" "$(version_or_na "${TRIMVALIDATE_BIN}")"
    printf 'remove_y_reads\t%s\trepo\tbinary\t%s\n' "${REMOVE_Y_READS_BIN}" "$(version_or_na "${REMOVE_Y_READS_BIN}")"
    printf 'tximport_compat\t%s\trepo\tbinary\t%s\n' "${TXIMPORT_BIN}" "$(version_or_na "${TXIMPORT_BIN}")"
    printf 'sample_fld\t%s\trepo\tbinary\t%s\n' "${SAMPLE_FLD_BIN}" "$(version_or_na "${SAMPLE_FLD_BIN}")"
    printf 'compute_expected_gc\t%s\trepo\tbinary\t%s\n' "${COMPUTE_EXPECTED_GC_BIN}" "$(version_or_na "${COMPUTE_EXPECTED_GC_BIN}")"
    printf 'compute_gc_bias\t%s\trepo\tbinary\t%s\n' "${COMPUTE_GC_BIAS_BIN}" "$(version_or_na "${COMPUTE_GC_BIAS_BIN}")"
    printf 'compare_salmon_star.py\t%s\trepo\tscript\tsha256=%s\n' "${COMPARE_SCRIPT}" "$(sha_or_na "${COMPARE_SCRIPT}")"
    printf 'make_gene_map_from_gtf.sh\t%s\trepo\tscript\tsha256=%s\n' "${MAKE_GENE_MAP_SCRIPT}" "$(sha_or_na "${MAKE_GENE_MAP_SCRIPT}")"
} > "${MANIFEST}"

{
    echo '#!/usr/bin/env bash'
    echo 'set -euo pipefail'
    printf 'export STAR_BIN=%q\n' "${STAR_BIN}"
    printf 'export SALMON_BIN=%q\n' "${SALMON_BIN}"
    printf 'export FASTQ_DUMP_BIN=%q\n' "${FASTQ_DUMP_BIN}"
    printf 'export FASTERQ_DUMP_BIN=%q\n' "${FASTERQ_DUMP_BIN}"
    printf 'export PREFETCH_BIN=%q\n' "${PREFETCH_BIN}"
    printf 'export PIGZ_BIN=%q\n' "${PIGZ_BIN}"
    printf 'export SEQTK_BIN=%q\n' "${SEQTK_BIN}"
    printf 'export SAMTOOLS_BIN=%q\n' "${SAMTOOLS_BIN}"
    printf 'export TRIMVALIDATE_BIN=%q\n' "${TRIMVALIDATE_BIN}"
    printf 'export REMOVE_Y_READS_BIN=%q\n' "${REMOVE_Y_READS_BIN}"
    printf 'export TXIMPORT_BIN=%q\n' "${TXIMPORT_BIN}"
    printf 'export SAMPLE_FLD_BIN=%q\n' "${SAMPLE_FLD_BIN}"
    printf 'export COMPUTE_EXPECTED_GC_BIN=%q\n' "${COMPUTE_EXPECTED_GC_BIN}"
    printf 'export COMPUTE_GC_BIAS_BIN=%q\n' "${COMPUTE_GC_BIAS_BIN}"
    printf 'export COMPARE_SALMON_STAR_SCRIPT=%q\n' "${COMPARE_SCRIPT}"
    printf 'export MAKE_GENE_MAP_SCRIPT=%q\n' "${MAKE_GENE_MAP_SCRIPT}"
} > "${ENV_FILE}"
chmod +x "${ENV_FILE}"

{
    echo "PE Bulk Benchmark Tool Summary"
    echo "Date (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "Repo root: ${REPO_ROOT}"
    echo "Manifest: ${MANIFEST}"
    echo "Env file: ${ENV_FILE}"
    echo
    awk 'NR==1 || $0 !~ /^$/' "${MANIFEST}"
} > "${SUMMARY}"

log "Wrote tool manifest: ${MANIFEST}"
log "Wrote tool env file: ${ENV_FILE}"
log "Wrote summary: ${SUMMARY}"
