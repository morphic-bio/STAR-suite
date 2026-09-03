#!/usr/bin/env bash
# JAX Flex benchmark matrix: three input paths x alignment on/off, on the
# post-Flex-fix build (fused alignQ deadlock fix + mate read-name validation).
#
# Protocol, per the recorded benchmark discipline:
#   * final code only  - STAR is built from the commit recorded in each
#     preflight file; the matrix refuses to start on a dirty working tree
#     unless ALLOW_DIRTY=1.
#   * SSD only         - FASTQ and CBQ inputs are staged on the /home NVMe and
#     verified byte-for-byte; outputs, STAR's temp directory, TMPDIR, and the
#     bucket spill directory all live on that same SSD. Nothing in a timed run
#     touches the disk array.
#   * cold page cache  - dropped before every timed run; a failure to drop is
#     fatal rather than a silent warm run.
#   * quiet box        - load and the top processes are recorded per run.
#   * one variable     - every arm shares one thread count, one bucket mode
#     (ram: JAX fits, so no spill I/O enters the comparison), one reference
#     set, and one hash cache. Only the input path and --flexNoAlign change.
#
# Arms (6): {gzip, bgzf, cbq} x {noalign, align}
#   gzip    delivered fastq.gz read through the serial zlib lane path
#   bgzf    the same files read in parallel as BGZF (--readFilesBgzfMode range)
#   cbq     the converted binary CBQ set, record-range partitioned
#   noalign hash-only counting, no genome (--flexNoAlign 1)
#   align   hash first, residual misses aligned (--flexNoAlign 0); the genome
#           loads and the fused help path drains alignQ
#
# Replicates: REPS (default 3), interleaved arm-by-arm, medians taken later.
#
# Usage: run_jax_matrix.sh [arm ...]     (default: all six)
#   REPS=1 run_jax_matrix.sh bgzf_align
set -uo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
STAR_BIN="${STAR_BIN:-${ROOT_DIR}/core/legacy/source/STAR}"
THREADS="${THREADS:-32}"
REPS="${REPS:-3}"
REP_START="${REP_START:-1}"   # resume later replicates without redoing rep 1
BUCKET_MODE="${BUCKET_MODE:-ram}"
BUCKET_COUNT="${BUCKET_COUNT:-256}"
ALLOW_DIRTY="${ALLOW_DIRTY:-0}"

# Array-resident originals (read only while staging, never during a timed run)
SRC_FASTQ="${JAX_SRC_FASTQ:-/mnt/pikachu/JAX_sequences/JAX_scRNAseq01}"
SRC_CBQ="${JAX_SRC_CBQ:-/mnt/pikachu/JAX_sequences/JAX_scRNAseq01_cbq_parallel_20260730_v1}"
FASTQ_GLOB='SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L*_R[12]_001.fastq.gz'

# Everything below is on the /home NVMe
STAGE_ROOT="${JAX_STAGE_ROOT:-/home/lhhung/jax_stage_20260903}"
STAGE_FASTQ="${STAGE_ROOT}/fastq"
STAGE_CBQ="${STAGE_ROOT}/cbq"
STAGE_REF="${STAGE_ROOT}/ref"
ARTIFACT_ROOT="${ARTIFACT_ROOT:-/home/lhhung/jax_matrix_20260903}"
export TMPDIR="${ARTIFACT_ROOT}/tmp"
# MATRIX_LOG_DIR lets diagnostic runs (profiled, or with experimental hooks)
# keep their logs and results rows away from the benchmark record.
LOG_DIR="${MATRIX_LOG_DIR:-${ROOT_DIR}/docs/benchmarks/jax_matrix_20260903}"
RESULTS="${LOG_DIR}/results.tsv"
SENTINEL="${ARTIFACT_ROOT}/matrix.done"

# References. Every one of these is read inside a timed run - the align arms
# load the 35 GiB genome index - and /storage/flex_filtered_reference_2024 is a
# symlink onto the disk array (offloaded 2026-07-29), so the whole reference
# set is staged onto the same SSD as the reads and read from there.
SRC_GENOME_DIR="${SRC_GENOME_DIR:-/storage/flex_filtered_reference_2024/star_index}"
SRC_CB_WHITELIST="${SRC_CB_WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}"
SRC_HASH_CACHE="${SRC_HASH_CACHE:-/storage/downsampled_100K/SC2300771/results/flex_h01_2024_20260320_081246/h01_cache.bin}"
SRC_SAMPLE_WHITELIST="${SRC_SAMPLE_WHITELIST:-/mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv}"
SRC_SAMPLE_PROBES="${SRC_SAMPLE_PROBES:-/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}"
GENOME_DIR="${STAGE_REF}/star_index"
PROBE_LIST="${GENOME_DIR}/flex_probe_artifacts/probe_list.txt"
CB_WHITELIST="${STAGE_REF}/$(basename "${SRC_CB_WHITELIST}")"
HASH_CACHE="${STAGE_REF}/$(basename "${SRC_HASH_CACHE}")"
SAMPLE_WHITELIST="${STAGE_REF}/$(basename "${SRC_SAMPLE_WHITELIST}")"
SAMPLE_PROBES="${STAGE_REF}/$(basename "${SRC_SAMPLE_PROBES}")"

ALL_ARMS=(gzip_noalign bgzf_noalign cbq_noalign gzip_align bgzf_align cbq_align)
if [[ $# -gt 0 ]]; then ARMS=("$@"); else ARMS=("${ALL_ARMS[@]}"); fi

say() { echo "[$(date -u +%FT%TZ)] $*"; }
fail() { echo "FAIL: $*" >&2; echo "failed=$*" >> "${SENTINEL}.failed"; }
die() { echo "FATAL: $*" >&2; printf 'status=fatal\nreason=%s\n' "$*" > "${SENTINEL}"; exit 1; }

mkdir -p "${ARTIFACT_ROOT}" "${TMPDIR}" "${LOG_DIR}" "${STAGE_REF}"
rm -f "${SENTINEL}" "${SENTINEL}.failed"

# ---------------------------------------------------------------- preflight --
for arm in "${ARMS[@]}"; do
    printf '%s\n' "${ALL_ARMS[@]}" | grep -qxF "${arm}" || die "unknown arm: ${arm}"
done
[[ -x "${STAR_BIN}" ]] || die "STAR binary is absent: ${STAR_BIN}"
COMMIT="$(git -C "${ROOT_DIR}" rev-parse HEAD)"
DIRTY="$(git -C "${ROOT_DIR}" status --porcelain --untracked-files=no -- core tests | head -5)"
if [[ -n "${DIRTY}" && "${ALLOW_DIRTY}" != 1 ]]; then
    die "working tree is dirty under core/ or tests/; timed runs require final code:"$'\n'"${DIRTY}"
fi
for path in "${SRC_GENOME_DIR}" "${SRC_CB_WHITELIST}" "${SRC_HASH_CACHE}" \
            "${SRC_SAMPLE_WHITELIST}" "${SRC_SAMPLE_PROBES}"; do
    [[ -e "${path}" ]] || die "required input is absent: ${path}"
done
# Outputs and temp must sit on the staging SSD. df resolves symlinks, so this
# also catches a reference tree that has been offloaded to the array.
stage_dev="$(df --output=source "$(dirname "${STAGE_ROOT}")" | tail -1 | tr -d ' ')"
case "${stage_dev}" in /dev/nvme*) ;; *) die "the staging root is not on an SSD: ${STAGE_ROOT} (${stage_dev})" ;; esac
for d in "${ARTIFACT_ROOT}" "${TMPDIR}"; do
    [[ "$(df --output=source "${d}" | tail -1 | tr -d ' ')" == "${stage_dev}" ]] \
        || die "${d} is not on the staging SSD (${stage_dev})"
done
sudo -n sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches' \
    || die "passwordless cache drop is unavailable; a warm run would be invalid"

# ------------------------------------------------------------------ staging --
cp -f "${SRC_SAMPLE_WHITELIST}" "${SAMPLE_WHITELIST}"
cp -f "${SRC_SAMPLE_PROBES}" "${SAMPLE_PROBES}"
cp -f "${SRC_CB_WHITELIST}" "${CB_WHITELIST}"
cp -f "${SRC_HASH_CACHE}" "${HASH_CACHE}"
for pair in "${SRC_CB_WHITELIST}:${CB_WHITELIST}" "${SRC_HASH_CACHE}:${HASH_CACHE}"; do
    src="${pair%%:*}"; dst="${pair##*:}"
    [[ "$(stat -c%s "${src}")" == "$(stat -c%s "${dst}" 2>/dev/null || echo 0)" ]] \
        || die "reference stage size mismatch: ${dst}"
done

# The genome index is only needed by the align arms, but stage it whenever one
# is requested so no timed run reads 35 GiB off the array.
need_genome=0
for arm in "${ARMS[@]}"; do [[ "${arm}" == *_align ]] && need_genome=1; done
if (( need_genome == 1 )); then
    src_bytes="$(du -sb "${SRC_GENOME_DIR}/" | cut -f1)"
    dst_bytes="$(du -sb "${GENOME_DIR}/" 2>/dev/null | cut -f1 || echo 0)"
    if [[ "${src_bytes}" != "${dst_bytes}" ]]; then
        avail="$(df --output=avail -B1 "${STAGE_REF}" | tail -1 | tr -d ' ')"
        (( avail > src_bytes + 21474836480 )) \
            || die "insufficient SSD space for the genome index: need $((src_bytes / 1073741824)) GiB"
        say "staging the genome index ($((src_bytes / 1073741824)) GiB) to ${GENOME_DIR}"
        mkdir -p "${GENOME_DIR}"
        rsync -a --delete "${SRC_GENOME_DIR}/" "${GENOME_DIR}/" || die "genome index staging failed"
        dst_bytes="$(du -sb "${GENOME_DIR}/" | cut -f1)"
    fi
    [[ "${src_bytes}" == "${dst_bytes}" ]] \
        || die "genome index stage size mismatch: src=${src_bytes} stage=${dst_bytes}"
    [[ -e "${GENOME_DIR}/Genome" && -e "${PROBE_LIST}" ]] \
        || die "staged genome index is incomplete (Genome or probe list absent)"
    say "genome index verified on the SSD: ${GENOME_DIR}"
fi

# Final locus assertion: nothing a timed run reads may resolve to the array.
for path in "${GENOME_DIR}" "${CB_WHITELIST}" "${HASH_CACHE}" "${PROBE_LIST}" \
            "${SAMPLE_WHITELIST}" "${SAMPLE_PROBES}"; do
    [[ -e "${path}" ]] || continue
    [[ "$(df --output=source "${path}" | tail -1 | tr -d ' ')" == "${stage_dev}" ]] \
        || die "reference is not on the staging SSD: ${path}"
done

stage_set() {  # stage_set <src> <dst> <find-name> <expected count>
    local src="$1" dst="$2" pattern="$3" expected="$4" missing=0
    mkdir -p "${dst}"
    mapfile -t sources < <(find "${src}" -maxdepth 1 -type f -name "${pattern}" | sort)
    [[ "${#sources[@]}" -eq "${expected}" ]] \
        || die "expected ${expected} files matching ${pattern} in ${src}, found ${#sources[@]}"
    local need=0
    for f in "${sources[@]}"; do
        local b s1 s2
        b="$(basename "${f}")"; s1="$(stat -c%s "${f}")"
        s2="$(stat -c%s "${dst}/${b}" 2>/dev/null || echo 0)"
        [[ "${s1}" == "${s2}" ]] || need=$((need + s1 - s2))
    done
    if (( need > 0 )); then
        local avail; avail="$(df --output=avail -B1 "${dst}" | tail -1 | tr -d ' ')"
        (( avail > need + 21474836480 )) \
            || die "insufficient SSD space: need $((need / 1073741824)) GiB plus headroom, have $((avail / 1073741824)) GiB"
        say "staging $((need / 1073741824)) GiB to ${dst}"
        rsync -a --info=progress2 "${sources[@]}" "${dst}/" || die "rsync to ${dst} failed"
    fi
    for f in "${sources[@]}"; do
        local b s1 s2
        b="$(basename "${f}")"; s1="$(stat -c%s "${f}")"
        s2="$(stat -c%s "${dst}/${b}" 2>/dev/null || echo 0)"
        if [[ "${s1}" != "${s2}" ]]; then echo "MISMATCH ${b} src=${s1} stage=${s2}" >&2; missing=1; fi
    done
    (( missing == 0 )) || die "stage verification failed for ${dst}"
    say "stage verified: ${dst} ($(du -sBG "${dst}" | cut -f1))"
}

need_fastq=0; need_cbq=0
for arm in "${ARMS[@]}"; do
    case "${arm}" in gzip_*|bgzf_*) need_fastq=1 ;; cbq_*) need_cbq=1 ;; esac
done
(( need_fastq == 1 )) && stage_set "${SRC_FASTQ}" "${STAGE_FASTQ}" "${FASTQ_GLOB}" 16
(( need_cbq == 1 )) && stage_set "${SRC_CBQ}" "${STAGE_CBQ}" 'lane_*.cbq' 8

mapfile -t reads_r2 < <(find "${STAGE_FASTQ}" -maxdepth 1 -type f -name '*_R2_001.fastq.gz' 2>/dev/null | sort)
mapfile -t reads_r1 < <(find "${STAGE_FASTQ}" -maxdepth 1 -type f -name '*_R1_001.fastq.gz' 2>/dev/null | sort)
mapfile -t cbq_files < <(find "${STAGE_CBQ}" -maxdepth 1 -type f -name 'lane_*.cbq' 2>/dev/null | sort)
join_comma() { local IFS=,; echo "$*"; }

[[ -f "${RESULTS}" ]] || printf 'arm\trep\tstamp\tinput\talign\twall\twall_s\tpeak_rss_gib\tcpu_pct\tinput_reads\thash_keep\thash_deny\thash_pass\tuniq_mapped\thelped\tsamples\tmanifest_sha\tstatus\n' > "${RESULTS}"

# ------------------------------------------------------------------ one arm --
run_arm() {
    local arm="$1" rep="$2"
    local input_mode="${arm%%_*}" align_mode="${arm##*_}"
    local stamp; stamp="$(date -u +%Y%m%dT%H%M%SZ)"
    local run_dir="${ARTIFACT_ROOT}/${arm}_rep${rep}_${stamp}"
    local prefix="${LOG_DIR}/${arm}_rep${rep}_${stamp}"
    local no_align=1; [[ "${align_mode}" == align ]] && no_align=0

    mkdir -p "${run_dir}"
    {
        date -u '+date_utc=%Y-%m-%dT%H:%M:%SZ'
        printf 'arm=%s\nrep=%s\ninput=%s\nalign=%s\nthreads=%s\nbucket_mode=%s\nbucket_count=%s\n' \
            "${arm}" "${rep}" "${input_mode}" "${align_mode}" "${THREADS}" "${BUCKET_MODE}" "${BUCKET_COUNT}"
        printf 'run_dir=%s\ncommit=%s\nstar_mtime=%s\n' "${run_dir}" "${COMMIT}" \
            "$(stat -c%y "${STAR_BIN}")"
        if [[ -n "${DIRTY}" ]]; then
            printf 'working_tree=DIRTY (ALLOW_DIRTY=1; not a benchmark of the recorded commit)\n'
            git -C "${ROOT_DIR}" diff --stat -- core tests | sed 's/^/  /'
        else
            printf 'working_tree=clean\n'
        fi
        printf 'env_STAR_FLEX_HASH_H0_ONLY=%s\n' "${STAR_FLEX_HASH_H0_ONLY:-unset}"
        printf 'stage_fastq=%s\nstage_cbq=%s\ntmpdir=%s\n' "${STAGE_FASTQ}" "${STAGE_CBQ}" "${TMPDIR}"
        uptime
        free -g | head -2
        ps -eo pid,comm,%cpu,%mem,etime --sort=-%cpu | sed -n '1,12p'
    } > "${prefix}.preflight.txt"

    local common=(
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
        --soloFlexOutputPrefix "${run_dir}/per_sample"
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
        --soloBucketMode "${BUCKET_MODE}"
        --soloBucketCount "${BUCKET_COUNT}"
        --flexPipeline yes
        --flexPipelineNTriage 0
        --flexPipelineNSolo 0
        --flexNoAlign "${no_align}"
        --dynamicThreadInterface 1
        --crAssignConsumerThreads -1
        --crAssignSearchThreads 1
        --outFileNamePrefix "${run_dir}/"
    )
    local input=()
    case "${input_mode}" in
        gzip) input=(--readFilesBgzfMode off
                     --readFilesIn "$(join_comma "${reads_r2[@]}")" "$(join_comma "${reads_r1[@]}")") ;;
        bgzf) input=(--readFilesBgzfMode range --bgzfReaderThreads "${THREADS}" --bgzfCrcCheck 1
                     --readFilesIn "$(join_comma "${reads_r2[@]}")" "$(join_comma "${reads_r1[@]}")") ;;
        cbq)  input=(--readFilesType Binseq PE --readFilesCbqRangeMode range
                     --readFilesIn "$(join_comma "${cbq_files[@]}")") ;;
    esac
    local cmd=("${common[@]}" "${input[@]}")
    { printf '%q ' "${cmd[@]}"; printf '\n'; } > "${prefix}.command.txt"

    sudo -n sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches' || { fail "${arm} rep${rep}: cache drop failed"; return 1; }
    say "${arm} rep${rep} starting (${input_mode}, align=${align_mode})"
    /usr/bin/time -v -o "${prefix}.time.txt" "${cmd[@]}" \
        > "${prefix}.stdout.txt" 2> "${prefix}.stderr.txt"
    local rc=$?
    for log in Log.out Log.final.out Log.progress.out; do
        [[ -f "${run_dir}/${log}" ]] && cp "${run_dir}/${log}" "${prefix}.${log}"
    done

    # ---- verification -----------------------------------------------------
    local status="ok"
    (( rc == 0 )) || status="exit_${rc}"
    local nogenome=0 helped=""
    grep -Fq 'Flex count-only no-genome: active' "${run_dir}/Log.out" 2>/dev/null && nogenome=1
    helped="$(sed -n 's/.*Fused producers aligned \([0-9][0-9]*\) queued reads.*/\1/p' "${run_dir}/Log.out" 2>/dev/null | tail -1)"
    if [[ "${align_mode}" == noalign ]]; then
        (( nogenome == 1 )) || status="${status};no_genome_path_inactive"
    else
        (( nogenome == 0 )) || status="${status};genome_unexpectedly_skipped"
        [[ -n "${helped}" ]] || status="${status};help_counter_absent"
    fi

    metric() { awk -F '|' -v l="$1" '$1 ~ l "[[:space:]]*$" {gsub(/[[:space:]]/,"",$2); print $2}' "${run_dir}/Log.final.out" 2>/dev/null; }
    local reads keep deny pass uniq
    reads="$(metric 'Number of input reads')"
    keep="$(metric 'Hash screen: KEEP')"
    deny="$(metric 'Hash screen: DENY')"
    pass="$(metric 'Hash screen: PASS')"
    uniq="$(metric 'Uniquely mapped reads number')"
    local wall wall_s rss cpu
    wall="$(awk -F': ' '/Elapsed \(wall clock\)/ {print $NF}' "${prefix}.time.txt" 2>/dev/null)"
    wall_s="$(awk -F': ' '/Elapsed \(wall clock\)/ {n=split($NF,t,":"); s=0; for(i=1;i<=n;i++) s=s*60+t[i]; printf "%.2f", s}' "${prefix}.time.txt" 2>/dev/null)"
    rss="$(awk -F': ' '/Maximum resident set size/ {printf "%.1f", $NF/1048576}' "${prefix}.time.txt" 2>/dev/null)"
    cpu="$(awk -F': ' '/Percent of CPU/ {print $NF}' "${prefix}.time.txt" 2>/dev/null)"

    # Canonical per-sample outputs, hashed so arms can be compared for identity.
    local samples manifest_sha
    samples="$(find "${run_dir}/per_sample" -name 'matrix.mtx*' 2>/dev/null | wc -l)"
    python3 - "${run_dir}" "${prefix}.manifest.tsv" <<'PY' 2>/dev/null
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
    manifest_sha="$(sha256sum "${prefix}.manifest.tsv" 2>/dev/null | cut -c1-16)"
    (( samples >= 16 )) || status="${status};samples_${samples}"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${arm}" "${rep}" "${stamp}" "${input_mode}" "${align_mode}" \
        "${wall:-NA}" "${wall_s:-NA}" "${rss:-NA}" "${cpu:-NA}" \
        "${reads:-NA}" "${keep:-NA}" "${deny:-NA}" "${pass:-NA}" "${uniq:-NA}" \
        "${helped:-NA}" "${samples}" "${manifest_sha:-NA}" "${status}" >> "${RESULTS}"

    if [[ "${status}" == ok ]]; then
        say "${arm} rep${rep} OK: wall=${wall} rss=${rss}GiB cpu=${cpu} reads=${reads}"
        # Keep rep 1 matrices for concordance work; later reps keep logs + hashes.
        (( rep > 1 )) && rm -rf "${run_dir}/per_sample"
        printf 'status=success\nrun_dir=%s\n' "${run_dir}" > "${prefix}.completion.txt"
        return 0
    fi
    fail "${arm} rep${rep}: ${status}"
    printf 'status=%s\nrun_dir=%s\n' "${status}" "${run_dir}" > "${prefix}.completion.txt"
    return 1
}

# --------------------------------------------------------------------- run --
say "matrix start: commit=${COMMIT} arms='${ARMS[*]}' reps=${REPS} threads=${THREADS} bucket=${BUCKET_MODE}"
n_ok=0; n_fail=0
for rep in $(seq "${REP_START}" "${REPS}"); do
    for arm in "${ARMS[@]}"; do
        if run_arm "${arm}" "${rep}"; then n_ok=$((n_ok + 1)); else n_fail=$((n_fail + 1)); fi
    done
done

{
    printf 'finished=%s\ncommit=%s\narms=%s\nreps=%s\nok=%s\nfailed=%s\n' \
        "$(date -u +%FT%TZ)" "${COMMIT}" "${ARMS[*]}" "${REPS}" "${n_ok}" "${n_fail}"
    printf 'results=%s\n' "${RESULTS}"
    column -t -s $'\t' "${RESULTS}" 2>/dev/null || cat "${RESULTS}"
} > "${SENTINEL}"
say "matrix complete: ${n_ok} ok, ${n_fail} failed; summary at ${SENTINEL}"
(( n_fail == 0 )) || exit 1
