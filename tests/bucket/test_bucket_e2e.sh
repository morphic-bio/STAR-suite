#!/usr/bin/env bash
set -euo pipefail

CASE="${1:-}"
OUT_ROOT="${2:-/tmp/star_suite_cb_bucket_tests}"
ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
STAR_BIN="${STAR_BIN:-${ROOT_DIR}/core/legacy/source/STAR}"

case "${CASE}" in B3|B4|B5|B6) ;; *) echo "usage: $0 B3|B4|B5|B6 [out-root]" >&2; exit 2 ;; esac

die() { echo "FAIL: $*" >&2; exit 1; }
[[ -x "${STAR_BIN}" ]] || die "STAR binary is absent: ${STAR_BIN}"
mkdir -p "${OUT_ROOT}"

canonical_manifest() {
    local run_dir="$1"
    local output="$2"
    python3 - "${run_dir}" "${output}" <<'PY'
import gzip
import hashlib
import sys
from pathlib import Path

root = Path(sys.argv[1])
output = Path(sys.argv[2])
selected = []
for base_name in ("Solo.out", "per_sample"):
    base = root / base_name
    if base.exists():
        selected.extend(path for path in base.rglob("*") if path.is_file())
if not selected:
    raise SystemExit(f"no canonical Flex outputs found under {root}")
lines = []
for path in sorted(selected):
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as handle:
            payload = handle.read()
    else:
        payload = path.read_bytes()
    lines.append(f"{path.relative_to(root)}\t{hashlib.sha256(payload).hexdigest()}")
output.write_text("\n".join(lines) + "\n", encoding="utf-8")
PY
}

compare_runs() {
    local expected="$1"
    local observed="$2"
    local label="$3"
    local manifests="${OUT_ROOT}/manifests"
    mkdir -p "${manifests}"
    local expected_manifest="${manifests}/$(printf '%s' "${label}" | tr '/ ' '__').expected"
    local observed_manifest="${manifests}/$(printf '%s' "${label}" | tr '/ ' '__').observed"
    canonical_manifest "${expected}" "${expected_manifest}"
    canonical_manifest "${observed}" "${observed_manifest}"
    diff -u "${expected_manifest}" "${observed_manifest}" \
        || die "canonical output mismatch: ${label}"
}

make_env_wrapper() {
    local wrapper="$1"
    local mode="$2"
    local bucket_count="$3"
    python3 - "${wrapper}" "${STAR_BIN}" "${mode}" "${bucket_count}" <<'PY'
import os
import shlex
import sys
from pathlib import Path

wrapper, star, mode, bucket_count = sys.argv[1:]
text = "#!/usr/bin/env bash\n"
if mode == "off":
    text += "unset STAR_SOLO_BUCKET_MODE STAR_SOLO_BUCKET_COUNT\n"
else:
    text += "export STAR_SOLO_BUCKET_MODE=" + shlex.quote(mode) + "\n"
    text += "export STAR_SOLO_BUCKET_COUNT=" + shlex.quote(bucket_count) + "\n"
text += "exec " + shlex.quote(star) + " \"$@\"\n"
Path(wrapper).write_text(text, encoding="utf-8")
os.chmod(wrapper, 0o755)
PY
}

run_gold_case() {
    local label="$1"
    local mode="$2"
    local threads="$3"
    local bucket_count="$4"
    local case_root="${OUT_ROOT}/gold_${label}"
    local wrapper="${OUT_ROOT}/STAR-${label}"
    if [[ -f "${case_root}/PASS" ]]; then
        return
    fi
    make_env_wrapper "${wrapper}" "${mode}" "${bucket_count}"
    if ! STAR_BIN="${wrapper}" BGZF_E2E_CASE=T4 \
        BGZF_E2E_THREADS="${threads}" BGZF_E2E_READ_LIMIT=2000 \
        BGZF_E2E_OUT_ROOT="${case_root}" \
        "${ROOT_DIR}/tests/bgzf/test_flex_e2e.sh" \
        >"${case_root}.log" 2>&1; then
        tail -100 "${case_root}.log" >&2 || true
        die "gold fixture failed for ${label}"
    fi
    if [[ "${mode}" != off ]]; then
        grep -F "Flex streaming CB buckets: active (ram, ${bucket_count} buckets)" \
            "${case_root}/runs/plain/Log.out" >/dev/null \
            || die "bucket path did not activate for ${label}"
    fi
    printf 'status=pass\n' > "${case_root}/PASS"
}

make_jax_800k_fixture() {
    local fixture="$1"
    local fastq_dir="${BUCKET_JAX_FASTQ_DIR:-/mnt/pikachu/JAX_sequences/JAX_scRNAseq01}"
    local stem="SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001"
    mkdir -p "${fixture}"
    for mate in R1 R2; do
        local source="${fastq_dir}/${stem}_${mate}_001.fastq.gz"
        local output="${fixture}/${mate}.fastq"
        [[ -f "${source}" ]] || return 1
        if [[ ! -f "${output}.complete" ]]; then
            local partial="${output}.partial"
            set +e
            gzip -cd "${source}" | head -n 3200000 > "${partial}"
            local head_status="${PIPESTATUS[1]}"
            set -e
            [[ "${head_status}" -eq 0 ]] || die "could not truncate ${source}"
            [[ "$(wc -l < "${partial}")" -eq 3200000 ]] \
                || die "${source} ended before 800000 FASTQ records"
            mv "${partial}" "${output}"
            printf 'records=800000\n' > "${output}.complete"
        fi
    done
    paste \
        <(awk 'NR%4==1 {sub(/^@/, ""); sub(/[ \/].*/, ""); print}' "${fixture}/R1.fastq") \
        <(awk 'NR%4==1 {sub(/^@/, ""); sub(/[ \/].*/, ""); print}' "${fixture}/R2.fastq") \
        | awk '$1!=$2 {print "mate mismatch at record " NR > "/dev/stderr"; exit 1}
               END {if (NR != 800000) exit 2}' \
        || die "JAX 800k fixture mate validation failed"
}

run_jax_case() {
    local label="$1"
    local mode="$2"
    local fixture="$3"
    local output="${OUT_ROOT}/jax800k_${label}"
    local genome_dir="${BUCKET_GENOME_DIR:-/storage/flex_filtered_reference_2024/star_index}"
    local cb_whitelist="${BUCKET_CB_WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}"
    local sample_whitelist="${BUCKET_SAMPLE_WHITELIST:-/mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv}"
    local probe_list="${BUCKET_PROBE_LIST:-${genome_dir}/flex_probe_artifacts/probe_list.txt}"
    local sample_probes="${BUCKET_SAMPLE_PROBES:-/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}"
    local hash_cache="${BUCKET_HASH_CACHE:-/storage/downsampled_100K/SC2300771/results/flex_h01_2024_20260320_081246/h01_cache.bin}"
    local required
    for required in "${genome_dir}" "${cb_whitelist}" "${sample_whitelist}" \
                    "${probe_list}" "${sample_probes}" "${hash_cache}"; do
        [[ -e "${required}" ]] || return 1
    done
    if [[ -f "${output}/PASS" ]]; then
        return
    fi
    mkdir -p "${output}"
    local bucket_env=(env -u STAR_SOLO_BUCKET_MODE -u STAR_SOLO_BUCKET_COUNT)
    if [[ "${mode}" != off ]]; then
        bucket_env=(env STAR_SOLO_BUCKET_MODE=ram STAR_SOLO_BUCKET_COUNT=256)
    fi
    "${bucket_env[@]}" "${STAR_BIN}" \
        --runThreadN 32 --genomeDir "${genome_dir}" \
        --soloType CB_UMI_Simple --soloCBstart 1 --soloUMIstart 17 \
        --soloCBlen 16 --soloUMIlen 12 --soloBarcodeReadLength 0 \
        --soloCBwhitelist "${cb_whitelist}" --flex yes \
        --soloFlexExpectedCellsPerTag 3000 \
        --soloSampleWhitelist "${sample_whitelist}" \
        --soloProbeList "${probe_list}" --soloSampleProbes "${sample_probes}" \
        --soloSampleProbeOffset 68 --soloFlexAllowedTags "${sample_whitelist}" \
        --soloFlexOutputPrefix "${output}/per_sample" \
        --limitIObufferSize 50000000 50000000 --outSJtype None \
        --outSAMtype None --outSAMattributes None --soloFeatures Gene \
        --soloCellFilter None --soloMultiMappers Rescue \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
        --soloStrand Unstranded --clipAdapterType CellRanger4 \
        --alignEndsType Local --chimSegmentMin 0 --soloKeysCompat cr \
        --soloSampleSearchNearby no --soloHashScreenFile "${hash_cache}" \
        --flexPipeline yes --flexPipelineNTriage 0 --flexPipelineNSolo 0 \
        --flexNoAlign 1 --dynamicThreadInterface 1 \
        --crAssignConsumerThreads -1 --crAssignSearchThreads 1 \
        --outFileNamePrefix "${output}/" --readFilesBgzfMode off \
        --readFilesIn "${fixture}/R2.fastq" "${fixture}/R1.fastq" \
        >"${output}/stdout.txt" 2>"${output}/stderr.txt"
    grep -F "Flex pipeline complete: total=800000" "${output}/Log.out" >/dev/null \
        || die "JAX ${label} did not consume exactly 800000 read pairs"
    if [[ "${mode}" != off ]]; then
        grep -F "Flex streaming CB buckets: active (ram, 256 buckets)" \
            "${output}/Log.out" >/dev/null \
            || die "JAX ${label} did not activate RAM buckets"
    fi
    printf 'status=pass\n' > "${output}/PASS"
}

case "${CASE}" in
B3)
    run_gold_case off off 4 256
    run_gold_case ram ram 4 256
    compare_runs "${OUT_ROOT}/gold_off/runs/plain" \
                 "${OUT_ROOT}/gold_ram/runs/plain" "B3_gold_off_vs_ram"

    jax_fixture="${OUT_ROOT}/jax800k_fixture"
    if make_jax_800k_fixture "${jax_fixture}" \
        && run_jax_case off off "${jax_fixture}" \
        && run_jax_case ram ram "${jax_fixture}"; then
        compare_runs "${OUT_ROOT}/jax800k_off" "${OUT_ROOT}/jax800k_ram" \
                     "B3_jax800k_off_vs_ram"
        off_counters="$(grep -F 'Flex pipeline complete:' "${OUT_ROOT}/jax800k_off/Log.out" | tail -1)"
        ram_counters="$(grep -F 'Flex pipeline complete:' "${OUT_ROOT}/jax800k_ram/Log.out" | tail -1)"
        [[ "${off_counters}" == "${ram_counters}" ]] \
            || die "JAX 800k pipeline counters differ"
        echo "PASS: B3 gold fixture and JAX 800k RAM-bucket equality"
    else
        echo "SKIP: B3 JAX 800k subcase (host fixture or reference assets absent)"
        echo "PASS: B3 gold fixture RAM-bucket equality"
    fi
    ;;
B4)
    run_gold_case off off 4 256
    baseline="${OUT_ROOT}/gold_off/runs/plain"
    for threads in 1 8 32; do
        for buckets in 64 256 1024; do
            label="ram_t${threads}_p${buckets}"
            run_gold_case "${label}" ram "${threads}" "${buckets}"
            compare_runs "${baseline}" "${OUT_ROOT}/gold_${label}/runs/plain" \
                         "B4_${label}"
        done
    done
    echo "PASS: B4 RAM-bucket equality at 1/8/32 threads and 64/256/1024 buckets"
    ;;
B5|B6)
    die "${CASE} is not implemented before Phase 3"
    ;;
esac
