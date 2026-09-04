#!/usr/bin/env bash
# Regression: a fully-fused Flex run with alignment enabled (--flexNoAlign 0)
# must finish when hash misses exceed the alignQ capacity (256).
#
# Before the fix every fused thread was a producer until it ran out of input
# and nobody drained alignQ, so a blocking push on the full queue waited for a
# consumer that did not exist yet: the run hung on the 257th miss. The fix
# makes a producer that finds the queue full align queued packets itself.
# This test therefore asserts three things: the run completes, the help path
# actually fired (queue pressure was real), and a single-thread run, which a
# reserved-consumer design cannot support, completes as well.
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
STAR_BIN="${STAR_BIN:-${ROOT_DIR}/core/legacy/source/STAR}"
WORKDIR="${BGZF_E2E_OUT_ROOT:-/tmp/star_suite_bgzf_flex_e2e}"
TIMEOUT_SECONDS="${BGZF_FUSED_ALIGN_TIMEOUT:-900}"

die() { echo "FAIL: $*" >&2; exit 1; }
[[ -x "${STAR_BIN}" ]] || die "STAR binary is absent: ${STAR_BIN}"

# Reuse the index, hash cache, and blocked FASTQ that test_flex_e2e.sh builds.
if [[ ! -f "${WORKDIR}/fused_inputs/hash_cache.bin" || ! -e "${WORKDIR}/star_index/Genome" ]]; then
    BGZF_E2E_CASE=T4 "${ROOT_DIR}/tests/bgzf/test_flex_e2e.sh"
fi
for required in "${WORKDIR}/star_index/Genome" "${WORKDIR}/fused_inputs/hash_cache.bin" \
                "${WORKDIR}/blocked.config.csv" \
                "${WORKDIR}/assets_base/whitelist.txt" \
                "${WORKDIR}/assets_base/sample_probe_catalog.tsv" \
                "${WORKDIR}/blocked_fastq/tinyflex_S1_L001_R1_001.fastq.gz" \
                "${WORKDIR}/blocked_fastq/tinyflex_S1_L001_R2_001.fastq.gz"; do
    [[ -e "${required}" ]] || die "required fixture is absent: ${required}"
done

# Every read must miss the hash cache so that every read is queued for
# alignment: scramble the probe region (the first 50 bases of R2) and keep the
# rest of the read, including the sample tag at offset 68, intact. Give every
# R1 a barcode one mismatch from two whitelist entries, so the aligned records
# also prove that the shared ambiguous store receives real observations rather
# than merely exercising an empty enable/drain lifecycle.
miss_dir="${WORKDIR}/miss_fastq"
rm -rf "${miss_dir}"
mkdir -p "${miss_dir}"
ambig_cb="$(python3 - "${WORKDIR}/assets_base/whitelist.txt" \
                          "${WORKDIR}/blocked_fastq/tinyflex_S1_L001_R1_001.fastq.gz" \
                          "${miss_dir}/r1_ambiguous.fastq" <<'PY'
import gzip
import sys

whitelist_path, source, dest = sys.argv[1:]
with open(whitelist_path, encoding="ascii") as handle:
    whitelist = [line.strip() for line in handle if line.strip()]
whitelist_set = set(whitelist)
candidate = None
for ii, left in enumerate(whitelist):
    for right in whitelist[ii + 1:]:
        differences = [jj for jj, pair in enumerate(zip(left, right)) if pair[0] != pair[1]]
        if len(differences) != 2:
            continue
        pos = differences[1]
        hybrid = left[:pos] + right[pos] + left[pos + 1:]
        if hybrid not in whitelist_set:
            candidate = hybrid
            break
    if candidate is not None:
        break
if candidate is None:
    raise SystemExit("could not construct a barcode one mismatch from two whitelist entries")

count = 0
with gzip.open(source, "rt") as fin, open(dest, "w", newline="\n") as fout:
    while True:
        name = fin.readline()
        if not name:
            break
        seq = fin.readline().rstrip("\n")
        plus = fin.readline()
        qual = fin.readline()
        if len(seq) < len(candidate):
            raise SystemExit("R1 is shorter than the generated ambiguous barcode")
        fout.write(name)
        fout.write(candidate + seq[len(candidate):] + "\n")
        fout.write(plus)
        fout.write(qual)
        count += 1
if count == 0:
    raise SystemExit("R1 fixture contained no reads")
print(candidate)
PY
)"
[[ -n "${ambig_cb}" ]] || die "failed to construct an ambiguous CB fixture"
"${ROOT_DIR}/tools/make_bgzf_fixture.sh" --block-bytes 1021 "${miss_dir}/r1_ambiguous.fastq" \
    "${miss_dir}/tinyflex_S1_L001_R1_001.fastq.gz"
rm -f "${miss_dir}/r1_ambiguous.fastq"
n_reads="$(python3 - "${WORKDIR}/blocked_fastq/tinyflex_S1_L001_R2_001.fastq.gz" \
                     "${miss_dir}/r2_scrambled.fastq" <<'PY'
import gzip
import random
import sys
src, dest = sys.argv[1:]
rng = random.Random(20260903)
count = 0
with gzip.open(src, "rt") as fin, open(dest, "w", newline="\n") as fout:
    while True:
        name = fin.readline()
        if not name:
            break
        seq = fin.readline().rstrip("\n")
        plus = fin.readline()
        qual = fin.readline()
        head = "".join(rng.choice("ACGT") for _ in range(min(50, len(seq))))
        fout.write(name)
        fout.write(head + seq[len(head):] + "\n")
        fout.write(plus)
        fout.write(qual)
        count += 1
print(count)
PY
)"
(( n_reads > 256 )) || die "fixture has ${n_reads} reads; need more than the alignQ capacity (256)"
"${ROOT_DIR}/tools/make_bgzf_fixture.sh" --block-bytes 1021 "${miss_dir}/r2_scrambled.fastq" \
    "${miss_dir}/tinyflex_S1_L001_R2_001.fastq.gz"
rm -f "${miss_dir}/r2_scrambled.fastq"
sed "s|${WORKDIR}/blocked_fastq|${miss_dir}|g" "${WORKDIR}/blocked.config.csv" \
    > "${WORKDIR}/miss.config.csv"
grep -qF "${miss_dir}" "${WORKDIR}/miss.config.csv" \
    || die "config rewrite did not point at ${miss_dir}"

# One wrapper per input path: "range" exercises the BGZF range reader,
# "off" forces the same BGZF files through the plain zlib lane reader, so the
# help path is covered at both producer sites.
make_wrapper() {
    local mode="$1"
    local wrapper="${WORKDIR}/STAR-fused-align-${mode}"
    cat > "${wrapper}" <<EOF
#!/usr/bin/env bash
exec "${STAR_BIN}" --readFilesBgzfMode ${mode} \\
    --soloHashScreenFile "${WORKDIR}/fused_inputs/hash_cache.bin" \\
    --flexPipeline yes --flexPipelineNTriage 0 --flexPipelineNSolo 0 --flexNoAlign 0 "\$@"
EOF
    chmod +x "${wrapper}"
}
make_wrapper range
make_wrapper off

metric() {
    awk -F '|' -v label="$2" '$1 ~ label "[[:space:]]*$" {gsub(/[[:space:]]/, "", $2); print $2}' "$1"
}

run_case() {
    local threads="$1"
    local mode="$2"
    local run_id="fused_align_${mode}_${threads}t"
    rm -rf "${WORKDIR}/runs/${run_id}"
    set +e
    STAR_BIN="${WORKDIR}/STAR-fused-align-${mode}" timeout "${TIMEOUT_SECONDS}" \
        "${ROOT_DIR}/scripts/run_flex_cr_config.sh" \
        --cr-config "${WORKDIR}/miss.config.csv" \
        --genome-dir "${WORKDIR}/star_index" \
        --cb-whitelist "${WORKDIR}/assets_base/whitelist.txt" \
        --solo-cb-start 1 --solo-cb-len 16 --solo-umi-start 17 --solo-umi-len 10 \
        --sample-probe-catalog "${WORKDIR}/assets_base/sample_probe_catalog.tsv" \
        --sample-probe-offset 68 --out-samtype none \
        --out-base "${WORKDIR}/runs" --run-id "${run_id}" --threads "${threads}" \
        > "${WORKDIR}/${run_id}.stdout" 2> "${WORKDIR}/${run_id}.stderr"
    local rc=$?
    set -e
    if (( rc == 124 )); then
        die "${run_id}: fully-fused alignment did not finish within ${TIMEOUT_SECONDS}s (alignQ deadlock regression)"
    fi
    (( rc == 0 )) || die "${run_id}: run exited ${rc} (see ${WORKDIR}/${run_id}.stderr)"

    local log="${WORKDIR}/runs/${run_id}/Log.out"
    local final="${WORKDIR}/runs/${run_id}/Log.final.out"
    if [[ "${mode}" == range ]]; then
        grep -F "BGZF parallel range readers: active" "${log}" >/dev/null \
            || die "${run_id}: range readers were not active"
    else
        grep -F "BGZF parallel range readers: active" "${log}" >/dev/null \
            && die "${run_id}: range readers were active with --readFilesBgzfMode off"
    fi
    grep -F "Flex shared ambiguous store: active for fully fused hash and alignment records" \
        "${log}" >/dev/null \
        || die "${run_id}: shared ambiguous store was not enabled for the alignment path"
    grep -F "Solo timing: sumThreads adopt shared ambiguous" "${log}" >/dev/null \
        || die "${run_id}: shared ambiguous store was not drained by Gene sumThreads"
    local shared_entries
    shared_entries="$(sed -n 's/.*sumThreads adopt shared ambiguous .* (\([0-9][0-9]*\) entries).*/\1/p' \
        "${log}" | tail -1)"
    [[ -n "${shared_entries}" ]] \
        || die "${run_id}: shared ambiguous drain did not report an entry count"
    (( shared_entries > 0 )) \
        || die "${run_id}: aligned ambiguous observations bypassed the shared store"
    local helped
    helped="$(sed -n 's/.*Fused producers aligned \([0-9][0-9]*\) queued reads while alignQ was full.*/\1/p' "${log}" | tail -1)"
    [[ -n "${helped}" ]] || die "${run_id}: help-path summary line is absent from Log.out"
    (( helped > 0 )) || die "${run_id}: producers never aligned while alignQ was full (queue pressure was not exercised)"
    local input_reads misses
    input_reads="$(metric "${final}" 'Number of input reads')"
    misses="$(metric "${final}" 'Hash screen: PASS')"
    (( input_reads >= 257 && input_reads <= n_reads )) \
        || die "${run_id}: ${input_reads} input reads counted, expected between 257 and ${n_reads}"
    (( misses >= 257 )) || die "${run_id}: only ${misses} reads missed the hash cache; the queue never filled"
    echo "PASS: ${run_id} finished; ${misses} misses queued, producers aligned ${helped} of them while alignQ was full; shared store drained ${shared_entries} ambiguous keys"
}

run_case 2 range
run_case 1 range
run_case 2 off
for run_id in fused_align_range_1t fused_align_off_2t; do
    for relative in \
        Solo.out/Barcodes.stats \
        Solo.out/Gene/raw/barcodes.tsv \
        Solo.out/Gene/raw/features.tsv \
        Solo.out/Gene/raw/matrix.mtx \
        per_sample/flex_gdna_library.json \
        per_sample/flex_gdna_summary.tsv \
        per_sample/flexfilter_summary.tsv; do
        cmp "${WORKDIR}/runs/fused_align_range_2t/${relative}" \
            "${WORKDIR}/runs/${run_id}/${relative}" \
            || die "${run_id}: ${relative} differs from the 2-thread range-reader result"
    done
done
echo "PASS: alignment outputs are byte-identical across input routes and thread counts"
echo "PASS: fully-fused alignment mode drains alignQ without a reserved consumer (${WORKDIR})"
