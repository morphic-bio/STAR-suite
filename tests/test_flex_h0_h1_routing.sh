#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
CBQ_ENCODER_BIN="${CBQ_ORDERED_ENCODER_BIN:-${REPO_ROOT}/core/legacy/source/cbq_ordered_encoder}"
GENOME_DIR="${FLEX_ROUTING_GENOME_DIR:-/home/lhhung/jax_stage_20260903/ref/star_index}"
PROBE_LIST="${FLEX_ROUTING_PROBE_LIST:-${GENOME_DIR}/flex_probe_artifacts/probe_list.txt}"

die() {
    echo "FAIL: $*" >&2
    exit 1
}

[[ -x "${STAR_BIN}" ]] || die "STAR binary is absent: ${STAR_BIN}"
[[ -x "${CBQ_ENCODER_BIN}" ]] || die "CBQ encoder is absent: ${CBQ_ENCODER_BIN}"
[[ -f "${GENOME_DIR}/Genome" ]] || die "STAR index is absent: ${GENOME_DIR}"
[[ -f "${PROBE_LIST}" ]] || die "Flex probe list is absent: ${PROBE_LIST}"

TEST_ROOT="$(mktemp -d /tmp/star_flex_h0_h1_routing.XXXXXX)"
cleanup() {
    local status=$?
    if [[ ${status} -eq 0 ]]; then
        rm -rf "${TEST_ROOT}"
    else
        echo "Routing-test artifacts preserved at ${TEST_ROOT}" >&2
    fi
}
trap cleanup EXIT

# Build a minimal v3 cache and paired FASTQs with one record for every routing
# row. The sample-miss read deliberately carries an H0 probe so the test proves
# sample eligibility is decided before the cache lookup.
python3 - "${TEST_ROOT}" <<'PY'
from pathlib import Path
import struct
import sys

root = Path(sys.argv[1])
root.mkdir(parents=True, exist_ok=True)

h0 = "A" * 50
h1_keep = "C" * 50
h1_deny = "G" * 50
missing = "T" * 50
n_window = "A" * 24 + "N" + "A" * 25
matched_tag = "ACTTTAGG"
unmatched_tag = "NNNNNNNN"

def read2(probe, tag):
    return probe + "A" * 18 + tag + "A" * 14

records = [
    ("h0_keep", read2(h0, matched_tag)),
    ("h1_keep", read2(h1_keep, matched_tag)),
    ("h1_deny", read2(h1_deny, matched_tag)),
    ("n_pass", read2(n_window, matched_tag)),
    ("miss_pass", read2(missing, matched_tag)),
    ("sample_deny", read2(h0, unmatched_tag)),
]
r1 = "ACGTACGTACGTACGT" + "AACCGGTTAACC"

with (root / "r2.fastq").open("w", encoding="ascii") as out2, \
     (root / "r1.fastq").open("w", encoding="ascii") as out1:
    for name, r2 in records:
        out2.write(f"@{name}\n{r2}\n+\n{'I' * len(r2)}\n")
        out1.write(f"@{name}\n{r1}\n+\n{'I' * len(r1)}\n")

def encode(seq):
    code = {"A": 0, "C": 1, "G": 2, "T": 3}
    lo = hi = 0
    mask = (1 << 64) - 1
    for base in seq:
        hi = ((hi << 2) | (lo >> 62)) & mask
        lo = ((lo << 2) | code[base]) & mask
    return lo, hi

# cacheClass: 0=H0 KEEP, 1=H1 KEEP, 2=certified H1 DENY.
cache_records = [
    (*encode(h0), 1, 0, 0, 1),
    (*encode(h1_keep), 2, 1, 0, 0),
    (*encode(h1_deny), 0, 2, 1, 0),
]
cache_records.sort(key=lambda rec: (rec[1], rec[0], rec[5]))
with (root / "routing_cache.bin").open("wb") as out:
    out.write(struct.pack("<8sHHIQ", b"FH01SEQ1", 3, 50, 24, len(cache_records)))
    for lo, hi, gene, cache_class, negative, sample in cache_records:
        out.write(struct.pack("<QQIBBH", lo, hi, gene, cache_class, negative, sample))

(root / "sample_whitelist.tsv").write_text("BC001\tACTTTAGG\n", encoding="ascii")
(root / "sample_probes.tsv").write_text(
    "ACTTTAGG\tACTTTAGG\tBC001\n", encoding="ascii")
PY

"${CBQ_ENCODER_BIN}" \
    --readFilesIn "${TEST_ROOT}/r2.fastq" "${TEST_ROOT}/r1.fastq" \
    --outFile "${TEST_ROOT}/routing.cbq"

metric() {
    local log_file="$1"
    local label="$2"
    awk -F '|' -v label="${label}" \
        '$1 ~ label "[[:space:]]*$" {gsub(/[[:space:]]/, "", $2); print $2}' \
        "${log_file}"
}

run_case() {
    local input_kind="$1"
    local no_align="$2"
    local out_dir="${TEST_ROOT}/${input_kind}_noalign${no_align}"
    local -a input_args
    mkdir -p "${out_dir}"

    if [[ "${input_kind}" == "fastq" ]]; then
        input_args=(--readFilesIn "${TEST_ROOT}/r2.fastq" "${TEST_ROOT}/r1.fastq")
    else
        input_args=(--readFilesType Binseq PE --readFilesCbqRangeMode range
                    --readFilesIn "${TEST_ROOT}/routing.cbq")
    fi

    "${STAR_BIN}" \
        --runThreadN 2 \
        --dynamicThreadInterface 1 \
        --genomeDir "${GENOME_DIR}" \
        "${input_args[@]}" \
        --soloType CB_UMI_Simple \
        --soloCBstart 1 --soloCBlen 16 \
        --soloUMIstart 17 --soloUMIlen 12 \
        --soloBarcodeReadLength 0 \
        --soloCBwhitelist None \
        --soloFeatures Gene \
        --soloCellFilter None \
        --soloProbeList "${PROBE_LIST}" \
        --soloSampleWhitelist "${TEST_ROOT}/sample_whitelist.tsv" \
        --soloSampleProbes "${TEST_ROOT}/sample_probes.tsv" \
        --soloSampleProbeOffset 68 \
        --soloSampleSearchNearby yes \
        --soloSampleStrictMatch no \
        --soloFlexAllowedTags "${TEST_ROOT}/sample_whitelist.tsv" \
        --soloHashScreenFile "${TEST_ROOT}/routing_cache.bin" \
        --soloInlineHashMode yes \
        --soloBucketMode ram --soloBucketCount 4 \
        --flex yes \
        --soloFlexExpectedCellsPerTag 1 \
        --flexPipeline yes --flexPipelineNTriage 0 --flexPipelineNSolo 0 \
        --flexNoAlign "${no_align}" \
        --outSAMtype None --outSAMattributes None --outSJtype None \
        --outTmpDir "${out_dir}/_STARtmp" \
        --soloFlexOutputPrefix "${out_dir}/per_sample" \
        --outFileNamePrefix "${out_dir}/" \
        >"${out_dir}/stdout.log" 2>"${out_dir}/stderr.log"

    grep -Fq \
        "Flex pipeline complete: total=6, triageKeep=2, triageDeny=1, sampleReject=1, triageMiss=2" \
        "${out_dir}/Log.out" || die "${input_kind} flexNoAlign=${no_align} routing counters differ"

    local expected_pass=0
    if [[ "${no_align}" == 0 ]]; then
        expected_pass=2
    fi
    [[ "$(metric "${out_dir}/Log.final.out" 'Hash screen: KEEP')" == 2 ]] \
        || die "${input_kind} flexNoAlign=${no_align}: H0/H1 KEEP count differs"
    [[ "$(metric "${out_dir}/Log.final.out" 'Hash screen: DENY')" == 1 ]] \
        || die "${input_kind} flexNoAlign=${no_align}: certified H1 DENY count differs"
    [[ "$(metric "${out_dir}/Log.final.out" 'Hash screen: unmatched sample tag')" == 1 ]] \
        || die "${input_kind} flexNoAlign=${no_align}: sample-tag DENY count differs"
    [[ "$(metric "${out_dir}/Log.final.out" 'Hash screen: PASS')" == "${expected_pass}" ]] \
        || die "${input_kind} flexNoAlign=${no_align}: only PASS records should depend on flexNoAlign"
}

run_case fastq 0
run_case fastq 1
run_case cbq 0
run_case cbq 1

echo "PASS: all six Flex H0/H1 routing rows agree for ASCII FASTQ and packed CBQ; only PASS depends on --flexNoAlign"
