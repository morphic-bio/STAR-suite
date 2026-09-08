#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
CBQ_ENCODER_BIN="${CBQ_ORDERED_ENCODER_BIN:-${REPO_ROOT}/core/legacy/source/cbq_ordered_encoder}"
DECISION_DUMP_BIN="${FLEX_DECISION_DUMP_BIN:-${REPO_ROOT}/flex/tools/molecule_first_resolver/flex_decision_sidecar_dump}"
GENOME_DIR="${FLEX_ROUTING_GENOME_DIR:-/home/lhhung/jax_stage_20260903/ref/star_index}"
PROBE_LIST="${FLEX_ROUTING_PROBE_LIST:-${GENOME_DIR}/flex_probe_artifacts/probe_list.txt}"

die() {
    echo "FAIL: $*" >&2
    exit 1
}

[[ -x "${STAR_BIN}" ]] || die "STAR binary is absent: ${STAR_BIN}"
[[ -x "${CBQ_ENCODER_BIN}" ]] || die "CBQ encoder is absent: ${CBQ_ENCODER_BIN}"
[[ -x "${DECISION_DUMP_BIN}" ]] || die "Flex decision sidecar dump tool is absent: ${DECISION_DUMP_BIN}"
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
h1x2_keep = "AC" * 25
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
    ("h1x2_keep", read2(h1x2_keep, matched_tag)),
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

# cacheClass: 0=H0 KEEP, 1=H1 KEEP, 2=certified DENY,
# 4=experimental H1X2 KEEP.
cache_records = [
    (*encode(h0), 1, 0, 0, 1),
    (*encode(h1_keep), 2, 1, 0, 0),
    (*encode(h1x2_keep), 2, 4, 0, 0),
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
(root / "cb_whitelist.txt").write_text("ACGTACGTACGTACGT\n", encoding="ascii")
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
    local probe_mismatch="$3"
    local sidecar_mode="${4:-on}"
    local out_dir="${TEST_ROOT}/${input_kind}_noalign${no_align}_probe${probe_mismatch}_sidecar${sidecar_mode}"
    local -a input_args
    local -a probe_mismatch_args=()
    local -a sidecar_args=()
    mkdir -p "${out_dir}"

    if [[ "${input_kind}" == "fastq" ]]; then
        input_args=(--readFilesIn "${TEST_ROOT}/r2.fastq" "${TEST_ROOT}/r1.fastq")
    else
        input_args=(--readFilesType Binseq PE --readFilesCbqRangeMode range
                    --readFilesIn "${TEST_ROOT}/routing.cbq")
    fi
    if [[ "${probe_mismatch}" != "default" ]]; then
        probe_mismatch_args=(--soloProbeMismatch "${probe_mismatch}")
    fi
    if [[ "${sidecar_mode}" == "on" ]]; then
        sidecar_args=(--soloFlexDecisionSidecar "${out_dir}/flex_decisions.bin")
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
        --soloCBwhitelist "${TEST_ROOT}/cb_whitelist.txt" \
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
        "${probe_mismatch_args[@]}" \
        "${sidecar_args[@]}" \
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

    local expected_keep=4
    local expected_miss=1
    if [[ "${probe_mismatch}" == "0" ]]; then
        expected_keep=3
        expected_miss=2
    fi
    local expected_deny=1
    if [[ "${no_align}" == 0 ]]; then
        # The only ordinary miss is T*50. With an H1X2 tier loaded, residual
        # alignment now requires a unique active-probe half anchor, so this
        # synthetic no-anchor read is denied before it reaches STAR.
        expected_deny=2
        expected_miss=0
    fi
    grep -Fq \
        "Flex pipeline complete: total=7, triageKeep=${expected_keep}, triageDeny=${expected_deny}, sampleReject=1, triageMiss=${expected_miss}" \
        "${out_dir}/Log.out" \
        || die "${input_kind} flexNoAlign=${no_align} soloProbeMismatch=${probe_mismatch}: routing counters differ"

    local expected_pass=0
    if [[ "${no_align}" == 0 ]]; then
        expected_pass="${expected_miss}"
    fi
    [[ "$(metric "${out_dir}/Log.final.out" 'Hash screen: KEEP')" == "${expected_keep}" ]] \
        || die "${input_kind} flexNoAlign=${no_align} soloProbeMismatch=${probe_mismatch}: H0/H1 KEEP count differs"
    [[ "$(metric "${out_dir}/Log.final.out" 'Hash screen: DENY')" == "${expected_deny}" ]] \
        || die "${input_kind} flexNoAlign=${no_align} soloProbeMismatch=${probe_mismatch}: cache/half-anchor DENY count differs"
    [[ "$(metric "${out_dir}/Log.final.out" 'Hash screen: unmatched sample tag')" == 1 ]] \
        || die "${input_kind} flexNoAlign=${no_align} soloProbeMismatch=${probe_mismatch}: sample-tag DENY count differs"
    [[ "$(metric "${out_dir}/Log.final.out" 'Hash screen: PASS')" == "${expected_pass}" ]] \
        || die "${input_kind} flexNoAlign=${no_align} soloProbeMismatch=${probe_mismatch}: only PASS records should depend on flexNoAlign"

    if [[ "${sidecar_mode}" == "off" ]]; then
        [[ ! -e "${out_dir}/flex_decisions.bin" ]] \
            || die "default-off run unexpectedly produced a decision sidecar"
        return
    fi

    "${DECISION_DUMP_BIN}" "${out_dir}/flex_decisions.bin" \
        >"${out_dir}/flex_decisions.tsv"
    [[ "$(($(wc -l <"${out_dir}/flex_decisions.tsv") - 1))" == 7 ]] \
        || die "${input_kind}: decision sidecar does not contain seven records"
    awk -F '\t' '
        NR == 2 && !($1 == 0 && $5 == "KEEP" && $6 == "H0" && $10 == 1 && $28 == "CACHE_KEEP") { exit 1 }
        NR == 3 && !($1 == 1 && $5 == "KEEP" && $6 == "H1" && $10 == 2 && $28 == "CACHE_KEEP") { exit 1 }
        NR == 4 && !($1 == 2 && $5 == "KEEP" && $6 == "H1X2" && $10 == 2 && $28 == "CACHE_KEEP") { exit 1 }
        NR == 5 && !($1 == 3 && $5 == "DENY" && $6 == "NEGATIVE" && $11 == 1 && $28 == "CACHE_DENY") { exit 1 }
        NR == 8 && !($1 == 6 && $5 == "DENY" && $6 == "." && $7 == "." && $15 == 1 && $28 == "SAMPLE_TAG_REJECT") { exit 1 }
    ' "${out_dir}/flex_decisions.tsv" \
        || die "${input_kind}: fixed H0/H1/H1X2/deny/sample records differ"
    if [[ "${probe_mismatch}" == "0" ]]; then
        awk -F '\t' 'NR == 6 { exit !($5 == "MISS" && $8 == 0 && $9 == 0) }' \
            "${out_dir}/flex_decisions.tsv" \
            || die "${input_kind}: disabled single-N retry was recorded incorrectly"
    else
        awk -F '\t' 'NR == 6 { exit !($5 == "KEEP" && $6 == "H1" && $7 == "H0" && $8 == 1 && $9 == 1) }' \
            "${out_dir}/flex_decisions.tsv" \
            || die "${input_kind}: single-N provenance differs"
    fi
    if [[ "${no_align}" == 1 ]]; then
        awk -F '\t' 'NR == 7 { exit !($5 == "MISS" && $17 == 0 && $18 == 0 && $27 == 1 && $28 == "CACHE_MISS_NO_ALIGN") }' \
            "${out_dir}/flex_decisions.tsv" \
            || die "${input_kind}: no-align miss provenance differs"
    else
        awk -F '\t' 'NR == 7 { exit !($5 == "DENY" && $17 == 0 && $18 == 0 && $22 == "ABSENT" && $28 == "RESIDUAL_NO_ANCHOR") }' \
            "${out_dir}/flex_decisions.tsv" \
            || die "${input_kind}: residual half-anchor denial provenance differs"
    fi
}

run_case fastq 0 default off
run_case fastq 0 default on
run_case fastq 1 default
run_case cbq 0 default
run_case cbq 1 default
run_case fastq 1 0
run_case cbq 1 0

run_bam_sidecar_case() {
    local input_kind="$1"
    local out_dir="${TEST_ROOT}/${input_kind}_bam_sidecar"
    local -a input_args
    if [[ "${input_kind}" == "fastq" ]]; then
        input_args=(--readFilesIn "${TEST_ROOT}/r2.fastq" "${TEST_ROOT}/r1.fastq")
    else
        input_args=(--readFilesType Binseq PE --readFilesCbqRangeMode off
                    --readFilesIn "${TEST_ROOT}/routing.cbq")
    fi
    mkdir -p "${out_dir}"
    "${STAR_BIN}" \
        --runThreadN 2 \
        --dynamicThreadInterface 1 \
        --genomeDir "${GENOME_DIR}" \
        "${input_args[@]}" \
        --soloType CB_UMI_Simple \
        --soloCBstart 1 --soloCBlen 16 \
        --soloUMIstart 17 --soloUMIlen 12 \
        --soloBarcodeReadLength 0 \
        --soloCBwhitelist "${TEST_ROOT}/cb_whitelist.txt" \
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
        --soloFlexDecisionSidecar "${out_dir}/flex_decisions.bin" \
        --flex yes --flexPipeline no \
        --soloFlexExpectedCellsPerTag 1 \
        --outSAMtype BAM Unsorted \
        --outSAMattributes NH HI AS nM NM CB UB \
        --outSJtype None \
        --outTmpDir "${out_dir}/_STARtmp" \
        --soloFlexOutputPrefix "${out_dir}/per_sample" \
        --outFileNamePrefix "${out_dir}/" \
        >"${out_dir}/stdout.log" 2>"${out_dir}/stderr.log"

    [[ -s "${out_dir}/Aligned.out.bam" ]] \
        || die "ordinary Flex BAM path did not produce BAM output"
    if command -v samtools >/dev/null 2>&1; then
        samtools quickcheck "${out_dir}/Aligned.out.bam" \
            || die "ordinary Flex BAM output failed samtools quickcheck"
    fi
    "${DECISION_DUMP_BIN}" "${out_dir}/flex_decisions.bin" \
        >"${out_dir}/flex_decisions.tsv"
    awk -F '\t' '
        NR == 2 && !($5 == "KEEP" && $6 == "H0") { exit 1 }
        NR == 3 && !($5 == "KEEP" && $6 == "H1") { exit 1 }
        NR == 4 && !($5 == "KEEP" && $6 == "H1X2") { exit 1 }
        NR == 5 && !($5 == "DENY" && $6 == "NEGATIVE") { exit 1 }
        NR == 6 && !($5 == "KEEP" && $6 == "H1" && $7 == "H0" && $8 == 1 && $9 == 1 && $17 == 0 && $28 == "CACHE_KEEP") { exit 1 }
        NR == 8 && !($5 == "DENY" && $6 == "." && $15 == 1 && $28 == "SAMPLE_TAG_REJECT") { exit 1 }
    ' "${out_dir}/flex_decisions.tsv" \
        || die "ordinary ${input_kind} Flex/BAM decision provenance differs"
}

run_bam_sidecar_case fastq
run_bam_sidecar_case cbq

echo "PASS: FASTQ and packed CBQ agree for H0/H1/H1X2 plus single-N routing; --soloProbeMismatch 0 restores exact-cache behavior"
