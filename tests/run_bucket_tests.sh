#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
PHASE_FILE="${ROOT_DIR}/tests/bucket/IMPLEMENTED_PHASE"
PHASES_FILE="${ROOT_DIR}/tests/bucket/PHASES.tsv"
HARNESS="${BUCKET_HARNESS_BIN:-${ROOT_DIR}/core/legacy/source/cb_bucket_harness}"
OUT_ROOT="${BUCKET_TEST_OUT_ROOT:-/tmp/star_suite_cb_bucket_tests}"
CURRENT_PHASE="${BUCKET_TEST_PHASE:-$(tr -d '[:space:]' < "${PHASE_FILE}")}"
RUN_REGRESSIONS="${BUCKET_RUN_REGRESSIONS:-1}"

die() { echo "FAIL: $*" >&2; exit 1; }
pass() { echo "PASS: $*"; }
skip() { echo "SKIP: $1 (pending Phase $2)"; }

[[ "${CURRENT_PHASE}" =~ ^[0-4]$ ]] || die "invalid implemented phase: ${CURRENT_PHASE}"
[[ -f "${PHASES_FILE}" ]] || die "missing phase manifest: ${PHASES_FILE}"

enabled() {
    local test_name="$1"
    local needed
    needed="$(awk -F '\t' -v name="${test_name}" '$1 == name { print $2 }' "${PHASES_FILE}")"
    [[ -n "${needed}" ]] || die "${test_name} is absent from ${PHASES_FILE}"
    if (( CURRENT_PHASE < needed )); then
        skip "${test_name}" "${needed}"
        return 1
    fi
    return 0
}

rm -rf "${OUT_ROOT}"
mkdir -p "${OUT_ROOT}/fixture" "${OUT_ROOT}/reference" "${OUT_ROOT}/results"
python3 "${ROOT_DIR}/tests/bucket/make_fixture.py" --outdir "${OUT_ROOT}/fixture"
for buckets in 64 256 1024; do
    python3 "${ROOT_DIR}/tests/bucket/reference_counter.py" \
        --r1 "${OUT_ROOT}/fixture/R1.fastq" \
        --r2 "${OUT_ROOT}/fixture/R2.fastq" \
        --barcodes "${OUT_ROOT}/fixture/barcodes.tsv" \
        --tags "${OUT_ROOT}/fixture/tags.tsv" \
        --probes "${OUT_ROOT}/fixture/probes.tsv" \
        --bucket-count "${buckets}" \
        > "${OUT_ROOT}/reference/reference_${buckets}.json"
done
observed_reference_sha="$(sha256sum "${OUT_ROOT}/reference/reference_256.json" | awk '{print $1}')"
expected_reference_sha="$(tr -d '[:space:]' < "${ROOT_DIR}/tests/bucket/reference_256.gold.sha256")"
[[ "${observed_reference_sha}" == "${expected_reference_sha}" ]] \
    || die "Phase 0 reference digest mismatch: ${observed_reference_sha} != ${expected_reference_sha}"
pass "Phase 0 deterministic multi-tag FASTQ and independent reference counter"

if enabled B1; then
    [[ -x "${HARNESS}" ]] || die "B1 requires ${HARNESS}"
    "${HARNESS}" --mode roundtrip --reference "${OUT_ROOT}/reference/reference_256.json"
    for producers in 1 8 32; do
        "${HARNESS}" --mode claims --producers "${producers}" \
            --records 8192 --bucket-count 256 \
            > "${OUT_ROOT}/results/claims_${producers}.json"
    done
    pass "B1 packed records, reference partition, and atomic claims"
fi

if enabled B2; then
    [[ -x "${HARNESS}" ]] || die "B2 requires ${HARNESS}"
    "${HARNESS}" --mode store --fixture "${OUT_ROOT}/fixture/truth.tsv" \
        --scratch "${OUT_ROOT}/spill"
    pass "B2 RAM/spill segment byte equality"
fi

if enabled B3; then
    "${ROOT_DIR}/tests/bucket/test_bucket_e2e.sh" B3 "${OUT_ROOT}"
    pass "B3 bucket tail end-to-end equality"
fi

if enabled B4; then
    "${ROOT_DIR}/tests/bucket/test_bucket_e2e.sh" B4 "${OUT_ROOT}"
    pass "B4 thread and bucket-count determinism"
fi

if enabled B5; then
    "${ROOT_DIR}/tests/bucket/test_bucket_e2e.sh" B5 "${OUT_ROOT}"
    pass "B5 spill and automatic-transition equality"
fi

if enabled B6; then
    "${ROOT_DIR}/tests/bucket/test_bucket_e2e.sh" B6 "${OUT_ROOT}"
    pass "B6 tag-parallel flexfilter equality"
fi

run_regressions() {
    local mode="$1"
    if [[ "${RUN_REGRESSIONS}" != 1 ]]; then
        echo "SKIP: B7-${mode} regression subprocesses (BUCKET_RUN_REGRESSIONS=${RUN_REGRESSIONS})"
        return
    fi
    local star_bin="${STAR_BIN:-${ROOT_DIR}/core/legacy/source/STAR}"
    local wrapper="${OUT_ROOT}/STAR-bucket-${mode}"
    if (( CURRENT_PHASE >= 3 )); then
        python3 - "${wrapper}" "${star_bin}" "${mode}" <<'PY'
import os
import shlex
import sys
from pathlib import Path
wrapper, star, mode = sys.argv[1:]
Path(wrapper).write_text(
    "#!/usr/bin/env bash\nexec " + shlex.quote(star) +
    " --soloBucketMode " + shlex.quote(mode) + " \"$@\"\n",
    encoding="utf-8")
os.chmod(wrapper, 0o755)
PY
        star_bin="${wrapper}"
    fi
    STAR_BIN="${star_bin}" \
        BGZF_TEST_OUT_ROOT="${OUT_ROOT}/bgzf_${mode}" \
        "${ROOT_DIR}/tests/run_bgzf_ingest_tests.sh"
}

if enabled B7-off; then
    run_regressions off
    pass "B7-off gzip/BGZF/CBQ regressions"
fi
if enabled B7-auto; then
    run_regressions auto
    pass "B7-auto gzip/BGZF/CBQ regressions"
fi
