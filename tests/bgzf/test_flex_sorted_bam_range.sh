#!/usr/bin/env bash
# Paper-scoped regression: paired BGZF Flex input may use the parallel reader
# while STAR's ordinary alignment path writes coordinate-sorted BAM. Other
# alignment-output and compression combinations must retain their old readers.
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
STAR_BIN="${STAR_BIN:-${ROOT_DIR}/core/legacy/source/STAR}"
WORKDIR="${BGZF_E2E_OUT_ROOT:-/tmp/star_suite_bgzf_flex_e2e}"
RUNS="${BGZF_SORTED_BAM_OUT_ROOT:-${WORKDIR}/sorted_bam_runs}"
THREADS="${BGZF_SORTED_BAM_THREADS:-4}"

die() { echo "FAIL: $*" >&2; exit 1; }
[[ -x "${STAR_BIN}" ]] || die "STAR binary is absent: ${STAR_BIN}"
command -v samtools >/dev/null 2>&1 || die "samtools is required"

# Reuse the deterministic tiny Flex fixture maintained by T4.
if [[ ! -e "${WORKDIR}/star_index/Genome" || ! -f "${WORKDIR}/blocked.config.csv" ]]; then
    BGZF_E2E_CASE=T4 BGZF_E2E_OUT_ROOT="${WORKDIR}" \
        "${ROOT_DIR}/tests/bgzf/test_flex_e2e.sh"
fi

for required in \
    "${WORKDIR}/star_index/Genome" \
    "${WORKDIR}/blocked.config.csv" \
    "${WORKDIR}/plain.config.csv" \
    "${WORKDIR}/assets_base/whitelist.txt" \
    "${WORKDIR}/assets_base/sample_probe_catalog.tsv"; do
    [[ -e "${required}" ]] || die "required fixture is absent: ${required}"
done

make_star_wrapper() {
    local name="$1"
    local mode="$2"
    local wrapper="${RUNS}/STAR-${name}"
    python3 - "${wrapper}" "${STAR_BIN}" "${mode}" <<'PY'
import os
import shlex
import sys
from pathlib import Path

wrapper, star, mode = sys.argv[1:]
Path(wrapper).write_text(
    "#!/usr/bin/env bash\nexec " + shlex.quote(star) +
    " --readFilesBgzfMode " + shlex.quote(mode) + " \"$@\"\n",
    encoding="utf-8",
)
os.chmod(wrapper, 0o755)
PY
}

run_case() {
    local name="$1"
    local layout="$2"
    local mode="$3"
    local output="$4"
    rm -rf "${RUNS}/runs/${name}"
    make_star_wrapper "${name}" "${mode}"
    local samtype="bam-${output}"
    STAR_BIN="${RUNS}/STAR-${name}" \
        "${ROOT_DIR}/scripts/run_flex_cr_config.sh" \
        --cr-config "${WORKDIR}/${layout}.config.csv" \
        --genome-dir "${WORKDIR}/star_index" \
        --cb-whitelist "${WORKDIR}/assets_base/whitelist.txt" \
        --solo-cb-start 1 --solo-cb-len 16 \
        --solo-umi-start 17 --solo-umi-len 10 \
        --sample-probe-catalog "${WORKDIR}/assets_base/sample_probe_catalog.tsv" \
        --sample-probe-offset 68 --out-samtype "${samtype}" \
        --out-base "${RUNS}/runs" --run-id "${name}" --threads "${THREADS}"
}

rm -rf "${RUNS}"
mkdir -p "${RUNS}/runs"

# Expected red before the core-range bridge exists.
run_case blocked_off_sorted blocked off sorted
run_case blocked_range_sorted blocked range sorted

off_run="${RUNS}/runs/blocked_off_sorted"
range_run="${RUNS}/runs/blocked_range_sorted"
grep -F "BGZF core range reader: active" "${range_run}/Log.out" >/dev/null \
    || die "paired BGZF Flex sorted-BAM path was not active"
grep -F "BGZF core range reader: active" "${off_run}/Log.out" >/dev/null \
    && die "--readFilesBgzfMode off activated the core range reader"
samtools quickcheck "${off_run}/Aligned.sortedByCoord.out.bam"
samtools quickcheck "${range_run}/Aligned.sortedByCoord.out.bam"
samtools view "${off_run}/Aligned.sortedByCoord.out.bam" > "${RUNS}/off.body.sam"
samtools view "${range_run}/Aligned.sortedByCoord.out.bam" > "${RUNS}/range.body.sam"
cmp "${RUNS}/off.body.sam" "${RUNS}/range.body.sam"
cmp "${off_run}/Solo.out/Barcodes.stats" "${range_run}/Solo.out/Barcodes.stats"

# Auto remains conservative outside the exact paper case.
run_case plain_auto_sorted plain auto sorted
grep -F "BGZF core range reader: active" "${RUNS}/runs/plain_auto_sorted/Log.out" >/dev/null \
    && die "plain gzip activated the BGZF core range reader"
run_case blocked_auto_unsorted blocked auto unsorted
grep -F "BGZF core range reader: active" "${RUNS}/runs/blocked_auto_unsorted/Log.out" >/dev/null \
    && die "unsorted BAM activated the BGZF core range reader"

# The bridge is not a general STAR/STARsolo BGZF reader. A coordinate-sorted
# BAM request without Flex must stay on the established Fastx/zlib path.
classic_run="${RUNS}/runs/blocked_auto_classic"
rm -rf "${classic_run}"
mkdir -p "${classic_run}"
"${STAR_BIN}" \
    --readFilesBgzfMode auto \
    --runThreadN "${THREADS}" \
    --genomeDir "${WORKDIR}/star_index" \
    --readFilesIn \
        "${WORKDIR}/blocked_fastq/tinyflex_S1_L001_R2_001.fastq.gz" \
        "${WORKDIR}/blocked_fastq/tinyflex_S1_L001_R1_001.fastq.gz" \
    --outFileNamePrefix "${classic_run}/" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes Standard
grep -F "BGZF core range reader: active" "${classic_run}/Log.out" >/dev/null \
    && die "classic STAR activated the paper-scoped Flex BGZF bridge"
samtools quickcheck "${classic_run}/Aligned.sortedByCoord.out.bam"

echo "PASS: Flex BGZF coordinate-sorted BAM core-range parity (${RUNS})"
