#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
PHASE_FILE="${ROOT_DIR}/tests/bgzf/IMPLEMENTED_PHASE"
PHASES_FILE="${ROOT_DIR}/tests/bgzf/PHASES.tsv"
HARNESS="${BGZF_INPUT_HARNESS_BIN:-${ROOT_DIR}/core/legacy/source/bgzf_input_harness}"
SCANNER="${ROOT_DIR}/tests/bgzf/reference_scanner.py"
MAKE_FIXTURE="${ROOT_DIR}/tools/make_bgzf_fixture.sh"
OUT_ROOT="${BGZF_TEST_OUT_ROOT:-/tmp/star_suite_bgzf_ingest_tests}"
CURRENT_PHASE="${BGZF_TEST_PHASE:-$(tr -d '[:space:]' < "${PHASE_FILE}")}"
BGZF_T7_FULL="${BGZF_T7_FULL:-0}"

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
mkdir -p "${OUT_ROOT}/inputs" "${OUT_ROOT}/results"

python3 - "${OUT_ROOT}/inputs/source.fastq" <<'PY'
import sys
from pathlib import Path

out = Path(sys.argv[1])
with out.open("w", encoding="ascii", newline="\n") as handle:
    for idx in range(257):
        length = 73 + (idx % 41)
        seq = ("ACGTN" * ((length + 4) // 5))[:length]
        qual = chr(33 + (idx % 40)) * length
        handle.write(f"@read{idx:06d}/1 instrument:run:{idx}\n{seq}\n+read{idx:06d}\n{qual}\n")
PY
gzip -n -c "${OUT_ROOT}/inputs/source.fastq" > "${OUT_ROOT}/inputs/plain.fastq.gz"
"${MAKE_FIXTURE}" --block-bytes 997 \
    "${OUT_ROOT}/inputs/plain.fastq.gz" "${OUT_ROOT}/inputs/blocked.fastq.gz"

python3 - "${OUT_ROOT}/inputs/overlong_line.fastq" <<'PY'
import random
import sys
from pathlib import Path

# Force one incompressible malformed sequence line across many BGZF members
# and multiple 64 KiB compressed work claims. The reader must reject it as
# soon as its fixed FASTQ capacity is exceeded, without growing its assembly
# buffer until a newline eventually appears.
length = 1 << 20
rng = random.Random(1903)
sequence = "".join(rng.choices("ACGT", k=length))
Path(sys.argv[1]).write_text(
    "@overlong/1\n" + sequence + "\n+\n" + "I" * length + "\n",
    encoding="ascii")
PY
"${MAKE_FIXTURE}" --block-bytes 4096 \
    "${OUT_ROOT}/inputs/overlong_line.fastq" \
    "${OUT_ROOT}/inputs/overlong_line.fastq.gz"
overlong_compressed_size=$(wc -c < "${OUT_ROOT}/inputs/overlong_line.fastq.gz")
(( overlong_compressed_size > 64 * 1024 )) \
    || die "overlong-line fixture does not cross multiple compressed work claims"

python3 - "${OUT_ROOT}/inputs/other_extra.fastq.gz" "${OUT_ROOT}/inputs/source.fastq" <<'PY'
import struct
import sys
import zlib
from pathlib import Path

dest = Path(sys.argv[1])
payload = Path(sys.argv[2]).read_bytes()
co = zlib.compressobj(6, zlib.DEFLATED, -15)
body = co.compress(payload) + co.flush()
extra = b"XY" + struct.pack("<H", 2) + b"zz"
header = b"\x1f\x8b\x08\x04" + b"\0" * 4 + b"\0\xff" + struct.pack("<H", len(extra)) + extra
trailer = struct.pack("<II", zlib.crc32(payload), len(payload) & 0xffffffff)
dest.write_bytes(header + body + trailer)
PY
python3 - "${OUT_ROOT}/inputs/bare_eof.fastq.gz" <<'PY'
import sys
from pathlib import Path
Path(sys.argv[1]).write_bytes(bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"))
PY
: > "${OUT_ROOT}/inputs/empty"
python3 - "${OUT_ROOT}/inputs/blocked.fastq.gz" "${OUT_ROOT}/inputs/missing_eof.fastq.gz" \
        "${OUT_ROOT}/inputs/mid_block.fastq.gz" <<'PY'
import json
import subprocess
import sys
from pathlib import Path

source, no_eof, mid_block = map(Path, sys.argv[1:])
data = source.read_bytes()
no_eof.write_bytes(data[:-28])
scanner = Path(__file__).resolve() if False else None
# The fixture uses many small members. Cutting ten bytes before the end of the
# second member guarantees a member-local truncation rather than an EOF-only loss.
first_size = int.from_bytes(data[16:18], "little") + 1
second_offset = first_size
second_size = int.from_bytes(data[second_offset + 16:second_offset + 18], "little") + 1
mid_block.write_bytes(data[:second_offset + second_size - 10])
PY

python3 "${SCANNER}" "${OUT_ROOT}/inputs/blocked.fastq.gz" \
    > "${OUT_ROOT}/results/reference.json"
python3 - "${OUT_ROOT}/results/reference.json" <<'PY'
import json
import sys
result = json.load(open(sys.argv[1], encoding="utf-8"))
assert result["bgzf"] and result["has_eof"]
assert result["record_count"] == 257
assert len(result["blocks"]) > 8
assert any(block["first_record_offset"] not in (None, 0) for block in result["blocks"])
PY
pass "Phase 0 fixture generator and independent reference scanner"

if enabled T1; then
    [[ -x "${HARNESS}" ]] || die "T1 requires ${HARNESS}"
    for case_name in blocked plain other_extra empty bare_eof missing_eof; do
        case "${case_name}" in
            blocked) input="${OUT_ROOT}/inputs/blocked.fastq.gz" ;;
            plain) input="${OUT_ROOT}/inputs/plain.fastq.gz" ;;
            other_extra) input="${OUT_ROOT}/inputs/other_extra.fastq.gz" ;;
            empty) input="${OUT_ROOT}/inputs/empty" ;;
            bare_eof) input="${OUT_ROOT}/inputs/bare_eof.fastq.gz" ;;
            missing_eof) input="${OUT_ROOT}/inputs/missing_eof.fastq.gz" ;;
        esac
        "${HARNESS}" --mode detect --input "${input}" \
            > "${OUT_ROOT}/results/detect_${case_name}.json"
    done
    python3 - "${OUT_ROOT}/results" <<'PY'
import json
import pathlib
import sys
root = pathlib.Path(sys.argv[1])
load = lambda name: json.loads((root / f"detect_{name}.json").read_text())
assert load("blocked") == {"bgzf": True, "has_eof": True}
assert load("plain") == {"bgzf": False, "has_eof": False}
assert load("other_extra") == {"bgzf": False, "has_eof": False}
assert load("empty") == {"bgzf": False, "has_eof": False}
assert load("bare_eof") == {"bgzf": True, "has_eof": True}
assert load("missing_eof") == {"bgzf": True, "has_eof": False}
PY
    pass "T1 BGZF detection matrix"
fi

if enabled T2; then
    [[ -x "${HARNESS}" ]] || die "T2 requires ${HARNESS}"
    "${HARNESS}" --mode blocks --input "${OUT_ROOT}/inputs/blocked.fastq.gz" \
        > "${OUT_ROOT}/results/blocks.json"
    python3 - "${OUT_ROOT}/results/reference.json" \
        "${OUT_ROOT}/results/blocks.json" <<'PY'
import json
import sys
ref, observed_blocks = [json.load(open(path, encoding="utf-8")) for path in sys.argv[1:]]
fields = ("compressed_offset", "compressed_size", "isize")
assert len(observed_blocks["blocks"]) == len(ref["blocks"])
for observed, expected in zip(observed_blocks["blocks"], ref["blocks"]):
    assert {k: observed[k] for k in fields} == {k: expected[k] for k in fields}
PY
    [[ ! -e "${OUT_ROOT}/inputs/blocked.fastq.gz.bgzi" ]] \
        || die "T2 runtime unexpectedly wrote a BGZF sidecar"
    pass "T2 inline BC/BSIZE traversal without a sidecar"
fi

if enabled T3; then
    [[ -x "${HARNESS}" ]] || die "T3 requires ${HARNESS}"
    for workers in 1 3 8; do
        "${HARNESS}" --mode records --input "${OUT_ROOT}/inputs/blocked.fastq.gz" \
            --workers "${workers}" --crc-check 1 \
            > "${OUT_ROOT}/results/records_${workers}.json"
        python3 - "${OUT_ROOT}/results/reference.json" "${OUT_ROOT}/results/records_${workers}.json" <<'PY'
import json
import sys
reference, observed = [json.load(open(path, encoding="utf-8")) for path in sys.argv[1:]]
assert observed["record_count"] == reference["record_count"]
assert observed["checksum64"] == reference["checksum64"]
PY
    done
    pass "T3 record equality at 1, 3, and 8 workers"
fi

if enabled T4; then
    BGZF_E2E_CASE=T4 "${ROOT_DIR}/tests/bgzf/test_flex_e2e.sh"
    pass "T4 Flex end-to-end equivalence"
fi

if enabled T5; then
    [[ -x "${HARNESS}" ]] || die "T5 requires ${HARNESS}"
    "${HARNESS}" --mode records --input "${OUT_ROOT}/inputs/missing_eof.fastq.gz" \
        --workers 3 --crc-check 1 > "${OUT_ROOT}/results/missing_eof_records.json"
    python3 - "${OUT_ROOT}/results/reference.json" \
        "${OUT_ROOT}/results/missing_eof_records.json" <<'PY'
import json
import sys
reference, observed = [json.load(open(path, encoding="utf-8")) for path in sys.argv[1:]]
assert observed["record_count"] == reference["record_count"]
assert observed["checksum64"] == reference["checksum64"]
PY
    if "${HARNESS}" --mode records --input "${OUT_ROOT}/inputs/mid_block.fastq.gz" \
        --workers 3 --crc-check 1 > /dev/null 2> "${OUT_ROOT}/results/mid_block.stderr"; then
        die "T5 mid-block truncation unexpectedly succeeded"
    fi
    grep -E "block offset [0-9]+" "${OUT_ROOT}/results/mid_block.stderr" >/dev/null \
        || die "T5 mid-block error did not include a block offset"
    if "${HARNESS}" --mode records --input "${OUT_ROOT}/inputs/overlong_line.fastq.gz" \
        --workers 3 --crc-check 1 > /dev/null \
        2> "${OUT_ROOT}/results/overlong_line.stderr"; then
        die "T5 overlong cross-member FASTQ line unexpectedly succeeded"
    fi
    grep -F "BGZF FASTQ line exceeds fixed Illumina capacity" \
        "${OUT_ROOT}/results/overlong_line.stderr" >/dev/null \
        || die "T5 overlong cross-member FASTQ line was not rejected during assembly"
    pass "T5 optional EOF and partial-member truncation"
fi

if enabled T8; then
    [[ -x "${HARNESS}" ]] || die "T8 requires ${HARNESS}"
    # A second mate with the same names ("/2"), different sequences, and a
    # different block size so the two files' block layouts are unrelated, plus
    # one variant with a single renamed record, one missing its last record,
    # and one whose optional FASTQ comments exceed the stored name-token bound.
    python3 - "${OUT_ROOT}/inputs/source.fastq" "${OUT_ROOT}/inputs" <<'PY'
import random
import sys
from pathlib import Path
lines = Path(sys.argv[1]).read_text(encoding="ascii").splitlines()
out = Path(sys.argv[2])
records = [lines[i:i + 4] for i in range(0, len(lines), 4)]
mate2 = [[rec[0].replace("/1 ", "/2 "), rec[1][::-1], rec[2], rec[3]] for rec in records]
renamed = [list(rec) for rec in mate2]
renamed[100][0] = "@read999999/2 instrument:run:100"
long_comment = [list(rec) for rec in mate2]
rng = random.Random(1904)
comment = "".join(rng.choices("abcdefghijklmnopqrstuvwxyz0123456789", k=1 << 20))
long_comment[100][0] = long_comment[100][0].split(maxsplit=1)[0] + " " + comment
variants = {
    "pair_r2": mate2,
    "pair_r2_renamed": renamed,
    "pair_r2_short": mate2[:-1],
    "pair_r2_long_comment": long_comment,
}
for name, recs in variants.items():
    (out / f"{name}.fastq").write_text(
        "\n".join(line for rec in recs for line in rec) + "\n", encoding="ascii")
PY
    for name in pair_r2 pair_r2_renamed pair_r2_short pair_r2_long_comment; do
        "${MAKE_FIXTURE}" --block-bytes 1013 "${OUT_ROOT}/inputs/${name}.fastq" \
            "${OUT_ROOT}/inputs/${name}.fastq.gz"
    done
    long_comment_compressed_size=$(wc -c < \
        "${OUT_ROOT}/inputs/pair_r2_long_comment.fastq.gz")
    (( long_comment_compressed_size > 64 * 1024 )) \
        || die "T8 long-comment fixture does not cross compressed work claims"
    "${HARNESS}" --mode pair --input "${OUT_ROOT}/inputs/blocked.fastq.gz" \
        --input2 "${OUT_ROOT}/inputs/pair_r2.fastq.gz" --workers 3 --crc-check 1 \
        > "${OUT_ROOT}/results/pair.json"
    python3 - "${OUT_ROOT}/results/pair.json" <<'PY'
import json
import sys
assert json.load(open(sys.argv[1], encoding="utf-8"))["record_count"] == 257
PY
    if "${HARNESS}" --mode pair --input "${OUT_ROOT}/inputs/blocked.fastq.gz" \
        --input2 "${OUT_ROOT}/inputs/pair_r2_renamed.fastq.gz" --workers 3 --crc-check 1 \
        > /dev/null 2> "${OUT_ROOT}/results/pair_renamed.stderr"; then
        die "T8 a renamed mate record unexpectedly paired"
    fi
    grep -F "read-name mismatch at record 100" "${OUT_ROOT}/results/pair_renamed.stderr" >/dev/null \
        || die "T8 read-name mismatch was not reported at record 100"
    "${HARNESS}" --mode pair --input "${OUT_ROOT}/inputs/blocked.fastq.gz" \
        --input2 "${OUT_ROOT}/inputs/pair_r2_renamed.fastq.gz" --workers 3 \
        --crc-check 1 --validate-read-names 0 \
        > "${OUT_ROOT}/results/pair_names_skipped.json"
    python3 - "${OUT_ROOT}/results/pair_names_skipped.json" <<'PY'
import json
import sys
result = json.load(open(sys.argv[1], encoding="utf-8"))
assert result["record_count"] == 257
assert result["mate1_name_bytes"] == 0
assert result["mate0_name_view_records"] > 240
assert result["sequence_view_records"] > 480
assert result["quality_view_records"] > 480
assert result["fastq_record_bytes"] <= 96
PY
    if "${HARNESS}" --mode pair --input "${OUT_ROOT}/inputs/blocked.fastq.gz" \
        --input2 "${OUT_ROOT}/inputs/pair_r2_short.fastq.gz" --workers 3 --crc-check 1 \
        > /dev/null 2> "${OUT_ROOT}/results/pair_short.stderr"; then
        die "T8 a short mate file unexpectedly paired"
    fi
    grep -F "record-count mismatch at record 256" "${OUT_ROOT}/results/pair_short.stderr" >/dev/null \
        || die "T8 record-count mismatch was not reported at record 256"
    "${HARNESS}" --mode pair --input "${OUT_ROOT}/inputs/blocked.fastq.gz" \
        --input2 "${OUT_ROOT}/inputs/pair_r2_long_comment.fastq.gz" \
        --workers 3 --crc-check 1 > "${OUT_ROOT}/results/pair_long_comment.json"
    python3 - "${OUT_ROOT}/results/pair_long_comment.json" <<'PY'
import json
import sys
assert json.load(open(sys.argv[1], encoding="utf-8"))["record_count"] == 257
PY
    pass "T8 paired-mate read-name and record-count validation"
fi

if enabled T6; then
    BGZF_E2E_CASE=T6 "${ROOT_DIR}/tests/bgzf/test_flex_e2e.sh"
    pass "T6 mixed BGZF/plain-gzip lanes"
fi

if enabled T9; then
    "${ROOT_DIR}/tests/bgzf/test_flex_fused_align.sh"
    pass "T9 fully-fused alignment mode drains alignQ without a reserved consumer"
fi

if enabled T10; then
    "${ROOT_DIR}/tests/bgzf/test_flex_sorted_bam_range.sh"
    pass "T10 paper-scoped Flex BGZF coordinate-sorted BAM parity and negative gates"
fi

run_t7() {
    local label="$1"
    local mode="$2"
    local star_bin="${STAR_BIN:-${ROOT_DIR}/core/legacy/source/STAR}"
    local wrapper="${OUT_ROOT}/STAR-${mode}"
    if (( CURRENT_PHASE >= 3 )); then
        python3 - "${wrapper}" "${star_bin}" "${mode}" <<'PY'
import os
import shlex
import sys
from pathlib import Path
wrapper, star, mode = sys.argv[1:]
Path(wrapper).write_text(
    "#!/usr/bin/env bash\nexec " + shlex.quote(star) +
    " --readFilesBgzfMode " + shlex.quote(mode) + " \"$@\"\n",
    encoding="utf-8")
os.chmod(wrapper, 0o755)
PY
        star_bin="${wrapper}"
    fi
    STAR_BIN="${star_bin}" OUT_ROOT="${OUT_ROOT}/t7_${mode}_fastx" \
        "${ROOT_DIR}/tests/run_fastx_contract_star_smoke.sh"
    if [[ "${BGZF_T7_FULL}" == 1 ]]; then
        STAR_BIN="${star_bin}" OUT_ROOT="${OUT_ROOT}/t7_${mode}_cbq" \
            "${ROOT_DIR}/tests/run_cbq_e2e_module_regression.sh"
    else
        BQTOOLS="${BQTOOLS:-bqtools}" OUT_ROOT="${OUT_ROOT}/t7_${mode}_cbq_reader" \
            "${ROOT_DIR}/tests/run_cbq_cpp_reader_smoke.sh"
    fi
    pass "${label} gzip/CBQ regressions"
}

if enabled T7-off; then
    run_t7 T7-off off
fi
if enabled T7-auto; then
    run_t7 T7-auto auto
fi

echo "PASS: BGZF ingest phase ${CURRENT_PHASE} gate completed (${OUT_ROOT})"
