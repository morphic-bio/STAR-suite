#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${FASTX_INPUT_HARNESS_BIN:-${REPO_ROOT}/core/legacy/source/fastx_input_harness}"
OUT_DIR="${FASTX_INPUT_HARNESS_SMOKE_OUTDIR:-/tmp/star_suite_fastx_input_harness_smoke}"

rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}/inputs" "${OUT_DIR}/dumps"

if [[ ! -x "${BIN}" ]]; then
  make -C "${REPO_ROOT}/core/legacy/source" fastx-input-harness
fi

cat > "${OUT_DIR}/inputs/lane1_R1.fastq" <<'FASTQ'
@readA/1 1:N:0:ACGT extraA
ACGTAC
+
IIIIII
@readB/1 1:Y:0:ACGT extraB
TTGGCA
+
HHHHHH
FASTQ

cat > "${OUT_DIR}/inputs/lane1_R2.fastq" <<'FASTQ'
@readA/2 2:N:0:ACGT extraA
TGCATG
+
IIIIII
@readB/2 2:Y:0:ACGT extraB
AACCGT
+
HHHHHH
FASTQ

cat > "${OUT_DIR}/inputs/lane1_R3.fastq" <<'FASTQ'
@readA/3 3:N:0:ACGT extraA
GGTTGG
+
IIIIII
@readB/3 3:Y:0:ACGT extraB
CCAACC
+
HHHHHH
FASTQ

cat > "${OUT_DIR}/inputs/lane2_R1.fastq" <<'FASTQ'
@readC/1 1:N:0:TGCA extraC
CCCCAA
+
FFFFFF
FASTQ

cat > "${OUT_DIR}/inputs/lane2_R2.fastq" <<'FASTQ'
@readC/2 2:N:0:TGCA extraC
GGGGTT
+
FFFFFF
FASTQ

gzip -c "${OUT_DIR}/inputs/lane1_R1.fastq" > "${OUT_DIR}/inputs/lane1_R1.fastq.gz"
gzip -c "${OUT_DIR}/inputs/lane1_R2.fastq" > "${OUT_DIR}/inputs/lane1_R2.fastq.gz"
gzip -c "${OUT_DIR}/inputs/lane1_R3.fastq" > "${OUT_DIR}/inputs/lane1_R3.fastq.gz"
gzip -c "${OUT_DIR}/inputs/lane2_R1.fastq" > "${OUT_DIR}/inputs/lane2_R1.fastq.gz"
gzip -c "${OUT_DIR}/inputs/lane2_R2.fastq" > "${OUT_DIR}/inputs/lane2_R2.fastq.gz"

cat > "${OUT_DIR}/inputs/truncated_R2.fastq" <<'FASTQ'
@readA/2 2:N:0:ACGT
TGCATG
+
FASTQ

# Enough records to fill both bounded lane queues when the consumer exits
# after its first record. This exercises cancellation of blocked producers.
: > "${OUT_DIR}/inputs/busy_R1.fastq"
: > "${OUT_DIR}/inputs/busy_R2.fastq"
for ((ii = 0; ii < 20; ++ii)); do
  printf '@busy%03d/1 1:N:0:ACGT\nACGTAC\n+\nIIIIII\n' "$ii" >> "${OUT_DIR}/inputs/busy_R1.fastq"
  printf '@busy%03d/2 2:N:0:ACGT\nTGCATG\n+\nIIIIII\n' "$ii" >> "${OUT_DIR}/inputs/busy_R2.fastq"
done

cat > "${OUT_DIR}/manifest.tsv" <<EOF_MANIFEST
${OUT_DIR}/inputs/lane1_R1.fastq.gz	${OUT_DIR}/inputs/lane1_R2.fastq.gz	ID:lane1
${OUT_DIR}/inputs/lane2_R1.fastq.gz	${OUT_DIR}/inputs/lane2_R2.fastq.gz	ID:lane2
EOF_MANIFEST

"${BIN}" --readNameSeparator / \
  --readFilesIn "${OUT_DIR}/inputs/lane1_R1.fastq" "${OUT_DIR}/inputs/lane1_R2.fastq" \
  > "${OUT_DIR}/dumps/plain_pair.tsv"

"${BIN}" --readNameSeparator / \
  --readFilesIn "${OUT_DIR}/inputs/lane1_R1.fastq.gz" "${OUT_DIR}/inputs/lane1_R2.fastq.gz" \
  > "${OUT_DIR}/dumps/gzip_internal_pair.tsv"

"${BIN}" --readNameSeparator / \
  --readFilesIn "${OUT_DIR}/inputs/lane1_R1.fastq.gz" "${OUT_DIR}/inputs/lane1_R2.fastq.gz" \
  --readFilesCommand zcat \
  > "${OUT_DIR}/dumps/gzip_zcat_pair.tsv"

"${BIN}" --readNameSeparator / \
  --readFilesIn "${OUT_DIR}/inputs/lane1_R1.fastq,${OUT_DIR}/inputs/lane2_R1.fastq" \
                "${OUT_DIR}/inputs/lane1_R2.fastq,${OUT_DIR}/inputs/lane2_R2.fastq" \
  > "${OUT_DIR}/dumps/comma_lanes_pair.tsv"

"${BIN}" --producer-pool --readNameSeparator / \
  --readFilesIn "${OUT_DIR}/inputs/lane1_R1.fastq,${OUT_DIR}/inputs/lane2_R1.fastq" \
                "${OUT_DIR}/inputs/lane1_R2.fastq,${OUT_DIR}/inputs/lane2_R2.fastq" \
  > "${OUT_DIR}/dumps/comma_lanes_pair_pool.tsv"

"${BIN}" --producer-pool --readNameSeparator / \
  --readFilesIn "${OUT_DIR}/inputs/lane1_R1.fastq.gz,${OUT_DIR}/inputs/lane2_R1.fastq.gz" \
                "${OUT_DIR}/inputs/lane1_R2.fastq.gz,${OUT_DIR}/inputs/lane2_R2.fastq.gz" \
  > "${OUT_DIR}/dumps/comma_lanes_gzip_pool.tsv"

"${BIN}" --readNameSeparator / \
  --readFilesManifest "${OUT_DIR}/manifest.tsv" \
  > "${OUT_DIR}/dumps/manifest_gzip_pair.tsv"

"${BIN}" --readNameSeparator / \
  --readFilesIn "${OUT_DIR}/inputs/lane1_R1.fastq" \
  > "${OUT_DIR}/dumps/plain_single.tsv"

"${BIN}" --readNameSeparator / \
  --readFilesIn "${OUT_DIR}/inputs/lane1_R1.fastq" \
                "${OUT_DIR}/inputs/lane1_R2.fastq" \
                "${OUT_DIR}/inputs/lane1_R3.fastq" \
  > "${OUT_DIR}/dumps/plain_three_mate.tsv"

"${BIN}" --readNameSeparator / --dump-fastq \
  --readFilesIn "${OUT_DIR}/inputs/lane1_R1.fastq.gz" "${OUT_DIR}/inputs/lane1_R2.fastq.gz" \
  > "${OUT_DIR}/dumps/debug.fastq"

diff -u "${OUT_DIR}/dumps/plain_pair.tsv" "${OUT_DIR}/dumps/gzip_internal_pair.tsv"
diff -u "${OUT_DIR}/dumps/plain_pair.tsv" "${OUT_DIR}/dumps/gzip_zcat_pair.tsv"
diff -u "${OUT_DIR}/dumps/comma_lanes_pair.tsv" "${OUT_DIR}/dumps/manifest_gzip_pair.tsv"
diff -u "${OUT_DIR}/dumps/comma_lanes_pair.tsv" "${OUT_DIR}/dumps/comma_lanes_pair_pool.tsv"
diff -u "${OUT_DIR}/dumps/comma_lanes_pair.tsv" "${OUT_DIR}/dumps/comma_lanes_gzip_pool.tsv"

# Closing a partially consumed pool must wake and join producers even when
# later-lane queues are full.
timeout 10 "${BIN}" --producer-pool --max-records 1 --readNameSeparator / \
  --readFilesIn "${OUT_DIR}/inputs/busy_R1.fastq,${OUT_DIR}/inputs/busy_R1.fastq" \
                "${OUT_DIR}/inputs/busy_R2.fastq,${OUT_DIR}/inputs/busy_R2.fastq" \
  > "${OUT_DIR}/dumps/pool_early_close.tsv"

# A producer-side parse error must propagate to the consumer without a hang.
if timeout 10 "${BIN}" --producer-pool --readNameSeparator / \
  --readFilesIn "${OUT_DIR}/inputs/lane1_R1.fastq" "${OUT_DIR}/inputs/truncated_R2.fastq" \
  > "${OUT_DIR}/dumps/pool_truncated.tsv" 2> "${OUT_DIR}/dumps/pool_truncated.err"; then
  echo "ERROR: truncated producer input unexpectedly succeeded" >&2
  exit 1
fi
grep -F "FASTQ record is missing quality line" "${OUT_DIR}/dumps/pool_truncated.err" >/dev/null

python3 - "${OUT_DIR}/dumps" <<'PY'
import pathlib
import sys

dump_dir = pathlib.Path(sys.argv[1])

def rows(name):
    with open(dump_dir / name, "r", encoding="utf-8") as handle:
        return [line.rstrip("\n").split("\t") for line in handle if line.strip()]

plain = rows("plain_pair.tsv")
assert len(plain) == 4, plain
assert {row[0] for row in plain} == {"0"}, plain
assert {row[3] for row in plain} == {"readA", "readB"}, plain
assert {row[4] for row in plain} == {"1", "2"}, plain
assert {row[2] for row in plain} == {"N", "Y"}, plain

lanes = rows("comma_lanes_pair.tsv")
assert len(lanes) == 6, lanes
assert [row[0] for row in lanes] == ["0", "0", "0", "0", "1", "1"], lanes
assert {row[3] for row in lanes} == {"readA", "readB", "readC"}, lanes

single = rows("plain_single.tsv")
assert len(single) == 2, single
assert {row[4] for row in single} == {"1"}, single

three_mate = rows("plain_three_mate.tsv")
assert len(three_mate) == 6, three_mate
assert {row[4] for row in three_mate} == {"1", "2", "3"}, three_mate

debug_fastq = dump_dir / "debug.fastq"
assert debug_fastq.stat().st_size > 0, debug_fastq
PY

echo "PASS: fastx input harness smoke (${OUT_DIR})"
