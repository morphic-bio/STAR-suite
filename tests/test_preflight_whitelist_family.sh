#!/usr/bin/env bash
# Smoke test for scripts/preflight_whitelist_family.py.

set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT="${ROOT_DIR}/scripts/preflight_whitelist_family.py"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "${TMPDIR}"' EXIT

python3 - <<'PY' "${TMPDIR}"
import gzip
import sys
from pathlib import Path

root = Path(sys.argv[1])
(root / "fastqs").mkdir()

(root / "family_a_TRU.txt").write_text("AAAAAAAAAAAAAAAA\nCCCCCCCCCCCCCCCC\n")
(root / "family_a_NXT.txt").write_text("AAAACCCCGGGGTTTT\n")
(root / "family_b_TRU.txt").write_text("TTTTTTTTTTTTTTTT\nGGGGGGGGGGGGGGGG\n")
(root / "family_b_NXT.txt").write_text("TGCATGCATGCATGCA\n")

with gzip.open(root / "fastqs" / "sample_L001_R1_001.fastq.gz", "wt") as handle:
    for i in range(100):
        bc = "TTTTTTTTTTTTTTTT" if i < 90 else "AAAAAAAAAAAAAAAA"
        seq = bc + "ACGTACGTACGT"
        handle.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")

manifest_header = (
    "provider_group\tstage\tlibrary\tchemistry\tfeature_ref\twhitelist\t"
    "fastq_root\tfastq_sample_ids\tlanes_by_sample_id\tfastq_roles\tfile_count\tnotes\n"
)
common_tail = (
    f"\t{root / 'fastqs'}\tsample\tsample:L001\tI1,I2,R1,R2\t1\ttest\n"
)
(root / "manifest_pass.tsv").write_text(
    manifest_header
    + f"SAMPLE\tS0\tGEX\tTRU\tNA\t{root / 'family_b_TRU.txt'}"
    + common_tail
)
(root / "manifest_fail.tsv").write_text(
    manifest_header
    + f"SAMPLE\tS0\tGEX\tTRU\tNA\t{root / 'family_a_TRU.txt'}"
    + common_tail
)
PY

COMMON_ARGS=(
  --whitelist "family_a:TRU:${TMPDIR}/family_a_TRU.txt"
  --whitelist "family_a:NXT:${TMPDIR}/family_a_NXT.txt"
  --whitelist "family_b:TRU:${TMPDIR}/family_b_TRU.txt"
  --whitelist "family_b:NXT:${TMPDIR}/family_b_NXT.txt"
  --sample-reads 100
  --max-fastqs-per-manifest-row 1
)

python3 "${SCRIPT}" \
  --manifest "${TMPDIR}/manifest_pass.tsv" \
  "${COMMON_ARGS[@]}" \
  --outdir "${TMPDIR}/pass" >/dev/null

grep -F $'SAMPLE:GEX\t100\tfamily_b:TRU' "${TMPDIR}/pass/whitelist_family_summary.tsv" >/dev/null
grep -F $'\tPASS\t' "${TMPDIR}/pass/whitelist_family_summary.tsv" >/dev/null

if python3 "${SCRIPT}" \
    --manifest "${TMPDIR}/manifest_fail.tsv" \
    "${COMMON_ARGS[@]}" \
    --outdir "${TMPDIR}/fail" >/dev/null 2>&1; then
  echo "Expected wrong-family manifest to fail" >&2
  exit 1
fi

grep -F "FAIL_EXPECTED_MISMATCH" "${TMPDIR}/fail/whitelist_family_summary.tsv" >/dev/null

echo "preflight_whitelist_family smoke test passed"
