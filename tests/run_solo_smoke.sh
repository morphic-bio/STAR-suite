#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-${SCRIPT_DIR}/../core/legacy/source/STAR}"
# Support injecting extra args (e.g., --defaultCoreScrna Yes) via STAR_EXTRA_ARGS
STAR_EXTRA_ARGS="${STAR_EXTRA_ARGS:-}"
TEST_DIR="${SOLO_SMOKE_OUTDIR:-${SCRIPT_DIR}/solo_smoke}"

# Ensure samtools exists
command -v samtools >/dev/null || { echo "ERROR: samtools not found"; exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 not found"; exit 1; }

# Ensure STAR binary exists
[ -f "$STAR_BIN" ] || { echo "ERROR: STAR binary not found at $STAR_BIN"; exit 1; }
REF_DIR="${TEST_DIR}/ref"
IDX_DIR="${REF_DIR}/star_index"
FASTQ_DIR="${TEST_DIR}/fastq"
OUT_DIR="${TEST_DIR}/output"
TMP_DIR="${TEST_DIR}/tmp"
WL="${TEST_DIR}/whitelist.txt"

if [[ -n "${SOLO_SMOKE_OUTDIR:-}" && -e "$TEST_DIR" ]]; then
  echo "ERROR: use a fresh SOLO_SMOKE_OUTDIR: $TEST_DIR" >&2
  exit 1
fi
rm -rf "$TEST_DIR"
mkdir -p "$REF_DIR" "$IDX_DIR" "$FASTQ_DIR" "$OUT_DIR"

# Unique mapping and known PCR duplicates make an empty MEX a real failure.
python3 - "$TEST_DIR" <<'PY'
from pathlib import Path
import random
import sys

root = Path(sys.argv[1])
rng = random.Random(20260927)
genome = "".join(rng.choice("ACGT") for _ in range(4000))
(root / "ref/chr1.fa").write_text(">chr1\n" + genome + "\n")
(root / "ref/genes.gtf").write_text(
    'chr1\ttest\texon\t1001\t2000\t.\t+\t.\tgene_id "GENE1"; transcript_id "T1";\n'
)
cb = "ACGTACGTACGTACGT"
(root / "whitelist.txt").write_text(cb + "\nTGCATGCATGCATGCA\n")
with (root / "fastq/R1.fastq").open("w") as r1, (root / "fastq/R2.fastq").open("w") as r2:
    for i, umi in enumerate(("ACGTACGTACGT", "ACGTACGTACGT", "TTAGCGATCGGA")):
        for mate, stream, seq in ((1, r1, cb + umi), (2, r2, genome[1100:1200])):
            stream.write(f"@read{i}/{mate}\n{seq}\n+\n{'I' * len(seq)}\n")
PY

# Build STAR index
"$STAR_BIN" \
  --runMode genomeGenerate \
  --genomeDir "$IDX_DIR" \
  --outFileNamePrefix "$REF_DIR/" \
  --genomeFastaFiles "$REF_DIR/chr1.fa" \
  --sjdbGTFfile "$REF_DIR/genes.gtf" \
  --sjdbOverhang 20 \
  --genomeSAindexNbases 4 >/dev/null

# Run STARsolo with samtools sorter
# shellcheck disable=SC2086
"$STAR_BIN" ${STAR_EXTRA_ARGS:-} \
  --runThreadN 2 \
  --genomeDir "$IDX_DIR" \
  --readFilesIn "$FASTQ_DIR/R2.fastq" "$FASTQ_DIR/R1.fastq" \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 --soloUMIlen 12 \
  --soloCBstart 1 --soloUMIstart 17 \
  --soloCBwhitelist "$WL" \
  --outSAMtype BAM SortedByCoordinate \
  --outBAMsortMethod samtools \
  --outFileNamePrefix "$OUT_DIR/" \
  --outTmpDir "$TMP_DIR"

# Checks
samtools quickcheck "$OUT_DIR/Aligned.sortedByCoord.out.bam"
samtools view -H "$OUT_DIR/Aligned.sortedByCoord.out.bam" | grep "SO:coordinate" >/dev/null
[ "$(samtools view -c -F 260 "$OUT_DIR/Aligned.sortedByCoord.out.bam")" -eq 3 ] || {
  echo "FAIL: expected three mapped primary reads" >&2
  exit 1
}

# Matrix presence alone does not prove UMI counts survived deduplication.
python3 - "$SCRIPT_DIR" "$OUT_DIR/Solo.out/Gene/raw" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
from compare_flex_hash_screen_mex import read_mex

_, _, counts = read_mex(sys.argv[2])
observed = {(gene, cb.removesuffix("-1")): count for (gene, cb), count in counts.items()}
expected = {("GENE1", "ACGTACGTACGTACGT"): 2}
if observed != expected:
    raise SystemExit(f"FAIL: default Solo counts {observed}; expected {expected}")
print("OK: three reads collapse to exactly two GEX molecules")
PY

# Plain STARsolo must not build the Flex-only one-mismatch barcode corrector.
# With a production 3M whitelist that precompute costs several GB and tens of
# seconds without contributing to the legacy matrix path.
grep -F "CbCorrector not required: inlineCBCorrection=0, inlineHashMode=0" \
  "$OUT_DIR/Log.out" >/dev/null || {
  echo "FAIL: plain STARsolo unexpectedly enabled the inline barcode corrector"
  exit 1
}

echo "OK: Solo smoke test passed"
