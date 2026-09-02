#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-${SCRIPT_DIR}/../core/legacy/source/STAR}"
# Support injecting extra args (e.g., --defaultCoreScrna Yes) via STAR_EXTRA_ARGS
STAR_EXTRA_ARGS="${STAR_EXTRA_ARGS:-}"
TEST_DIR="${SCRIPT_DIR}/solo_smoke"

# Ensure samtools exists
command -v samtools >/dev/null || { echo "ERROR: samtools not found"; exit 1; }

# Ensure STAR binary exists
[ -f "$STAR_BIN" ] || { echo "ERROR: STAR binary not found at $STAR_BIN"; exit 1; }
REF_DIR="${TEST_DIR}/ref"
IDX_DIR="${REF_DIR}/star_index"
FASTQ_DIR="${TEST_DIR}/fastq"
OUT_DIR="${TEST_DIR}/output"
TMP_DIR="${TEST_DIR}/tmp"
WL="${TEST_DIR}/whitelist.txt"

rm -rf "$TEST_DIR"
mkdir -p "$REF_DIR" "$IDX_DIR" "$FASTQ_DIR" "$OUT_DIR"

# Minimal reference (single chr)
cat > "$REF_DIR/chr1.fa" << 'EOF'
>chr1
AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT
EOF

cat > "$REF_DIR/genes.gtf" << 'EOF'
chr1	.	gene	1	60	.	+	.	gene_id "GENE1";
chr1	.	exon	1	60	.	+	.	gene_id "GENE1"; transcript_id "T1";
EOF

# Whitelist (two barcodes)
cat > "$WL" << 'EOF'
ACGTACGTACGTACGT
TGCATGCATGCATGCA
EOF

# Generate R1 (CB+UMI) and R2 (cDNA)
CB="ACGTACGTACGTACGT"
UMI="AAAAAAAAAAAA"
READ1="${CB}${UMI}"
READ2="AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT"

q1=$(printf 'I%.0s' $(seq 1 ${#READ1}))
q2=$(printf 'I%.0s' $(seq 1 ${#READ2}))

cat > "$FASTQ_DIR/R1.fastq" << EOF
@read1/1
$READ1
+
$q1
EOF

cat > "$FASTQ_DIR/R2.fastq" << EOF
@read1/2
$READ2
+
$q2
EOF

# Build STAR index
"$STAR_BIN" \
  --runMode genomeGenerate \
  --genomeDir "$IDX_DIR" \
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

# Minimal Solo output check
[ -f "$OUT_DIR/Solo.out/Gene/raw/matrix.mtx" ] || {
  echo "FAIL: Solo matrix.mtx missing"
  exit 1
}

# Plain STARsolo must not build the Flex-only one-mismatch barcode corrector.
# With a production 3M whitelist that precompute costs several GB and tens of
# seconds without contributing to the legacy matrix path.
grep -F "CbCorrector not required: inlineCBCorrection=0, inlineHashMode=0" \
  "$OUT_DIR/Log.out" >/dev/null || {
  echo "FAIL: plain STARsolo unexpectedly enabled the inline barcode corrector"
  exit 1
}

echo "OK: Solo smoke test passed"
