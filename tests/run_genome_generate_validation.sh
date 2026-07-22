#!/usr/bin/env bash
set -euo pipefail
# Validate STAR genomeGenerate produces a correct and complete index.
#
# Builds a genome index from a synthetic reference (two small chromosomes
# plus a minimal GTF) and validates the resulting index structure.
#
# Tier A: self-contained, no external data, container-safe.
#
# Usage:
#   tests/run_genome_generate_validation.sh [options]
#
# Options:
#   --threads N    STAR threads (default: 2)
#   --outdir DIR   Output directory
#   -h, --help     Show help

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
THREADS="${THREADS:-2}"
OUTDIR=""

die() { echo "FAIL: $*" >&2; exit 1; }
log() { printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*"; }
usage() { sed -n '2,/^$/s/^# \?//p' "$0"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads) THREADS="$2"; shift 2 ;;
    --outdir)  OUTDIR="$2";  shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

TIMESTAMP="$(date +'%Y%m%d_%H%M%S')"
OUTDIR="${OUTDIR:-/tmp/genome_generate_validation_${TIMESTAMP}_$$}"
mkdir -p "${OUTDIR}"

[[ -x "${STAR_BIN}" ]] || die "STAR binary not found: ${STAR_BIN}"

log "=== genomeGenerate Validation ==="
log "STAR:    ${STAR_BIN}"
log "Threads: ${THREADS}"
log "Outdir:  ${OUTDIR}"

PASS=0
FAIL=0

check() {
  local desc="$1" path="$2"
  if [[ -e "${path}" ]]; then
    log "  PASS: ${desc}"
    ((++PASS))
  else
    log "  FAIL: ${desc} — not found: ${path}"
    ((++FAIL))
  fi
}

check_nonempty() {
  local desc="$1" path="$2"
  if [[ -s "${path}" ]]; then
    log "  PASS: ${desc}"
    ((++PASS))
  else
    log "  FAIL: ${desc} — missing or empty: ${path}"
    ((++FAIL))
  fi
}

# ── Generate synthetic reference ──────────────────────────────────────
REF_DIR="${OUTDIR}/ref"
mkdir -p "${REF_DIR}"

# Two chromosomes: chr1 (500bp) and chrY (300bp)
python3 - "${REF_DIR}" << 'PYGEN'
import sys, random, os

ref_dir = sys.argv[1]
random.seed(42)
bases = "ACGT"

# chr1: 500bp
chr1_seq = "".join(random.choice(bases) for _ in range(500))
# chrY: 300bp
chrY_seq = "".join(random.choice(bases) for _ in range(300))

fasta_path = os.path.join(ref_dir, "genome.fa")
with open(fasta_path, "w") as f:
    f.write(">chr1\n")
    for i in range(0, len(chr1_seq), 80):
        f.write(chr1_seq[i:i+80] + "\n")
    f.write(">chrY\n")
    for i in range(0, len(chrY_seq), 80):
        f.write(chrY_seq[i:i+80] + "\n")

# Minimal GTF with 2 genes (one per chrom)
gtf_path = os.path.join(ref_dir, "genes.gtf")
with open(gtf_path, "w") as f:
    f.write('chr1\tsynthetic\tgene\t50\t200\t.\t+\t.\tgene_id "GENE1"; gene_name "Gene1";\n')
    f.write('chr1\tsynthetic\ttranscript\t50\t200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "Gene1";\n')
    f.write('chr1\tsynthetic\texon\t50\t200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "Gene1";\n')
    f.write('chrY\tsynthetic\tgene\t30\t150\t.\t+\t.\tgene_id "GENEY"; gene_name "GeneY";\n')
    f.write('chrY\tsynthetic\ttranscript\t30\t150\t.\t+\t.\tgene_id "GENEY"; transcript_id "TXY"; gene_name "GeneY";\n')
    f.write('chrY\tsynthetic\texon\t30\t150\t.\t+\t.\tgene_id "GENEY"; transcript_id "TXY"; gene_name "GeneY";\n')

print(f"Reference: {fasta_path} (chr1=500bp, chrY=300bp)")
print(f"GTF: {gtf_path} (2 genes)")
PYGEN

# ── Run genomeGenerate ────────────────────────────────────────────────
INDEX_DIR="${OUTDIR}/star_index"
mkdir -p "${INDEX_DIR}"

log "Running genomeGenerate..."
"${STAR_BIN}" \
  --runMode genomeGenerate \
  --genomeDir "${INDEX_DIR}" \
  --runThreadN "${THREADS}" \
  --genomeFastaFiles "${REF_DIR}/genome.fa" \
  --sjdbGTFfile "${REF_DIR}/genes.gtf" \
  --genomeSAindexNbases 3 \
  > "${OUTDIR}/genomeGenerate.log" 2>&1 \
  || die "genomeGenerate failed (log: ${OUTDIR}/genomeGenerate.log)"

log "genomeGenerate completed."

# ── Validate index structure ──────���───────────────────────────────────
log "--- Validating index files ---"

check_nonempty "Genome"              "${INDEX_DIR}/Genome"
check_nonempty "SA"                  "${INDEX_DIR}/SA"
check_nonempty "SAindex"             "${INDEX_DIR}/SAindex"
check_nonempty "chrName.txt"         "${INDEX_DIR}/chrName.txt"
check_nonempty "chrLength.txt"       "${INDEX_DIR}/chrLength.txt"
check_nonempty "chrStart.txt"        "${INDEX_DIR}/chrStart.txt"
check_nonempty "chrNameLength.txt"   "${INDEX_DIR}/chrNameLength.txt"
check          "genomeParameters.txt" "${INDEX_DIR}/genomeParameters.txt"
# sjdbList.fromGTF.out.tab may be empty on tiny synthetic genomes with no
# introns.  Check existence only.
check          "sjdbList.fromGTF.out.tab" "${INDEX_DIR}/sjdbList.fromGTF.out.tab"
check_nonempty "sjdbInfo.txt"        "${INDEX_DIR}/sjdbInfo.txt"
check_nonempty "geneInfo.tab"        "${INDEX_DIR}/geneInfo.tab"

# ── Validate chrName.txt content ──────────────────────────────────────
log "--- Validating chromosome names ---"

if grep -qx "chr1" "${INDEX_DIR}/chrName.txt" 2>/dev/null; then
  log "  PASS: chrName.txt contains chr1"
  ((++PASS))
else
  log "  FAIL: chrName.txt missing chr1"
  ((++FAIL))
fi

if grep -qx "chrY" "${INDEX_DIR}/chrName.txt" 2>/dev/null; then
  log "  PASS: chrName.txt contains chrY"
  ((++PASS))
else
  log "  FAIL: chrName.txt missing chrY"
  ((++FAIL))
fi

NCHR=$(wc -l < "${INDEX_DIR}/chrName.txt")
if [[ "${NCHR}" -eq 2 ]]; then
  log "  PASS: chrName.txt has exactly 2 chromosomes"
  ((++PASS))
else
  log "  FAIL: chrName.txt has ${NCHR} chromosomes (expected 2)"
  ((++FAIL))
fi

# ── Validate chrLength.txt content ────────────────────────────────────
log "--- Validating chromosome lengths ---"

CHR1_LEN=$(sed -n '1p' "${INDEX_DIR}/chrLength.txt")
CHRY_LEN=$(sed -n '2p' "${INDEX_DIR}/chrLength.txt")

if [[ "${CHR1_LEN}" -eq 500 ]]; then
  log "  PASS: chr1 length = 500"
  ((++PASS))
else
  log "  FAIL: chr1 length = ${CHR1_LEN} (expected 500)"
  ((++FAIL))
fi

if [[ "${CHRY_LEN}" -eq 300 ]]; then
  log "  PASS: chrY length = 300"
  ((++PASS))
else
  log "  FAIL: chrY length = ${CHRY_LEN} (expected 300)"
  ((++FAIL))
fi

# ── Validate geneInfo.tab content ───────��─────────────────────────────
log "--- Validating gene annotation ---"

NGENES=$(tail -n +2 "${INDEX_DIR}/geneInfo.tab" | wc -l)
if [[ "${NGENES}" -eq 2 ]]; then
  log "  PASS: geneInfo.tab has 2 genes"
  ((++PASS))
else
  log "  FAIL: geneInfo.tab has ${NGENES} genes (expected 2)"
  ((++FAIL))
fi

if grep -q "GENE1" "${INDEX_DIR}/geneInfo.tab" 2>/dev/null; then
  log "  PASS: geneInfo.tab contains GENE1"
  ((++PASS))
else
  log "  FAIL: geneInfo.tab missing GENE1"
  ((++FAIL))
fi

# ── Validate genomeParameters.txt ─���─────────────────────────────��─────
log "--- Validating genome parameters ---"

if grep -q "genomeSAindexNbases.*3" "${INDEX_DIR}/genomeParameters.txt" 2>/dev/null; then
  log "  PASS: genomeSAindexNbases = 3 in parameters"
  ((++PASS))
else
  log "  FAIL: genomeSAindexNbases not found or wrong in parameters"
  ((++FAIL))
fi

GSIZE=$(grep "^versionGenome" "${INDEX_DIR}/genomeParameters.txt" 2>/dev/null || true)
if [[ -n "${GSIZE}" ]]; then
  log "  PASS: versionGenome present in parameters"
  ((++PASS))
else
  log "  FAIL: versionGenome missing from parameters"
  ((++FAIL))
fi

# ── Quick alignment sanity check ───────���──────────────────────────────
# Generate a few reads that should map, verify STAR can load the index
log "--- Quick alignment test ---"

ALIGN_DIR="${OUTDIR}/align_test"
mkdir -p "${ALIGN_DIR}"

python3 - "${REF_DIR}/genome.fa" "${ALIGN_DIR}" << 'PYGEN'
import sys, os

ref_fa = sys.argv[1]
align_dir = sys.argv[2]

# Read chr1 sequence
seq = ""
reading = False
with open(ref_fa) as f:
    for line in f:
        if line.startswith(">chr1"):
            reading = True
            continue
        elif line.startswith(">"):
            break
        if reading:
            seq += line.strip()

# Write a few reads from known positions
with open(os.path.join(align_dir, "test.fastq"), "w") as f:
    for i in range(5):
        pos = 50 + i * 20
        read_seq = seq[pos:pos+50]
        f.write(f"@read{i}\n{read_seq}\n+\n{'I'*50}\n")
PYGEN

"${STAR_BIN}" \
  --runMode alignReads \
  --genomeDir "${INDEX_DIR}" \
  --readFilesIn "${ALIGN_DIR}/test.fastq" \
  --outFileNamePrefix "${ALIGN_DIR}/" \
  --outSAMtype None \
  --runThreadN 1 \
  > "${ALIGN_DIR}/align.log" 2>&1 \
  || die "alignment test failed (log: ${ALIGN_DIR}/align.log)"

check "Log.final.out from alignment" "${ALIGN_DIR}/Log.final.out"

# Check that some reads mapped
if grep -q "Uniquely mapped reads number" "${ALIGN_DIR}/Log.final.out" 2>/dev/null; then
  MAPPED=$(grep "Uniquely mapped reads number" "${ALIGN_DIR}/Log.final.out" | awk '{print $NF}')
  if [[ "${MAPPED}" -gt 0 ]]; then
    log "  PASS: ${MAPPED} reads mapped to synthetic index"
    ((++PASS))
  else
    log "  FAIL: 0 reads mapped (expected >0)"
    ((++FAIL))
  fi
else
  log "  FAIL: could not parse mapping stats"
  ((++FAIL))
fi

# ── Summary ────��──────────────────────────────────────────────────────
log "--- Results: ${PASS} passed, ${FAIL} failed ---"

if [[ "${FAIL}" -gt 0 ]]; then
  die "${FAIL} validation(s) failed. Outputs at: ${OUTDIR}"
fi

log "PASS: genomeGenerate validation"
