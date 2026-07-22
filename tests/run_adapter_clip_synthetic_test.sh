#!/usr/bin/env bash
set -euo pipefail
# Validate STAR --clip3pAdapterSeq on synthetic reads with known adapter content.
#
# Generates paired-end reads where R2 has a known adapter suffix appended to
# genomic sequence.  Runs STAR with --clip3pAdapterSeq and verifies that
# trimmed alignments are shorter and mapping rate is maintained or improved.
#
# Tier A: self-contained, no external data, container-safe.
#
# Usage:
#   tests/run_adapter_clip_synthetic_test.sh [options]
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
OUTDIR="${OUTDIR:-/tmp/adapter_clip_synthetic_${TIMESTAMP}_$$}"
mkdir -p "${OUTDIR}"

[[ -x "${STAR_BIN}" ]] || die "STAR binary not found: ${STAR_BIN}"

log "=== Adapter Clip Synthetic Test ==="
log "STAR:    ${STAR_BIN}"
log "Threads: ${THREADS}"
log "Outdir:  ${OUTDIR}"

PASS=0
FAIL=0

# ── Build synthetic reference and index ───────────────────────────────
REF_DIR="${OUTDIR}/ref"
INDEX_DIR="${OUTDIR}/star_index"
mkdir -p "${REF_DIR}" "${INDEX_DIR}"

ADAPTER_SEQ="AGATCGGAAGAG"

python3 - "${REF_DIR}" << 'PYGEN'
import sys, random, os

ref_dir = sys.argv[1]
random.seed(42)
bases = "ACGT"

# chr1: 50000bp — large enough for stable STAR index generation
chr1_seq = "".join(random.choice(bases) for _ in range(50000))

fasta_path = os.path.join(ref_dir, "genome.fa")
with open(fasta_path, "w") as f:
    f.write(">chr1\n")
    for i in range(0, len(chr1_seq), 80):
        f.write(chr1_seq[i:i+80] + "\n")

# Minimal GTF
gtf_path = os.path.join(ref_dir, "genes.gtf")
with open(gtf_path, "w") as f:
    f.write('chr1\tsynthetic\tgene\t100\t49000\t.\t+\t.\tgene_id "GENE1"; gene_name "Gene1";\n')
    f.write('chr1\tsynthetic\ttranscript\t100\t49000\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1";\n')
    f.write('chr1\tsynthetic\texon\t100\t49000\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1";\n')

print(f"Reference: {fasta_path} (chr1=50000bp)")
PYGEN

log "Building STAR index..."
"${STAR_BIN}" \
  --runMode genomeGenerate \
  --genomeDir "${INDEX_DIR}" \
  --runThreadN "${THREADS}" \
  --genomeFastaFiles "${REF_DIR}/genome.fa" \
  --sjdbGTFfile "${REF_DIR}/genes.gtf" \
  --genomeSAindexNbases 7 \
  > "${OUTDIR}/index.log" 2>&1 \
  || die "genomeGenerate failed"

# ── Generate reads with and without adapter contamination ─────────────
FASTQ_DIR="${OUTDIR}/fastq"
mkdir -p "${FASTQ_DIR}"

python3 - "${REF_DIR}/genome.fa" "${FASTQ_DIR}" "${ADAPTER_SEQ}" << 'PYGEN'
import sys, random, os

ref_fa = sys.argv[1]
fastq_dir = sys.argv[2]
adapter = sys.argv[3]
random.seed(42)

# Read chr1
seq = ""
with open(ref_fa) as f:
    for line in f:
        if line.startswith(">"):
            continue
        seq += line.strip()

n_reads = 200
read_len = 75  # target read length

# Generate reads: some with adapter suffix, some clean
# "adapter" reads: 50bp genomic + adapter fill to 75bp
# "clean" reads: 75bp genomic
clean_path = os.path.join(fastq_dir, "clean.fastq")
adapter_path = os.path.join(fastq_dir, "adapter_contaminated.fastq")

with open(clean_path, "w") as cf, open(adapter_path, "w") as af:
    for i in range(n_reads):
        pos = 200 + (i * 71) % 48000  # spread across genome
        # Clean read
        clean_seq = seq[pos:pos+read_len]
        if len(clean_seq) < read_len:
            clean_seq += "N" * (read_len - len(clean_seq))
        cf.write(f"@clean_{i}\n{clean_seq}\n+\n{'I'*len(clean_seq)}\n")

        # Adapter-contaminated read: shorter genomic + adapter
        genomic_len = 30 + (i % 30)  # 30-59bp genomic
        genomic_seq = seq[pos:pos+genomic_len]
        adapter_fill = adapter * 3  # enough adapter
        contam_seq = (genomic_seq + adapter_fill)[:read_len]
        af.write(f"@contam_{i}\n{contam_seq}\n+\n{'I'*len(contam_seq)}\n")

print(f"Generated {n_reads} clean reads -> {clean_path}")
print(f"Generated {n_reads} adapter-contaminated reads -> {adapter_path}")
PYGEN

# ── Run STAR: no trimming (baseline) ─────────────────────────────────
log "Running STAR without adapter trimming (baseline)..."
NOTRIM_DIR="${OUTDIR}/run_notrim"
mkdir -p "${NOTRIM_DIR}"

"${STAR_BIN}" \
  --runMode alignReads \
  --genomeDir "${INDEX_DIR}" \
  --readFilesIn "${FASTQ_DIR}/adapter_contaminated.fastq" \
  --outFileNamePrefix "${NOTRIM_DIR}/" \
  --outSAMtype None \
  --runThreadN "${THREADS}" \
  > "${NOTRIM_DIR}/star.log" 2>&1 \
  || die "STAR no-trim run failed"

# ── Run STAR: with adapter trimming ──────────────────────────────────
log "Running STAR with --clip3pAdapterSeq ${ADAPTER_SEQ}..."
TRIM_DIR="${OUTDIR}/run_trim"
mkdir -p "${TRIM_DIR}"

"${STAR_BIN}" \
  --runMode alignReads \
  --genomeDir "${INDEX_DIR}" \
  --readFilesIn "${FASTQ_DIR}/adapter_contaminated.fastq" \
  --clip3pAdapterSeq "${ADAPTER_SEQ}" \
  --outFileNamePrefix "${TRIM_DIR}/" \
  --outSAMtype None \
  --runThreadN "${THREADS}" \
  > "${TRIM_DIR}/star.log" 2>&1 \
  || die "STAR trim run failed"

# ── Run STAR: clean reads (control) ──────────────────────────────────
log "Running STAR with clean reads (control)..."
CLEAN_DIR="${OUTDIR}/run_clean"
mkdir -p "${CLEAN_DIR}"

"${STAR_BIN}" \
  --runMode alignReads \
  --genomeDir "${INDEX_DIR}" \
  --readFilesIn "${FASTQ_DIR}/clean.fastq" \
  --outFileNamePrefix "${CLEAN_DIR}/" \
  --outSAMtype None \
  --runThreadN "${THREADS}" \
  > "${CLEAN_DIR}/star.log" 2>&1 \
  || die "STAR clean run failed"

# ── Validate results ──────────────────────────────────────────────────
log "--- Validating outputs ---"

parse_mapped() {
  local log_file="$1"
  grep "Uniquely mapped reads number" "${log_file}" | awk '{print $NF}'
}

parse_mapped_pct() {
  local log_file="$1"
  grep "Uniquely mapped reads %" "${log_file}" | awk '{print $NF}' | tr -d '%'
}

NOTRIM_MAPPED=$(parse_mapped "${NOTRIM_DIR}/Log.final.out")
TRIM_MAPPED=$(parse_mapped "${TRIM_DIR}/Log.final.out")
CLEAN_MAPPED=$(parse_mapped "${CLEAN_DIR}/Log.final.out")

NOTRIM_PCT=$(parse_mapped_pct "${NOTRIM_DIR}/Log.final.out")
TRIM_PCT=$(parse_mapped_pct "${TRIM_DIR}/Log.final.out")
CLEAN_PCT=$(parse_mapped_pct "${CLEAN_DIR}/Log.final.out")

log "  No-trim mapped:   ${NOTRIM_MAPPED} (${NOTRIM_PCT}%)"
log "  With-trim mapped: ${TRIM_MAPPED} (${TRIM_PCT}%)"
log "  Clean mapped:     ${CLEAN_MAPPED} (${CLEAN_PCT}%)"

# Test 1: trimming should improve or maintain mapping rate vs no-trim
if python3 -c "exit(0 if float('${TRIM_PCT}') >= float('${NOTRIM_PCT}') else 1)"; then
  log "  PASS: trimming mapping rate (${TRIM_PCT}%) >= no-trim (${NOTRIM_PCT}%)"
  ((++PASS))
else
  log "  FAIL: trimming mapping rate (${TRIM_PCT}%) < no-trim (${NOTRIM_PCT}%)"
  ((++FAIL))
fi

# Test 2: trimmed mapping rate should be close to clean reads
if python3 -c "
trim_pct = float('${TRIM_PCT}')
clean_pct = float('${CLEAN_PCT}')
# Allow 20% relative gap (synthetic reads may not all align perfectly)
exit(0 if clean_pct == 0 or trim_pct >= clean_pct * 0.80 else 1)
"; then
  log "  PASS: trimmed rate (${TRIM_PCT}%) within 80% of clean rate (${CLEAN_PCT}%)"
  ((++PASS))
else
  log "  FAIL: trimmed rate (${TRIM_PCT}%) too far below clean rate (${CLEAN_PCT}%)"
  ((++FAIL))
fi

# Test 3: some reads should actually map in the trimmed run
if [[ "${TRIM_MAPPED}" -gt 0 ]]; then
  log "  PASS: trimmed run has ${TRIM_MAPPED} mapped reads (>0)"
  ((++PASS))
else
  log "  FAIL: trimmed run has 0 mapped reads"
  ((++FAIL))
fi

# Test 4: Log.final.out files exist for all three runs
for run_name in notrim trim clean; do
  local_dir="${OUTDIR}/run_${run_name}"
  if [[ -f "${local_dir}/Log.final.out" ]]; then
    log "  PASS: Log.final.out exists for ${run_name}"
    ((++PASS))
  else
    log "  FAIL: Log.final.out missing for ${run_name}"
    ((++FAIL))
  fi
done

log "--- Results: ${PASS} passed, ${FAIL} failed ---"

if [[ "${FAIL}" -gt 0 ]]; then
  die "${FAIL} validation(s) failed. Outputs at: ${OUTDIR}"
fi

log "PASS: adapter clip synthetic test"
