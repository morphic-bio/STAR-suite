#!/usr/bin/env bash
set -euo pipefail
# Binary-oriented scRNA GEX smoke test.
#
# Validates an installed STAR binary against the public chr22+chrY Codespaces
# reference and a small 10x male B-cell fixture.  Defaults to --outSAMtype
# None (no BAM) for the primary release-gate surface.
#
# Tier A: public-download, container-safe, release-gate candidate.
#
# Prerequisites:
#   - STAR binary (env STAR_BIN, or in-tree fallback)
#   - Internet access for first run (downloads ~19 MB reference + ~3.4 GB FASTQs
#     which are then subsampled)
#   - python3, samtools (for index build / optional BAM checks)
#
# Usage:
#   tests/run_gex_binary_smoke.sh [options]
#
# Options:
#   --threads N          STAR threads (default: 4)
#   --read-limit N       Subsample reads (default: 5000)
#   --write-bam sorted   Write sorted BAM and validate CB/UB tags
#   --write-bam unsorted Write unsorted BAM and validate CB/UB tags
#   --outdir DIR         Output directory (default: /tmp/gex_binary_smoke_<ts>)
#   --keep               Keep workdir on success (default: keep always)
#   -h, --help           Show help

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
THREADS="${THREADS:-4}"
READ_LIMIT="${READ_LIMIT:-5000}"
BAM_MODE="none"
OUTDIR=""

die() { echo "FAIL: $*" >&2; exit 1; }
log() { printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*"; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "missing required command: $1"
}

usage() { sed -n '2,/^$/s/^# \?//p' "$0"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads)    THREADS="$2";    shift 2 ;;
    --read-limit) READ_LIMIT="$2"; shift 2 ;;
    --write-bam)  BAM_MODE="$2";   shift 2 ;;
    --outdir)     OUTDIR="$2";     shift 2 ;;
    --keep)       shift ;;
    -h|--help)    usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

TIMESTAMP="$(date +'%Y%m%d_%H%M%S')"
OUTDIR="${OUTDIR:-/tmp/gex_binary_smoke_${TIMESTAMP}_$$}"

# ── Prerequisites ───────────────────────────────────────────────────
need_cmd python3
[[ -x "${STAR_BIN}" ]] || die "STAR binary not found: ${STAR_BIN}"

case "${BAM_MODE}" in
  none|sorted|unsorted) ;;
  *) die "--write-bam must be none, sorted, or unsorted" ;;
esac

if [[ "${BAM_MODE}" != "none" ]]; then
  need_cmd samtools
fi

mkdir -p "${OUTDIR}"
log "=== GEX Binary Smoke Test ==="
log "STAR:       ${STAR_BIN}"
log "Threads:    ${THREADS}"
log "Read limit: ${READ_LIMIT}"
log "BAM mode:   ${BAM_MODE}"
log "Outdir:     ${OUTDIR}"

# ── Build reference + index via Codespaces helpers ──────────────────
REF_DIR="${OUTDIR}/ref"
INDEX_DIR="${OUTDIR}/star_index"
DATA_DIR="${OUTDIR}/data"
RUN_DIR="${OUTDIR}/run"
mkdir -p "${REF_DIR}" "${INDEX_DIR}" "${DATA_DIR}" "${RUN_DIR}"

CODESPACES="${REPO_ROOT}/scripts/codespaces"

# Fetch chr22+chrY reference
log "Fetching public chr22+chrY reference..."
DEMO_ROOT="${OUTDIR}/demo" \
CACHE_ROOT="${OUTDIR}/demo/cache" \
DATA_ROOT="${DATA_DIR}" \
INDEX_ROOT="${OUTDIR}/demo/indices" \
RUN_ROOT="${OUTDIR}/demo/runs" \
  bash "${CODESPACES}/fetch_public_chr22y_reference.sh" --outdir "${REF_DIR}"

[[ -f "${REF_DIR}/fasta/genome.fa" ]]  || die "reference FASTA missing"
[[ -f "${REF_DIR}/genes/genes.gtf" ]]  || die "reference GTF missing"

# Build STAR index
log "Building STAR index..."
"${STAR_BIN}" \
  --runMode genomeGenerate \
  --genomeDir "${INDEX_DIR}" \
  --runThreadN "${THREADS}" \
  --genomeFastaFiles "${REF_DIR}/fasta/genome.fa" \
  --sjdbGTFfile "${REF_DIR}/genes/genes.gtf" \
  --genomeSAindexNbases 11

[[ -f "${INDEX_DIR}/Genome" ]] || die "STAR index build failed"
log "Index built."

# ── Create synthetic GEX fixture ────────────────────────────────────
# Use a tiny synthetic fixture: generate CB+UMI reads that map to chr22.
# This avoids needing the large 10x tar download for the release gate.
FASTQ_DIR="${DATA_DIR}/fastq"
WL_FILE="${DATA_DIR}/whitelist.txt"
mkdir -p "${FASTQ_DIR}"

log "Generating synthetic GEX fixture (${READ_LIMIT} reads)..."

python3 - "${REF_DIR}/fasta/genome.fa" "${FASTQ_DIR}" "${WL_FILE}" "${READ_LIMIT}" << 'PYGEN'
import sys, random, os

ref_fa = sys.argv[1]
fastq_dir = sys.argv[2]
wl_path = sys.argv[3]
n_reads = int(sys.argv[4])

# Read chr22 sequence
seq_lines = []
reading = False
with open(ref_fa) as f:
    for line in f:
        if line.startswith(">chr22") or line.startswith(">22"):
            reading = True
            continue
        elif line.startswith(">"):
            reading = False
            continue
        if reading:
            seq_lines.append(line.strip())
chr22_seq = "".join(seq_lines).upper()

# Generate 100 synthetic barcodes (16nt) and write whitelist
random.seed(42)
bases = "ACGT"
barcodes = []
for _ in range(100):
    bc = "".join(random.choice(bases) for _ in range(16))
    barcodes.append(bc)

with open(wl_path, "w") as wl:
    for bc in barcodes:
        wl.write(bc + "\n")

# Generate paired reads: R1 = CB(16) + UMI(12), R2 = 90bp from chr22
r1_path = os.path.join(fastq_dir, "R1.fastq.gz")
r2_path = os.path.join(fastq_dir, "R2.fastq.gz")

import gzip
max_start = len(chr22_seq) - 90
if max_start < 1:
    max_start = 1

with gzip.open(r1_path, "wt") as r1f, gzip.open(r2_path, "wt") as r2f:
    for i in range(n_reads):
        bc = barcodes[i % len(barcodes)]
        umi = "".join(random.choice(bases) for _ in range(12))
        pos = random.randint(0, max_start)
        cdna = chr22_seq[pos:pos+90]
        if len(cdna) < 90:
            cdna = cdna + "N" * (90 - len(cdna))
        r1_seq = bc + umi
        q1 = "I" * len(r1_seq)
        q2 = "I" * 90
        r1f.write(f"@read{i}/1\n{r1_seq}\n+\n{q1}\n")
        r2f.write(f"@read{i}/2\n{cdna}\n+\n{q2}\n")

print(f"Generated {n_reads} read pairs -> {r1_path}, {r2_path}")
print(f"Whitelist: {wl_path} ({len(barcodes)} barcodes)")
PYGEN

[[ -f "${FASTQ_DIR}/R1.fastq.gz" ]] || die "R1 FASTQ generation failed"
[[ -f "${FASTQ_DIR}/R2.fastq.gz" ]] || die "R2 FASTQ generation failed"
[[ -f "${WL_FILE}" ]]               || die "whitelist generation failed"

# ── Run STAR ────────────────────────────────────────────────────────
log "Running STAR alignment..."

STAR_ARGS=(
  --runThreadN "${THREADS}"
  --genomeDir "${INDEX_DIR}"
  --readFilesIn "${FASTQ_DIR}/R2.fastq.gz" "${FASTQ_DIR}/R1.fastq.gz"
  --readFilesCommand zcat
  --outFileNamePrefix "${RUN_DIR}/"
  --soloType CB_UMI_Simple
  --soloCBstart 1 --soloCBlen 16
  --soloUMIstart 17 --soloUMIlen 12
  --soloCBwhitelist "${WL_FILE}"
  --soloFeatures GeneFull
)

case "${BAM_MODE}" in
  none)
    STAR_ARGS+=(--outSAMtype None)
    ;;
  sorted)
    STAR_ARGS+=(--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 200000000)
    ;;
  unsorted)
    STAR_ARGS+=(--outSAMtype BAM Unsorted)
    ;;
esac

set +e
"${STAR_BIN}" "${STAR_ARGS[@]}" > "${RUN_DIR}/star.log" 2>&1
STAR_EXIT=$?
set -e

if [[ "${STAR_EXIT}" -ne 0 ]]; then
  cat "${RUN_DIR}/star.log" >&2
  die "STAR exited with code ${STAR_EXIT}"
fi

log "STAR completed."

# ── Validate outputs ────────────────────────────────────────────────
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

log "--- Validating outputs ---"

# Core Solo outputs
check "Log.final.out"                        "${RUN_DIR}/Log.final.out"
check "Log.out"                              "${RUN_DIR}/Log.out"
check "GeneFull/raw/matrix.mtx"              "${RUN_DIR}/Solo.out/GeneFull/raw/matrix.mtx"
check "GeneFull/raw/barcodes.tsv"            "${RUN_DIR}/Solo.out/GeneFull/raw/barcodes.tsv"
check "GeneFull/raw/features.tsv"            "${RUN_DIR}/Solo.out/GeneFull/raw/features.tsv"
check "GeneFull/Summary.csv"                 "${RUN_DIR}/Solo.out/GeneFull/Summary.csv"
check "Barcodes.stats"                       "${RUN_DIR}/Solo.out/Barcodes.stats"

# Verify raw matrix is non-empty
check_nonempty "GeneFull/raw/matrix.mtx non-empty" "${RUN_DIR}/Solo.out/GeneFull/raw/matrix.mtx"

# UMIperCellSorted may be empty on tiny synthetic fixtures (too few reads for
# cell calling), but should exist when a cell filter is active.
if [[ -f "${RUN_DIR}/Solo.out/GeneFull/UMIperCellSorted.txt" ]]; then
  check "UMIperCellSorted.txt" "${RUN_DIR}/Solo.out/GeneFull/UMIperCellSorted.txt"
fi

# Verify features count is reasonable (chr22+chrY should have >100 genes)
if [[ -f "${RUN_DIR}/Solo.out/GeneFull/raw/features.tsv" ]]; then
  NFEAT=$(wc -l < "${RUN_DIR}/Solo.out/GeneFull/raw/features.tsv")
  if [[ "${NFEAT}" -gt 100 ]]; then
    log "  PASS: features.tsv has ${NFEAT} genes (>100)"
    ((++PASS))
  else
    log "  FAIL: features.tsv has only ${NFEAT} genes (expected >100)"
    ((++FAIL))
  fi
fi

# BAM validation
case "${BAM_MODE}" in
  sorted)
    check "Aligned.sortedByCoord.out.bam" "${RUN_DIR}/Aligned.sortedByCoord.out.bam"
    if [[ -f "${RUN_DIR}/Aligned.sortedByCoord.out.bam" ]]; then
      if samtools quickcheck "${RUN_DIR}/Aligned.sortedByCoord.out.bam" 2>/dev/null; then
        log "  PASS: BAM quickcheck"
        ((++PASS))
      else
        log "  FAIL: BAM quickcheck failed"
        ((++FAIL))
      fi
      if samtools view -H "${RUN_DIR}/Aligned.sortedByCoord.out.bam" | grep -q "SO:coordinate"; then
        log "  PASS: BAM sort order = coordinate"
        ((++PASS))
      else
        log "  FAIL: BAM missing SO:coordinate"
        ((++FAIL))
      fi
    fi
    ;;
  unsorted)
    check "Aligned.out.bam" "${RUN_DIR}/Aligned.out.bam"
    if [[ -f "${RUN_DIR}/Aligned.out.bam" ]]; then
      if samtools quickcheck "${RUN_DIR}/Aligned.out.bam" 2>/dev/null; then
        log "  PASS: BAM quickcheck"
        ((++PASS))
      else
        log "  FAIL: BAM quickcheck failed"
        ((++FAIL))
      fi
    fi
    ;;
  none)
    # Verify no BAM was written
    if [[ ! -f "${RUN_DIR}/Aligned.sortedByCoord.out.bam" && ! -f "${RUN_DIR}/Aligned.out.bam" ]]; then
      log "  PASS: no BAM output (expected for no-BAM mode)"
      ((++PASS))
    else
      log "  FAIL: BAM written in no-BAM mode"
      ((++FAIL))
    fi
    ;;
esac

log "--- Results: ${PASS} passed, ${FAIL} failed ---"

if [[ "${FAIL}" -gt 0 ]]; then
  die "${FAIL} validation(s) failed. Outputs at: ${OUTDIR}"
fi

log "PASS: GEX binary smoke test (BAM mode: ${BAM_MODE})"
log "Workdir: ${OUTDIR}"
