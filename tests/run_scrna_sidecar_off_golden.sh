#!/usr/bin/env bash
set -euo pipefail

# Self-contained normal scRNA regression for the default-off spatial sidecar.
# The golden manifest was captured from the feature branch base commit using
# this exact fixture and compares stable count outputs byte-for-byte by SHA-256.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
GOLDEN="${SCRNA_SIDECAR_OFF_GOLDEN:-${SCRIPT_DIR}/fixtures/scrna_sidecar_off_golden.json}"
OUTDIR=""
EMIT_MANIFEST=""

die() { echo "FAIL: $*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="${2:?}"; shift 2 ;;
    --emit-manifest) EMIT_MANIFEST="${2:?}"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 [--outdir DIR] [--emit-manifest FILE]"
      exit 0
      ;;
    *) die "unknown option: $1" ;;
  esac
done

[[ -x "${STAR_BIN}" ]] || die "STAR binary not found: ${STAR_BIN}"
command -v python3 >/dev/null || die "python3 is required"

if [[ -z "${OUTDIR}" ]]; then
  OUTDIR="$(mktemp -d /tmp/star-scrna-sidecar-off-golden-XXXXXX)"
else
  [[ ! -e "${OUTDIR}" ]] || die "output path already exists: ${OUTDIR}"
  mkdir -p "${OUTDIR}"
fi

REF="${OUTDIR}/reference"
INDEX="${OUTDIR}/index"
FASTQ="${OUTDIR}/fastq"
RUN="${OUTDIR}/run"
mkdir -p "${REF}" "${INDEX}" "${FASTQ}" "${RUN}"

python3 - "${REF}" "${FASTQ}" <<'PY'
import gzip
import random
import sys
from pathlib import Path

reference = Path(sys.argv[1])
fastq = Path(sys.argv[2])
randomizer = random.Random(20260722)
bases = "ACGT"
genome = "".join(randomizer.choice(bases) for _ in range(12000))
genes = [
    ("GENE_A", 1001, 1300),
    ("GENE_B", 3501, 3800),
    ("GENE_C", 6001, 6300),
    ("GENE_D", 8501, 8800),
]
(reference / "genome.fa").write_text(">chrSynthetic\n" + genome + "\n")
with (reference / "genes.gtf").open("w") as handle:
    for gene, start, end in genes:
        attributes = f'gene_id "{gene}"; gene_name "{gene}";'
        handle.write(f"chrSynthetic\tfixture\tgene\t{start}\t{end}\t.\t+\t.\t{attributes}\n")
        handle.write(
            f"chrSynthetic\tfixture\ttranscript\t{start}\t{end}\t.\t+\t.\t"
            f'{attributes} transcript_id "{gene}.T1";\n'
        )
        handle.write(
            f"chrSynthetic\tfixture\texon\t{start}\t{end}\t.\t+\t.\t"
            f'{attributes} transcript_id "{gene}.T1"; exon_number "1";\n'
        )

barcode = "ACGTACGTACGTACGT"
(fastq / "whitelist.txt").write_text(barcode + "\n")
records = []

def add(gene, umi, count):
    gene_index = [item[0] for item in genes].index(gene)
    start = genes[gene_index][1] - 1
    cdna = genome[start + 31 : start + 106]
    for _ in range(count):
        records.append((gene, umi, cdna))

# Multi-gene corrected-UMI group with unequal read support.
add("GENE_A", "CCCCCCCCCCCC", 5)
add("GENE_B", "CCCCCCCCCCCC", 3)
# Multi-gene corrected-UMI group with equal read support.
add("GENE_C", "GGGGGGGGGGGG", 2)
add("GENE_D", "GGGGGGGGGGGG", 2)
# Asymmetric one-mismatch UMI neighborhoods exercise the original-root support
# part of MultiGeneUMI_CR after per-gene correction.
add("GENE_A", "AAAAAAAAAAAA", 6)
add("GENE_A", "CAAAAAAAAAAA", 3)
add("GENE_A", "GAAAAAAAAAAA", 3)
add("GENE_B", "AAAAAAAAAAAA", 7)
add("GENE_B", "TAAAAAAAAAAA", 4)
# Independent low-depth molecule.
add("GENE_D", "TTTTTTTTTTTT", 1)

with gzip.open(fastq / "R1.fastq.gz", "wt") as r1, gzip.open(
    fastq / "R2.fastq.gz", "wt"
) as r2:
    for index, (gene, umi, cdna) in enumerate(records):
        name = f"normal_scrna_{index}_{gene}"
        read1 = barcode + umi
        r1.write(f"@{name}/1\n{read1}\n+\n{'I' * len(read1)}\n")
        r2.write(f"@{name}/2\n{cdna}\n+\n{'I' * len(cdna)}\n")
print(f"generated {len(records)} normal scRNA read pairs")
PY

"${STAR_BIN}" \
  --runMode genomeGenerate \
  --runThreadN 1 \
  --genomeDir "${INDEX}" \
  --genomeFastaFiles "${REF}/genome.fa" \
  --sjdbGTFfile "${REF}/genes.gtf" \
  --sjdbOverhang 74 \
  --genomeSAindexNbases 4 \
  >"${OUTDIR}/genome.stdout.log" 2>"${OUTDIR}/genome.stderr.log"

# Deliberately omit --soloSpatialFeatureSidecar. This is a normal Solo run.
"${STAR_BIN}" \
  --runThreadN 1 \
  --genomeDir "${INDEX}" \
  --readFilesIn "${FASTQ}/R2.fastq.gz" "${FASTQ}/R1.fastq.gz" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${RUN}/" \
  --outSAMtype None \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloCBwhitelist "${FASTQ}/whitelist.txt" \
  --soloCBmatchWLtype Exact \
  --soloFeatures GeneFull \
  --soloStrand Unstranded \
  --soloUMIdedup 1MM_CR \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloMultiMappers Unique \
  --soloCellFilter None \
  >"${OUTDIR}/star.stdout.log" 2>"${OUTDIR}/star.stderr.log"

FILES=(
  "Solo.out/GeneFull/raw/matrix.mtx"
  "Solo.out/GeneFull/raw/barcodes.tsv"
  "Solo.out/GeneFull/raw/features.tsv"
)
for relative in "${FILES[@]}"; do
  [[ -f "${RUN}/${relative}" ]] || die "missing count output: ${relative}"
done

if find "${RUN}" -maxdepth 1 -name '*.bin' -print -quit | grep -q .; then
  die "normal scRNA run unexpectedly emitted a sidecar binary"
fi

if [[ -n "${EMIT_MANIFEST}" ]]; then
  python3 - "${RUN}" "${EMIT_MANIFEST}" "${FILES[@]}" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
output = Path(sys.argv[2])
files = {}
for relative in sys.argv[3:]:
    path = root / relative
    files[relative] = {"bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
manifest = {
    "schema": "star_suite.scrna_sidecar_off_golden.v1",
    "base_commit": "a996107e271c013e39f9398151deda0017da35d6",
    "files": files,
}
output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
print(output)
PY
else
  [[ -f "${GOLDEN}" ]] || die "golden manifest not found: ${GOLDEN}"
  python3 - "${RUN}" "${GOLDEN}" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
manifest = json.loads(Path(sys.argv[2]).read_text())
if manifest.get("schema") != "star_suite.scrna_sidecar_off_golden.v1":
    raise SystemExit("invalid golden manifest schema")
for relative, expected in manifest["files"].items():
    path = root / relative
    actual = {"bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    if actual != expected:
        raise SystemExit(f"normal scRNA count output diverged from branch base: {relative}: {actual} != {expected}")
print(f"normal scRNA sidecar-off count outputs match {manifest['base_commit']}")
PY
fi

echo "PASS: normal scRNA sidecar-off golden parity"
