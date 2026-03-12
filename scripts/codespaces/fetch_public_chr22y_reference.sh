#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"
# shellcheck source=public_demo_sources.sh
source "${SCRIPT_DIR}/public_demo_sources.sh"

OUTDIR="${OUTDIR:-${DATA_ROOT}/public_human_chr22y_ref}"

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Fetch a tiny public human reference root containing only chr22 and chrY.
The output layout matches the STAR/CR-style reference surface used by the
Codespaces demo scripts.

Options:
  --outdir DIR   Output directory. Default: ${DATA_ROOT}/public_human_chr22y_ref
  -h, --help     Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

need_cmd curl
need_cmd gzip

OUTDIR="$(mkdir -p "${OUTDIR}" && cd "${OUTDIR}" && pwd)"
RAW_DIR="${OUTDIR}/raw"
FASTA_DIR="${OUTDIR}/fasta"
GENES_DIR="${OUTDIR}/genes"
mkdir -p "${RAW_DIR}" "${FASTA_DIR}" "${GENES_DIR}"

CHR22_GZ="${RAW_DIR}/chr22.fa.gz"
CHRY_GZ="${RAW_DIR}/chrY.fa.gz"
GTF_GZ="${RAW_DIR}/genes.gtf.gz"
GENOME_FA="${FASTA_DIR}/genome.fa"
GENES_GTF="${GENES_DIR}/genes.gtf"
MANIFEST="${OUTDIR}/MANIFEST.txt"

[[ -f "${CHR22_GZ}" ]] || curl -L --fail -o "${CHR22_GZ}" "${CODESPACES_PUBLIC_HUMAN_CHR22_FASTA_URL}"
[[ -f "${CHRY_GZ}" ]] || curl -L --fail -o "${CHRY_GZ}" "${CODESPACES_PUBLIC_HUMAN_CHRY_FASTA_URL}"
[[ -f "${GTF_GZ}" ]] || curl -L --fail -o "${GTF_GZ}" "${CODESPACES_PUBLIC_HUMAN_GTF_URL}"

if [[ ! -f "${GENOME_FA}" ]]; then
  {
    gzip -cd "${CHR22_GZ}"
    gzip -cd "${CHRY_GZ}"
  } > "${GENOME_FA}"
fi

if [[ ! -f "${GENES_GTF}" ]]; then
  python3 - <<'PY' "${GTF_GZ}" "${GENES_GTF}"
import gzip
import sys
from pathlib import Path
src = Path(sys.argv[1])
out = Path(sys.argv[2])
keep = {'chr22', 'chrY', '22', 'Y'}
with gzip.open(src, 'rt') as inp, out.open('w') as handle:
    for raw in inp:
        if raw.startswith('#'):
            handle.write(raw)
            continue
        fields = raw.split('\t', 1)
        if fields and fields[0] in keep:
            handle.write(raw)
PY
fi

python3 - <<'PY' "${GENOME_FA}" "${GENES_GTF}" "${MANIFEST}"
import sys
from pathlib import Path
fa = Path(sys.argv[1])
gtf = Path(sys.argv[2])
man = Path(sys.argv[3])
chroms = []
with fa.open() as handle:
    for line in handle:
        if line.startswith('>'):
            chroms.append(line[1:].strip().split()[0])
records = 0
with gtf.open() as handle:
    for line in handle:
        if line.startswith('#'):
            continue
        records += 1
man.write_text(
    'Public human chr22+chrY mini-reference\n'
    f'Ensembl release: 110\n'
    f'Chromosomes: {",".join(chroms)}\n'
    f'GTF records: {records}\n',
    encoding='utf-8',
)
PY

log "Prepared ${OUTDIR}"
