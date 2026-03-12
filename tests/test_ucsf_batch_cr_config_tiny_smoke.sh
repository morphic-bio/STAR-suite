#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
RUN_SCRIPT="${REPO_ROOT}/scripts/run_ucsf_perturb_yremove_batch.sh"
STAR_BIN="${REPO_ROOT}/core/legacy/source/STAR"
TINY_REF="${REPO_ROOT}/tests/solo_smoke/ref/star_index"
WL="${REPO_ROOT}/tests/solo_smoke/whitelist.txt"

tmpdir="$(mktemp -d /tmp/ucsf_batch_cr_tiny_XXXXXX)"
trap 'rm -rf "${tmpdir}"' EXIT

mkdir -p "${tmpdir}/dataset" "${tmpdir}/gex/sample1" "${tmpdir}/guides/sample1"

python3 - <<'PY' "${tmpdir}"
import gzip
import sys
from pathlib import Path

root = Path(sys.argv[1])

gex_r1 = root / "gex/sample1/sample1_S1_L001_R1_001.fastq.gz"
gex_r2 = root / "gex/sample1/sample1_S1_L001_R2_001.fastq.gz"
guide_r1 = root / "guides/sample1/sample1_S1_L001_R1_001.fastq.gz"
guide_r2 = root / "guides/sample1/sample1_S1_L001_R2_001.fastq.gz"

cb = "ACGTACGTACGTACGT"
gex_seq = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT"
guide_seq = "ACGTACGTACGTACGTACGT"

gex_reads = []
guide_reads = []
for i in range(5):
    umi = f"{i:012d}"[-12:]
    gex_reads.append((f"@gex{i}", cb + umi, gex_seq))
    guide_reads.append((f"@guide{i}", cb + umi, guide_seq))

for path, reads in ((gex_r1, gex_reads), (guide_r1, guide_reads)):
    with gzip.open(path, "wt") as handle:
        for name, r1, _ in reads:
            handle.write(f"{name}/1\n{r1}\n+\n{'I'*len(r1)}\n")

for path, reads in ((gex_r2, gex_reads), (guide_r2, guide_reads)):
    with gzip.open(path, "wt") as handle:
        for name, _, r2 in reads:
            handle.write(f"{name}/2\n{r2}\n+\n{'I'*len(r2)}\n")
PY

cat > "${tmpdir}/feature_ref.csv" <<EOF
id,name,sequence,feature_type
g1,G1,ACGTACGTACGTACGTACGT,CRISPR Guide Capture
EOF

cat > "${tmpdir}/config.csv" <<EOF
[gene-expression]
reference,/storage/A375-CR-9.01/refdata-gex-GRCh38-2024-A
create-bam,false
chemistry,auto

[feature]
reference,${tmpdir}/feature_ref.csv

[libraries]
fastq_id,fastqs,feature_types
sample1,${tmpdir}/gex/sample1,Gene Expression
sample1,${tmpdir}/guides/sample1,CRISPR Guide Capture
EOF

out="${tmpdir}/out"

"${RUN_SCRIPT}" \
  --star-bin "${STAR_BIN}" \
  --dataset-root "${tmpdir}/dataset" \
  --genome-dir "${TINY_REF}" \
  --solo-cb-whitelist "${WL}" \
  --cr-whitelist "${WL}" \
  --threads 2 \
  --cr-config "${tmpdir}/config.csv" \
  --out-root "${out}"

sample_root="${out}/samples/sample1"

[[ -f "${sample_root}/RUN_COMPLETE.ok" ]]
[[ -f "${sample_root}/run/Aligned.out_Y.bam" ]]
[[ -f "${sample_root}/run/Aligned.out_noY.bam" ]]
[[ -d "${sample_root}/run/y_separated" ]]
[[ "$(find "${sample_root}/run/y_separated" -maxdepth 1 -type f -name '*.fastq.gz' | wc -l)" == "4" ]]

grep -F "cr_config=${tmpdir}/config.csv" "${sample_root}/RUN_MANIFEST.txt" >/dev/null
grep -F $'sample1\tdone\t'"${sample_root}" "${out}/RUNS.tsv" >/dev/null

echo "PASS"
