#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
RUN_SCRIPT="${REPO_ROOT}/scripts/run_ucsf_perturb_yremove_batch.sh"
STAR_BIN="${REPO_ROOT}/core/legacy/source/STAR"
TINY_REF="${REPO_ROOT}/tests/solo_smoke/ref/star_index"
WL="${REPO_ROOT}/tests/solo_smoke/whitelist.txt"

tmpdir="$(mktemp -d /tmp/ucsf_batch_cr_multifeature_tiny_XXXXXX)"
trap 'rm -rf "${tmpdir}"' EXIT

mkdir -p "${tmpdir}/dataset" "${tmpdir}/gex/sample1" "${tmpdir}/guides/sample1" "${tmpdir}/custom/sample1"

python3 - <<'PY' "${tmpdir}"
import gzip
import sys
from pathlib import Path

root = Path(sys.argv[1])
cb = 'ACGTACGTACGTACGT'
gex_seq = 'AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT'
guide_seq = 'ACGTACGTACGTACGTACGT'
custom_seq = 'TGCATGCATGCATGCATGCA'

for lib, seq in [('gex', gex_seq), ('guides', guide_seq), ('custom', custom_seq)]:
    r1 = root / lib / 'sample1' / 'sample1_S1_L001_R1_001.fastq.gz'
    r2 = root / lib / 'sample1' / 'sample1_S1_L001_R2_001.fastq.gz'
    with gzip.open(r1, 'wt') as h1, gzip.open(r2, 'wt') as h2:
        for i in range(5):
            umi = f'{i:012d}'[-12:]
            name = f'@{lib}{i}'
            r1_seq = cb + umi
            r1_qual = 'I' * len(r1_seq)
            r2_qual = 'I' * len(seq)
            h1.write(f'{name}/1\n{r1_seq}\n+\n{r1_qual}\n')
            h2.write(f'{name}/2\n{seq}\n+\n{r2_qual}\n')
PY

cat > "${tmpdir}/guide_ref.csv" <<EOF
id,name,sequence,feature_type
g1,G1,ACGTACGTACGTACGTACGT,CRISPR Guide Capture
EOF

cat > "${tmpdir}/custom_ref.csv" <<EOF
id,name,sequence,feature_type
c1,C1,TGCATGCATGCATGCATGCA,Custom
EOF

cat > "${tmpdir}/config.csv" <<EOF
[gene-expression]
reference,/storage/A375-CR-9.01/refdata-gex-GRCh38-2024-A
create-bam,false
chemistry,auto

[star]
sample-id,DE_30KO

[feature]
reference,${tmpdir}/guide_ref.csv

[libraries]
fastq_id,fastqs,feature_types,star_chemistry,star_feature_ref,star_library_id
mRNA_DE_30KO,${tmpdir}/gex/sample1,Gene Expression,,,gex_main
PolyIII_DE_30KO,${tmpdir}/guides/sample1,CRISPR Guide Capture,NXT,,guide_main
LARRY_DE_30KO,${tmpdir}/custom/sample1,Custom,TRU,${tmpdir}/custom_ref.csv,custom_main
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

sample_root="${out}/samples/DE_30KO"
[[ -f "${sample_root}/RUN_COMPLETE.ok" ]]
[[ -f "${sample_root}/run/Aligned.out_Y.bam" ]]
[[ -f "${sample_root}/run/Aligned.out_noY.bam" ]]
[[ -d "${sample_root}/run/y_separated" ]]
[[ "$(find "${sample_root}/run/y_separated" -maxdepth 1 -type f -name '*.fastq.gz' | wc -l)" == "4" ]]
grep -F "feature_dirs=${tmpdir}/guides/sample1,${tmpdir}/custom/sample1" "${sample_root}/RUN_MANIFEST.txt" >/dev/null
grep -F "cr_sample_id=DE_30KO" "${sample_root}/RUN_MANIFEST.txt" >/dev/null
grep -F "${tmpdir}/guides/sample1,DE_30KO,CRISPR Guide Capture,CRISPR Guide Capture,NXT,,guide_main" "${sample_root}/pf_multi_config.csv" >/dev/null
grep -F "${tmpdir}/custom/sample1,DE_30KO,Custom,Custom,TRU,${tmpdir}/custom_ref.csv,custom_main" "${sample_root}/pf_multi_config.csv" >/dev/null

echo "PASS"
