#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
RUN_SCRIPT="${REPO_ROOT}/scripts/run_ucsf_perturb_yremove_batch.sh"

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

mkdir -p "${tmpdir}/dataset" "${tmpdir}/gex/sample1" "${tmpdir}/guides/sample1" "${tmpdir}/genome"

make_fastq() {
  local path="$1"
  python3 - <<'PY' "${path}"
import gzip, sys
with gzip.open(sys.argv[1], "wt") as handle:
    handle.write("@r1\nACGT\n+\nFFFF\n")
PY
}

make_fastq "${tmpdir}/gex/sample1/sample1_S1_L001_R1_001.fastq.gz"
make_fastq "${tmpdir}/gex/sample1/sample1_S1_L001_R2_001.fastq.gz"
make_fastq "${tmpdir}/guides/sample1/sample1_S1_L001_R1_001.fastq.gz"
make_fastq "${tmpdir}/guides/sample1/sample1_S1_L001_R2_001.fastq.gz"

printf 'AAAA\n' > "${tmpdir}/solo_wl.txt"
printf 'AAAA\n' > "${tmpdir}/cr_wl.txt"

cat > "${tmpdir}/feature_ref.csv" <<EOF
id,name,sequence,feature_type
g1,G1,ACGT,CRISPR Guide Capture
EOF

config="${tmpdir}/config.csv"
cat > "${config}" <<EOF
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

STAR_BIN=/bin/true \
  "${RUN_SCRIPT}" \
  --dataset-root "${tmpdir}/dataset" \
  --genome-dir "${tmpdir}/genome" \
  --solo-cb-whitelist "${tmpdir}/solo_wl.txt" \
  --cr-whitelist "${tmpdir}/cr_wl.txt" \
  --cr-config "${config}" \
  --out-root "${out}" \
  --dry-run

manifest="${out}/samples/sample1/RUN_MANIFEST.txt"
command="${out}/samples/sample1/RUN_COMMAND.sh"
env_file="${out}/samples/sample1/cr_config_inputs.env"
pf_multi="${out}/samples/sample1/pf_multi_config.csv"

[[ -f "${manifest}" ]]
[[ -f "${command}" ]]
[[ -f "${env_file}" ]]
[[ -f "${pf_multi}" ]]

grep -F "cr_config=${config}" "${manifest}" >/dev/null
grep -F "cr_gene_expression_chemistry=auto" "${manifest}" >/dev/null
grep -F "${tmpdir}/gex/sample1" "${pf_multi}" >/dev/null
grep -F "${tmpdir}/guides/sample1" "${pf_multi}" >/dev/null
grep -F -- "--pfMultiConfig ${pf_multi}" "${command}" >/dev/null
grep -F -- "${tmpdir}/gex/sample1/sample1_S1_L001_R2_001.fastq.gz" "${command}" >/dev/null
grep -F -- "${tmpdir}/guides/sample1/sample1_S1_L001_R2_001.fastq.gz" "${command}" >/dev/null
grep -F $'sample\tstatus\tsample_root' "${out}/RUNS.tsv" >/dev/null
grep -F $'sample1\tdone\t'"${out}/samples/sample1" "${out}/RUNS.tsv" >/dev/null

echo "PASS"
