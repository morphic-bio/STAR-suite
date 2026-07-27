#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
STAR_SOURCE_DIR="${ROOT_DIR}/core/legacy/source"
TMP_DIR="$(mktemp -d /tmp/star-flex-gdna-test.XXXXXX)"
trap 'rm -rf -- "${TMP_DIR}"' EXIT

make -C "${STAR_SOURCE_DIR}" test_FlexGdna

STAR_BIN="${STAR_BIN:-${STAR_SOURCE_DIR}/STAR}"
if [[ -x "${STAR_BIN}" ]]; then
    help_text="$("${STAR_BIN}" --help 2>&1)"
    grep -q '^soloFlexGdna[[:space:]]' <<<"${help_text}"
    grep -q '^soloFlexGdnaProbeSet[[:space:]]' <<<"${help_text}"

    # The required diagnostic flag must remain inert without --flex yes.
    # A deliberately missing genome should fail at genome loading, not at
    # Flex probe-metadata validation.
    set +e
    "${STAR_BIN}" \
        --runMode alignReads \
        --genomeDir "${TMP_DIR}/missing-index" \
        --readFilesIn /dev/null \
        --soloFlexGdna yes \
        --outFileNamePrefix "${TMP_DIR}/nonflex." \
        >"${TMP_DIR}/nonflex.stdout" 2>"${TMP_DIR}/nonflex.stderr"
    nonflex_status=$?
    set -e
    [[ "${nonflex_status}" -ne 0 ]]
    ! grep -q 'requires complete probe-region metadata' \
        "${TMP_DIR}/nonflex.stderr" "${TMP_DIR}/nonflex.Log.out"

    # A required diagnostic cannot silently run without the final FlexFilter
    # cell set that defines its molecule denominator.
    set +e
    "${STAR_BIN}" \
        --genomeDir /nonexistent \
        --readFilesIn \
            "${ROOT_DIR}/tests/fixtures/trim_qc_fastq_tiny.fastq" \
            "${ROOT_DIR}/tests/fixtures/trim_qc_fastq_tiny.fastq" \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist None \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 12 \
        --soloFeatures Gene \
        --soloProbeList \
            "${ROOT_DIR}/tests/fixtures/flex_probe_gene_list_tiny.txt" \
        --flex yes \
        --soloRunFlexFilter no \
        --soloFlexGdna yes \
        --outFileNamePrefix "${TMP_DIR}/no-filter." \
        >"${TMP_DIR}/no-filter.stdout" 2>"${TMP_DIR}/no-filter.stderr"
    no_filter_status=$?
    set -e
    [[ "${no_filter_status}" -ne 0 ]]
    grep -q 'requires final FlexFilter cell calls' \
        "${TMP_DIR}/no-filter.stderr"
fi

# The Python cache pilot is an independent FH01SEQ1 reader used by the H2
# utilities. Confirm that v3 region bits do not contaminate its 15-bit gene.
ROOT_DIR="${ROOT_DIR}" CACHE_PATH="${TMP_DIR}/cache_v3.bin" python3 - <<'PY'
import os
import struct
import sys

sys.path.insert(0, os.environ["ROOT_DIR"])
from scripts.flex_h01_pilot import load_binary_hash_cache

path = os.environ["CACHE_PATH"]
with open(path, "wb") as handle:
    handle.write(struct.pack("<8sHHIQ", b"FH01SEQ1", 3, 50, 24, 1))
    handle.write(struct.pack("<QQIBBH", 11, 22, 77 | (2 << 30), 1, 0, 9))

rows = load_binary_hash_cache(path)
assert len(rows) == 1
assert rows[0]["resolved_gene_idx15"] == 77
assert rows[0]["probe_region"] == 2
assert rows[0]["sample_idx"] == 9
PY

echo "Flex gDNA focused suite passed"
