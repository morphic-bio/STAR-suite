#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR="${OUTDIR:-$(mktemp -d /tmp/star_suite_transcriptvb_gc_effective_length_XXXXXX)}"
CXX_BIN="${CXX:-g++}"

mkdir -p "${OUTDIR}"

cat > "${OUTDIR}/transcripts.fa" <<'FASTA'
>low_gc
AAAAAAAAAAAAAAAAAAAA
>high_gc
GCGCGCGCGCGCGCGCGCGC
FASTA

cat > "${OUTDIR}/gc_effective_length_test.cpp" <<'CPP'
#include "alignment_model.h"
#include "effective_length_wrapper.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

static bool near(double a, double b, double eps = 1e-9) {
    return std::fabs(a - b) <= eps;
}

static int expect_near(const std::string& label, double got, double want) {
    if (!near(got, want)) {
        std::cerr << "FAIL " << label << ": got " << got << ", want " << want << "\n";
        return 1;
    }
    std::cout << "PASS " << label << ": " << got << "\n";
    return 0;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: gc_effective_length_test transcripts.fa\n";
        return 2;
    }

    libem::Transcriptome transcriptome;
    if (!transcriptome.loadFromFasta(argv[1])) {
        return 2;
    }
    if (!transcriptome.reorderByNames({"low_gc", "high_gc"})) {
        std::cerr << "failed to reorder transcriptome\n";
        return 2;
    }

    std::vector<int32_t> raw_lengths = {20, 20};
    std::vector<double> fld(8, 0.0);
    fld[6] = 1.0;

    std::vector<double> neutral = computeEffectiveLengthsFromPMFWrapper(fld, raw_lengths);

    std::vector<double> unit_gc(101, 1.0);
    std::vector<double> neutral_gc =
        computeGCBiasedEffectiveLengthsWrapper(transcriptome, fld, raw_lengths, unit_gc);

    int failed = 0;
    failed += expect_near("neutral_low", neutral[0], 14.0);
    failed += expect_near("neutral_high", neutral[1], 14.0);
    failed += expect_near("unit_gc_low_matches_neutral", neutral_gc[0], neutral[0]);
    failed += expect_near("unit_gc_high_matches_neutral", neutral_gc[1], neutral[1]);

    std::vector<double> biased_gc(101, 1.0);
    biased_gc[0] = 0.5;
    biased_gc[100] = 2.0;
    std::vector<double> biased =
        computeGCBiasedEffectiveLengthsWrapper(transcriptome, fld, raw_lengths, biased_gc);

    failed += expect_near("low_gc_biased", biased[0], 7.0);
    failed += expect_near("high_gc_biased", biased[1], 28.0);

    return failed == 0 ? 0 : 1;
}
CPP

"${CXX_BIN}" -std=c++11 \
    -I"${ROOT_DIR}/core/features/vbem/source/libem" \
    "${OUTDIR}/gc_effective_length_test.cpp" \
    "${ROOT_DIR}/core/features/vbem/source/libem/effective_length_wrapper.cpp" \
    "${ROOT_DIR}/core/features/vbem/source/libem/effective_length.cpp" \
    "${ROOT_DIR}/core/features/vbem/source/libem/alignment_model.cpp" \
    -o "${OUTDIR}/gc_effective_length_test"

"${OUTDIR}/gc_effective_length_test" "${OUTDIR}/transcripts.fa" \
    | tee "${OUTDIR}/test.log"

echo "Output: ${OUTDIR}"
