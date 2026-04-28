// Thin CLI wrapper around libscrna::atac::RunAtacEvidenceFromPeaks.
// Same on-disk output as the original standalone tool; the implementation
// now lives in libscrna so STAR's chromap orchestration can call it
// directly without spawning a child process.

#include <cstring>
#include <iostream>
#include <string>

#include "libscrna/AtacEvidenceFromPeaks.h"

namespace {

void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " --fragments fragments.tsv[.gz]"
        << " --peaks peaks.narrowPeak"
        << " --out atac_evidence.tsv"
        << " [--barcode-whitelist whitelist.txt]"
        << " [--barcode-suffix -1]\n"
        << "Reads a Chromap/ARC 5-col fragments TSV (chrom start end barcode count)\n"
        << "and a narrowPeak file (BED-like; first 3 columns used) and emits per-barcode\n"
        << "ATAC evidence with the following columns:\n"
        << "  barcode, atac_peak_region_cutsites, atac_peak_region_fragments,\n"
        << "  atac_fragments, atac_peak_fraction\n";
}

bool parse_args(int argc, char** argv,
                libscrna::atac::AtacEvidenceFromPeaksOptions* opts) {
    if (opts == nullptr) {
        return false;
    }
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--fragments" && i + 1 < argc) {
            opts->fragments_path = argv[++i];
        } else if (arg == "--peaks" && i + 1 < argc) {
            opts->peaks_path = argv[++i];
        } else if (arg == "--out" && i + 1 < argc) {
            opts->out_path = argv[++i];
        } else if (arg == "--barcode-whitelist" && i + 1 < argc) {
            opts->whitelist_path = argv[++i];
        } else if (arg == "--barcode-suffix" && i + 1 < argc) {
            opts->barcode_suffix = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            return false;
        } else {
            std::cerr << "Unknown or incomplete argument: " << arg << "\n";
            usage(argv[0]);
            return false;
        }
    }
    if (opts->fragments_path.empty() || opts->peaks_path.empty() ||
        opts->out_path.empty()) {
        usage(argv[0]);
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    libscrna::atac::AtacEvidenceFromPeaksOptions opts;
    if (!parse_args(argc, argv, &opts)) {
        return 1;
    }
    return libscrna::atac::RunAtacEvidenceFromPeaks(opts, &std::cerr);
}
