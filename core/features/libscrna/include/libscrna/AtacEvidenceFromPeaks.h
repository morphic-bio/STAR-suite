#ifndef LIBSCRNA_ATAC_EVIDENCE_FROM_PEAKS_H
#define LIBSCRNA_ATAC_EVIDENCE_FROM_PEAKS_H

#include <ostream>
#include <string>

namespace libscrna {
namespace atac {

struct AtacEvidenceFromPeaksOptions {
    std::string fragments_path;       // Chromap/ARC 5-col TSV (.tsv or .tsv.gz)
    std::string peaks_path;           // narrowPeak (BED-like; first 3 cols used)
    std::string out_path;             // per-barcode evidence TSV out
    std::string whitelist_path;       // optional CB whitelist (one per line)
    std::string barcode_suffix;       // optional suffix appended to whitelist barcodes
};

// Read fragments (column 5 = per-row count, ARC-style) and peaks; emit
// per-barcode evidence with columns:
//   barcode, atac_peak_region_cutsites, atac_peak_region_fragments,
//   atac_fragments, atac_peak_fraction
//
// Returns 0 on success, non-zero on failure. Diagnostic / progress messages
// are written to `*err` (default std::cerr).
int RunAtacEvidenceFromPeaks(const AtacEvidenceFromPeaksOptions& opts,
                             std::ostream* err = nullptr);

}  // namespace atac
}  // namespace libscrna

#endif  // LIBSCRNA_ATAC_EVIDENCE_FROM_PEAKS_H
