#ifndef LIBSCRNA_ATAC_EVIDENCE_FROM_PEAKS_H
#define LIBSCRNA_ATAC_EVIDENCE_FROM_PEAKS_H

#include <cstdint>
#include <ostream>
#include <string>

namespace libscrna {
namespace atac {

typedef bool (*AtacEvidenceBarcodeDecoderFn)(uint64_t key,
                                             uint32_t barcode_length,
                                             void* ctx,
                                             std::string* out);

int RunAtacEvidenceFromBinary(
    const std::string& sidecar_path,
    const std::string& peaks_path,
    const std::string& out_path,
    AtacEvidenceBarcodeDecoderFn decoder,
    void* decoder_ctx,
    std::ostream* err = nullptr);

}  // namespace atac
}  // namespace libscrna

#endif  // LIBSCRNA_ATAC_EVIDENCE_FROM_PEAKS_H
