#ifndef MULTIOME_ATAC_PEAK_MEX_H_
#define MULTIOME_ATAC_PEAK_MEX_H_

#include <cstdint>
#include <string>

namespace star {
namespace multiome {

enum class MultiomeAtacPeakMexThresholdMode {
  P_VALUE,
  Q_VALUE
};

struct MultiomeAtacPeakMexArgs {
  std::string fragments;
  std::string sidecar;
  std::string barcode_translate;
  bool barcode_translate_from_first = false;

  std::string peaks;
  std::string summits_out;
  std::string out_dir;
  std::string metrics_tsv;
  std::string temp_dir;
  std::string keep_intermediates_dir;

  int threads = 1;
  bool call_peaks_from_sidecar = false;
  bool force = false;
  uint64_t max_barcodes = 0;

  MultiomeAtacPeakMexThresholdMode macs3_threshold_mode =
      MultiomeAtacPeakMexThresholdMode::P_VALUE;
  double macs3_pvalue = 1e-5;
  double macs3_qvalue = 0.0;
  int macs3_min_length = 200;
  int macs3_max_gap = 30;
  bool macs3_uint8_counts = true;
};

// Returns 0 on success and non-zero on validation or materialization errors.
// This entry point is safe to call after the alignment phase has fully closed
// the AEV1 sidecar and associated chromosome metadata file.
int RunMultiomeAtacPeakMex(const MultiomeAtacPeakMexArgs &args);

}  // namespace multiome
}  // namespace star

#endif  // MULTIOME_ATAC_PEAK_MEX_H_
