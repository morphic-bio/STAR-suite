#ifndef STAR_CHROMAP_CONTRACT_H_
#define STAR_CHROMAP_CONTRACT_H_

#include <cstdint>
#include <string>
#include <vector>

namespace star {
namespace multiome {

enum class ChromapOutputFormat {
  BED,
  TAGALIGN,
  PAF,
  SAM,
  BAM,
  CRAM,
  PAIRS
};

enum class Tn5ShiftMode {
  CLASSICAL,
  SYMMETRIC
};

enum class ChromapMacs3FragPeaksSource {
  FILE,
  MEMORY
};

enum class ChromapMacs3FragThresholdMode {
  P_VALUE,
  Q_VALUE
};

enum class ChromapInputFormat {
  FASTQ,
  CBQ
};

enum class ChromapContractStatus {
  OK,
  INVALID_CONFIG,
  CHROMAP_FAILED
};

// Acquire returns waitNs (time blocked while obtaining the permit). Release
// reports waitNs from the matching acquire plus per-batch telemetry. Signatures
// mirror process_features' pf_permit_{acquire,release}_fn so the same STAR-side
// shim shape can serve both domains.
typedef uint64_t (*ChromapPermitAcquireFn)(void *hook_ctx);
typedef void (*ChromapPermitReleaseFn)(void *hook_ctx,
                                       uint64_t wait_ns,
                                       uint64_t work_units,
                                       uint64_t work_bytes,
                                       uint64_t work_ns);

struct ChromapPermitHooks {
  ChromapPermitAcquireFn acquire = nullptr;
  ChromapPermitReleaseFn release = nullptr;
  void *hook_ctx = nullptr;
};

struct ChromapAtacConfig {
  std::string reference_fasta;
  std::string chromap_index;
  ChromapInputFormat input_format = ChromapInputFormat::FASTQ;
  std::vector<std::string> read1_fastqs;
  std::vector<std::string> read2_fastqs;
  std::vector<std::string> barcode_fastqs;
  std::vector<std::string> read_pair_cbqs;
  std::vector<std::string> barcode_cbqs;
  std::string barcode_whitelist;
  std::string barcode_translate_table;
  // If true, the barcode translate table is read in natural <from_bc>\t<to_bc>
  // order (col1 = source / hash key). Default false preserves the historical
  // Chromap convention where col1 is the destination and col2 is the key.
  bool barcode_translate_from_first_column = false;
  std::string output_path;
  // Optional: with BAM/CRAM output_path, emit scATAC fragments to this path
  // (one Chromap pass; same semantics as Chromap --atac-fragments).
  std::string fragment_output_path;
  std::string summary_path;
  std::string temp_dir;
  std::string read_format;
  std::string chromosome_order;
  std::string pairs_flipping_chromosome_order;

  int threads = 1;
  int max_insert_size = 2000;
  int mapq_threshold = 30;
  int min_read_length = 30;
  int barcode_correction_error_threshold = 1;
  double barcode_correction_probability_threshold = 0.9;
  int error_threshold = 8;
  int min_num_seeds_required_for_mapping = 2;
  int drop_repetitive_reads = 500000;
  int hts_threads = 0;

  bool trim_adapters = true;
  bool remove_pcr_duplicates = true;
  bool remove_pcr_duplicates_at_cell_level = true;
  bool tn5_shift = true;
  Tn5ShiftMode tn5_shift_mode = Tn5ShiftMode::CLASSICAL;
  bool output_mappings_not_in_whitelist = false;
  bool allocate_multi_mappings = false;
  bool low_memory_mode = false;
  uint64_t low_memory_ram_limit = 0;
  bool skip_barcode_check = false;
  bool sort_bam = false;
  bool write_index = false;
  uint64_t sort_bam_ram_limit = 8ULL * 1024 * 1024 * 1024;
  bool emit_no_y_bam = false;
  bool emit_y_bam = false;
  std::string no_y_output_path;
  std::string y_output_path;

  bool call_macs3_frag_peaks = false;
  std::string macs3_frag_peaks_output;
  std::string macs3_frag_summits_output;
  std::string macs3_frag_keep_intermediates_dir;
  ChromapMacs3FragThresholdMode macs3_frag_threshold_mode =
      ChromapMacs3FragThresholdMode::P_VALUE;
  double macs3_frag_pvalue = 1e-5;
  double macs3_frag_qvalue = 0.0;
  int macs3_frag_min_length = 200;
  int macs3_frag_max_gap = 30;
  bool macs3_frag_uint8_counts = true;
  ChromapMacs3FragPeaksSource macs3_frag_peaks_source =
      ChromapMacs3FragPeaksSource::MEMORY;
  bool macs3_frag_low_mem = false;

  // Optional per-barcode ATAC evidence-from-peaks output. When set with
  // MACS3 FRAG peaks from the memory source, libchromap writes
  // fragment_output_path as the fixed-record binary sidecar and the
  // contract calls libscrna's sidecar reader after peak calling. Empty
  // disables.
  std::string atac_evidence_from_peaks_output;

  ChromapOutputFormat output_format = ChromapOutputFormat::BED;
  ChromapPermitHooks permit_hooks;
};

struct ChromapAtacResult {
  ChromapContractStatus status = ChromapContractStatus::INVALID_CONFIG;
  int exit_code = 1;
  std::string message;
  std::string output_path;
  std::string summary_path;
};

ChromapAtacResult runChromapAtac(const ChromapAtacConfig &config);

const char *chromapContractStatusName(ChromapContractStatus status);

}  // namespace multiome
}  // namespace star

#endif  // STAR_CHROMAP_CONTRACT_H_
