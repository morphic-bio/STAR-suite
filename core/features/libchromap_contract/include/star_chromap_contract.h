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

enum class ChromapContractStatus {
  OK,
  INVALID_CONFIG,
  CHROMAP_FAILED
};

struct ChromapPermitWork {
  uint64_t work_units = 0;
  uint64_t work_bytes = 0;
  uint64_t elapsed_ns = 0;
};

typedef bool (*ChromapPermitAcquireFn)(void *context,
                                       uint32_t requested_permits,
                                       ChromapPermitWork *work);
typedef void (*ChromapPermitReleaseFn)(void *context,
                                       uint32_t released_permits,
                                       const ChromapPermitWork *work);

struct ChromapPermitHooks {
  ChromapPermitAcquireFn acquire = nullptr;
  ChromapPermitReleaseFn release = nullptr;
  void *context = nullptr;
};

struct ChromapAtacConfig {
  std::string reference_fasta;
  std::string chromap_index;
  std::vector<std::string> read1_fastqs;
  std::vector<std::string> read2_fastqs;
  std::vector<std::string> barcode_fastqs;
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

  bool call_macs3_frag_peaks = false;
  std::string macs3_frag_peaks_output;
  std::string macs3_frag_summits_output;
  std::string macs3_frag_keep_intermediates_dir;
  double macs3_frag_pvalue = 1e-5;
  int macs3_frag_min_length = 200;
  int macs3_frag_max_gap = 30;
  bool macs3_frag_uint8_counts = true;
  ChromapMacs3FragPeaksSource macs3_frag_peaks_source =
      ChromapMacs3FragPeaksSource::FILE;
  bool macs3_frag_low_mem = false;

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
