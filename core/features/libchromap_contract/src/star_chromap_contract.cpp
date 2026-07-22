#include "star_chromap_contract.h"

#include <cerrno>
#include <cctype>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

#include "barcode_translator.h"
#include "libchromap.h"
#include "libscrna/AtacEvidenceFromPeaks.h"

namespace star {
namespace multiome {
namespace {

ChromapAtacResult makeResult(ChromapContractStatus status,
                             int exit_code,
                             const std::string &message,
                             const ChromapAtacConfig &config) {
  ChromapAtacResult result;
  result.status = status;
  result.exit_code = exit_code;
  result.message = message;
  result.output_path = config.output_path;
  result.summary_path = config.summary_path;
  return result;
}

std::string trimCopy(const std::string &input) {
  size_t first = input.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) {
    return "";
  }
  size_t last = input.find_last_not_of(" \t\r\n");
  return input.substr(first, last - first + 1);
}

std::string lowerCopy(std::string s) {
  for (size_t i = 0; i < s.size(); ++i) {
    s[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(s[i])));
  }
  return s;
}

bool isUnsetToken(const std::string &input) {
  const std::string value = lowerCopy(trimCopy(input));
  return value.empty() || value == "-" || value == "none";
}

bool hasUnsetPath(const std::vector<std::string> &paths) {
  for (size_t i = 0; i < paths.size(); ++i) {
    if (isUnsetToken(paths[i])) {
      return true;
    }
  }
  return false;
}

bool isDirectory(const std::string &path) {
  struct stat st;
  return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

bool ensureDirectory(const std::string &path, std::string *error) {
  const std::string trimmed = trimCopy(path);
  if (trimmed.empty() || trimmed == "/") {
    return true;
  }
  if (isDirectory(trimmed)) {
    return true;
  }

  std::string current;
  size_t pos = 0;
  if (!trimmed.empty() && trimmed[0] == '/') {
    current = "/";
    pos = 1;
  }

  while (pos <= trimmed.size()) {
    size_t next = trimmed.find('/', pos);
    std::string part = trimmed.substr(
        pos, next == std::string::npos ? std::string::npos : next - pos);
    pos = (next == std::string::npos) ? trimmed.size() + 1 : next + 1;
    if (part.empty()) {
      continue;
    }
    if (!current.empty() && current[current.size() - 1] != '/') {
      current += '/';
    }
    current += part;
    if (isDirectory(current)) {
      continue;
    }
    if (mkdir(current.c_str(), 0775) != 0 && errno != EEXIST) {
      if (error != nullptr) {
        *error = "could not create directory " + current + ": " +
                 std::strerror(errno);
      }
      return false;
    }
    if (!isDirectory(current)) {
      if (error != nullptr) {
        *error = "path exists but is not a directory: " + current;
      }
      return false;
    }
  }

  return true;
}

std::vector<std::string> trimPaths(const std::vector<std::string> &paths) {
  std::vector<std::string> trimmed;
  trimmed.reserve(paths.size());
  for (size_t i = 0; i < paths.size(); ++i) {
    trimmed.push_back(trimCopy(paths[i]));
  }
  return trimmed;
}

std::string validateConfig(const ChromapAtacConfig &config) {
  if (isUnsetToken(config.reference_fasta)) {
    return "reference FASTA is required";
  }
  if (isUnsetToken(config.chromap_index)) {
    return "Chromap index is required";
  }
  if (config.input_format == ChromapInputFormat::CBQ) {
    if (!config.read1_fastqs.empty() || !config.read2_fastqs.empty() ||
        !config.barcode_fastqs.empty()) {
      return "FASTQ inputs cannot be mixed with CBQ input";
    }
    if (config.read_pair_cbqs.empty() || hasUnsetPath(config.read_pair_cbqs)) {
      return "ATAC read-pair CBQs are required";
    }
    if (config.barcode_cbqs.empty() || hasUnsetPath(config.barcode_cbqs)) {
      return "ATAC barcode CBQs are required";
    }
    if (config.read_pair_cbqs.size() != config.barcode_cbqs.size()) {
      return "read-pair CBQ and barcode CBQ counts differ";
    }
  } else {
    if (!config.read_pair_cbqs.empty() || !config.barcode_cbqs.empty()) {
      return "CBQ inputs require CBQ input format";
    }
    if (config.read1_fastqs.empty() || hasUnsetPath(config.read1_fastqs)) {
      return "ATAC read1 FASTQs are required";
    }
    if (config.read2_fastqs.empty() || hasUnsetPath(config.read2_fastqs)) {
      return "ATAC read2 FASTQs are required";
    }
    if (config.barcode_fastqs.empty() || hasUnsetPath(config.barcode_fastqs)) {
      return "ATAC barcode FASTQs are required";
    }
    if (config.read1_fastqs.size() != config.read2_fastqs.size()) {
      return "read1 and read2 FASTQ counts differ";
    }
    if (config.read1_fastqs.size() != config.barcode_fastqs.size()) {
      return "read1 and barcode FASTQ counts differ";
    }
  }
  if (isUnsetToken(config.barcode_whitelist)) {
    return "ATAC barcode whitelist is required";
  }
  if (isUnsetToken(config.output_path)) {
    return "Chromap output path is required";
  }
  if (config.threads < 1) {
    return "thread count must be >= 1";
  }
  if (config.mapq_threshold < 0 || config.mapq_threshold > 255) {
    return "MAPQ threshold must be in [0, 255]";
  }
  if (config.hts_threads < 0) {
    return "HTS thread count must be >= 0";
  }
  if (config.write_index && !config.sort_bam) {
    return "BAM/CRAM index output requires coordinate sorting";
  }
  if (config.sort_bam &&
      config.output_format != ChromapOutputFormat::BAM &&
      config.output_format != ChromapOutputFormat::CRAM) {
    return "coordinate sorting requires BAM or CRAM output";
  }
  if (config.write_index &&
      config.output_format != ChromapOutputFormat::BAM &&
      config.output_format != ChromapOutputFormat::CRAM) {
    return "index output requires BAM or CRAM output";
  }
  if ((config.emit_no_y_bam || config.emit_y_bam) &&
      config.output_format != ChromapOutputFormat::SAM &&
      config.output_format != ChromapOutputFormat::BAM &&
      config.output_format != ChromapOutputFormat::CRAM) {
    return "Y/noY stream output requires SAM, BAM, or CRAM output";
  }
  if (config.emit_no_y_bam && isUnsetToken(config.no_y_output_path)) {
    return "noY output path is required when noY stream output is enabled";
  }
  if (config.emit_y_bam && isUnsetToken(config.y_output_path)) {
    return "Y output path is required when Y stream output is enabled";
  }
  if (config.emit_no_y_bam &&
      trimCopy(config.no_y_output_path) == trimCopy(config.output_path)) {
    return "noY output path must differ from primary output path";
  }
  if (config.emit_y_bam &&
      trimCopy(config.y_output_path) == trimCopy(config.output_path)) {
    return "Y output path must differ from primary output path";
  }
  if (config.emit_no_y_bam && config.emit_y_bam &&
      trimCopy(config.no_y_output_path) == trimCopy(config.y_output_path)) {
    return "Y and noY output paths must differ";
  }
  if (!isUnsetToken(config.fragment_output_path) &&
      config.output_format != ChromapOutputFormat::BAM &&
      config.output_format != ChromapOutputFormat::CRAM) {
    return "fragment output path requires BAM or CRAM Chromap output";
  }
  if (!isUnsetToken(config.fragment_output_path) &&
      !isUnsetToken(config.output_path) &&
      trimCopy(config.fragment_output_path) == trimCopy(config.output_path)) {
    return "fragment output path must differ from primary output path";
  }
  if (config.call_macs3_frag_peaks) {
    if (isUnsetToken(config.macs3_frag_peaks_output)) {
      return "MACS3 FRAG narrowPeak output is required when peak calling is enabled";
    }
    if (isUnsetToken(config.macs3_frag_summits_output)) {
      return "MACS3 FRAG summits output is required when peak calling is enabled";
    }
    if (config.macs3_frag_threshold_mode ==
        ChromapMacs3FragThresholdMode::Q_VALUE) {
      if (!(config.macs3_frag_qvalue > 0.0 &&
            config.macs3_frag_qvalue <= 1.0)) {
        return "MACS3 FRAG q-value must be in (0, 1]";
      }
    } else if (!(config.macs3_frag_pvalue > 0.0 &&
                 config.macs3_frag_pvalue <= 1.0)) {
      return "MACS3 FRAG p-value must be in (0, 1]";
    }
  }
  if (!isUnsetToken(config.atac_evidence_from_peaks_output) &&
      isUnsetToken(config.fragment_output_path)) {
    return "in-process ATAC evidence requires --chromapAtacSecondaryFragments to point at the binary sidecar output";
  }
  if (!isUnsetToken(config.atac_evidence_from_peaks_output) &&
      (!config.call_macs3_frag_peaks ||
       config.macs3_frag_peaks_source != ChromapMacs3FragPeaksSource::MEMORY)) {
    return "in-process ATAC evidence requires MACS3 FRAG peak calling with memory source";
  }
  if ((config.permit_hooks.acquire == nullptr) !=
      (config.permit_hooks.release == nullptr)) {
    return "permit acquire and release hooks must be supplied together";
  }
  return "";
}

chromap::MappingOutputFormat toChromapOutputFormat(
    ChromapOutputFormat format) {
  switch (format) {
    case ChromapOutputFormat::BED:
      return chromap::MAPPINGFORMAT_BED;
    case ChromapOutputFormat::TAGALIGN:
      return chromap::MAPPINGFORMAT_TAGALIGN;
    case ChromapOutputFormat::PAF:
      return chromap::MAPPINGFORMAT_PAF;
    case ChromapOutputFormat::SAM:
      return chromap::MAPPINGFORMAT_SAM;
    case ChromapOutputFormat::BAM:
      return chromap::MAPPINGFORMAT_BAM;
    case ChromapOutputFormat::CRAM:
      return chromap::MAPPINGFORMAT_CRAM;
    case ChromapOutputFormat::PAIRS:
      return chromap::MAPPINGFORMAT_PAIRS;
  }
  return chromap::MAPPINGFORMAT_UNKNOWN;
}

chromap::MappingParameters toChromapParameters(
    const ChromapAtacConfig &config) {
  chromap::MappingParameters parameters;

  parameters.reference_file_path = trimCopy(config.reference_fasta);
  parameters.index_file_path = trimCopy(config.chromap_index);
  if (config.input_format == ChromapInputFormat::CBQ) {
    parameters.read_input_format = chromap::ReadInputFormat::kCbq;
    parameters.read_pair_cbq_paths = trimPaths(config.read_pair_cbqs);
    parameters.barcode_cbq_paths = trimPaths(config.barcode_cbqs);
  } else {
    parameters.read_input_format = chromap::ReadInputFormat::kFastq;
    parameters.read_file1_paths = trimPaths(config.read1_fastqs);
    parameters.read_file2_paths = trimPaths(config.read2_fastqs);
    parameters.barcode_file_paths = trimPaths(config.barcode_fastqs);
  }
  parameters.barcode_whitelist_file_path = trimCopy(config.barcode_whitelist);
  if (!isUnsetToken(config.barcode_translate_table)) {
    parameters.barcode_translate_table_file_path =
        trimCopy(config.barcode_translate_table);
  }
  parameters.barcode_translate_from_first_column =
      config.barcode_translate_from_first_column;
  parameters.mapping_output_file_path = trimCopy(config.output_path);
  if (!isUnsetToken(config.fragment_output_path)) {
    parameters.atac_fragment_output_file_path =
        trimCopy(config.fragment_output_path);
    parameters.atac_fragment_binary_output_file_path =
        trimCopy(config.fragment_output_path);
  }
  if (!isUnsetToken(config.summary_path)) {
    parameters.summary_metadata_file_path = trimCopy(config.summary_path);
  }
  if (!isUnsetToken(config.temp_dir)) {
    parameters.temp_directory_path = trimCopy(config.temp_dir);
  }
  if (!isUnsetToken(config.read_format)) {
    parameters.read_format = trimCopy(config.read_format);
  }
  if (!isUnsetToken(config.chromosome_order)) {
    parameters.custom_rid_order_file_path = trimCopy(config.chromosome_order);
  }
  if (!isUnsetToken(config.pairs_flipping_chromosome_order)) {
    parameters.pairs_flipping_custom_rid_order_file_path =
        trimCopy(config.pairs_flipping_chromosome_order);
  }

  parameters.num_threads = config.threads;
  parameters.max_insert_size = config.max_insert_size;
  parameters.mapq_threshold = static_cast<uint8_t>(config.mapq_threshold);
  parameters.min_read_length = config.min_read_length;
  parameters.barcode_correction_error_threshold =
      config.barcode_correction_error_threshold;
  parameters.barcode_correction_probability_threshold =
      config.barcode_correction_probability_threshold;
  parameters.error_threshold = config.error_threshold;
  parameters.min_num_seeds_required_for_mapping =
      config.min_num_seeds_required_for_mapping;
  parameters.drop_repetitive_reads = config.drop_repetitive_reads;
  parameters.hts_threads = config.hts_threads;

  parameters.trim_adapters = config.trim_adapters;
  parameters.remove_pcr_duplicates = config.remove_pcr_duplicates;
  parameters.remove_pcr_duplicates_at_bulk_level =
      !config.remove_pcr_duplicates_at_cell_level;
  parameters.is_bulk_data = false;
  parameters.Tn5_shift = config.tn5_shift;
  if (config.tn5_shift_mode == Tn5ShiftMode::CLASSICAL) {
    parameters.Tn5_forward_shift = 4;
    parameters.Tn5_reverse_shift = -5;
  } else {
    parameters.Tn5_forward_shift = 4;
    parameters.Tn5_reverse_shift = -4;
  }
  parameters.output_mappings_not_in_whitelist =
      config.output_mappings_not_in_whitelist;
  parameters.allocate_multi_mappings = config.allocate_multi_mappings;
  parameters.only_output_unique_mappings = !config.allocate_multi_mappings;
  parameters.low_memory_mode = config.low_memory_mode;
  parameters.low_mem_ram_limit = config.low_memory_ram_limit;
  parameters.skip_barcode_check = config.skip_barcode_check;
  parameters.mapping_output_format = toChromapOutputFormat(config.output_format);
  parameters.sort_bam = config.sort_bam;
  parameters.write_index = config.write_index;
  parameters.sort_bam_ram_limit = config.sort_bam_ram_limit;
  parameters.emit_noY_stream = config.emit_no_y_bam;
  parameters.emit_Y_stream = config.emit_y_bam;
  if (config.emit_no_y_bam) {
    parameters.noY_output_path = trimCopy(config.no_y_output_path);
  }
  if (config.emit_y_bam) {
    parameters.Y_output_path = trimCopy(config.y_output_path);
  }
  parameters.call_macs3_frag_peaks = config.call_macs3_frag_peaks;
  if (!isUnsetToken(config.macs3_frag_peaks_output)) {
    parameters.macs3_frag_peaks_narrowpeak_path =
        trimCopy(config.macs3_frag_peaks_output);
  }
  if (!isUnsetToken(config.macs3_frag_summits_output)) {
    parameters.macs3_frag_peaks_summits_path =
        trimCopy(config.macs3_frag_summits_output);
  }
  if (!isUnsetToken(config.macs3_frag_keep_intermediates_dir)) {
    parameters.macs3_frag_keep_intermediates_dir =
        trimCopy(config.macs3_frag_keep_intermediates_dir);
  }
  parameters.macs3_frag_threshold_mode =
      (config.macs3_frag_threshold_mode ==
       ChromapMacs3FragThresholdMode::Q_VALUE)
          ? chromap::peaks::Macs3FragThresholdMode::kQValue
          : chromap::peaks::Macs3FragThresholdMode::kPValue;
  parameters.macs3_frag_pvalue = config.macs3_frag_pvalue;
  parameters.macs3_frag_qvalue = config.macs3_frag_qvalue;
  parameters.macs3_frag_min_length = config.macs3_frag_min_length;
  parameters.macs3_frag_max_gap = config.macs3_frag_max_gap;
  parameters.macs3_frag_uint8_counts = config.macs3_frag_uint8_counts;
  parameters.macs3_frag_peaks_source =
      (config.macs3_frag_peaks_source == ChromapMacs3FragPeaksSource::MEMORY)
          ? chromap::Macs3FragPeaksSource::kMemory
          : chromap::Macs3FragPeaksSource::kFile;
  parameters.macs3_frag_low_mem = config.macs3_frag_low_mem;

  return parameters;
}

struct EvidenceBarcodeDecoderCtx {
  chromap::BarcodeTranslator *translator = nullptr;
};

bool decodeEvidenceBarcode(uint64_t barcode_key,
                           uint32_t barcode_length,
                           void *ctx,
                           std::string *out) {
  if (ctx == nullptr || out == nullptr) {
    return false;
  }
  EvidenceBarcodeDecoderCtx *decoder =
      static_cast<EvidenceBarcodeDecoderCtx *>(ctx);
  if (decoder->translator == nullptr) {
    return false;
  }
  *out = decoder->translator->Translate(barcode_key, barcode_length);
  return true;
}

}  // namespace

ChromapAtacResult runChromapAtac(const ChromapAtacConfig &config) {
  const std::string validation_error = validateConfig(config);
  if (!validation_error.empty()) {
    return makeResult(ChromapContractStatus::INVALID_CONFIG, 2,
                      validation_error, config);
  }

  if (!isUnsetToken(config.temp_dir)) {
    std::string mkdir_error;
    if (!ensureDirectory(config.temp_dir, &mkdir_error)) {
      return makeResult(ChromapContractStatus::INVALID_CONFIG, 2,
                        "Chromap temp directory: " + mkdir_error, config);
    }
  }

  chromap::MappingParameters parameters = toChromapParameters(config);
  parameters.permit_acquire_hook = config.permit_hooks.acquire;
  parameters.permit_release_hook = config.permit_hooks.release;
  parameters.permit_hook_ctx = config.permit_hooks.hook_ctx;

  const std::string evidence_out =
      trimCopy(config.atac_evidence_from_peaks_output);
  const bool run_binary_evidence =
      !evidence_out.empty() && evidence_out != "-" &&
      config.call_macs3_frag_peaks &&
      config.macs3_frag_peaks_source == ChromapMacs3FragPeaksSource::MEMORY;

  const chromap::ChromapRunResult chromap_result =
      chromap::RunAtacMapping(parameters);
  if (!chromap_result.ok) {
    return makeResult(ChromapContractStatus::CHROMAP_FAILED,
                      chromap_result.exit_code == 0 ? 1
                                                    : chromap_result.exit_code,
                      chromap_result.message, config);
  }

  if (run_binary_evidence) {
    const std::string sidecar_path = trimCopy(config.fragment_output_path);
    chromap::BarcodeTranslator evidence_barcode_translator;
    if (!isUnsetToken(config.barcode_translate_table)) {
      evidence_barcode_translator.SetTranslateTable(
          trimCopy(config.barcode_translate_table),
          config.barcode_translate_from_first_column);
    }
    EvidenceBarcodeDecoderCtx decoder_ctx;
    decoder_ctx.translator = &evidence_barcode_translator;
    std::ostringstream evidence_log;
    const int rc = libscrna::atac::RunAtacEvidenceFromBinary(
        sidecar_path, parameters.macs3_frag_peaks_narrowpeak_path,
        evidence_out, decodeEvidenceBarcode, &decoder_ctx, &evidence_log);
    if (rc != 0) {
      return makeResult(ChromapContractStatus::CHROMAP_FAILED, rc,
                        "atac_evidence_from_peaks (binary sidecar): " +
                            evidence_log.str(),
                        config);
    }
    return makeResult(ChromapContractStatus::OK, 0,
                      chromap_result.message + "\n" + evidence_log.str(),
                      config);
  }

  return makeResult(ChromapContractStatus::OK, 0, chromap_result.message,
                    config);
}

const char *chromapContractStatusName(ChromapContractStatus status) {
  switch (status) {
    case ChromapContractStatus::OK:
      return "OK";
    case ChromapContractStatus::INVALID_CONFIG:
      return "INVALID_CONFIG";
    case ChromapContractStatus::CHROMAP_FAILED:
      return "CHROMAP_FAILED";
  }
  return "UNKNOWN";
}

}  // namespace multiome
}  // namespace star
