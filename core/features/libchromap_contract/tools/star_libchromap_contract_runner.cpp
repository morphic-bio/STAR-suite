#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "star_chromap_contract.h"

namespace {

void usage() {
  std::cerr
      << "Usage: star_libchromap_contract_runner --ref FASTA --index INDEX "
         "[--read1 CSV --read2 CSV --barcode CSV | "
         "--input-format cbq --read-pair-cbq CSV --barcode-cbq CSV] "
         "--barcode-whitelist FILE "
         "--output FILE [options]\n"
      << "\nOptions:\n"
      << "  --input-format fastq|cbq (default fastq)\n"
      << "  --read-pair-cbq CSV paired-read CBQ file(s), comma separated\n"
      << "  --barcode-cbq CSV   barcode CBQ file(s), comma separated\n"
      << "  --output-format BED|BAM|CRAM|SAM|TagAlign|pairs (default BED)\n"
      << "  --sort-bam            coordinate-sort BAM/CRAM (needs --output-format BAM|CRAM)\n"
      << "  --atac-fragments FILE secondary scATAC fragments path (dual mode; BAM|CRAM only)\n"
      << "  --emit-noY-bam        emit a noY SAM/BAM/CRAM stream\n"
      << "  --noY-output FILE     explicit noY stream output path\n"
      << "  --emit-Y-bam          emit a Y-only SAM/BAM/CRAM stream\n"
      << "  --Y-output FILE       explicit Y stream output path\n"
      << "  --summary FILE\n"
      << "  --call-macs3-frag-peaks\n"
      << "  --macs3-frag-peaks-output FILE\n"
      << "  --macs3-frag-summits-output FILE\n"
      << "  --macs3-frag-pvalue FLOAT\n"
      << "  --macs3-frag-qvalue FLOAT\n"
      << "  --macs3-frag-min-length N\n"
      << "  --macs3-frag-max-gap N\n"
      << "  --macs3-frag-peaks-source file|memory\n"
      << "  --macs3-frag-keep-intermediates DIR\n"
      << "  --macs3-frag-no-uint8-counts\n"
      << "  --macs3-frag-low-mem\n"
      << "  --barcode-translate FILE\n"
      << "  --read-format FORMAT\n"
      << "        Chromap read-format string, e.g. bc:8:23:- to extract and\n"
      << "        reverse-complement ATAC barcode-read bases 9-24.\n"
      << "  --barcode-translate-from-first\n"
      << "        Read translation table as <from_bc>\\t<to_bc> (col1 is the\n"
      << "        hash key / source). Default is the historical Chromap\n"
      << "        convention <to_bc>\\t<from_bc> (col2 is the hash key).\n"
      << "  --temp-dir DIR\n"
      << "  --threads N\n"
      << "  --tn5-shift-mode classical|symmetric\n"
      << "  --tagalign\n"
      << "  --low-mem\n"
      << "  --low-mem-ram BYTES\n"
      << "  --skip-barcode-check\n";
}

std::vector<std::string> splitCsv(const std::string &value) {
  std::vector<std::string> fields;
  std::stringstream stream(value);
  std::string field;
  while (std::getline(stream, field, ',')) {
    if (!field.empty()) {
      fields.push_back(field);
    }
  }
  return fields;
}

bool requireValue(int argc, char **argv, int index) {
  return index + 1 < argc && argv[index + 1][0] != '\0';
}

uint64_t parseUint64(const std::string &value, const std::string &name) {
  char *end = nullptr;
  const unsigned long long parsed = std::strtoull(value.c_str(), &end, 10);
  if (end == value.c_str() || *end != '\0') {
    std::cerr << "Invalid integer for " << name << ": " << value << "\n";
    std::exit(2);
  }
  return static_cast<uint64_t>(parsed);
}

int parseInt(const std::string &value, const std::string &name) {
  char *end = nullptr;
  const long parsed = std::strtol(value.c_str(), &end, 10);
  if (end == value.c_str() || *end != '\0') {
    std::cerr << "Invalid integer for " << name << ": " << value << "\n";
    std::exit(2);
  }
  return static_cast<int>(parsed);
}

double parseDouble(const std::string &value, const std::string &name) {
  char *end = nullptr;
  const double parsed = std::strtod(value.c_str(), &end);
  if (end == value.c_str() || *end != '\0') {
    std::cerr << "Invalid number for " << name << ": " << value << "\n";
    std::exit(2);
  }
  return parsed;
}

}  // namespace

int main(int argc, char **argv) {
  star::multiome::ChromapAtacConfig config;
  bool saw_macs3_frag_pvalue = false;
  bool saw_macs3_frag_qvalue = false;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      usage();
      return 0;
    } else if (arg == "--ref" && requireValue(argc, argv, i)) {
      config.reference_fasta = argv[++i];
    } else if (arg == "--index" && requireValue(argc, argv, i)) {
      config.chromap_index = argv[++i];
    } else if (arg == "--input-format" && requireValue(argc, argv, i)) {
      const std::string f = argv[++i];
      if (f == "fastq" || f == "FASTQ") {
        config.input_format = star::multiome::ChromapInputFormat::FASTQ;
      } else if (f == "cbq" || f == "CBQ" || f == "binseq" ||
                 f == "BINSEQ") {
        config.input_format = star::multiome::ChromapInputFormat::CBQ;
      } else {
        std::cerr << "Invalid --input-format: " << f << "\n";
        return 2;
      }
    } else if (arg == "--read1" && requireValue(argc, argv, i)) {
      config.read1_fastqs = splitCsv(argv[++i]);
    } else if (arg == "--read2" && requireValue(argc, argv, i)) {
      config.read2_fastqs = splitCsv(argv[++i]);
    } else if (arg == "--barcode" && requireValue(argc, argv, i)) {
      config.barcode_fastqs = splitCsv(argv[++i]);
    } else if (arg == "--read-pair-cbq" && requireValue(argc, argv, i)) {
      config.read_pair_cbqs = splitCsv(argv[++i]);
    } else if (arg == "--barcode-cbq" && requireValue(argc, argv, i)) {
      config.barcode_cbqs = splitCsv(argv[++i]);
    } else if (arg == "--barcode-whitelist" && requireValue(argc, argv, i)) {
      config.barcode_whitelist = argv[++i];
    } else if (arg == "--barcode-translate" && requireValue(argc, argv, i)) {
      config.barcode_translate_table = argv[++i];
    } else if (arg == "--read-format" && requireValue(argc, argv, i)) {
      config.read_format = argv[++i];
    } else if (arg == "--barcode-translate-from-first") {
      config.barcode_translate_from_first_column = true;
    } else if (arg == "--output" && requireValue(argc, argv, i)) {
      config.output_path = argv[++i];
    } else if (arg == "--output-format" && requireValue(argc, argv, i)) {
      const std::string f = argv[++i];
      if (f == "BED" || f == "bed" || f == "fragments") {
        config.output_format = star::multiome::ChromapOutputFormat::BED;
      } else if (f == "TagAlign" || f == "tagalign") {
        config.output_format = star::multiome::ChromapOutputFormat::TAGALIGN;
      } else if (f == "SAM" || f == "sam") {
        config.output_format = star::multiome::ChromapOutputFormat::SAM;
      } else if (f == "BAM" || f == "bam") {
        config.output_format = star::multiome::ChromapOutputFormat::BAM;
      } else if (f == "CRAM" || f == "cram") {
        config.output_format = star::multiome::ChromapOutputFormat::CRAM;
      } else if (f == "pairs") {
        config.output_format = star::multiome::ChromapOutputFormat::PAIRS;
      } else {
        std::cerr << "Invalid --output-format: " << f << "\n";
        return 2;
      }
    } else if (arg == "--sort-bam") {
      config.sort_bam = true;
    } else if (arg == "--atac-fragments" && requireValue(argc, argv, i)) {
      config.fragment_output_path = argv[++i];
    } else if (arg == "--emit-noY-bam") {
      config.emit_no_y_bam = true;
    } else if (arg == "--noY-output" && requireValue(argc, argv, i)) {
      config.no_y_output_path = argv[++i];
    } else if (arg == "--emit-Y-bam") {
      config.emit_y_bam = true;
    } else if (arg == "--Y-output" && requireValue(argc, argv, i)) {
      config.y_output_path = argv[++i];
    } else if (arg == "--summary" && requireValue(argc, argv, i)) {
      config.summary_path = argv[++i];
    } else if (arg == "--call-macs3-frag-peaks") {
      config.call_macs3_frag_peaks = true;
    } else if (arg == "--macs3-frag-peaks-output" && requireValue(argc, argv, i)) {
      config.macs3_frag_peaks_output = argv[++i];
    } else if (arg == "--macs3-frag-summits-output" && requireValue(argc, argv, i)) {
      config.macs3_frag_summits_output = argv[++i];
    } else if (arg == "--macs3-frag-pvalue" && requireValue(argc, argv, i)) {
      saw_macs3_frag_pvalue = true;
      config.macs3_frag_pvalue = parseDouble(argv[++i], "--macs3-frag-pvalue");
    } else if (arg == "--macs3-frag-qvalue" && requireValue(argc, argv, i)) {
      saw_macs3_frag_qvalue = true;
      config.macs3_frag_threshold_mode =
          star::multiome::ChromapMacs3FragThresholdMode::Q_VALUE;
      config.macs3_frag_qvalue = parseDouble(argv[++i], "--macs3-frag-qvalue");
    } else if (arg == "--macs3-frag-min-length" && requireValue(argc, argv, i)) {
      config.macs3_frag_min_length = parseInt(argv[++i], "--macs3-frag-min-length");
    } else if (arg == "--macs3-frag-max-gap" && requireValue(argc, argv, i)) {
      config.macs3_frag_max_gap = parseInt(argv[++i], "--macs3-frag-max-gap");
    } else if (arg == "--macs3-frag-peaks-source" && requireValue(argc, argv, i)) {
      const std::string source = argv[++i];
      if (source == "file") {
        config.macs3_frag_peaks_source =
            star::multiome::ChromapMacs3FragPeaksSource::FILE;
      } else if (source == "memory") {
        config.macs3_frag_peaks_source =
            star::multiome::ChromapMacs3FragPeaksSource::MEMORY;
      } else {
        std::cerr << "Invalid --macs3-frag-peaks-source: " << source << "\n";
        return 2;
      }
    } else if (arg == "--macs3-frag-keep-intermediates" && requireValue(argc, argv, i)) {
      config.macs3_frag_keep_intermediates_dir = argv[++i];
    } else if (arg == "--macs3-frag-no-uint8-counts") {
      config.macs3_frag_uint8_counts = false;
    } else if (arg == "--macs3-frag-low-mem") {
      config.macs3_frag_low_mem = true;
    } else if (arg == "--temp-dir" && requireValue(argc, argv, i)) {
      config.temp_dir = argv[++i];
    } else if (arg == "--threads" && requireValue(argc, argv, i)) {
      config.threads = parseInt(argv[++i], "--threads");
    } else if (arg == "--tn5-shift-mode" && requireValue(argc, argv, i)) {
      const std::string mode = argv[++i];
      if (mode == "classical") {
        config.tn5_shift_mode = star::multiome::Tn5ShiftMode::CLASSICAL;
      } else if (mode == "symmetric") {
        config.tn5_shift_mode = star::multiome::Tn5ShiftMode::SYMMETRIC;
      } else {
        std::cerr << "Invalid --tn5-shift-mode: " << mode << "\n";
        return 2;
      }
    } else if (arg == "--tagalign") {
      config.output_format = star::multiome::ChromapOutputFormat::TAGALIGN;
    } else if (arg == "--low-mem") {
      config.low_memory_mode = true;
    } else if (arg == "--low-mem-ram" && requireValue(argc, argv, i)) {
      config.low_memory_ram_limit = parseUint64(argv[++i], "--low-mem-ram");
    } else if (arg == "--skip-barcode-check") {
      config.skip_barcode_check = true;
    } else {
      std::cerr << "Unknown or incomplete option: " << arg << "\n";
      usage();
      return 2;
    }
  }

  if (saw_macs3_frag_pvalue && saw_macs3_frag_qvalue) {
    std::cerr << "Config error: --macs3-frag-pvalue and --macs3-frag-qvalue "
              << "are mutually exclusive\n";
    return 2;
  }
  if (saw_macs3_frag_pvalue &&
      !(config.macs3_frag_pvalue > 0.0 &&
        config.macs3_frag_pvalue <= 1.0)) {
    std::cerr << "Config error: --macs3-frag-pvalue must be in (0, 1]\n";
    return 2;
  }
  if (saw_macs3_frag_qvalue &&
      !(config.macs3_frag_qvalue > 0.0 &&
        config.macs3_frag_qvalue <= 1.0)) {
    std::cerr << "Config error: --macs3-frag-qvalue must be in (0, 1]\n";
    return 2;
  }

  if (!config.fragment_output_path.empty()) {
    if (config.output_format != star::multiome::ChromapOutputFormat::BAM &&
        config.output_format != star::multiome::ChromapOutputFormat::CRAM) {
      std::cerr
          << "Config error: --atac-fragments requires --output-format BAM or CRAM\n";
      return 2;
    }
    if (config.fragment_output_path == config.output_path) {
      std::cerr << "Config error: --atac-fragments must differ from --output\n";
      return 2;
    }
  }

  const star::multiome::ChromapAtacResult result =
      star::multiome::runChromapAtac(config);

  std::cerr << "status="
            << star::multiome::chromapContractStatusName(result.status)
            << "\n";
  if (!result.message.empty()) {
    std::cerr << "message=" << result.message << "\n";
  }
  if (!result.output_path.empty()) {
    std::cerr << "output=" << result.output_path << "\n";
  }
  if (!result.summary_path.empty()) {
    std::cerr << "summary=" << result.summary_path << "\n";
  }

  return result.exit_code;
}
