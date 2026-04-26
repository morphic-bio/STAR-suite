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
         "--read1 CSV --read2 CSV --barcode CSV --barcode-whitelist FILE "
         "--output FILE [options]\n"
      << "\nOptions:\n"
      << "  --output-format BED|BAM|CRAM|SAM|TagAlign|pairs (default BED)\n"
      << "  --sort-bam            coordinate-sort BAM/CRAM (needs --output-format BAM|CRAM)\n"
      << "  --atac-fragments FILE secondary scATAC fragments path (dual mode; BAM|CRAM only)\n"
      << "  --summary FILE\n"
      << "  --barcode-translate FILE\n"
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

}  // namespace

int main(int argc, char **argv) {
  star::multiome::ChromapAtacConfig config;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      usage();
      return 0;
    } else if (arg == "--ref" && requireValue(argc, argv, i)) {
      config.reference_fasta = argv[++i];
    } else if (arg == "--index" && requireValue(argc, argv, i)) {
      config.chromap_index = argv[++i];
    } else if (arg == "--read1" && requireValue(argc, argv, i)) {
      config.read1_fastqs = splitCsv(argv[++i]);
    } else if (arg == "--read2" && requireValue(argc, argv, i)) {
      config.read2_fastqs = splitCsv(argv[++i]);
    } else if (arg == "--barcode" && requireValue(argc, argv, i)) {
      config.barcode_fastqs = splitCsv(argv[++i]);
    } else if (arg == "--barcode-whitelist" && requireValue(argc, argv, i)) {
      config.barcode_whitelist = argv[++i];
    } else if (arg == "--barcode-translate" && requireValue(argc, argv, i)) {
      config.barcode_translate_table = argv[++i];
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
    } else if (arg == "--summary" && requireValue(argc, argv, i)) {
      config.summary_path = argv[++i];
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
