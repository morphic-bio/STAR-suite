#include "multiome_atac_peak_mex.h"

#include <cstdint>
#include <cstdlib>
#include <cctype>
#include <iostream>
#include <string>

namespace {

void usage(const char *prog) {
  std::cerr
      << "Usage: " << prog
      << " (--fragments fragments.tsv[.gz] | --sidecar atac_fragments.bin) "
         "--peaks peaks.narrowPeak --out-dir peak_mex "
         "--metrics-tsv atac_metrics.tsv [options]\n"
      << "\nOptions:\n"
      << "  --summits-out PATH              Summits BED output when calling peaks\n"
      << "  --call-peaks-from-sidecar       Call MACS3-compatible FRAG peaks from "
         "the binary sidecar before building MEX\n"
      << "  --peak-call-mode MODE           Peak caller: frag (default) or "
         "macs-bed\n"
      << "  --macs-profile NAME             MACS BED profile, e.g. signac-atac; "
         "implies --peak-call-mode macs-bed\n"
      << "  --barcode-translate PATH        Optional source-to-destination barcode "
         "translation table for sidecar input\n"
      << "  --barcode-translate-from-first  Treat column 1 as source and column 2 "
         "as destination\n"
      << "  --threads N                     Peak-caller threads (default 1)\n"
      << "  --macs3-frag-pvalue P           MACS3 FRAG p-value (default 1e-5)\n"
      << "  --macs3-frag-qvalue Q           MACS3 FRAG q-value/FDR threshold\n"
      << "  --macs3-frag-min-length N       Peak min length (default 200)\n"
      << "  --macs3-frag-max-gap N          Peak max gap (default 30)\n"
      << "  --macs3-frag-no-uint8-counts    Disable uint8 count workspace\n"
      << "  --temp-dir PATH                 MACS3 temporary parent directory\n"
      << "  --keep-intermediates-dir PATH   Keep MACS3 intermediate bedGraphs\n"
      << "  --max-barcodes N                Keep only top N barcodes by total "
         "fragments; 0 keeps all (default 0)\n"
      << "  --force                         Rebuild peaks even if output paths exist\n";
}

bool take_value(const std::string &arg, int argc, char **argv, int *index,
                std::string *out) {
  if (*index + 1 >= argc) {
    std::cerr << "Missing value for " << arg << "\n";
    return false;
  }
  *out = argv[++(*index)];
  return true;
}

bool parse_uint64(const std::string &value, uint64_t *out) {
  char *end = nullptr;
  const unsigned long long parsed = std::strtoull(value.c_str(), &end, 10);
  if (value.empty() || end == value.c_str() || *end != '\0') {
    return false;
  }
  *out = static_cast<uint64_t>(parsed);
  return true;
}

bool parse_int(const std::string &value, int *out) {
  char *end = nullptr;
  const long parsed = std::strtol(value.c_str(), &end, 10);
  if (value.empty() || end == value.c_str() || *end != '\0') {
    return false;
  }
  *out = static_cast<int>(parsed);
  return true;
}

bool parse_double(const std::string &value, double *out) {
  char *end = nullptr;
  const double parsed = std::strtod(value.c_str(), &end);
  if (value.empty() || end == value.c_str() || *end != '\0') {
    return false;
  }
  *out = parsed;
  return true;
}

std::string lower_copy(std::string value) {
  for (size_t i = 0; i < value.size(); ++i) {
    value[i] =
        static_cast<char>(std::tolower(static_cast<unsigned char>(value[i])));
  }
  return value;
}

bool parse_peak_call_mode(const std::string &value,
                          star::multiome::MultiomeAtacPeakCallMode *out) {
  if (out == nullptr) {
    return false;
  }
  const std::string mode = lower_copy(value);
  if (mode == "frag" || mode == "macs3-frag" || mode == "macs3_frag") {
    *out = star::multiome::MultiomeAtacPeakCallMode::FRAG;
    return true;
  }
  if (mode == "macs-bed" || mode == "macs_bed" || mode == "bed") {
    *out = star::multiome::MultiomeAtacPeakCallMode::MACS_BED;
    return true;
  }
  return false;
}

bool parse_args(int argc, char **argv,
                star::multiome::MultiomeAtacPeakMexArgs *args, bool *help) {
  if (args == nullptr || help == nullptr) {
    return false;
  }
  bool saw_macs3_pvalue = false;
  bool saw_macs3_qvalue = false;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      usage(argv[0]);
      *help = true;
      return true;
    } else if (arg == "--fragments") {
      if (!take_value(arg, argc, argv, &i, &args->fragments)) return false;
    } else if (arg == "--sidecar") {
      if (!take_value(arg, argc, argv, &i, &args->sidecar)) return false;
    } else if (arg == "--peaks") {
      if (!take_value(arg, argc, argv, &i, &args->peaks)) return false;
    } else if (arg == "--summits-out") {
      if (!take_value(arg, argc, argv, &i, &args->summits_out)) return false;
    } else if (arg == "--out-dir") {
      if (!take_value(arg, argc, argv, &i, &args->out_dir)) return false;
    } else if (arg == "--metrics-tsv") {
      if (!take_value(arg, argc, argv, &i, &args->metrics_tsv)) return false;
    } else if (arg == "--barcode-translate") {
      if (!take_value(arg, argc, argv, &i, &args->barcode_translate)) {
        return false;
      }
    } else if (arg == "--temp-dir") {
      if (!take_value(arg, argc, argv, &i, &args->temp_dir)) return false;
    } else if (arg == "--keep-intermediates-dir") {
      if (!take_value(arg, argc, argv, &i, &args->keep_intermediates_dir)) {
        return false;
      }
    } else if (arg == "--barcode-translate-from-first") {
      args->barcode_translate_from_first = true;
    } else if (arg == "--call-peaks-from-sidecar") {
      args->call_peaks_from_sidecar = true;
    } else if (arg == "--peak-call-mode") {
      std::string value;
      if (!take_value(arg, argc, argv, &i, &value) ||
          !parse_peak_call_mode(value, &args->peak_call_mode)) {
        std::cerr << "Invalid --peak-call-mode value\n";
        return false;
      }
    } else if (arg == "--macs-profile") {
      if (!take_value(arg, argc, argv, &i, &args->macs_profile)) {
        return false;
      }
      args->peak_call_mode =
          star::multiome::MultiomeAtacPeakCallMode::MACS_BED;
    } else if (arg == "--force") {
      args->force = true;
    } else if (arg == "--macs3-frag-no-uint8-counts") {
      args->macs3_uint8_counts = false;
    } else if (arg == "--max-barcodes") {
      std::string value;
      if (!take_value(arg, argc, argv, &i, &value) ||
          !parse_uint64(value, &args->max_barcodes)) {
        std::cerr << "Invalid --max-barcodes value\n";
        return false;
      }
    } else if (arg == "--threads") {
      std::string value;
      if (!take_value(arg, argc, argv, &i, &value) ||
          !parse_int(value, &args->threads) || args->threads < 1) {
        std::cerr << "Invalid --threads value\n";
        return false;
      }
    } else if (arg == "--macs3-frag-pvalue") {
      std::string value;
      saw_macs3_pvalue = true;
      if (!take_value(arg, argc, argv, &i, &value) ||
          !parse_double(value, &args->macs3_pvalue) ||
          !(args->macs3_pvalue > 0.0 && args->macs3_pvalue <= 1.0)) {
        std::cerr << "Invalid --macs3-frag-pvalue value\n";
        return false;
      }
      args->macs3_threshold_mode =
          star::multiome::MultiomeAtacPeakMexThresholdMode::P_VALUE;
      args->macs3_threshold_explicit = true;
    } else if (arg == "--macs3-frag-qvalue") {
      std::string value;
      saw_macs3_qvalue = true;
      if (!take_value(arg, argc, argv, &i, &value) ||
          !parse_double(value, &args->macs3_qvalue) ||
          !(args->macs3_qvalue > 0.0 && args->macs3_qvalue <= 1.0)) {
        std::cerr << "Invalid --macs3-frag-qvalue value\n";
        return false;
      }
      args->macs3_threshold_mode =
          star::multiome::MultiomeAtacPeakMexThresholdMode::Q_VALUE;
      args->macs3_threshold_explicit = true;
    } else if (arg == "--macs3-frag-min-length") {
      std::string value;
      if (!take_value(arg, argc, argv, &i, &value) ||
          !parse_int(value, &args->macs3_min_length) ||
          args->macs3_min_length < 1) {
        std::cerr << "Invalid --macs3-frag-min-length value\n";
        return false;
      }
      args->macs3_min_length_explicit = true;
    } else if (arg == "--macs3-frag-max-gap") {
      std::string value;
      if (!take_value(arg, argc, argv, &i, &value) ||
          !parse_int(value, &args->macs3_max_gap) ||
          args->macs3_max_gap < 0) {
        std::cerr << "Invalid --macs3-frag-max-gap value\n";
        return false;
      }
      args->macs3_max_gap_explicit = true;
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      usage(argv[0]);
      return false;
    }
  }
  if (saw_macs3_pvalue && saw_macs3_qvalue) {
    std::cerr << "--macs3-frag-pvalue and --macs3-frag-qvalue are mutually "
              << "exclusive\n";
    return false;
  }
  return true;
}

}  // namespace

int main(int argc, char **argv) {
  star::multiome::MultiomeAtacPeakMexArgs args;
  bool help = false;
  if (!parse_args(argc, argv, &args, &help)) {
    return 2;
  }
  if (help) {
    return 0;
  }
  return star::multiome::RunMultiomeAtacPeakMex(args);
}
