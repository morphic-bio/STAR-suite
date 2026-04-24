#include "star_chromap_orchestration.h"

#if !WITH_CHROMAP
#error "star_chromap_orchestration.cpp must be compiled with WITH_CHROMAP=1"
#endif

#include "IncludeDefine.h"
#include "Parameters.h"
#include "TimeFunctions.h"
#include "star_chromap_contract.h"

#include <cctype>
#include <exception>
#include <thread>

namespace {

std::string lowerCopy(std::string s) {
  for (size_t i = 0; i < s.size(); ++i) {
    s[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(s[i])));
  }
  return s;
}

std::string trimCopy(const std::string &input) {
  size_t first = input.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) {
    return "";
  }
  size_t last = input.find_last_not_of(" \t\r\n");
  return input.substr(first, last - first + 1);
}

std::vector<std::string> splitCsvPaths(const std::string &value) {
  std::vector<std::string> fields;
  std::stringstream stream(value);
  std::string field;
  while (std::getline(stream, field, ',')) {
    fields.push_back(trimCopy(field));
  }
  if (!value.empty() && value[value.size() - 1] == ',') {
    fields.push_back("");
  }
  return fields;
}

bool isUnsetToken(const std::string &input) {
  std::string t = lowerCopy(trimCopy(input));
  if (t.empty()) {
    return true;
  }
  return t == "-" || t == "none";
}

bool isConcurrentStartMode(const std::string &input) {
  return lowerCopy(trimCopy(input)) == "concurrent";
}

bool isPostMappingStartMode(const std::string &input) {
  const std::string mode = lowerCopy(trimCopy(input));
  return mode == "postmapping" || mode == "post_mapping" ||
         mode == "sequential";
}

bool hasUnsetPath(const std::vector<std::string> &paths) {
  for (size_t i = 0; i < paths.size(); ++i) {
    if (isUnsetToken(paths[i])) {
      return true;
    }
  }
  return false;
}

bool validateAndBuildConfig(Parameters &P,
                            bool batchModeActive,
                            star::multiome::ChromapAtacConfig *cfg) {
  if (P.chromapAtac.enabled == 0) {
    return true;
  }
  if (batchModeActive) {
    P.inOut->logMain
        << "ERROR: --chromapAtacEnable 1 is not supported in batch mode. "
           "Run Chromap separately or disable batch mode for this integration gate.\n";
    return false;
  }

  if (P.chromapAtac.threads < 1) {
    P.inOut->logMain << "ERROR: --chromapAtacThreads must be >= 1 (got "
                     << P.chromapAtac.threads << ")\n";
    return false;
  }
  if (P.chromapAtac.htsThreads < 0) {
    P.inOut->logMain << "ERROR: --chromapAtacHtsThreads must be >= 0 (got "
                     << P.chromapAtac.htsThreads << ")\n";
    return false;
  }

  if (!isConcurrentStartMode(P.chromapAtac.startMode) &&
      !isPostMappingStartMode(P.chromapAtac.startMode)) {
    P.inOut->logMain
        << "ERROR: --chromapAtacStartMode must be postMapping or concurrent (got \""
        << P.chromapAtac.startMode << "\")\n";
    return false;
  }

  const std::string tn5Raw = trimCopy(P.chromapAtac.tn5ShiftMode);
  const std::string tn5 = lowerCopy(tn5Raw);
  star::multiome::Tn5ShiftMode tn5Mode;
  if (tn5 == "classical") {
    tn5Mode = star::multiome::Tn5ShiftMode::CLASSICAL;
  } else if (tn5 == "symmetric") {
    tn5Mode = star::multiome::Tn5ShiftMode::SYMMETRIC;
  } else {
    P.inOut->logMain
        << "ERROR: --chromapAtacTn5ShiftMode must be classical or symmetric (got \""
        << tn5Raw << "\")\n";
    return false;
  }

  const std::string formatRaw = trimCopy(P.chromapAtac.outputFormat);
  const std::string format = lowerCopy(formatRaw);
  star::multiome::ChromapOutputFormat outputFormat;
  if (format == "bed" || format == "fragments") {
    outputFormat = star::multiome::ChromapOutputFormat::BED;
  } else if (format == "tagalign") {
    outputFormat = star::multiome::ChromapOutputFormat::TAGALIGN;
  } else if (format == "sam") {
    outputFormat = star::multiome::ChromapOutputFormat::SAM;
  } else if (format == "bam") {
    outputFormat = star::multiome::ChromapOutputFormat::BAM;
  } else if (format == "cram") {
    outputFormat = star::multiome::ChromapOutputFormat::CRAM;
  } else if (format == "pairs") {
    outputFormat = star::multiome::ChromapOutputFormat::PAIRS;
  } else {
    P.inOut->logMain
        << "ERROR: --chromapAtacOutputFormat must be BED, fragments, TagAlign, SAM, BAM, CRAM, or pairs (got \""
        << formatRaw << "\")\n";
    return false;
  }
  if (P.chromapAtac.sortBam != 0 &&
      outputFormat != star::multiome::ChromapOutputFormat::BAM &&
      outputFormat != star::multiome::ChromapOutputFormat::CRAM) {
    P.inOut->logMain
        << "ERROR: --chromapAtacSortBam 1 requires --chromapAtacOutputFormat BAM or CRAM\n";
    return false;
  }
  if (P.chromapAtac.writeIndex != 0 && P.chromapAtac.sortBam == 0) {
    P.inOut->logMain
        << "ERROR: --chromapAtacWriteIndex 1 requires --chromapAtacSortBam 1\n";
    return false;
  }

  if (!isUnsetToken(P.chromapAtac.secondaryFragments)) {
    if (outputFormat != star::multiome::ChromapOutputFormat::BAM &&
        outputFormat != star::multiome::ChromapOutputFormat::CRAM) {
      P.inOut->logMain
          << "ERROR: --chromapAtacSecondaryFragments requires --chromapAtacOutputFormat BAM or CRAM\n";
      return false;
    }
    const std::string primaryOut = trimCopy(P.chromapAtac.outputFragments);
    const std::string secondaryOut = trimCopy(P.chromapAtac.secondaryFragments);
    if (primaryOut == secondaryOut) {
      P.inOut->logMain
          << "ERROR: --chromapAtacSecondaryFragments must differ from --chromapAtacOutputFragments\n";
      return false;
    }
  }

  if (isUnsetToken(P.chromapAtac.referenceFasta)) {
    P.inOut->logMain
        << "ERROR: --chromapAtacReferenceFasta is required when --chromapAtacEnable 1\n";
    return false;
  }
  if (isUnsetToken(P.chromapAtac.chromapIndex)) {
    P.inOut->logMain << "ERROR: --chromapAtacIndex is required when --chromapAtacEnable 1\n";
    return false;
  }
  if (isUnsetToken(P.chromapAtac.read1Csv)) {
    P.inOut->logMain << "ERROR: --chromapAtacRead1 is required when --chromapAtacEnable 1\n";
    return false;
  }
  if (isUnsetToken(P.chromapAtac.read2Csv)) {
    P.inOut->logMain << "ERROR: --chromapAtacRead2 is required when --chromapAtacEnable 1\n";
    return false;
  }
  if (isUnsetToken(P.chromapAtac.barcodeCsv)) {
    P.inOut->logMain << "ERROR: --chromapAtacBarcode is required when --chromapAtacEnable 1\n";
    return false;
  }
  if (isUnsetToken(P.chromapAtac.barcodeWhitelist)) {
    P.inOut->logMain
        << "ERROR: --chromapAtacBarcodeWhitelist is required when --chromapAtacEnable 1\n";
    return false;
  }
  if (isUnsetToken(P.chromapAtac.outputFragments)) {
    P.inOut->logMain
        << "ERROR: --chromapAtacOutputFragments is required when --chromapAtacEnable 1\n";
    return false;
  }

  const std::vector<std::string> r1 = splitCsvPaths(P.chromapAtac.read1Csv);
  const std::vector<std::string> r2 = splitCsvPaths(P.chromapAtac.read2Csv);
  const std::vector<std::string> bc = splitCsvPaths(P.chromapAtac.barcodeCsv);
  if (r1.empty() || r2.empty() || bc.empty() ||
      hasUnsetPath(r1) || hasUnsetPath(r2) || hasUnsetPath(bc)) {
    P.inOut->logMain
        << "ERROR: --chromapAtacRead1, --chromapAtacRead2, and --chromapAtacBarcode must "
           "list at least one non-empty FASTQ path each\n";
    return false;
  }
  if (r1.size() != r2.size() || r1.size() != bc.size()) {
    P.inOut->logMain
        << "ERROR: --chromapAtacRead1, --chromapAtacRead2, and --chromapAtacBarcode "
           "must list the same number of FASTQ paths\n";
    return false;
  }

  if (cfg == nullptr) {
    return true;
  }

  cfg->reference_fasta = trimCopy(P.chromapAtac.referenceFasta);
  cfg->chromap_index = trimCopy(P.chromapAtac.chromapIndex);
  cfg->read1_fastqs = r1;
  cfg->read2_fastqs = r2;
  cfg->barcode_fastqs = bc;
  cfg->barcode_whitelist = trimCopy(P.chromapAtac.barcodeWhitelist);
  if (!isUnsetToken(P.chromapAtac.barcodeTranslate)) {
    cfg->barcode_translate_table = trimCopy(P.chromapAtac.barcodeTranslate);
  }
  cfg->output_path = trimCopy(P.chromapAtac.outputFragments);
  cfg->fragment_output_path.clear();
  if (!isUnsetToken(P.chromapAtac.secondaryFragments)) {
    cfg->fragment_output_path = trimCopy(P.chromapAtac.secondaryFragments);
  }
  if (!isUnsetToken(P.chromapAtac.summary)) {
    cfg->summary_path = trimCopy(P.chromapAtac.summary);
  }
  if (!isUnsetToken(P.chromapAtac.tempDir)) {
    cfg->temp_dir = trimCopy(P.chromapAtac.tempDir);
  }
  cfg->threads = P.chromapAtac.threads;
  cfg->hts_threads = P.chromapAtac.htsThreads;
  cfg->tn5_shift_mode = tn5Mode;
  cfg->output_format = outputFormat;
  cfg->sort_bam = P.chromapAtac.sortBam != 0;
  cfg->write_index = P.chromapAtac.writeIndex != 0;
  cfg->sort_bam_ram_limit = P.chromapAtac.sortBamRam;
  return true;
}

}  // namespace

struct StarChromapAtacAsyncRun::Impl {
  explicit Impl(const star::multiome::ChromapAtacConfig &config)
      : cfg(config),
        worker([this]() {
          try {
            result = star::multiome::runChromapAtac(cfg);
          } catch (...) {
            exception = std::current_exception();
          }
        }) {}

  star::multiome::ChromapAtacConfig cfg;
  star::multiome::ChromapAtacResult result;
  std::exception_ptr exception;
  std::thread worker;
};

StarChromapAtacAsyncRun::StarChromapAtacAsyncRun() : impl_(nullptr) {}

StarChromapAtacAsyncRun::~StarChromapAtacAsyncRun() {
  if (impl_ != nullptr) {
    if (impl_->worker.joinable()) {
      impl_->worker.join();
    }
    delete impl_;
    impl_ = nullptr;
  }
}

bool preflightStarChromapAtacIfEnabled(Parameters &P, bool batchModeActive) {
  return validateAndBuildConfig(P, batchModeActive, nullptr);
}

bool startStarChromapAtacIfEnabled(Parameters &P,
                                   bool batchModeActive,
                                   StarChromapAtacAsyncRun &run) {
  star::multiome::ChromapAtacConfig cfg;
  if (!validateAndBuildConfig(P, batchModeActive, &cfg)) {
    return false;
  }
  if (P.chromapAtac.enabled == 0) {
    return true;
  }
  if (!isConcurrentStartMode(P.chromapAtac.startMode)) {
    return true;
  }
  if (run.impl_ != nullptr) {
    P.inOut->logMain
        << "ERROR: Chromap ATAC concurrent run was already started\n";
    return false;
  }

  P.inOut->logMain << timeMonthDayTime()
                   << " ..... starting concurrent in-process Chromap ATAC "
                      "(libchromap contract)\n"
                   << flush;

  try {
    run.impl_ = new StarChromapAtacAsyncRun::Impl(cfg);
  } catch (const std::exception &error) {
    P.inOut->logMain << "ERROR: failed to start concurrent Chromap ATAC: "
                     << error.what() << "\n";
    return false;
  } catch (...) {
    P.inOut->logMain
        << "ERROR: failed to start concurrent Chromap ATAC: unknown exception\n";
    return false;
  }

  return true;
}

bool runStarChromapAtacIfEnabled(Parameters &P,
                                 bool batchModeActive,
                                 StarChromapAtacAsyncRun &run) {
  if (run.impl_ != nullptr) {
    P.inOut->logMain << timeMonthDayTime()
                     << " ..... waiting for concurrent in-process Chromap ATAC\n"
                     << flush;

    if (run.impl_->worker.joinable()) {
      run.impl_->worker.join();
    }

    if (run.impl_->exception) {
      try {
        std::rethrow_exception(run.impl_->exception);
      } catch (const std::exception &error) {
        P.inOut->logMain
            << "ERROR: concurrent Chromap ATAC threw exception: "
            << error.what() << "\n";
      } catch (...) {
        P.inOut->logMain
            << "ERROR: concurrent Chromap ATAC threw unknown exception\n";
      }
      delete run.impl_;
      run.impl_ = nullptr;
      return false;
    }

    const star::multiome::ChromapAtacResult result = run.impl_->result;
    delete run.impl_;
    run.impl_ = nullptr;

    if (result.status != star::multiome::ChromapContractStatus::OK) {
      P.inOut->logMain << "ERROR: Chromap ATAC contract failed: status="
                       << star::multiome::chromapContractStatusName(result.status)
                       << " exit_code=" << result.exit_code;
      if (!result.message.empty()) {
        P.inOut->logMain << " message=" << result.message;
      }
      P.inOut->logMain << "\n";
      return false;
    }

    P.inOut->logMain << timeMonthDayTime()
                     << " ..... finished concurrent in-process Chromap ATAC successfully\n"
                     << flush;
    return true;
  }

  star::multiome::ChromapAtacConfig cfg;
  if (!validateAndBuildConfig(P, batchModeActive, &cfg)) {
    return false;
  }
  if (P.chromapAtac.enabled == 0) {
    return true;
  }
  if (isConcurrentStartMode(P.chromapAtac.startMode)) {
    P.inOut->logMain
        << "ERROR: --chromapAtacStartMode concurrent was requested, but Chromap "
           "was not started before STAR mapping\n";
    return false;
  }

  P.inOut->logMain << timeMonthDayTime()
                   << " ..... starting in-process Chromap ATAC (libchromap contract)\n"
                   << flush;

  const star::multiome::ChromapAtacResult result =
      star::multiome::runChromapAtac(cfg);

  if (result.status != star::multiome::ChromapContractStatus::OK) {
    P.inOut->logMain << "ERROR: Chromap ATAC contract failed: status="
                     << star::multiome::chromapContractStatusName(result.status)
                     << " exit_code=" << result.exit_code;
    if (!result.message.empty()) {
      P.inOut->logMain << " message=" << result.message;
    }
    P.inOut->logMain << "\n";
    return false;
  }

  P.inOut->logMain << timeMonthDayTime()
                   << " ..... finished in-process Chromap ATAC successfully\n"
                   << flush;
  return true;
}
