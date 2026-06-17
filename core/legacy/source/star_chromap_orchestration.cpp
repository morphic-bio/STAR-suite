#include "star_chromap_orchestration.h"

#if !WITH_CHROMAP
#error "star_chromap_orchestration.cpp must be compiled with WITH_CHROMAP=1"
#endif

#include "GlobalVariables.h"
#include "IncludeDefine.h"
#include "Parameters.h"
#include "ThreadControl.h"
#include "TimeFunctions.h"
#include "multiome_atac_peak_mex.h"
#include "star_chromap_contract.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <thread>

namespace {

// Permit-hook shims for the chromap libchromap_contract that mirror the
// process_features pf_api shape (PfMultiAssign.cpp). Each chromap PE worker
// thread acquires one permit per mini-batch (~64 mapped pairs) via these
// functions, which forward to ThreadControl's PermitDomain::ATAC counters.
ThreadControl::PermitHookContext kAtacPermitHookContext{
    ThreadControl::PermitDomain::ATAC};

extern "C" uint64_t chromapStarDynamicPermitAcquire(void *hookCtx) {
  const ThreadControl::PermitHookContext *permitCtx =
      static_cast<const ThreadControl::PermitHookContext *>(hookCtx);
  const ThreadControl::PermitDomain domain =
      (permitCtx == nullptr) ? ThreadControl::PermitDomain::ATAC
                             : permitCtx->domain;
  return g_threadChunks.mapPermitAcquireForDomain(domain);
}

extern "C" void chromapStarDynamicPermitRelease(
    void *hookCtx,
    uint64_t waitNs,
    uint64_t workUnits,
    uint64_t workBytes,
    uint64_t workNs) {
  const ThreadControl::PermitHookContext *permitCtx =
      static_cast<const ThreadControl::PermitHookContext *>(hookCtx);
  const ThreadControl::PermitDomain domain =
      (permitCtx == nullptr) ? ThreadControl::PermitDomain::ATAC
                             : permitCtx->domain;
  g_threadChunks.mapPermitReleaseForDomain(domain, waitNs, workUnits, workBytes,
                                           workNs);
}

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

int parameterInputLevel(const Parameters &P, const std::string &name) {
  for (size_t i = 0; i < P.parArray.size(); ++i) {
    if (P.parArray[i] != nullptr && P.parArray[i]->nameString == name) {
      return P.parArray[i]->inputLevel;
    }
  }
  return -1;
}

bool parseYesNo(const std::string &input, bool *out) {
  if (out == nullptr) {
    return false;
  }
  const std::string t = lowerCopy(trimCopy(input));
  if (t == "yes" || t == "y" || t == "true" || t == "1") {
    *out = true;
    return true;
  }
  if (t == "no" || t == "n" || t == "false" || t == "0") {
    *out = false;
    return true;
  }
  return false;
}

bool parsePeakCallMode(const std::string &input,
                       star::multiome::MultiomeAtacPeakCallMode *out) {
  if (out == nullptr) {
    return false;
  }
  const std::string mode = lowerCopy(trimCopy(input));
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

const char *peakCallModeName(star::multiome::MultiomeAtacPeakCallMode mode) {
  return mode == star::multiome::MultiomeAtacPeakCallMode::MACS_BED
             ? "macs-bed"
             : "frag";
}

bool resolveInlinePeakCallMode(
    const Parameters &P,
    star::multiome::MultiomeAtacPeakCallMode *mode,
    std::string *profile,
    std::ostream &log) {
  if (mode == nullptr || profile == nullptr) {
    return false;
  }
  if (!parsePeakCallMode(P.multiomeAtacPeakMex.peakCallMode, mode)) {
    log << "ERROR: --multiomeAtacPeakCallMode must be frag or macs-bed (got \""
        << P.multiomeAtacPeakMex.peakCallMode << "\")\n";
    return false;
  }

  const bool modeExplicit =
      parameterInputLevel(P, "multiomeAtacPeakCallMode") > 0;
  const bool profileSet = !isUnsetToken(P.multiomeAtacPeakMex.macsProfile);
  profile->clear();
  if (profileSet) {
    *profile = trimCopy(P.multiomeAtacPeakMex.macsProfile);
    if (*mode == star::multiome::MultiomeAtacPeakCallMode::FRAG) {
      if (modeExplicit) {
        log << "ERROR: --multiomeAtacPeakMacsProfile requires "
               "--multiomeAtacPeakCallMode macs-bed, not frag\n";
        return false;
      }
      *mode = star::multiome::MultiomeAtacPeakCallMode::MACS_BED;
    }
  } else if (*mode == star::multiome::MultiomeAtacPeakCallMode::MACS_BED) {
    *profile = "signac-atac";
  }
  return true;
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

std::string outPrefixPath(const Parameters &P, const std::string &relative) {
  std::string base = trimCopy(P.outFileNamePrefix);
  if (base.empty()) {
    base = "./";
  }
  if (!base.empty() && base[base.size() - 1] != '/') {
    base += "/";
  }
  return base + relative;
}

std::string chooseOutputPath(const std::string &explicitPath,
                             const std::string &fallbackPath,
                             const Parameters &P,
                             const std::string &defaultRelative) {
  if (!isUnsetToken(explicitPath)) {
    return trimCopy(explicitPath);
  }
  if (!isUnsetToken(fallbackPath)) {
    return trimCopy(fallbackPath);
  }
  return outPrefixPath(P, defaultRelative);
}

std::string chooseOutputPath(const std::string &explicitPath,
                             const Parameters &P,
                             const std::string &defaultRelative) {
  if (!isUnsetToken(explicitPath)) {
    return trimCopy(explicitPath);
  }
  return outPrefixPath(P, defaultRelative);
}

std::string deriveSecondaryOutputPath(const std::string &primaryPath,
                                      const std::string &suffix) {
  const std::string primary = trimCopy(primaryPath);
  if (primary == "-" || primary == "/dev/stdout" || primary == "/dev/stderr") {
    return "chromap_output" + suffix + ".sam";
  }

  const size_t slashPos = primary.find_last_of("/\\");
  const size_t searchStart = slashPos == std::string::npos ? 0 : slashPos + 1;
  const std::string lower = lowerCopy(primary);

  if (lower.size() > 7 &&
      lower.compare(lower.size() - 7, 7, ".sam.gz") == 0) {
    return primary.substr(0, primary.size() - 7) + suffix + ".sam.gz";
  }
  if (lower.size() > 4 &&
      lower.compare(lower.size() - 4, 4, ".bam") == 0) {
    return primary.substr(0, primary.size() - 4) + suffix + ".bam";
  }
  if (lower.size() > 5 &&
      lower.compare(lower.size() - 5, 5, ".cram") == 0) {
    return primary.substr(0, primary.size() - 5) + suffix + ".cram";
  }
  if (lower.size() > 4 &&
      lower.compare(lower.size() - 4, 4, ".sam") == 0) {
    return primary.substr(0, primary.size() - 4) + suffix + ".sam";
  }

  const size_t dotPos = primary.rfind('.');
  if (dotPos != std::string::npos && dotPos > searchStart) {
    return primary.substr(0, dotPos) + suffix + primary.substr(dotPos);
  }
  return primary + suffix + ".sam";
}

bool validateAndBuildConfig(Parameters &P,
                            bool batchModeActive,
                            star::multiome::ChromapAtacConfig *cfg) {
  bool inlinePeakMex = false;
  if (!parseYesNo(P.multiomeAtacPeakMex.inlineMode, &inlinePeakMex)) {
    P.inOut->logMain
        << "ERROR: --multiomeAtacPeakMexInline must be yes or no (got \""
        << P.multiomeAtacPeakMex.inlineMode << "\")\n";
    return false;
  }
  bool peakTranslateFromFirst = false;
  if (!parseYesNo(P.multiomeAtacPeakMex.barcodeTranslateFromFirst,
                  &peakTranslateFromFirst)) {
    P.inOut->logMain
        << "ERROR: --multiomeAtacPeakBarcodeTranslateFromFirst must be yes or no (got \""
        << P.multiomeAtacPeakMex.barcodeTranslateFromFirst << "\")\n";
    return false;
  }
  (void)peakTranslateFromFirst;
  if (P.multiomeAtacPeakMex.threads < 0) {
    P.inOut->logMain
        << "ERROR: --multiomeAtacPeakThreads must be >= 0 (0 inherits --runThreadN)\n";
    return false;
  }
  if (inlinePeakMex) {
    star::multiome::MultiomeAtacPeakCallMode peakCallMode;
    std::string macsProfile;
    if (!resolveInlinePeakCallMode(P, &peakCallMode, &macsProfile,
                                   P.inOut->logMain)) {
      return false;
    }
  }
  if (inlinePeakMex && P.chromapAtac.enabled == 0) {
    P.inOut->logMain
        << "ERROR: --multiomeAtacPeakMexInline yes requires --chromapAtacEnable 1\n";
    return false;
  }
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

  const std::string inputFormatRaw = trimCopy(P.chromapAtac.inputFormat);
  const std::string inputFormat = lowerCopy(inputFormatRaw);
  star::multiome::ChromapInputFormat chromapInputFormat;
  if (inputFormat == "fastq") {
    chromapInputFormat = star::multiome::ChromapInputFormat::FASTQ;
  } else if (inputFormat == "cbq" || inputFormat == "binseq") {
    chromapInputFormat = star::multiome::ChromapInputFormat::CBQ;
  } else {
    P.inOut->logMain
        << "ERROR: --chromapAtacInputFormat must be fastq or cbq (got \""
        << inputFormatRaw << "\")\n";
    return false;
  }
  const bool cbqInput =
      chromapInputFormat == star::multiome::ChromapInputFormat::CBQ;

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
  if ((P.chromapAtac.emitNoYBam != 0 || P.chromapAtac.emitYBam != 0) &&
      outputFormat != star::multiome::ChromapOutputFormat::SAM &&
      outputFormat != star::multiome::ChromapOutputFormat::BAM &&
      outputFormat != star::multiome::ChromapOutputFormat::CRAM) {
    P.inOut->logMain
        << "ERROR: chromapAtacEmitNoYBam/chromapAtacEmitYBam require "
           "--chromapAtacOutputFormat SAM, BAM, or CRAM\n";
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
  if (inlinePeakMex && isUnsetToken(P.chromapAtac.secondaryFragments)) {
    P.inOut->logMain
        << "ERROR: --multiomeAtacPeakMexInline yes requires "
           "--chromapAtacSecondaryFragments to point at the AEV1 binary sidecar\n";
    return false;
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
  if (P.chromapAtac.callMacs3FragPeaks != 0) {
    if (isUnsetToken(P.chromapAtac.macs3FragPeaksOutput)) {
      P.inOut->logMain
          << "ERROR: --chromapAtacMacs3FragPeaksOutput is required when "
             "--chromapAtacCallMacs3FragPeaks 1\n";
      return false;
    }
    if (isUnsetToken(P.chromapAtac.macs3FragSummitsOutput)) {
      P.inOut->logMain
          << "ERROR: --chromapAtacMacs3FragSummitsOutput is required when "
             "--chromapAtacCallMacs3FragPeaks 1\n";
      return false;
    }
  }

  const int macs3PvalueInputLevel =
      parameterInputLevel(P, "chromapAtacMacs3FragPvalue");
  const int macs3QvalueInputLevel =
      parameterInputLevel(P, "chromapAtacMacs3FragQvalue");
  const bool macs3QvalueExplicit = macs3QvalueInputLevel > 0;
  const bool macs3UseQvalue =
      macs3QvalueExplicit && P.chromapAtac.macs3FragQvalue > 0.0;
  const bool needsMacs3Threshold =
      P.chromapAtac.callMacs3FragPeaks != 0 || inlinePeakMex;
  if (needsMacs3Threshold && macs3UseQvalue && macs3PvalueInputLevel > 0) {
    P.inOut->logMain
        << "ERROR: --chromapAtacMacs3FragPvalue and "
           "--chromapAtacMacs3FragQvalue are mutually exclusive\n";
    return false;
  }
  if (needsMacs3Threshold && macs3QvalueExplicit) {
    if (!(P.chromapAtac.macs3FragQvalue >= 0.0 &&
          P.chromapAtac.macs3FragQvalue <= 1.0)) {
      P.inOut->logMain
          << "ERROR: --chromapAtacMacs3FragQvalue must be 0 or in (0, 1]\n";
      return false;
    }
  }
  if (needsMacs3Threshold && !macs3UseQvalue &&
      !(P.chromapAtac.macs3FragPvalue > 0.0 &&
        P.chromapAtac.macs3FragPvalue <= 1.0)) {
    P.inOut->logMain
        << "ERROR: --chromapAtacMacs3FragPvalue must be in (0, 1]\n";
    return false;
  }

  std::vector<std::string> r1;
  std::vector<std::string> r2;
  std::vector<std::string> bc;
  std::vector<std::string> readPairCbq;
  std::vector<std::string> barcodeCbq;
  if (cbqInput) {
    if (!isUnsetToken(P.chromapAtac.read1Csv) ||
        !isUnsetToken(P.chromapAtac.read2Csv) ||
        !isUnsetToken(P.chromapAtac.barcodeCsv)) {
      P.inOut->logMain
          << "ERROR: --chromapAtacInputFormat cbq cannot be mixed with "
             "--chromapAtacRead1, --chromapAtacRead2, or --chromapAtacBarcode\n";
      return false;
    }
    if (isUnsetToken(P.chromapAtac.readPairCbqCsv)) {
      P.inOut->logMain
          << "ERROR: --chromapAtacReadPairCbq is required when "
             "--chromapAtacInputFormat cbq\n";
      return false;
    }
    if (isUnsetToken(P.chromapAtac.barcodeCbqCsv)) {
      P.inOut->logMain
          << "ERROR: --chromapAtacBarcodeCbq is required when "
             "--chromapAtacInputFormat cbq\n";
      return false;
    }
    readPairCbq = splitCsvPaths(P.chromapAtac.readPairCbqCsv);
    barcodeCbq = splitCsvPaths(P.chromapAtac.barcodeCbqCsv);
    if (readPairCbq.empty() || barcodeCbq.empty() ||
        hasUnsetPath(readPairCbq) || hasUnsetPath(barcodeCbq)) {
      P.inOut->logMain
          << "ERROR: --chromapAtacReadPairCbq and --chromapAtacBarcodeCbq must "
             "list at least one non-empty CBQ path each\n";
      return false;
    }
    if (readPairCbq.size() != barcodeCbq.size()) {
      P.inOut->logMain
          << "ERROR: --chromapAtacReadPairCbq and --chromapAtacBarcodeCbq "
             "must list the same number of CBQ paths\n";
      return false;
    }
  } else {
    if (!isUnsetToken(P.chromapAtac.readPairCbqCsv) ||
        !isUnsetToken(P.chromapAtac.barcodeCbqCsv)) {
      P.inOut->logMain
          << "ERROR: --chromapAtacReadPairCbq and --chromapAtacBarcodeCbq "
             "require --chromapAtacInputFormat cbq\n";
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
    r1 = splitCsvPaths(P.chromapAtac.read1Csv);
    r2 = splitCsvPaths(P.chromapAtac.read2Csv);
    bc = splitCsvPaths(P.chromapAtac.barcodeCsv);
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
  }

  if (cfg == nullptr) {
    return true;
  }

  cfg->reference_fasta = trimCopy(P.chromapAtac.referenceFasta);
  cfg->chromap_index = trimCopy(P.chromapAtac.chromapIndex);
  cfg->input_format = chromapInputFormat;
  if (cbqInput) {
    cfg->read_pair_cbqs = readPairCbq;
    cfg->barcode_cbqs = barcodeCbq;
  } else {
    cfg->read1_fastqs = r1;
    cfg->read2_fastqs = r2;
    cfg->barcode_fastqs = bc;
  }
  if (!isUnsetToken(P.chromapAtac.readFormat)) {
    cfg->read_format = trimCopy(P.chromapAtac.readFormat);
  }
  cfg->barcode_whitelist = trimCopy(P.chromapAtac.barcodeWhitelist);
  if (!isUnsetToken(P.chromapAtac.barcodeTranslate)) {
    cfg->barcode_translate_table = trimCopy(P.chromapAtac.barcodeTranslate);
  }
  cfg->barcode_translate_from_first_column =
      P.chromapAtac.barcodeTranslateFromFirst != 0;
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
  cfg->emit_no_y_bam = P.chromapAtac.emitNoYBam != 0;
  cfg->emit_y_bam = P.chromapAtac.emitYBam != 0;
  if (cfg->emit_no_y_bam) {
    cfg->no_y_output_path =
        isUnsetToken(P.chromapAtac.noYOutput)
            ? deriveSecondaryOutputPath(cfg->output_path, ".noY")
            : trimCopy(P.chromapAtac.noYOutput);
  }
  if (cfg->emit_y_bam) {
    cfg->y_output_path =
        isUnsetToken(P.chromapAtac.YOutput)
            ? deriveSecondaryOutputPath(cfg->output_path, ".Y")
            : trimCopy(P.chromapAtac.YOutput);
  }
  cfg->low_memory_mode = P.chromapAtac.lowMem != 0;
  cfg->low_memory_ram_limit = P.chromapAtac.lowMemRam;
  cfg->call_macs3_frag_peaks = P.chromapAtac.callMacs3FragPeaks != 0;
  if (!isUnsetToken(P.chromapAtac.macs3FragPeaksOutput)) {
    cfg->macs3_frag_peaks_output = trimCopy(P.chromapAtac.macs3FragPeaksOutput);
  }
  if (!isUnsetToken(P.chromapAtac.macs3FragSummitsOutput)) {
    cfg->macs3_frag_summits_output = trimCopy(P.chromapAtac.macs3FragSummitsOutput);
  }
  if (!isUnsetToken(P.chromapAtac.macs3FragKeepIntermediates)) {
    cfg->macs3_frag_keep_intermediates_dir =
        trimCopy(P.chromapAtac.macs3FragKeepIntermediates);
  }
  if (!isUnsetToken(P.chromapAtac.evidenceFromPeaksOutput)) {
    cfg->atac_evidence_from_peaks_output =
        trimCopy(P.chromapAtac.evidenceFromPeaksOutput);
  }
  cfg->macs3_frag_pvalue = P.chromapAtac.macs3FragPvalue;
  if (macs3UseQvalue) {
    cfg->macs3_frag_threshold_mode =
        star::multiome::ChromapMacs3FragThresholdMode::Q_VALUE;
    cfg->macs3_frag_qvalue = P.chromapAtac.macs3FragQvalue;
  } else {
    cfg->macs3_frag_threshold_mode =
        star::multiome::ChromapMacs3FragThresholdMode::P_VALUE;
  }
  cfg->macs3_frag_min_length = P.chromapAtac.macs3FragMinLength;
  cfg->macs3_frag_max_gap = P.chromapAtac.macs3FragMaxGap;
  cfg->macs3_frag_uint8_counts = P.chromapAtac.macs3FragUint8Counts != 0;
  cfg->macs3_frag_peaks_source =
      star::multiome::ChromapMacs3FragPeaksSource::MEMORY;
  cfg->macs3_frag_low_mem = P.chromapAtac.macs3FragLowMem != 0;

  // Wire the ATAC permit shims into the contract whenever STAR's dynamic
  // thread interface is enabled, mirroring how PfMultiProcess sets
  // assignOpts.enableStarDynamicPermitHooks = (P.dynamicThreadInterface == 1).
  // chromap's MapPairedEndReads then issues acquire/release per per-thread
  // mini-batch (~64 read pairs).
  if (P.dynamicThreadInterface == 1) {
    cfg->permit_hooks.acquire = chromapStarDynamicPermitAcquire;
    cfg->permit_hooks.release = chromapStarDynamicPermitRelease;
    cfg->permit_hooks.hook_ctx = &kAtacPermitHookContext;
  } else {
    cfg->permit_hooks.acquire = nullptr;
    cfg->permit_hooks.release = nullptr;
    cfg->permit_hooks.hook_ctx = nullptr;
  }
  return true;
}

// In-process per-barcode ATAC evidence is produced inside the
// libchromap_contract from the binary secondary-fragments sidecar. The
// orchestration only routes --chromapAtacEvidenceFromPeaksOutput into the
// contract config; this compatibility stub remains a no-op at old call sites.
bool runEvidenceFromPeaksIfEnabled(Parameters & /*P*/) { return true; }

bool runInlinePeakMexIfEnabled(Parameters &P) {
  bool inlinePeakMex = false;
  if (!parseYesNo(P.multiomeAtacPeakMex.inlineMode, &inlinePeakMex)) {
    P.inOut->logMain
        << "ERROR: --multiomeAtacPeakMexInline must be yes or no (got \""
        << P.multiomeAtacPeakMex.inlineMode << "\")\n";
    return false;
  }
  if (!inlinePeakMex) {
    return true;
  }

  bool translateFromFirst = true;
  if (!parseYesNo(P.multiomeAtacPeakMex.barcodeTranslateFromFirst,
                  &translateFromFirst)) {
    P.inOut->logMain
        << "ERROR: --multiomeAtacPeakBarcodeTranslateFromFirst must be yes or no (got \""
        << P.multiomeAtacPeakMex.barcodeTranslateFromFirst << "\")\n";
    return false;
  }

  star::multiome::MultiomeAtacPeakMexArgs args;
  args.sidecar = trimCopy(P.chromapAtac.secondaryFragments);
  if (!isUnsetToken(P.multiomeAtacPeakMex.barcodeTranslate)) {
    args.barcode_translate = trimCopy(P.multiomeAtacPeakMex.barcodeTranslate);
  } else if (!isUnsetToken(P.chromapAtac.barcodeTranslate)) {
    args.barcode_translate = trimCopy(P.chromapAtac.barcodeTranslate);
  }
  args.barcode_translate_from_first = translateFromFirst;
  args.peaks = chooseOutputPath(P.multiomeAtacPeakMex.narrowPeak,
                                P.chromapAtac.macs3FragPeaksOutput,
                                P, "atac/atac_peaks.narrowPeak");
  args.summits_out = chooseOutputPath(P.multiomeAtacPeakMex.summits,
                                      P.chromapAtac.macs3FragSummitsOutput,
                                      P, "atac/atac_summits.bed");
  args.out_dir = chooseOutputPath(P.multiomeAtacPeakMex.mexOutDir,
                                  P, "atac/peak_mex");
  args.metrics_tsv = chooseOutputPath(P.multiomeAtacPeakMex.metricsTsv,
                                      P, "atac/atac_metrics.tsv");
  if (!isUnsetToken(P.chromapAtac.tempDir)) {
    args.temp_dir = trimCopy(P.chromapAtac.tempDir);
  }
  if (!isUnsetToken(P.chromapAtac.macs3FragKeepIntermediates)) {
    args.keep_intermediates_dir = trimCopy(P.chromapAtac.macs3FragKeepIntermediates);
  }
  args.threads = (P.multiomeAtacPeakMex.threads > 0)
                     ? P.multiomeAtacPeakMex.threads
                     : P.runThreadN;
  args.call_peaks_from_sidecar = true;
  args.max_barcodes = P.multiomeAtacPeakMex.maxBarcodes;
  if (!resolveInlinePeakCallMode(P, &args.peak_call_mode, &args.macs_profile,
                                 P.inOut->logMain)) {
    return false;
  }
  const int macs3PvalueInputLevel =
      parameterInputLevel(P, "chromapAtacMacs3FragPvalue");
  const int macs3QvalueInputLevel =
      parameterInputLevel(P, "chromapAtacMacs3FragQvalue");
  if (macs3QvalueInputLevel > 0 &&
      P.chromapAtac.macs3FragQvalue > 0.0) {
    args.macs3_threshold_mode =
        star::multiome::MultiomeAtacPeakMexThresholdMode::Q_VALUE;
    args.macs3_qvalue = P.chromapAtac.macs3FragQvalue;
    args.macs3_threshold_explicit = true;
  } else {
    args.macs3_threshold_mode =
        star::multiome::MultiomeAtacPeakMexThresholdMode::P_VALUE;
    args.macs3_threshold_explicit = macs3PvalueInputLevel > 0;
  }
  args.macs3_pvalue = P.chromapAtac.macs3FragPvalue;
  args.macs3_min_length = P.chromapAtac.macs3FragMinLength;
  args.macs3_max_gap = P.chromapAtac.macs3FragMaxGap;
  args.macs3_uint8_counts = P.chromapAtac.macs3FragUint8Counts != 0;
  args.macs3_min_length_explicit =
      parameterInputLevel(P, "chromapAtacMacs3FragMinLength") > 0;
  args.macs3_max_gap_explicit =
      parameterInputLevel(P, "chromapAtacMacs3FragMaxGap") > 0;

  P.inOut->logMain << timeMonthDayTime()
                   << " ..... starting inline Multiome ATAC peak/MEX materialization\n"
                   << "       sidecar=" << args.sidecar << "\n"
                   << "       peak_call_mode="
                   << peakCallModeName(args.peak_call_mode) << "\n"
                   << "       macs_profile="
                   << (args.macs_profile.empty() ? "-" : args.macs_profile)
                   << "\n"
                   << "       peaks=" << args.peaks << "\n"
                   << "       summits=" << args.summits_out << "\n"
                   << "       peak_mex=" << args.out_dir << "\n"
                   << "       metrics=" << args.metrics_tsv << "\n"
                   << flush;

  const int rc = star::multiome::RunMultiomeAtacPeakMex(args);
  if (rc != 0) {
    P.inOut->logMain
        << "ERROR: inline Multiome ATAC peak/MEX materialization failed "
        << "with exit_code=" << rc << "\n";
    return false;
  }
  P.inOut->logMain << timeMonthDayTime()
                   << " ..... finished inline Multiome ATAC peak/MEX materialization\n"
                   << flush;
  return true;
}

}  // namespace

struct StarChromapAtacAsyncRun::Impl {
  explicit Impl(const star::multiome::ChromapAtacConfig &config,
                Parameters *p_for_sampler,
                int telemetryIntervalSec)
      : cfg(config),
        worker([this]() {
          try {
            result = star::multiome::runChromapAtac(cfg);
          } catch (...) {
            exception = std::current_exception();
          }
        }),
        samplerStop_(false),
        samplerStartTime_(std::chrono::steady_clock::now()) {
    // Periodic mapPermitSnapshot sampler: locked column schema (see
    // multiomic-atac-scrna plans/2026-04-27-atac-permits-controller-followups.md
    // v4). Header + per-emit flush so SIGTERM-kill keeps trajectory data.
    // Only runs when the dynamic-thread interface AND telemetry are on,
    // and the user-set interval is positive.
    if (p_for_sampler != nullptr && telemetryIntervalSec > 0 &&
        p_for_sampler->dynamicThreadInterface == 1 &&
        p_for_sampler->dynamicThreadTelemetry == 1) {
      P_for_sampler_ = p_for_sampler;
      samplerThread_ = std::thread([this, telemetryIntervalSec]() {
        runSampler(telemetryIntervalSec);
      });
    }
  }

  ~Impl() {
    stopSampler();
    if (worker.joinable()) {
      worker.join();
    }
  }

  void stopSampler() {
    if (samplerThread_.joinable()) {
      {
        std::lock_guard<std::mutex> lk(samplerMutex_);
        samplerStop_.store(true);
      }
      samplerCv_.notify_all();
      samplerThread_.join();
    }
  }

  void runSampler(int intervalSec) {
    if (P_for_sampler_ == nullptr) {
      return;
    }
    // Header (column names) + samples (values). Tab-separated, prefixed
    // with "[ATAC permit telemetry]" so they grep cleanly out of Log.out.
    // Header uses the same prefix + column names as values, so awk-style
    // parsing of "header row vs data rows" works naturally.
    {
      pthread_mutex_lock(&g_threadChunks.mutexLogMain);
      P_for_sampler_->inOut->logMain
          << "[ATAC permit telemetry]"
             "\ttimeSec\tconfigured\ttarget"
             "\tavailable\tinUse\twaiters"
             "\tmapAcquire\tmapWaitMaxMs\tmapWorkAvgMs"
             "\tatacAcquire\tatacWaitMaxMs\tatacWorkAvgMs"
             "\tretunes\n";
      P_for_sampler_->inOut->logMain.flush();
      pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
    }

    auto avgMs = [](uint64_t totalNs, uint64_t calls) -> double {
      return calls > 0 ? (static_cast<double>(totalNs) / 1.0e6 /
                          static_cast<double>(calls))
                       : 0.0;
    };

    // ATAC drain-time controller state (Step 6). Activated when
    // P.dynamicThreadAtacController == 1 AND at least one of MAP/ATAC
    // floor is set initially. The controller observes per-domain
    // workUnitsTotal deltas, smooths them via EWMA, and adjusts the
    // MAP/ATAC floors so both domains finish their mapping phase at
    // roughly the same rate. Mirrors the pf eta/chunked load-balancing
    // mechanism that has been the load-balancing core of pf/STAR since
    // the start.
    const bool controllerEnabled =
        (P_for_sampler_->dynamicThreadAtacController == 1) &&
        (P_for_sampler_->dynamicThreadInterface == 1) &&
        (P_for_sampler_->chromapAtac.enabled == 1) &&
        (P_for_sampler_->dynamicThreadMapFloor > 0 ||
         P_for_sampler_->dynamicThreadAtacFloor > 0);
    int curMapFloor = P_for_sampler_->dynamicThreadMapFloor;
    int curAtacFloor = P_for_sampler_->dynamicThreadAtacFloor;
    const int curFeatureFloor = P_for_sampler_->dynamicThreadFeatureFloor;
    uint64_t prevMapDone = 0;
    uint64_t prevAtacDone = 0;
    auto prevTick = std::chrono::steady_clock::now();
    double mapRateEwma = 0.0;
    double atacRateEwma = 0.0;
    constexpr double kEwmaAlpha = 0.30;
    constexpr double kRateGapThreshold = 0.20; // 20% imbalance triggers a retune
    bool primed = false; // first tick: just sample baseline, don't retune
    if (controllerEnabled) {
      pthread_mutex_lock(&g_threadChunks.mutexLogMain);
      P_for_sampler_->inOut->logMain
          << "[ATAC drain-time controller] enabled: mapFloor=" << curMapFloor
          << " atacFloor=" << curAtacFloor
          << " featureFloor=" << curFeatureFloor
          << " intervalSec=" << intervalSec
          << " ewmaAlpha=" << kEwmaAlpha
          << " gapThreshold=" << kRateGapThreshold << "\n";
      P_for_sampler_->inOut->logMain.flush();
      pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
    }

    while (true) {
      std::unique_lock<std::mutex> lk(samplerMutex_);
      samplerCv_.wait_for(lk, std::chrono::seconds(intervalSec),
                          [this]() { return samplerStop_.load(); });
      const bool stop = samplerStop_.load();
      lk.unlock();
      if (stop) {
        break;
      }

      const auto snap = g_threadChunks.mapPermitSnapshot();
      const auto now = std::chrono::steady_clock::now();
      const double timeSec = std::chrono::duration<double>(
                                 now - samplerStartTime_)
                                 .count();

      std::ostringstream line;
      line << "[ATAC permit telemetry]"
           << "\t" << std::fixed << std::setprecision(1) << timeSec
           << "\t" << snap.configuredPermits
           << "\t" << snap.targetPermits
           << "\t" << snap.availablePermits
           << "\t" << snap.inUsePermits
           << "\t" << snap.currentWaiters
           << "\t" << snap.mapDomain.acquireCalls
           << "\t" << std::fixed << std::setprecision(3)
           << (snap.mapDomain.waitNsMax / 1.0e6)
           << "\t" << avgMs(snap.mapDomain.workNsTotal,
                            snap.mapDomain.acquireCalls)
           << "\t" << snap.atacDomain.acquireCalls
           << "\t" << (snap.atacDomain.waitNsMax / 1.0e6)
           << "\t" << avgMs(snap.atacDomain.workNsTotal,
                            snap.atacDomain.acquireCalls)
           << "\t" << snap.retuneCalls
           << "\n";

      pthread_mutex_lock(&g_threadChunks.mutexLogMain);
      P_for_sampler_->inOut->logMain << line.str();
      P_for_sampler_->inOut->logMain.flush();
      pthread_mutex_unlock(&g_threadChunks.mutexLogMain);

      // ATAC drain-time controller: rate-balance retune of MAP/ATAC floors.
      if (controllerEnabled) {
        const uint64_t mapDone = snap.mapDomain.workUnitsTotal;
        const uint64_t atacDone = snap.atacDomain.workUnitsTotal;
        const double dt = std::chrono::duration<double>(now - prevTick).count();
        prevTick = now;

        // Raw deltas this tick; used both for the EWMA update and for the
        // per-tick idle guards below. Computing these before the EWMA
        // smoothing means the smoothed rate can stay nonzero from prior
        // history while we still see "this domain did literally no work
        // this tick" — and refuse to grow its floor in that case.
        const uint64_t mapDelta =
            (mapDone >= prevMapDone) ? (mapDone - prevMapDone) : 0;
        const uint64_t atacDelta =
            (atacDone >= prevAtacDone) ? (atacDone - prevAtacDone) : 0;
        const double mapInst = (dt > 0)
            ? static_cast<double>(mapDelta) / dt : 0.0;
        const double atacInst = (dt > 0)
            ? static_cast<double>(atacDelta) / dt : 0.0;
        prevMapDone = mapDone;
        prevAtacDone = atacDone;

        if (!primed) {
          // Seed the EWMA from the first interval rather than skipping
          // entirely — at high ratios we already know the right direction
          // and a 10s primer-skip costs us a full step of correction.
          mapRateEwma = mapInst;
          atacRateEwma = atacInst;
          primed = true;
          // Fall through and possibly retune on this same tick when the
          // imbalance is already large.
        } else {
          mapRateEwma = kEwmaAlpha * mapInst + (1.0 - kEwmaAlpha) * mapRateEwma;
          atacRateEwma = kEwmaAlpha * atacInst + (1.0 - kEwmaAlpha) * atacRateEwma;
        }

        // Skip retune when either domain is idle (no work this tick AND
        // no smoothed history) or the pool is barely contended (waiters=0
        // means no domain is being starved).
        if (mapRateEwma <= 0.0 || atacRateEwma <= 0.0 ||
            snap.currentWaiters == 0) {
          continue;
        }

        // Imbalance: ratio of fast/slow rates. Step size scales with the
        // log of the imbalance so a 100x ratio gets a chunky correction
        // instead of +1.
        const double ratio = (atacRateEwma > mapRateEwma)
            ? (atacRateEwma / mapRateEwma)
            : (mapRateEwma / atacRateEwma);
        if (ratio < (1.0 + kRateGapThreshold)) {
          continue; // balanced enough
        }

        // Step grows with imbalance: log2(ratio) clamped to [1, headroom].
        // ratio=2 -> step=1, ratio=4 -> 2, ratio=16 -> 4, ratio=100 -> 6+,
        // ratio=1024 -> 10. The clamp against headroom is applied below.
        const int chunkyStep = std::max<int>(
            1, static_cast<int>(std::ceil(std::log2(std::max(2.0, ratio)))));

        // Cap: let the slow domain reserve nearly the whole pool, leaving
        // only a one-permit floor for the fast domain so it can still make
        // forward progress (and so the FIFO grant logic doesn't deadlock
        // on a fully-claimed pool). Was poolHalf; now poolFull-1.
        const int slowDomainCap = std::max(1, snap.configuredPermits - 1);
        int newMap = curMapFloor;
        int newAtac = curAtacFloor;
        const char *direction = nullptr;

        if (mapRateEwma > atacRateEwma) {
          // MAP is draining faster → ATAC is being starved → grow ATAC floor.
          // Raw-delta idle guard: if ATAC didn't actually advance this tick,
          // the EWMA may still be nonzero from history but ATAC is already
          // done (or stuck on something other than permits). Don't grow its
          // floor — that would just leave permits idle.
          if (atacDelta == 0) {
            // skip atac++ this tick
          } else if (curAtacFloor < slowDomainCap) {
            const int step = std::min(chunkyStep, slowDomainCap - curAtacFloor);
            newAtac = curAtacFloor + step;
            // Drop MAP floor symmetrically but keep at least 0 (fast
            // domain doesn't need a floor to make progress; its threads
            // grab permits via FIFO when slow domain isn't using them).
            newMap = std::max(0, curMapFloor - step);
            direction = "atac++";
          }
        } else {
          // ATAC is draining faster → MAP is being starved (rare) →
          // grow MAP floor. Same raw-delta idle guard for MAP.
          if (mapDelta == 0) {
            // skip map++ this tick
          } else if (curMapFloor < slowDomainCap) {
            const int step = std::min(chunkyStep, slowDomainCap - curMapFloor);
            newMap = curMapFloor + step;
            newAtac = std::max(0, curAtacFloor - step);
            direction = "map++";
          }
        }

        if (direction != nullptr &&
            (newMap != curMapFloor || newAtac != curAtacFloor)) {
          std::vector<int> floors(3, 0);
          floors[0] = newMap;          // MAP
          floors[1] = curFeatureFloor; // FEATURE (unchanged)
          floors[2] = newAtac;         // ATAC
          g_threadChunks.mapPermitConfigureDomainFloors(floors);
          curMapFloor = newMap;
          curAtacFloor = newAtac;

          std::ostringstream ctrLine;
          ctrLine << "[ATAC drain-time controller] retune"
                  << "\ttimeSec=" << std::fixed << std::setprecision(1) << timeSec
                  << "\tmapRate=" << std::fixed << std::setprecision(1)
                  << mapRateEwma
                  << "\tatacRate=" << atacRateEwma
                  << "\tratio=" << std::setprecision(3) << ratio
                  << "\tdirection=" << direction
                  << "\tnewMapFloor=" << newMap
                  << "\tnewAtacFloor=" << newAtac
                  << "\n";
          pthread_mutex_lock(&g_threadChunks.mutexLogMain);
          P_for_sampler_->inOut->logMain << ctrLine.str();
          P_for_sampler_->inOut->logMain.flush();
          pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
        }
      }
    }
  }

  star::multiome::ChromapAtacConfig cfg;
  star::multiome::ChromapAtacResult result;
  std::exception_ptr exception;
  std::thread worker;

  // Periodic permit-snapshot sampler. Started in the ctor when the
  // dynamic-thread interface and telemetry are both on; joined in
  // ~Impl or via stopSampler() called from runStarChromapAtacIfEnabled
  // before the result is consumed.
  Parameters *P_for_sampler_ = nullptr;
  std::atomic<bool> samplerStop_;
  std::mutex samplerMutex_;
  std::condition_variable samplerCv_;
  std::thread samplerThread_;
  std::chrono::steady_clock::time_point samplerStartTime_;
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
    run.impl_ = new StarChromapAtacAsyncRun::Impl(
        cfg, &P, P.dynamicThreadTelemetryIntervalSec);
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
    // Stop the periodic snapshot sampler before reading the result so the
    // end-of-chromap summary line lands cleanly after the last sample.
    run.impl_->stopSampler();

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

    // ATAC permit telemetry summary (only when dynamicThreadInterface=1).
    if (P.dynamicThreadInterface == 1) {
      const auto snap = g_threadChunks.mapPermitSnapshot();
      P.inOut->logMain
          << "ATAC permit telemetry: enabled="
          << (snap.telemetryEnabled ? "yes" : "no")
          << " acquireCalls=" << snap.atacDomain.acquireCalls
          << " waitNsTotal=" << snap.atacDomain.waitNsTotal
          << " waitNsMax=" << snap.atacDomain.waitNsMax
          << " workUnitsTotal=" << snap.atacDomain.workUnitsTotal
          << " workNsTotal=" << snap.atacDomain.workNsTotal
          << " workNsMax=" << snap.atacDomain.workNsMax << "\n"
          << flush;
      // Integration guard: acquireCalls are telemetry counters, so they are
      // intentionally zero when --dynamicThreadTelemetry=0. When telemetry is
      // enabled, a successful Chromap ATAC mapping must exercise the hook at
      // least once; otherwise the linked libchromap/contract combination is
      // non-diagnostic for permit-sharing benchmarks.
      if (snap.telemetryEnabled && snap.atacDomain.acquireCalls == 0) {
        P.inOut->logMain
            << "ERROR: ATAC permit integration check failed: acquireCalls=0 "
               "with --dynamicThreadInterface 1, --dynamicThreadTelemetry 1, "
               "and --chromapAtacEnable 1.\n"
               "SOLUTION: build against a hook-enabled Chromap-suite tree "
               "whose mapping_parameters.h defines permit_acquire_hook / "
               "permit_release_hook / permit_hook_ctx, and point both the "
               "libchromap_contract and STAR builds at it via "
               "CHROMAP_SUITE_DIR=/path/to/that/tree (and "
               "CHROMAP_DIR=/path/to/that/tree for the contract sub-Makefile).\n"
            << flush;
        return false;
      }
      if (!snap.telemetryEnabled) {
        P.inOut->logMain
            << "NOTICE: ATAC permit integration check skipped because "
               "--dynamicThreadTelemetry=0; acquire counters are disabled.\n"
            << flush;
      }
    }
    if (!runEvidenceFromPeaksIfEnabled(P)) {
      return false;
    }
    if (!runInlinePeakMexIfEnabled(P)) {
      return false;
    }
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
  if (!runEvidenceFromPeaksIfEnabled(P)) {
    return false;
  }
  if (!runInlinePeakMexIfEnabled(P)) {
    return false;
  }
  return true;
}
