#include "star_chromap_orchestration.h"

#if !WITH_CHROMAP
#error "star_chromap_orchestration.cpp must be compiled with WITH_CHROMAP=1"
#endif

#include "GlobalVariables.h"
#include "IncludeDefine.h"
#include "Parameters.h"
#include "ThreadControl.h"
#include "TimeFunctions.h"
#include "star_chromap_contract.h"

#include <atomic>
#include <cctype>
#include <chrono>
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

  const std::string peaksSourceRaw = trimCopy(P.chromapAtac.macs3FragPeaksSource);
  const std::string peaksSource = lowerCopy(peaksSourceRaw);
  star::multiome::ChromapMacs3FragPeaksSource macs3FragPeaksSource;
  if (peaksSource == "file") {
    macs3FragPeaksSource = star::multiome::ChromapMacs3FragPeaksSource::FILE;
  } else if (peaksSource == "memory") {
    macs3FragPeaksSource = star::multiome::ChromapMacs3FragPeaksSource::MEMORY;
  } else {
    P.inOut->logMain
        << "ERROR: --chromapAtacMacs3FragPeaksSource must be file or memory (got \""
        << peaksSourceRaw << "\")\n";
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
  cfg->macs3_frag_pvalue = P.chromapAtac.macs3FragPvalue;
  cfg->macs3_frag_min_length = P.chromapAtac.macs3FragMinLength;
  cfg->macs3_frag_max_gap = P.chromapAtac.macs3FragMaxGap;
  cfg->macs3_frag_uint8_counts = P.chromapAtac.macs3FragUint8Counts != 0;
  cfg->macs3_frag_peaks_source = macs3FragPeaksSource;
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

        if (!primed) {
          prevMapDone = mapDone;
          prevAtacDone = atacDone;
          primed = true;
          continue; // baseline only, no retune on first tick
        }

        const double mapInst = (dt > 0 && mapDone >= prevMapDone)
            ? static_cast<double>(mapDone - prevMapDone) / dt
            : 0.0;
        const double atacInst = (dt > 0 && atacDone >= prevAtacDone)
            ? static_cast<double>(atacDone - prevAtacDone) / dt
            : 0.0;
        prevMapDone = mapDone;
        prevAtacDone = atacDone;

        mapRateEwma = (mapRateEwma <= 0.0)
            ? mapInst
            : (kEwmaAlpha * mapInst + (1.0 - kEwmaAlpha) * mapRateEwma);
        atacRateEwma = (atacRateEwma <= 0.0)
            ? atacInst
            : (kEwmaAlpha * atacInst + (1.0 - kEwmaAlpha) * atacRateEwma);

        // Skip retune when either domain is idle (no work this tick AND
        // no smoothed history) or the pool is barely contended (waiters=0
        // means no domain is being starved).
        if (mapRateEwma <= 0.0 || atacRateEwma <= 0.0 ||
            snap.currentWaiters == 0) {
          continue;
        }

        // Imbalance: ratio of fast/slow rates. When imbalance > threshold,
        // shift one floor up by 1 and the other floor down by 1 to push
        // more permits toward the slow domain.
        const double ratio = (atacRateEwma > mapRateEwma)
            ? (atacRateEwma / mapRateEwma)
            : (mapRateEwma / atacRateEwma);
        if (ratio < (1.0 + kRateGapThreshold)) {
          continue; // balanced enough
        }

        const int poolHalf = std::max(1, snap.configuredPermits / 2);
        int newMap = curMapFloor;
        int newAtac = curAtacFloor;
        const char *direction = nullptr;

        if (mapRateEwma > atacRateEwma) {
          // MAP is draining faster → ATAC is being starved → grow ATAC floor
          if (curAtacFloor < poolHalf) {
            newAtac = curAtacFloor + 1;
            if (curMapFloor > 1) newMap = curMapFloor - 1;
            direction = "atac++";
          }
        } else {
          // ATAC is draining faster → MAP is being starved (rare) →
          // grow MAP floor
          if (curMapFloor < poolHalf) {
            newMap = curMapFloor + 1;
            if (curAtacFloor > 1) newAtac = curAtacFloor - 1;
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
          << "ATAC permit telemetry: acquireCalls=" << snap.atacDomain.acquireCalls
          << " waitNsTotal=" << snap.atacDomain.waitNsTotal
          << " waitNsMax=" << snap.atacDomain.waitNsMax
          << " workUnitsTotal=" << snap.atacDomain.workUnitsTotal
          << " workNsTotal=" << snap.atacDomain.workNsTotal
          << " workNsMax=" << snap.atacDomain.workNsMax << "\n"
          << flush;
      // Integration guard: when STAR has wired the ATAC permit hooks into
      // the contract (dynamicThreadInterface == 1) and Chromap actually ran
      // an ATAC mapping, the linked libchromap *must* invoke the hook per
      // mini-batch. If acquireCalls is still 0 after the run, the linked
      // libchromap.a does not honor the hooks (header/lib mismatch) and any
      // fairness benchmark over this binary is non-diagnostic. Fail loudly
      // rather than producing misleading numbers.
      if (snap.atacDomain.acquireCalls == 0) {
        P.inOut->logMain
            << "ERROR: ATAC permit integration check failed: acquireCalls=0 "
               "with --dynamicThreadInterface 1 and --chromapAtacEnable 1.\n"
               "SOLUTION: build against a hook-enabled Chromap-suite tree "
               "whose mapping_parameters.h defines permit_acquire_hook / "
               "permit_release_hook / permit_hook_ctx, and point both the "
               "libchromap_contract and STAR builds at it via "
               "CHROMAP_SUITE_DIR=/path/to/that/tree (and "
               "CHROMAP_DIR=/path/to/that/tree for the contract sub-Makefile).\n"
            << flush;
        return false;
      }
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
  return true;
}
