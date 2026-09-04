#include "mapThreadsSpawn.h"
#include "FlexPipeline.h"
#include "FlexHashScreen.h"
#include "SoloReadFeature.h"
#include "SoloReadBarcode.h"
#include "SoloFeature.h"
#include "SoloFeatureTypes.h"
#include "Transcriptome.h"
#include "ThreadControl.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"
#include "InlineCBCorrection.h"
#include "SaturationPermitController.h"
#include "TimeFunctions.h"
#include "streamFuns.h"
#include "systemFunctions.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <zlib.h>

namespace {
std::string serializeRetuneTrace(const std::vector<int> &traceTargets) {
    if (traceTargets.empty()) {
        return "-";
    }

    std::ostringstream traceStream;
    for (size_t i = 0; i < traceTargets.size(); ++i) {
        if (i > 0) {
            traceStream << "|";
        }
        traceStream << traceTargets[i];
    }
    return traceStream.str();
}

std::string lowerCopy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

bool unsetToken(const std::string &value) {
    std::string v = lowerCopy(value);
    return v.empty() || v == "-" || v == "none" || v == "no";
}

bool readableFile(const std::string &path) {
    if (path.empty() || path == "-") {
        return false;
    }
    std::ifstream in(path.c_str(), std::ios::binary);
    return in.good();
}

[[noreturn]] void fatalBgzfRangeMode(Parameters &P, const std::string &reason) {
    std::ostringstream errOut;
    errOut << "EXITING because of fatal input ERROR: --readFilesBgzfMode range could not be activated.\n"
           << reason << "\n"
           << "SOLUTION: use paired BGZF FASTQ in a supported fully-fused Flex run, or set --readFilesBgzfMode auto/off.\n";
    exitWithError(errOut.str(), std::cerr, P.inOut->logMain,
                  EXIT_CODE_PARAMETER, P);
}

[[noreturn]] void fatalBgzfInput(Parameters &P, const std::string &reason) {
    std::ostringstream errOut;
    errOut << "EXITING because of fatal BGZF input ERROR:\n"
           << reason << "\n";
    exitWithError(errOut.str(), std::cerr, P.inOut->logMain,
                  EXIT_CODE_INPUT_FILES, P);
}

void finalizeFlexDynamicPermits(Parameters &P,
                                const FlexPipelineState &state) {
    if (!state.dynamicPermitsEnabled) {
        return;
    }

    // All fused/staged alignment workers have joined at the call sites. A
    // consumed BGZF stream likewise implies that every claimed inflate item
    // has completed, so both application domains are durably finished.
    g_threadChunks.mapPermitMarkDomainComplete(
        ThreadControl::PermitDomain::MAP);
    if (state.bgzfRangeActive) {
        g_threadChunks.mapPermitMarkDomainComplete(
            ThreadControl::PermitDomain::FEATURE);
    }

    if (P.dynamicThreadTelemetry == 1) {
        const ThreadControl::MapPermitSnapshot snapshot =
            g_threadChunks.mapPermitSnapshot();
        P.inOut->logMain
            << "Flex dynamic permit telemetry: configured="
            << snapshot.configuredPermits
            << ", map(acquires/blocked/maxInUse/workReads/waitMs/workMs)="
            << snapshot.mapDomain.acquireCalls << "/"
            << snapshot.mapDomain.blockedAcquireCalls << "/"
            << snapshot.mapDomain.maxInUse << "/"
            << snapshot.mapDomain.workUnitsTotal << "/"
            << snapshot.mapDomain.waitNsTotal / 1.0e6 << "/"
            << snapshot.mapDomain.workNsTotal / 1.0e6
            << ", bgzf(acquires/blocked/maxInUse/workBlocks/workBytes/waitMs/workMs)="
            << snapshot.featureDomain.acquireCalls << "/"
            << snapshot.featureDomain.blockedAcquireCalls << "/"
            << snapshot.featureDomain.maxInUse << "/"
            << snapshot.featureDomain.workUnitsTotal << "/"
            << snapshot.featureDomain.workBytesTotal << "/"
            << snapshot.featureDomain.waitNsTotal / 1.0e6 << "/"
            << snapshot.featureDomain.workNsTotal / 1.0e6
            << "\n" << std::flush;
    }
}
}

bool flexPipelineActivationGuard(Parameters &P, std::string *reason, bool logMessages) {
    const auto &ps = P.pSolo;
    auto reject = [&](const std::string &message) {
        if (reason != nullptr) {
            *reason = message;
        }
        if (logMessages) {
            P.inOut->logMain << "Flex pipeline: not active (" << message << ")\n" << std::flush;
        }
        return false;
    };
    if (ps.flexPipelineStr == "no") {
        if (reason != nullptr) {
            *reason = "disabled by --flexPipeline no";
        }
        if (logMessages) {
            P.inOut->logMain << "Flex pipeline: disabled by --flexPipeline no\n" << std::flush;
        }
        return false;
    }
    if (!ps.flexMode) {
        return reject("flexMode=false");
    }
    int nLanes = static_cast<int>(P.readFilesN);
    int nTriage = ps.flexPipelineNTriage;
    int nSolo = ps.flexPipelineNSolo;
    bool fullyFused = (nTriage == 0 && nSolo == 0);
    if (P.readFilesTypeN == 20 && P.outSAMtype.at(0) != "None") {
        return reject("CBQ/Binseq input uses the standard STAR CBQ adapter path");
    }
    if (!ps.hashScreenEnabled) {
        return reject("hash screen not enabled");
    }
    if (P.outSAMtype.at(0) != "None") {
        return reject("outSAMtype is not None");
    }
    if (P.readFilesTypeN == 20 && !fullyFused) {
        return reject("CBQ/Binseq input currently requires fully-fused mode: --flexPipelineNTriage 0 --flexPipelineNSolo 0");
    }
    int minThreads = fullyFused ? 1 : (nLanes + nTriage + nSolo + 1);
    if (P.runThreadN < minThreads) {
        std::ostringstream msg;
        msg << "runThreadN=" << P.runThreadN
            << " < minimum " << minThreads << " for " << nLanes << " lanes + "
            << nTriage << " triage + "
            << nSolo << " solo + 1 worker";
        return reject(msg.str());
    }
    if (reason != nullptr) {
        reason->clear();
    }
    return true;
}

bool flexNoGenomeCountOnlyActivationGuard(Parameters &P, std::string *reason) {
    auto reject = [&](const std::string &message) {
        if (reason != nullptr) {
            *reason = message;
        }
        return false;
    };

    if (std::getenv("STAR_DISABLE_FLEX_NO_GENOME") != nullptr) {
        return reject("disabled by STAR_DISABLE_FLEX_NO_GENOME");
    }

    std::string baseReason;
    if (!flexPipelineActivationGuard(P, &baseReason, false)) {
        if (baseReason == "hash screen not enabled") {
            return reject("prebuilt Flex hash cache is missing or unreadable");
        }
        return reject(baseReason.empty() ? "Flex pipeline guard rejected command" : baseReason);
    }

    const auto &ps = P.pSolo;
    if (P.runMode != "alignReads") {
        return reject("runMode is not alignReads");
    }
    if (ps.flexNoAlign == 0) {
        return reject("--flexNoAlign is not enabled");
    }
    if (ps.flexPipelineNTriage != 0 || ps.flexPipelineNSolo != 0) {
        return reject("requires fully fused mode: --flexPipelineNTriage 0 --flexPipelineNSolo 0");
    }
    if (P.readFilesN == 0) {
        return reject("no input read files");
    }
    if (P.runThreadN < 1) {
        return reject("runThreadN < 1");
    }
    if (!(P.readFilesTypeN == 1 || (P.readFilesTypeN == 20 && P.cbqInputActive))) {
        return reject("only FASTQ and CBQ fused inputs are supported");
    }
    if (P.readFilesTypeN == 1) {
        if (P.readFilesNames.size() < 2 ||
            P.readFilesNames[0].size() < static_cast<size_t>(P.readFilesN) ||
            P.readFilesNames[1].size() < static_cast<size_t>(P.readFilesN)) {
            return reject("FASTQ Flex input requires two complete read-file mate lists");
        }
    }
    if (P.readFilesTypeN == 20) {
        if (!P.cbqInputActive) {
            return reject("CBQ input is not active");
        }
        if (P.readFilesNames.empty() ||
            P.readFilesNames[0].size() < static_cast<size_t>(P.readFilesN)) {
            return reject("CBQ Flex input requires one complete read-file list");
        }
    }
    if (P.outSAMtype.empty() || P.outSAMtype[0] != "None" || P.outSAMbool ||
        P.outBAMcoord || P.outBAMunsorted || P.quant.trSAM.bamYes) {
        return reject("SAM/BAM output is enabled");
    }
    if (P.outSAMunmapped.yes || P.outReadsUnmapped != "None") {
        return reject("unmapped-read output is enabled");
    }
    if (P.outSJ.yes || P.outFilterBySJoutStage != 0 || P.outFilterType == "BySJout") {
        return reject("splice-junction output/filtering is enabled");
    }
    if (P.twoPass.yes || P.sjdbInsert.yes) {
        return reject("two-pass or run-time SJ insertion is enabled");
    }
    if (P.pGe.transform.outYes) {
        return reject("genome transform output is enabled");
    }
    if (P.pCh.segmentMin > 0) {
        return reject("chimeric detection/output is enabled");
    }
    if (P.outWigFlags.yes || P.emitNoYBAMyes || P.emitYReadNamesyes || P.emitYNoYFastqyes || P.emitYNoYCbqyes) {
        return reject("genome-backed auxiliary output is enabled");
    }
    if (P.quant.geCount.yes || P.quant.trSAM.yes || P.quant.transcriptVB.yes || P.quant.slam.yes) {
        return reject("genome/transcriptome quantification output is enabled");
    }
    if (P.quant.geneFull.yes || P.quant.geneFull_Ex50pAS.yes || P.quant.geneFull_ExonOverIntron.yes) {
        return reject("GeneFull-style Solo feature requires genome-backed annotation");
    }
    if (P.var.yes || P.wasp.yes || P.trimQcEnabled) {
        return reject("variant/WASP/trim-QC mode is enabled");
    }
    if (P.chromapAtac.enabled != 0 || !unsetToken(P.multiomeAtacPeakMex.inlineMode)) {
        return reject("Chromap/multiome ATAC output is enabled");
    }
    if (!unsetToken(P.pfMulti.pfMultiConfig) || !unsetToken(P.pfMulti.ocmMultiEnable)) {
        return reject("pf-multi/OCM post-processing is enabled");
    }
    if (P.batchMode || P.batchModeInt != 0 || P.quant.slam.batchMode) {
        return reject("batch mode is enabled");
    }
    if (P.runRestart.type != 0) {
        return reject("restart mode is enabled");
    }
    if (ps.nFeatures != 1 || !ps.featureYes[SoloFeatureTypes::Gene]) {
        return reject("requires exactly --soloFeatures Gene");
    }
    if (ps.type == ParametersSolo::SoloTypes::None || ps.type == ParametersSolo::SoloTypes::CB_samTagOut) {
        return reject("requires count-producing Solo CB/UMI mode");
    }
    if (ps.trackReadIdsForTags) {
        return reject("BAM CB/UB tag replay state is enabled");
    }
    if (ps.hashScreenFile.empty() || !readableFile(ps.hashScreenFile)) {
        return reject("prebuilt Flex hash cache is missing or unreadable");
    }
    if (ps.probeListPath.empty() || !readableFile(ps.probeListPath)) {
        return reject("Flex probe list is missing or unreadable");
    }

    if (reason != nullptr) {
        reason->clear();
    }
    return true;
}

static void mapThreadsSpawnFlexPipeline(Parameters &P, ReadAlignChunk** RAchunk) {
    // Kill legacy FIFO decompression children (zcat/gzip) — the pipeline
    // opens its own gzFile handles and does not use the FIFO streams.
    P.closeReadsFiles();

    const int nLanes = static_cast<int>(P.readFilesN);
    const int nTriage = P.pSolo.flexPipelineNTriage;
    const int nSolo = P.pSolo.flexPipelineNSolo;

    // Mode selection: nTriage=0 + nSolo=0 → fully fused (read+route+solo inline)
    //                 nTriage=0 + nSolo>0 → fused reader+router, separate solo consumers
    //                 nTriage>0            → separate reader, triage, solo threads
    const bool fullyFused = (nTriage == 0 && nSolo == 0);
    const bool fusedReaderRouter = (nTriage == 0 && nSolo > 0);
    const int actualNTriage = (nTriage == 0) ? 0 : nTriage;
    const int actualNSolo = fullyFused ? 0 : nSolo;
    const bool noAlign = (P.pSolo.flexNoAlign != 0);

    // Fully fused: every thread reads and hashes, and with alignment enabled
    // aligns queued misses itself whenever alignQ is full (enqueueForAlign in
    // FlexPipeline.cpp), so no thread is reserved as a consumer and
    // --runThreadN 1 is valid.
    const int nWorkers = fullyFused ? 0 : (P.runThreadN - nLanes - actualNTriage - actualNSolo);
    const int nFusedThreads = fullyFused ? P.runThreadN : 0;

    if (P.readFilesTypeN == 1 && !fullyFused &&
        P.readFilesBgzfMode == "range") {
        fatalBgzfRangeMode(P,
            "BGZF range readers currently require --flexPipelineNTriage 0 --flexPipelineNSolo 0");
    }

    P.inOut->logMain << "Flex pipeline: runThreadN=" << P.runThreadN
                     << ", nLanes=" << nLanes
                     << ", triage=" << actualNTriage
                     << (fullyFused ? " (fully fused + lane steal + role switch)" :
                         fusedReaderRouter ? " (fused reader+router)" : "")
                     << ", soloConsumers=" << actualNSolo
                     << ", fusedThreads=" << nFusedThreads
                     << ", dedicatedWorkers=" << nWorkers
                     << (noAlign ? ", noAlign=ON (alignment skipped)" : "")
                     << "\n" << std::flush;

    FlexPipelineState state;
    state.init(nLanes, actualNSolo, actualNTriage);
    state.dynamicPermitsEnabled =
        P.dynamicThreadInterface == 1 && g_threadChunks.mapPermitEnabled();

    // Store lane file paths for dynamic claiming (fully-fused mode)
    if (fullyFused) {
        state.laneFiles.resize(nLanes);
        for (int lane = 0; lane < nLanes; ++lane) {
            if (P.readFilesTypeN == 20 && P.cbqInputActive) {
                state.laneFiles[lane].r2path = P.readFilesNames[0][lane];
                state.laneFiles[lane].r1path.clear();
            } else {
                state.laneFiles[lane].r2path = P.readFilesNames[0][lane];
                state.laneFiles[lane].r1path = P.readFilesNames[1][lane];
            }
        }
        state.nFusedThreads = nFusedThreads;
        if (P.readFilesTypeN == 20 && P.cbqInputActive) {
            std::string cbqRangeReason;
            if (flexPrepareCbqRangeTasks(&state, P, nFusedThreads, &cbqRangeReason)) {
                P.inOut->logMain << "Flex CBQ range: active ("
                                 << cbqRangeReason << ")\n" << std::flush;
            } else {
                P.inOut->logMain << "Flex CBQ range: not active ("
                                 << cbqRangeReason
                                 << "); using whole-lane CBQ readers\n" << std::flush;
            }
        } else if (P.readFilesTypeN == 1 && P.readFilesBgzfMode != "off") {
            std::string bgzfRangeReason;
            bool fatalError = false;
            if (flexPrepareBgzfRangeTasks(&state, P, nFusedThreads,
                                          &bgzfRangeReason, &fatalError)) {
                P.inOut->logMain << "BGZF parallel range readers: active ("
                                 << bgzfRangeReason << ")\n";
            } else if (fatalError) {
                if (P.readFilesBgzfMode == "range") {
                    fatalBgzfRangeMode(P, bgzfRangeReason);
                }
                fatalBgzfInput(P, bgzfRangeReason);
            } else {
                P.inOut->logMain << "BGZF parallel range readers: not active ("
                                 << bgzfRangeReason << "); using zlib\n";
            }
            for (int lane = 0; lane < nLanes; ++lane) {
                P.inOut->logMain << "BGZF input lane " << lane << ": "
                                 << (state.laneUsesBgzfRange(lane) ? "range" : "zlib")
                                 << "\n";
            }
            P.inOut->logMain << std::flush;
        }
    }

    // Per-thread SoloReadFeature + Stats for fully-fused mode
    std::vector<SoloReadFeature *> fusedFeats(nFusedThreads);
    std::vector<Stats> fusedStats(nFusedThreads);

    if (fullyFused) {
        // Hash hits and residual-alignment misses both contribute Flex gene
        // observations in this mode. Accumulate their ambiguous barcodes in
        // the same striped store, so the alignment-enabled path avoids the
        // per-thread field-wise fold just like the no-genome count-only path.
        // The store remains active until Gene sumThreads drains it.
        const bool useSharedAmbig =
            P.pSolo.inlineHashMode && P.pSolo.flexMode;
        SoloReadFeature::sharedAmbigEnable(useSharedAmbig);
        P.inOut->logMain << "Flex shared ambiguous store: "
                         << (useSharedAmbig ? "active" : "inactive")
                         << " for fully fused hash and alignment records\n"
                         << std::flush;

        // All runThreadN threads are fused: they steal lanes, align queued
        // misses whenever alignQ fills, then switch to alignment
        std::vector<pthread_t> fusedThreads(nFusedThreads);
        std::vector<FlexLaneReaderArgs> fusedArgs(nFusedThreads);
        std::vector<std::unique_ptr<SoloReadBarcode>> fusedBars;
        fusedBars.reserve(nFusedThreads);
        for (int i = 0; i < nFusedThreads; ++i) {
            fusedBars.emplace_back(new SoloReadBarcode(P));
        }

        for (int i = 0; i < nFusedThreads; ++i) {
            fusedFeats[i] = new SoloReadFeature(SoloFeatureTypes::Gene, P, -(200 + i));
            fusedStats[i] = Stats();
            fusedArgs[i].state = &state;
            fusedArgs[i].P = &P;
            fusedArgs[i].gzR2 = nullptr;
            fusedArgs[i].gzR1 = nullptr;
            fusedArgs[i].laneId = -1;
            fusedArgs[i].readFeat = fusedFeats[i];
            fusedArgs[i].stats = &fusedStats[i];
            fusedArgs[i].RA = RAchunk[i]->RA;
            fusedArgs[i].threadId = i;
            fusedArgs[i].readBar = fusedBars[i].get();
            pthread_create(&fusedThreads[i], nullptr, flexLaneReaderFullThread, &fusedArgs[i]);
        }
        P.inOut->logMain << "  " << nFusedThreads << " fused threads started (lane-steal + role-switch)\n" << std::flush;

        // Stats reporter
        pthread_t reporterThread;
        FlexStatsReporterArgs reporterArgs;
        reporterArgs.state = &state;
        reporterArgs.P = &P;
        pthread_create(&reporterThread, nullptr, flexStatsReporterThread, &reporterArgs);

        // Join all fused threads (they self-transition from reader to aligner)
        for (int i = 0; i < nFusedThreads; ++i) pthread_join(fusedThreads[i], nullptr);
        P.inOut->logMain << "  All fused threads joined\n" << std::flush;

        finalizeFlexDynamicPermits(P, state);

        P.inOut->logMain << "  Fused producers aligned "
                         << state.counters.alignHelped.load(std::memory_order_relaxed)
                         << " queued reads while alignQ was full\n" << std::flush;

        if (state.inputFailed.load(std::memory_order_relaxed)) {
            fatalBgzfInput(P, state.inputError.empty()
                ? "a fused input reader failed" : state.inputError);
        }

        state.pipelineDone.store(true, std::memory_order_relaxed);
        pthread_join(reporterThread, nullptr);
        P.inOut->logMain << "  Stats reporter joined\n" << std::flush;

        // Merge stats from all fused threads
        for (int i = 0; i < nFusedThreads; ++i) {
            // Screened KEEP/DENY reads used fusedBars; alignment misses used
            // the ReadAlign-owned barcode object. Merge the former before
            // Solo constructs its global exact-CB prior.
            RAchunk[i]->RA->soloRead->readBar->addCounts(*fusedBars[i]);
            RAchunk[i]->RA->soloRead->readBar->addStats(*fusedBars[i]);
            g_statsAll.addStats(fusedStats[i]);
            if (!noAlign)
                g_statsAll.addStats(RAchunk[i]->RA->statsRA);
        }
        g_statsAll.readN = state.counters.readsTotal.load();

        // Merge Solo hashes
        if (P.pSolo.inlineHashMode && P.pSolo.featureYes[SoloFeatureTypes::Gene]) {
            SoloReadFeature *workerGeneFeat = RAchunk[0]->RA->soloRead->readFeat[P.pSolo.featureInd[SoloFeatureTypes::Gene]];
            for (int i = 0; i < nFusedThreads; ++i) {
                workerGeneFeat->mergeInlineHash(*fusedFeats[i]);
            }
        }

        for (int i = 0; i < nFusedThreads; ++i) delete fusedFeats[i];

        P.inOut->logMain << "Flex pipeline complete: total=" << state.counters.readsTotal.load()
                         << ", triageKeep=" << state.counters.triageKeep.load()
                         << ", triageDeny=" << state.counters.triageDeny.load()
                         << ", sampleReject=" << state.counters.triageSampleReject.load()
                         << ", triageMiss=" << state.counters.triageMiss.load() << "\n";
        for (int i = 0; i < nLanes; ++i) {
            P.inOut->logMain << "  Lane " << i << ": " << state.counters.perLaneReads[i] << " reads\n";
        }
        P.inOut->logMain << std::flush;
        return;
    }

    // --- Non-fused paths (legacy pipeline modes) ---

    // Per-lane SoloReadFeature + Stats (not used in non-fused)
    std::vector<SoloReadFeature *> laneFeats;
    std::vector<Stats> laneStats;

    std::vector<pthread_t> readerThreads(nLanes);
    std::vector<FlexLaneReaderArgs> readerArgs(nLanes);
    for (int lane = 0; lane < nLanes; ++lane) {
        const std::string &r2path = P.readFilesNames[0][lane];
        const std::string &r1path = P.readFilesNames[1][lane];
        gzFile gzR2 = gzopen(r2path.c_str(), "rb");
        gzFile gzR1 = gzopen(r1path.c_str(), "rb");
        if (!gzR2 || !gzR1) {
            std::ostringstream errOut;
            errOut << "EXITING: Flex pipeline cannot open lane " << lane
                   << " files: " << r2path << " / " << r1path << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, 1, P);
        }
        gzbuffer(gzR2, 1 << 20);
        gzbuffer(gzR1, 1 << 20);

        readerArgs[lane].state = &state;
        readerArgs[lane].P = &P;
        readerArgs[lane].gzR2 = gzR2;
        readerArgs[lane].gzR1 = gzR1;
        readerArgs[lane].laneId = lane;
        readerArgs[lane].readFeat = nullptr;
        readerArgs[lane].stats = nullptr;
        readerArgs[lane].RA = nullptr;
        readerArgs[lane].threadId = lane;
        readerArgs[lane].readBar = nullptr;

        if (fusedReaderRouter) {
            pthread_create(&readerThreads[lane], nullptr, flexLaneReaderRouterThread, &readerArgs[lane]);
            P.inOut->logMain << "  Flex reader+router lane " << lane << " started\n" << std::flush;
        } else {
            pthread_create(&readerThreads[lane], nullptr, flexLaneReaderThread, &readerArgs[lane]);
            P.inOut->logMain << "  Flex reader lane " << lane << " started\n" << std::flush;
        }
    }

    // --- Stage 2: Triage (skipped when any fused mode) ---
    std::vector<pthread_t> triageThreads(actualNTriage);
    std::vector<FlexTriageArgs> triageArgs(actualNTriage);
    for (int i = 0; i < actualNTriage; ++i) {
        triageArgs[i].state = &state;
        triageArgs[i].P = &P;
        pthread_create(&triageThreads[i], nullptr, flexTriageThread, &triageArgs[i]);
        P.inOut->logMain << "  Flex triage " << i << " started\n" << std::flush;
    }

    // --- Stage 3a: Solo consumers ---
    std::vector<pthread_t> soloThreads(actualNSolo);
    std::vector<FlexSoloConsumerArgs> soloArgs(actualNSolo);
    std::vector<SoloReadFeature *> soloFeats(actualNSolo);
    std::vector<Stats> soloStats(actualNSolo);
    std::vector<std::unique_ptr<SoloReadBarcode>> soloBars;
    soloBars.reserve(actualNSolo);
    for (int i = 0; i < actualNSolo; ++i) {
        soloBars.emplace_back(new SoloReadBarcode(P));
    }

    for (int i = 0; i < actualNSolo; ++i) {
        soloFeats[i] = new SoloReadFeature(SoloFeatureTypes::Gene, P, -(100 + i));
        soloArgs[i].state = &state;
        soloArgs[i].P = &P;
        soloArgs[i].consumerId = i;
        soloArgs[i].readFeat = soloFeats[i];
        soloArgs[i].stats = &soloStats[i];
        soloArgs[i].readBar = soloBars[i].get();
        pthread_create(&soloThreads[i], nullptr, flexSoloConsumerThread, &soloArgs[i]);
        P.inOut->logMain << "  Flex solo consumer " << i << " started\n" << std::flush;
    }

    // --- Stage 3b: Alignment workers ---
    std::vector<pthread_t> workerThreads(nWorkers);
    std::vector<FlexAlignWorkerArgs> workerArgs(nWorkers);
    for (int i = 0; i < nWorkers; ++i) {
        workerArgs[i].state = &state;
        workerArgs[i].RA = RAchunk[i]->RA;
        pthread_create(&workerThreads[i], nullptr, flexAlignWorkerThread, &workerArgs[i]);
        P.inOut->logMain << "  Flex align worker " << i << " started\n" << std::flush;
    }

    // --- Stats reporter (does not count against runThreadN) ---
    pthread_t reporterThread;
    FlexStatsReporterArgs reporterArgs;
    reporterArgs.state = &state;
    reporterArgs.P = &P;
    pthread_create(&reporterThread, nullptr, flexStatsReporterThread, &reporterArgs);
    P.inOut->logMain << "  Flex stats reporter started\n" << std::flush;

    // --- Join all threads ---
    for (int i = 0; i < nLanes; ++i) pthread_join(readerThreads[i], nullptr);
    P.inOut->logMain << "  All readers joined\n" << std::flush;

    for (int i = 0; i < actualNTriage; ++i) pthread_join(triageThreads[i], nullptr);
    if (actualNTriage > 0)
        P.inOut->logMain << "  All triage threads joined\n" << std::flush;

    for (int i = 0; i < actualNSolo; ++i) pthread_join(soloThreads[i], nullptr);
    if (actualNSolo > 0)
        P.inOut->logMain << "  All solo consumers joined\n" << std::flush;

    for (int i = 0; i < nWorkers; ++i) pthread_join(workerThreads[i], nullptr);
    P.inOut->logMain << "  All alignment workers joined\n" << std::flush;

    finalizeFlexDynamicPermits(P, state);

    // Stop and join stats reporter
    state.pipelineDone.store(true, std::memory_order_relaxed);
    pthread_join(reporterThread, nullptr);
    P.inOut->logMain << "  Stats reporter joined\n" << std::flush;

    // --- Merge stats ---
    for (int i = 0; i < actualNSolo; ++i) {
        // Alignment workers already retain their own CB evidence. Put the
        // screened KEEP/DENY evidence on one worker so Solo's ordinary
        // readBarSum reduction sees each consumer exactly once.
        RAchunk[0]->RA->soloRead->readBar->addCounts(*soloBars[i]);
        RAchunk[0]->RA->soloRead->readBar->addStats(*soloBars[i]);
        g_statsAll.addStats(soloStats[i]);
    }
    for (int i = 0; i < nWorkers; ++i) {
        g_statsAll.addStats(RAchunk[i]->RA->statsRA);
    }

    g_statsAll.readN = state.counters.readsTotal.load();

    // --- Merge Solo hashes into per-feature sum ---
    if (P.pSolo.inlineHashMode && P.pSolo.featureYes[SoloFeatureTypes::Gene]) {
        SoloReadFeature *workerGeneFeat = RAchunk[0]->RA->soloRead->readFeat[P.pSolo.featureInd[SoloFeatureTypes::Gene]];
        for (int i = 0; i < actualNSolo; ++i) {
            workerGeneFeat->mergeInlineHash(*soloFeats[i]);
        }
    }

    // Cleanup
    for (int i = 0; i < actualNSolo; ++i) delete soloFeats[i];

    // Log pipeline counters
    P.inOut->logMain << "Flex pipeline complete: total=" << state.counters.readsTotal.load()
                     << ", triageKeep=" << state.counters.triageKeep.load()
                     << ", triageDeny=" << state.counters.triageDeny.load()
                     << ", sampleReject=" << state.counters.triageSampleReject.load()
                     << ", triageMiss=" << state.counters.triageMiss.load() << "\n";
    for (int i = 0; i < nLanes; ++i) {
        P.inOut->logMain << "  Lane " << i << ": " << state.counters.perLaneReads[i] << " reads\n";
    }
    P.inOut->logMain << std::flush;
}

void runFlexNoGenomeCountOnly(Parameters &P) {
    P.closeReadsFiles();

    g_statsAll.resetN();
    time(&g_statsAll.timeStartMap);
    *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeStartMap) << " ..... started mapping\n"
                        << std::flush;
    g_statsAll.timeLastReport = g_statsAll.timeStartMap;
    g_statsAll.progressReportHeader(P.inOut->logProgress);

    P.inOut->logMain << "Flex count-only no-genome: hash cache " << P.pSolo.hashScreenFile << "\n"
                     << std::flush;

    if (P.pSolo.inlineCBCorrection && P.pSolo.cbWLyes && !P.pSolo.cbWLstr.empty()) {
        InlineCBCorrection::initializeWhitelist(P.pSolo);
        P.inOut->logMain << "[INLINE-CB-INIT] size=" << P.pSolo.cbWLstr.size()
                         << " exact_map=" << InlineCBCorrection::exactMapSize()
                         << " variant_map=" << InlineCBCorrection::variantMapSize()
                         << " variant_collisions=" << InlineCBCorrection::variantCollisionSize()
                         << " collision_max_fanout=" << InlineCBCorrection::variantCollisionMaxFanout()
                         << "\n";
    }

    const int nLanes = static_cast<int>(P.readFilesN);
    const int nFusedThreads = P.runThreadN;

    P.inOut->logMain << "Flex pipeline: runThreadN=" << P.runThreadN
                     << ", nLanes=" << nLanes
                     << ", triage=0 (fully fused + lane steal + role switch)"
                     << ", soloConsumers=0"
                     << ", fusedThreads=" << nFusedThreads
                     << ", dedicatedWorkers=0"
                     << ", noAlign=ON (alignment skipped)"
                     << ", noGenome=ON"
                     << "\n" << std::flush;

    FlexPipelineState state;
    state.init(nLanes, 0, 0);
    state.laneFiles.resize(nLanes);
    for (int lane = 0; lane < nLanes; ++lane) {
        if (P.readFilesTypeN == 20 && P.cbqInputActive) {
            state.laneFiles[lane].r2path = P.readFilesNames[0][lane];
            state.laneFiles[lane].r1path.clear();
        } else {
            state.laneFiles[lane].r2path = P.readFilesNames[0][lane];
            state.laneFiles[lane].r1path = P.readFilesNames[1][lane];
        }
    }
    state.nFusedThreads = nFusedThreads;
    if (P.readFilesTypeN == 20 && P.cbqInputActive) {
        std::string cbqRangeReason;
        if (flexPrepareCbqRangeTasks(&state, P, nFusedThreads, &cbqRangeReason)) {
            P.inOut->logMain << "Flex CBQ range: active ("
                             << cbqRangeReason << ")\n" << std::flush;
        } else {
            P.inOut->logMain << "Flex CBQ range: not active ("
                             << cbqRangeReason
                             << "); using whole-lane CBQ readers\n" << std::flush;
        }
    } else if (P.readFilesTypeN == 1 && P.readFilesBgzfMode != "off") {
        std::string bgzfRangeReason;
        bool fatalError = false;
        if (flexPrepareBgzfRangeTasks(&state, P, nFusedThreads,
                                      &bgzfRangeReason, &fatalError)) {
            P.inOut->logMain << "BGZF parallel range readers: active ("
                             << bgzfRangeReason << ")\n";
        } else if (fatalError) {
            if (P.readFilesBgzfMode == "range") {
                fatalBgzfRangeMode(P, bgzfRangeReason);
            }
            fatalBgzfInput(P, bgzfRangeReason);
        } else {
            P.inOut->logMain << "BGZF parallel range readers: not active ("
                             << bgzfRangeReason << "); using zlib\n";
        }
        for (int lane = 0; lane < nLanes; ++lane) {
            P.inOut->logMain << "BGZF input lane " << lane << ": "
                             << (state.laneUsesBgzfRange(lane) ? "range" : "zlib")
                             << "\n";
        }
        P.inOut->logMain << std::flush;
    }

    // Accumulate ambiguous barcodes into one striped store rather than per
    // thread. Each barcode is then combined in a single place during mapping,
    // instead of being reconciled field by field afterwards, which was the
    // largest serial block in the Flex tail. Only for this fused path; the
    // staged pipeline keeps per-thread maps and its existing merge.
    SoloReadFeature::sharedAmbigEnable(P.pSolo.inlineHashMode && P.pSolo.flexMode);
    P.inOut->logMain << "Flex shared ambiguous store: active for fully fused "
                        "no-genome hash records\n" << std::flush;

    std::vector<SoloReadFeature *> fusedFeats(nFusedThreads, nullptr);
    std::vector<Stats> fusedStats(nFusedThreads);
    std::vector<pthread_t> fusedThreads(nFusedThreads);
    std::vector<FlexLaneReaderArgs> fusedArgs(nFusedThreads);
    std::vector<std::unique_ptr<SoloReadBarcode>> fusedBars;
    fusedBars.reserve(nFusedThreads);
    for (int i = 0; i < nFusedThreads; ++i) {
        fusedBars.emplace_back(new SoloReadBarcode(P));
    }

    for (int i = 0; i < nFusedThreads; ++i) {
        fusedFeats[i] = new SoloReadFeature(SoloFeatureTypes::Gene, P, -(200 + i));
        fusedStats[i] = Stats();
        fusedArgs[i].state = &state;
        fusedArgs[i].P = &P;
        fusedArgs[i].gzR2 = nullptr;
        fusedArgs[i].gzR1 = nullptr;
        fusedArgs[i].laneId = -1;
        fusedArgs[i].readFeat = fusedFeats[i];
        fusedArgs[i].stats = &fusedStats[i];
        fusedArgs[i].RA = nullptr;
        fusedArgs[i].threadId = i;
        fusedArgs[i].readBar = fusedBars[i].get();
        int threadStatus = pthread_create(&fusedThreads[i], nullptr, flexLaneReaderFullThread, &fusedArgs[i]);
        if (threadStatus > 0) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR: pthread error while creating Flex no-genome thread # "
                   << i << ", error code: " << threadStatus;
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_RUNTIME, P);
        }
    }
    P.inOut->logMain << "  " << nFusedThreads << " fused no-genome threads started\n" << std::flush;

    pthread_t reporterThread;
    FlexStatsReporterArgs reporterArgs;
    reporterArgs.state = &state;
    reporterArgs.P = &P;
    int reporterStatus = pthread_create(&reporterThread, nullptr, flexStatsReporterThread, &reporterArgs);
    if (reporterStatus > 0) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: pthread error while creating Flex stats reporter, error code: "
               << reporterStatus;
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_RUNTIME, P);
    }

    for (int i = 0; i < nFusedThreads; ++i) {
        pthread_join(fusedThreads[i], nullptr);
    }
    P.inOut->logMain << "  All fused no-genome threads joined\n" << std::flush;

    if (state.inputFailed.load(std::memory_order_relaxed)) {
        fatalBgzfInput(P, state.inputError.empty()
            ? "a fused input reader failed" : state.inputError);
    }

    state.pipelineDone.store(true, std::memory_order_relaxed);
    pthread_join(reporterThread, nullptr);
    P.inOut->logMain << "  Stats reporter joined\n" << std::flush;

    for (int i = 0; i < nFusedThreads; ++i) {
        g_statsAll.addStats(fusedStats[i]);
    }
    g_statsAll.readN = state.counters.readsTotal.load();

    P.inOut->logMain << "Flex pipeline complete: total=" << state.counters.readsTotal.load()
                     << ", triageKeep=" << state.counters.triageKeep.load()
                     << ", triageDeny=" << state.counters.triageDeny.load()
                     << ", sampleReject=" << state.counters.triageSampleReject.load()
                     << ", triageMiss=" << state.counters.triageMiss.load() << "\n";
    for (int i = 0; i < nLanes; ++i) {
        P.inOut->logMain << "  Lane " << i << ": " << state.counters.perLaneReads[i] << " reads\n";
    }
    P.inOut->logMain << std::flush;

    time(&g_statsAll.timeFinishMap);
    *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeFinishMap) << " ..... finished mapping\n"
                        << std::flush;
    P.inOut->logMain << timeMonthDayTime(g_statsAll.timeFinishMap) << " ..... finished mapping\n"
                     << "RAM after mapping:\n"
                     << linuxProcMemory() << std::flush;

    SoloReadBarcode readBarSum(P);
    {
        for (int i = 0; i < nFusedThreads; ++i) {
            readBarSum.addCounts(*fusedBars[i]);
            readBarSum.addStats(*fusedBars[i]);
        }
        ofstream *statsStream = &ofstrOpen(P.outFileNamePrefix + P.pSolo.outFileNames[0] + "Barcodes.stats",
                                           ERROR_OUT, P);
        readBarSum.statsOut(*statsStream);
        statsStream->close();

        if (P.pSolo.CBmatchWL.mm1_multi_pc) {
            for (uint32 ii = 0; ii < P.pSolo.cbWLsize; ii++) {
                readBarSum.cbReadCountExact[ii]++;
            }
        }
    }

    bool quantYesOriginal = P.quant.yes;
    P.quant.yes = false;
    Transcriptome transcriptomeNoGenome(P);
    P.quant.yes = quantYesOriginal;

    std::vector<SoloFeature *> soloFeatAll(P.pSolo.nFeatures, nullptr);
    SoloFeature geneFeature(P, nullptr, transcriptomeNoGenome, SoloFeatureTypes::Gene,
                            &readBarSum, soloFeatAll.data());
    uint32 geneFeatureIndex = P.pSolo.featureInd[SoloFeatureTypes::Gene];
    soloFeatAll[geneFeatureIndex] = &geneFeature;
    for (int i = 0; i < nFusedThreads; ++i) {
        geneFeature.readFeatAll[i] = fusedFeats[i];
    }

    *P.inOut->logStdOut << timeMonthDayTime() << " ..... started Solo counting\n" << std::flush;
    P.inOut->logMain << timeMonthDayTime() << " ..... started Solo counting\n" << std::flush;
    geneFeature.processRecords();
    *P.inOut->logStdOut << timeMonthDayTime() << " ..... finished Solo counting\n" << std::flush;
    P.inOut->logMain << timeMonthDayTime() << " ..... finished Solo counting\n" << std::flush;

    if (P.pSolo.inlineCBCorrection) {
        std::unordered_map<uint64_t, InlineCBCorrection::MergedAmbigEntry> mergedAmbig;
        InlineCBCorrection::mergeAmbiguousShards(mergedAmbig);
        size_t mergedUnique = mergedAmbig.size();
        uint64_t mergedTotal = 0;
        size_t mergedMaxParents = 0;
        for (const auto &kv : mergedAmbig) {
            mergedTotal += kv.second.count;
            if (kv.second.parents.size() > mergedMaxParents) {
                mergedMaxParents = kv.second.parents.size();
            }
        }

        static const bool disableAmbigResolve =
            (std::getenv("STAR_DISABLE_AMBIG_CB_RESOLVE") != nullptr);
        InlineCBCorrection::AmbigResolveStats ambigStats;
        std::vector<uint32_t> resolvedIdx;
        if (!disableAmbigResolve) {
            ambigStats = InlineCBCorrection::resolveAmbiguousMerged(mergedAmbig, P.pSolo, resolvedIdx);
        }

        P.inOut->logMain << "[INLINE-CB] exact=" << getInlineCbExactCount()
                         << " corrected=" << getInlineCbCorrectedCount()
                         << " n_rescued=" << getInlineCbNRescuedCount()
                         << " with_N=" << getInlineCbWithNCount()
                         << " rejected=" << getInlineCbRejectedCount()
                         << " parent_evidence=" << InlineCBCorrection::parentEvidenceTotal()
                         << " ambig_captured=" << InlineCBCorrection::ambigCapturedTotal()
                         << " ambig_unique=" << InlineCBCorrection::ambigUniqueVariants()
                         << " ambig_max_parents=" << InlineCBCorrection::ambigMaxParents()
                         << " ambig_inline_resolved=" << InlineCBCorrection::getResolvedAmbigCount()
                         << " ambig_merged_unique=" << mergedUnique
                         << " ambig_merged_total=" << mergedTotal
                         << " ambig_merged_max_parents=" << mergedMaxParents
                         << " ambig_resolved=" << ambigStats.resolved
                         << " ambig_ambiguous=" << ambigStats.ambiguous
                         << " ambig_unresolved=" << ambigStats.unresolved
                         << "\n";
        InlineCBCorrection::clearWhitelist();
        InlineCBCorrection::clearAmbiguous();
    }

    for (int i = 0; i < nFusedThreads; ++i) {
        delete fusedFeats[i];
        fusedFeats[i] = nullptr;
        geneFeature.readFeatAll[i] = nullptr;
    }
}

void mapThreadsSpawn (Parameters &P, ReadAlignChunk** RAchunk) {
    // Check activation guard for Flex pipeline mode
    std::string flexActivationReason;
    const bool flexPipelineActive =
        flexPipelineActivationGuard(P, &flexActivationReason, true);
    const bool interfaceEnabled = (P.dynamicThreadInterface == 1);
    const bool telemetryEnabled = (P.dynamicThreadTelemetry == 1);
    const bool variableThreadsEnabled = (P.variableThreads == 1);
    // Permit-pool budget. With chromapAtac concurrent the pool spans STAR's
    // GEX MAP/FEATURE workers AND chromap's ATAC workers, so the budget is
    // runThreadN + chromapAtac.threads (a separate thread budget than runThreadN
    // alone). When chromapAtac is off, pool stays at runThreadN. The user can
    // override via --dynamicThreadConstMapPermits, which is now honored as-is
    // (no clamp to runThreadN); pass 0 for the auto-sized default.
    const int permitTotalThreads = (P.chromapAtac.enabled == 1)
        ? (P.runThreadN + std::max(0, P.chromapAtac.threads))
        : P.runThreadN;
    const int configuredPermits = (P.dynamicThreadConstMapPermits > 0)
        ? P.dynamicThreadConstMapPermits
        : permitTotalThreads;
    g_threadChunks.mapPermitConfigure(interfaceEnabled, permitTotalThreads, configuredPermits, telemetryEnabled, variableThreadsEnabled);
    g_threadChunks.mapPermitConfigureCpuAware(
        interfaceEnabled && P.dynamicThreadPfControllerCpuAware == 1,
        P.dynamicThreadPfControllerCpuSampleMs,
        P.dynamicThreadPfControllerCpuEmaAlpha
    );
    g_threadChunks.mapPermitConfigureRetunePlan(P.variableThreadsPermitSequence, P.variableThreadsRetuneEveryAcquires);

    // Per-domain borrowable floors (Step 5a). Index order must match
    // ThreadControl::permitDomainIndex(): MAP=0, FEATURE=1, ATAC=2.
    int configuredMapFloor = std::max(0, P.dynamicThreadMapFloor);
    int configuredFeatureFloor = std::max(0, P.dynamicThreadFeatureFloor);
    int configuredAtacFloor = std::max(0, P.dynamicThreadAtacFloor);
    if (interfaceEnabled && P.chromapAtac.enabled == 1 &&
        P.dynamicThreadAtacController == 2) {
        const bool featureActive = P.dynamicThreadFeatureWorkEstimate > 0;
        if (featureActive && configuredPermits < 3) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR: three-domain saturation "
                   << "control requires at least 3 configured permits, got "
                   << configuredPermits;
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, 1, P);
        }
        star::multiome::SaturationPermitController::Config controllerConfig;
        controllerConfig.configuredPermits = configuredPermits;
        controllerConfig.fixedFeatureFloor = configuredFeatureFloor;
        controllerConfig.featureActive = featureActive;
        controllerConfig.workEstimates.map = P.dynamicThreadMapWorkEstimate;
        controllerConfig.workEstimates.feature = P.dynamicThreadFeatureWorkEstimate;
        controllerConfig.workEstimates.atac = P.dynamicThreadAtacWorkEstimate;
        const star::multiome::SaturationPermitController controller(controllerConfig);
        const auto initial = controller.initialDecision();
        configuredMapFloor = initial.mapFloor;
        configuredFeatureFloor = initial.featureFloor;
        configuredAtacFloor = initial.atacFloor;
    }
    {
        std::vector<int> domainFloors(3, 0);
        domainFloors[0] = configuredMapFloor;
        domainFloors[1] = configuredFeatureFloor;
        domainFloors[2] = configuredAtacFloor;
        g_threadChunks.mapPermitConfigureDomainFloors(domainFloors);
    }
    // FIFO waiter queue (Step 7). When enabled, ThreadControl serves
    // queued waiters in strict arrival order; new arrivals cannot
    // fast-path past existing waiters.
    g_threadChunks.mapPermitConfigureFifoWaiters(P.dynamicThreadFifoWaiters == 1);

    if (interfaceEnabled) {
        pthread_mutex_lock(&g_threadChunks.mutexLogMain);
        P.inOut->logMain << "Dynamic thread interface enabled: map permits=" << configuredPermits
                         << " (runThreadN=" << P.runThreadN
                         << ", telemetry=" << (telemetryEnabled ? "on" : "off")
                         << ", variableThreads=" << (variableThreadsEnabled ? "on" : "off")
                         << ", cpuAware=" << ((P.dynamicThreadPfControllerCpuAware == 1) ? "on" : "off")
                         << ", cpuSampleMs=" << P.dynamicThreadPfControllerCpuSampleMs
                         << ", cpuEmaAlpha=" << P.dynamicThreadPfControllerCpuEmaAlpha
                         << ", retuneEveryAcquires=" << P.variableThreadsRetuneEveryAcquires
                         << ", retuneSequenceLength=" << P.variableThreadsPermitSequence.size()
                         << ", floors(map/feature/atac)="
                         << configuredMapFloor << "/"
                         << configuredFeatureFloor << "/"
                         << configuredAtacFloor
                         << ", atacController=" << P.dynamicThreadAtacController
                         << ", fifo=" << ((P.dynamicThreadFifoWaiters == 1) ? "on" : "off")
                         << ")\n" << flush;
        pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
    }

    // Flex shares the same controller as ordinary mapping. This check must be
    // after controller initialization so fused alignment and BGZF inflater
    // workers can borrow from one runThreadN-sized pool.
    if (flexPipelineActive) {
        mapThreadsSpawnFlexPipeline(P, RAchunk);
        return;
    }
    if (P.readFilesBgzfMode == "range") {
        fatalBgzfRangeMode(P, flexActivationReason.empty()
            ? "command is not a supported fused Flex run"
            : flexActivationReason);
    }

    for (int ithread=1;ithread<P.runThreadN;ithread++) {//spawn threads
        int threadStatus=pthread_create(&g_threadChunks.threadArray[ithread], NULL, &g_threadChunks.threadRAprocessChunks, (void *) RAchunk[ithread]);
        if (threadStatus>0) {//something went wrong with one of threads
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: phtread error while creating thread # " << ithread <<", error code: "<<threadStatus ;
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, 1, P);
        };
        pthread_mutex_lock(&g_threadChunks.mutexLogMain);
        P.inOut->logMain << "Created thread # " <<ithread <<"\n"<<flush;
        pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
    };

    RAchunk[0]->processChunks(); //start main thread

    for (int ithread=1;ithread<P.runThreadN;ithread++) {//wait for all threads to complete
        int threadStatus = pthread_join(g_threadChunks.threadArray[ithread], NULL);
        if (threadStatus>0) {//something went wrong with one of threads
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: phtread error while joining thread # " << ithread <<", error code: "<<threadStatus ;
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, 1, P);
        };
        pthread_mutex_lock(&g_threadChunks.mutexLogMain);
        P.inOut->logMain << "Joined thread # " <<ithread <<"\n"<<flush;
        pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
    };

    if (interfaceEnabled) {
        g_threadChunks.mapPermitMarkDomainComplete(
            ThreadControl::PermitDomain::MAP);
    }

    if (interfaceEnabled && telemetryEnabled) {
        ThreadControl::MapPermitSnapshot snapshot = g_threadChunks.mapPermitSnapshot();
        const double waitMsTotal = snapshot.waitNsTotal / 1.0e6;
        const double waitMsMax = snapshot.waitNsMax / 1.0e6;
        const double workMsTotal = snapshot.workNsTotal / 1.0e6;
        const double workMsMax = snapshot.workNsMax / 1.0e6;
        const double lastReleaseAgoMs = snapshot.lastReleaseAgoNs / 1.0e6;
        const double avgWaitMs = snapshot.acquireCalls > 0 ? waitMsTotal / static_cast<double>(snapshot.acquireCalls) : 0.0;
        const double avgWorkMs = snapshot.acquireCalls > 0 ? workMsTotal / static_cast<double>(snapshot.acquireCalls) : 0.0;
        const double elapsedNs = static_cast<double>(snapshot.telemetryElapsedNs);
        const double mapOccupancy = elapsedNs > 0.0
            ? snapshot.mapDomain.inUsePermitNs / elapsedNs : 0.0;
        const double featureOccupancy = elapsedNs > 0.0
            ? snapshot.featureDomain.inUsePermitNs / elapsedNs : 0.0;
        const double atacOccupancy = elapsedNs > 0.0
            ? snapshot.atacDomain.inUsePermitNs / elapsedNs : 0.0;
        const double idlePermitAverage = elapsedNs > 0.0
            ? snapshot.availablePermitNs / elapsedNs : 0.0;

        pthread_mutex_lock(&g_threadChunks.mutexLogMain);
        P.inOut->logMain << "Dynamic thread telemetry: acquires=" << snapshot.acquireCalls
                         << ", retunes=" << snapshot.retuneCalls
                         << ", blockedAcquires=" << snapshot.blockedAcquireCalls
                         << ", waitTimeouts=" << snapshot.waitTimeoutEvents
                         << ", stallWarnings=" << snapshot.stallWarnEvents
                         << ", waiters(current/max)=" << snapshot.currentWaiters << "/" << snapshot.maxWaiters
                         << ", retuneEveryAcquires=" << snapshot.retuneEveryAcquires
                         << ", retuneSequenceLength=" << snapshot.sequenceLength
                         << ", targetPermits=" << snapshot.targetPermits
                         << ", configuredPermits=" << snapshot.configuredPermits
                         << ", inUsePermits=" << snapshot.inUsePermits
                         << ", cpuAware=" << (snapshot.cpuAwareEnabled ? "on" : "off")
                         << ", cpuInitialized=" << (snapshot.cpuInitialized ? "yes" : "no")
                         << ", cpuSampleCount=" << snapshot.cpuSampleCount
                         << ", cpuBusyInstant=" << snapshot.cpuBusyInstant
                         << ", cpuBusyEma=" << snapshot.cpuBusyEma
                         << ", cpuIdleEma=" << snapshot.cpuIdleEma
                         << ", retuneTrace=" << serializeRetuneTrace(snapshot.retuneTraceTargets)
                         << ", retuneTraceDropped=" << snapshot.retuneTraceDropped
                         << ", lastReleaseAgoMs=" << lastReleaseAgoMs
                         << ", workUnits=" << snapshot.workUnitsTotal
                         << ", workBytes=" << snapshot.workBytesTotal
                         << ", waitMs(total/avg/max)=" << waitMsTotal << "/" << avgWaitMs << "/" << waitMsMax
                         << ", workMs(total/avg/max)=" << workMsTotal << "/" << avgWorkMs << "/" << workMsMax
                         << ", floorsActive=" << (snapshot.floorsActive ? "yes" : "no")
                         << ", fifo=" << (snapshot.fifoEnabled ? "on" : "off")
                         << ", fifoDepth=" << snapshot.fifoQueueDepth
                         << ", occupancyAvg(map/feature/atac/idle)="
                         << mapOccupancy << "/" << featureOccupancy << "/"
                         << atacOccupancy << "/" << idlePermitAverage
                         << ", contendedIdlePermitMs="
                         << (snapshot.contendedIdlePermitNs / 1.0e6)
                         << ", noAdmissibleGrants="
                         << snapshot.noAdmissibleGrantEvents
                         << ", floorChanges=" << snapshot.floorChangeCalls
                         << ", domainWork(mapUnits,mapBytes,atacUnits,atacBytes)="
                         << snapshot.mapDomain.workUnitsTotal << ","
                         << snapshot.mapDomain.workBytesTotal << ","
                         << snapshot.atacDomain.workUnitsTotal << ","
                         << snapshot.atacDomain.workBytesTotal
                         << ", mapState(floor,inUse,maxInUse,waiters,maxWaiters,blocked,fast,queued,releases)="
                         << snapshot.mapDomain.floor << ","
                         << snapshot.mapDomain.inUse << ","
                         << snapshot.mapDomain.maxInUse << ","
                         << snapshot.mapDomain.currentWaiters << ","
                         << snapshot.mapDomain.maxWaiters << ","
                         << snapshot.mapDomain.blockedAcquireCalls << ","
                         << snapshot.mapDomain.fastAcquireCalls << ","
                         << snapshot.mapDomain.queuedGrantCalls << ","
                         << snapshot.mapDomain.releaseCalls
                         << ", atacState(floor,inUse,maxInUse,waiters,maxWaiters,blocked,fast,queued,releases)="
                         << snapshot.atacDomain.floor << ","
                         << snapshot.atacDomain.inUse << ","
                         << snapshot.atacDomain.maxInUse << ","
                         << snapshot.atacDomain.currentWaiters << ","
                         << snapshot.atacDomain.maxWaiters << ","
                         << snapshot.atacDomain.blockedAcquireCalls << ","
                         << snapshot.atacDomain.fastAcquireCalls << ","
                         << snapshot.atacDomain.queuedGrantCalls << ","
                         << snapshot.atacDomain.releaseCalls
                         << "\n" << flush;
        pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
    }
};
