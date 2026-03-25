#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "systemFunctions.h"
#include "SoloReadInfoLoader.h"
#include "SoloReadInfoSink.h"
#include "hash_shims_cpp_compat.h"  // For unpackReadIdCbUmi
#include "ErrorWarning.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>

namespace {
double soloElapsedSeconds(const std::chrono::steady_clock::time_point &start)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

bool nonFlexHashBridgeEnabled()
{
    return std::getenv("STAR_SOLO_NONFLEX_HASH_BRIDGE") != nullptr;
}

bool nonFlexHashBridgeApplies(const SoloFeature& feat)
{
    if (!feat.pSolo.inlineHashMode || feat.pSolo.flexMode || !nonFlexHashBridgeEnabled()) {
        return false;
    }
    switch (feat.featureType) {
        case SoloFeatureTypes::Gene:
        case SoloFeatureTypes::GeneFull:
        case SoloFeatureTypes::GeneFull_Ex50pAS:
        case SoloFeatureTypes::GeneFull_ExonOverIntron:
            return true;
        default:
            return false;
    }
}

void populateBridgeReadAccounting(SoloFeature &feat,
                                  std::vector<uint32_t> &nReadPerCBunique1,
                                  std::vector<uint32_t> &nReadPerCBmulti1)
{
    feat.readFlagCounts.flagCounts.clear();
    feat.readFlagCounts.flagCounts.swap(feat.readFeatSum->readFlag.flagCounts);
    feat.readFlagCounts.flagCountsNoCB = feat.readFeatSum->readFlag.flagCountsNoCB;
    feat.readFeatSum->readFlag.flagCountsNoCB = {};

    for (const auto &kv : feat.readFeatSum->bridgeImmediateReadCounts_) {
        if (kv.first >= nReadPerCBunique1.size()) {
            continue;
        }
        nReadPerCBunique1[kv.first] += static_cast<uint32_t>(kv.second & 0xFFFFFFFFu);
        nReadPerCBmulti1[kv.first] += static_cast<uint32_t>(kv.second >> 32);
    }

    const size_t pinN = feat.readFeatSum->bridgePinNreadUnique_.size();
    if (pinN > 0) {
        const size_t lim = std::min(nReadPerCBunique1.size(), pinN);
        for (size_t i = 0; i < lim; ++i) {
            nReadPerCBunique1[i] += feat.readFeatSum->bridgePinNreadUnique_[i];
            nReadPerCBmulti1[i] += feat.readFeatSum->bridgePinNreadMulti_[i];
        }
        feat.readFeatSum->bridgePinNreadUnique_.clear();
        feat.readFeatSum->bridgePinNreadMulti_.clear();
    }

    decltype(feat.readFeatSum->bridgeImmediateReadCounts_)().swap(feat.readFeatSum->bridgeImmediateReadCounts_);
}
}

void SoloFeature::countCBgeneUMI()
{    
    const auto countStart = std::chrono::steady_clock::now();

    // Skip legacy Solo counting when inline CB correction is active, but continue to inline-hash flow
    if (pSolo.inlineCBCorrection) {
        // Assert that Solo structures are unused
        assert(packedReadInfo.data.empty());
        assert(rGeneUMI == nullptr);
        P.inOut->logMain << "Skipping legacy Solo counting (inline CB correction active); continuing with inline hash collapse" << endl;
    }
    
    time_t rawTime;
    
    rguStride=2;
    if (pSolo.readIndexYes[featureType])
        rguStride=3; //to keep readI column

#ifdef DEBUG_CB_UB_PARITY
    // Skip parity validation when minimal memory flag is on (parity requires packed storage)
    if (pSolo.soloFlexMinimalMemory && pSolo.inlineHashMode) {
        parityEnabled = false;
    } else {
        bool parityEnv = (std::getenv("STAR_DEBUG_CB_UB_PARITY") != nullptr);
        if (nonFlexBridgePath) {
            parityEnv = false;
        }
        if (parityEnv && packedReadInfo.data.empty()) {
            resetPackedStorage(nReadsInput);
        }
        resetDebugStatusCounters();
        parityEnabled = parityEnv;
    }
#endif

    const bool nonFlexBridgePath = nonFlexHashBridgeApplies(*this);

    // Allocate packedReadInfo if:
    // 1. readInfoYes is set for this feature type, OR
    // 2. trackReadIdsForTags is enabled (for sorted BAM CB/UB tag injection)
    // Non-Flex direct bridge mode does not use packed read info for this benchmark path.
    // Skip only if soloFlexMinimalMemory is on AND inlineHashMode is on AND trackReadIdsForTags is off
    bool needPackedReadInfo = pSolo.trackReadIdsForTags
        || (pSolo.readInfoYes[featureType] && !nonFlexBridgePath);
    bool skipForMinimalMemory = pSolo.soloFlexMinimalMemory && pSolo.inlineHashMode && !pSolo.trackReadIdsForTags;
    if (needPackedReadInfo && !skipForMinimalMemory) {
        resetPackedStorage(nReadsInput);
        time(&rawTime);
#ifdef SOLO_USE_PACKED_READINFO
        P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Allocated and initialized packed readInfo array, nReadsInput = " << nReadsInput <<endl;
#else
        P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Allocated and initialized readInfo array, nReadsInput = " << nReadsInput <<endl;
#endif
    };
    
    // Inline hash path: resolve/correct, then walk the hash directly (no materialization)
    if (pSolo.inlineHashMode) {
        // Resolve ambiguous CBs (before collapse)
        resolveAmbiguousCBs();
        
        // Run clique correction if enabled (operates on hash)
        if (pSolo.umiCorrectionMode > 0) {
            runCliqueCorrection();
        }

        if (nonFlexBridgePath) {
            time(&rawTime);
            P.inOut->logMain << timeMonthDayTime(rawTime)
                             << " ... Experimental non-Flex inline-hash bridge: direct hash collapse"
                             << " (no materializeRGUFromHash / no legacy collapseUMIall)"
                             << endl;

            const auto bridgeCollapseStart = std::chrono::steady_clock::now();
            collapseUMIall_fromBridgeHash();
            P.inOut->logMain << "Solo timing: collapseUMIall_fromBridgeHash "
                             << soloElapsedSeconds(bridgeCollapseStart) << " s" << endl;

            if (nCB == 0) {
                P.inOut->logMain << "WARNING: non-Flex inline-hash bridge produced no cells after direct collapse"
                                 << endl;
                P.inOut->logMain << "Solo timing: countCBgeneUMI " << soloElapsedSeconds(countStart) << " s" << endl;
                return;
            }

            std::vector<uint32> nReadPerCBunique1(pSolo.cbWLsize), nReadPerCBmulti1(pSolo.cbWLsize);
            populateBridgeReadAccounting(*this, nReadPerCBunique1, nReadPerCBmulti1);

            nReadPerCBunique.resize(nCB);
            nReadPerCBtotal.resize(nCB);
            for (uint32 icb = 0; icb < nCB; icb++) {
                uint32 wlIndex = indCB[icb];
                nReadPerCBunique[icb] = nReadPerCBunique1[wlIndex];
                nReadPerCBtotal[icb] = nReadPerCBunique1[wlIndex] + nReadPerCBmulti1[wlIndex];
            }

            for (uint32 icb = 0; icb < nCB; icb++) {
                readFeatSum->stats.V[readFeatSum->stats.yesUMIs] += nUMIperCB[icb];
                if (nGenePerCB[icb] > 0)
                    ++readFeatSum->stats.V[readFeatSum->stats.yesCellBarcodes];
                readFeatSum->stats.V[readFeatSum->stats.yesWLmatch] += nReadPerCBtotal[icb];
                readFeatSum->stats.V[readFeatSum->stats.yessubWLmatch_UniqueFeature] += nReadPerCBunique[icb];
            }

            time(&rawTime);
            P.inOut->logMain << timeMonthDayTime(rawTime)
                             << " ... Finished collapsing UMIs (non-Flex inline-hash bridge, direct hash)" << endl;
            P.inOut->logMain << "Solo timing: countCBgeneUMI " << soloElapsedSeconds(countStart) << " s" << endl;
            return;
        }
        
        // Direct hash consumption: no materialization/legacy collapse
        const auto hashCollapseStart = std::chrono::steady_clock::now();
        collapseUMIall_fromHash();
        P.inOut->logMain << "Solo timing: collapseUMIall_fromHash " << soloElapsedSeconds(hashCollapseStart) << " s" << endl;
        
        // Populate packedReadInfo from readIdTracker_ for sorted BAM CB/UB tag injection
        if (pSolo.trackReadIdsForTags && readFeatSum && readFeatSum->readIdTracker_) {
            time(&rawTime);
            P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Populating packedReadInfo from readIdTracker for sorted BAM CB/UB tags" << endl;
            
            size_t trackerSize = kh_size(readFeatSum->readIdTracker_);
            size_t populated = 0;
            
            for (khiter_t iter = kh_begin(readFeatSum->readIdTracker_); iter != kh_end(readFeatSum->readIdTracker_); ++iter) {
                if (!kh_exist(readFeatSum->readIdTracker_, iter)) continue;
                
                uint32_t readId = kh_key(readFeatSum->readIdTracker_, iter);
                uint64_t val = kh_val(readFeatSum->readIdTracker_, iter);
                
                uint32_t cbIdx, umi24;
                uint8_t status;
                unpackReadIdCbUmi(val, &cbIdx, &umi24, &status);
                
                // Record into packedReadInfo
                if (readId < nReadsInput) {
                    recordReadInfo(readId, cbIdx, umi24, status);
                    populated++;
                }
            }
            
            P.inOut->logMain << "  Populated " << populated << " readIds from tracker (tracker size: " << trackerSize << ")" << endl;
        }
        
        // Export readId/CB/UB/status TSV table (env var gated, after CB/UB finalized)
        writeReadIdTagTable();
        
        time(&rawTime);
        P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished collapsing UMIs (direct hash mode)" << endl;
        P.inOut->logMain << "Solo timing: countCBgeneUMI " << soloElapsedSeconds(countStart) << " s" << endl;
        return;
    }
    
    // Early return when minimal memory flag is on (prevents accidental fallthrough if code is reordered)
    if (pSolo.soloFlexMinimalMemory && pSolo.inlineHashMode) {
        return; // Skip legacy path entirely
    }
    
    // Packed-only path: parse with loader, buffer via CountingSink, and collapse.
    {
        readFlagCounts.flagCounts.reserve((pSolo.cbWLsize ? pSolo.cbWLsize : 1)*3/2);
        readFlagCounts.flagCountsNoCB = {};
        vector<uint32> nReadPerCBunique1(pSolo.cbWLsize), nReadPerCBmulti1(pSolo.cbWLsize);

        SoloReadInfoLoader loader;
        CountingSink sink;
        for (int ii=0; ii<P.runThreadN; ii++) {
            // Defensive check: verify readFeatAll[ii] and its streamReads are valid
            // This guards against wiring issues in SoloFeature_sumThreads or featureInd mismatch
            if (readFeatAll[ii] == nullptr) {
                ostringstream errOut;
                errOut << "EXITING because of fatal ERROR: readFeatAll[" << ii << "] is null in countCBgeneUMI\n"
                       << "This indicates a wiring issue in SoloFeature_sumThreads.cpp\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            }
            if (!pSolo.inlineHashMode && !readFeatAll[ii]->binarySpoolInMemory && readFeatAll[ii]->streamReads == nullptr) {
                ostringstream errOut;
                errOut << "EXITING because of fatal ERROR: readFeatAll[" << ii << "]->streamReads is null in non-inline-hash mode\n"
                       << "featureType=" << featureType << " thread=" << ii << "\n"
                       << "This indicates streamReads was not opened during mapping phase\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            }
            loader.load(*readFeatAll[ii], SoloReadInfoMode::Counting,
                        [&](const ReadInfoRecord &rec){ sink.onRecord(*this, rec); },
                        readBarSum->cbReadCountExact,
                        readFlagCounts,
                        nReadPerCBunique1, nReadPerCBmulti1);
            readFeatSum->addStats(*readFeatAll[ii]);
            if (readFeatAll[ii]->binarySpoolInMemory) {
                readFeatAll[ii]->binarySpoolBuffer.clear();
                readFeatAll[ii]->binarySpoolReadPos = 0;
            }
        }
        readFlagCounts.countsAddNoCBarray(readFeatSum->readFlag.flagCountsNoCB);

        sink.finalize(*this);

#ifdef DEBUG_CB_UB_PARITY
        // Skip parity validation when minimal memory flag is on (parity requires packed storage)
        if (pSolo.soloFlexMinimalMemory && pSolo.inlineHashMode) {
            // parityEnabled already set to false in first DEBUG_CB_UB_PARITY block
        } else {
            // Initialize parity validation if enabled
            parityEnabled = (std::getenv("STAR_DEBUG_CB_UB_PARITY") != nullptr);
            if (parityEnabled) {
                std::string mismatchPath = P.outFileNamePrefix + "Aligned.out.cbub_mismatches.txt";
                mismatchLog.open(mismatchPath, std::ios::trunc);
                if (mismatchLog.is_open()) {
                mismatchLog << "iCB\tiRead\tlegacyCB\tinlineCB\tlegacyUB\tinlineUB\tlegacyGene\tinlineGene\tumi\tcliqueSize\tinlineStatus\tinlineCbIdx\treason\n";
                }
                parityReadsLegacy = parityReadsInline = parityMatches = parityMismatches = 0;
            }
        }
#endif

        // After finalize(), rGeneUMI/rCBp are filled for collapseUMIall().
        // Compute per-CB sizes and initialize matrices to mirror legacy path.
        nReadPerCB.resize(nCB);
        nReadPerCBmax=0;
        for (uint32 iCB=0; iCB<nCB; iCB++) {
            nReadPerCB[iCB] = rCBn[iCB];
            if (nReadPerCB[iCB] > nReadPerCBmax) nReadPerCBmax = nReadPerCB[iCB];
        }

        // Populate per-CB unique/total from loader-accumulated vectors
        nReadPerCBunique.resize(nCB);
        nReadPerCBtotal.resize(nCB);
        for (uint32 icb=0; icb<nCB; icb++) {
            uint32 wlIndex = indCB[icb];
            nReadPerCBunique[icb] = nReadPerCBunique1[wlIndex];
            nReadPerCBtotal[icb] = nReadPerCBunique1[wlIndex] + nReadPerCBmulti1[wlIndex];
        }

        // Initialize count matrices similar to legacy path
        nUMIperCB.resize(nCB);
        nGenePerCB.resize(nCB);
        countMatStride = pSolo.umiDedup.yes.N + 1;
        countCellGeneUMI.resize(nReadsMapped*countMatStride/5+16);
        countCellGeneUMIindex.resize(nCB+1, 0);
        if (pSolo.multiMap.yes.multi) {
            countMatMult.s = 1 + pSolo.multiMap.yes.N * pSolo.umiDedup.yes.N;
            countMatMult.m.resize(nReadsMapped*countMatMult.s/5+16);
            countMatMult.i.resize(nCB+1, 0);
        }

        // Collapse UMIs once here; CountingSink no longer calls collapseUMIall.
        const auto collapseStart = std::chrono::steady_clock::now();
        collapseUMIall();
        P.inOut->logMain << "Solo timing: collapseUMIall (countCBgeneUMI wrapper) " << soloElapsedSeconds(collapseStart) << " s" << endl;

        // Free temporary arrays allocated via CountingSink::finalize
        if (rGeneUMI) { delete[] rGeneUMI; rGeneUMI=nullptr; }
        if (rCBp) { delete[] rCBp; rCBp=nullptr; }
        if (rCBn) { delete[] rCBn; rCBn=nullptr; }

        P.inOut->logMain << "RAM for solo feature "<< SoloFeatureTypes::Names[featureType] <<"\n" <<  linuxProcMemory() << flush;

        // rGeneUMI/rCBp/rCBn already freed above

        time(&rawTime);
        P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished collapsing UMIs" <<endl;
#ifdef DEBUG_CB_UB_PARITY
        logDebugStatusCounters();
        logDebugStageCounters();
#endif
        P.inOut->logMain << "Solo timing: countCBgeneUMI " << soloElapsedSeconds(countStart) << " s" << endl;
        return;
    }
};
