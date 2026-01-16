#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "systemFunctions.h"
#include "SoloReadInfoLoader.h"
#include "SoloReadInfoSink.h"

void SoloFeature::countCBgeneUMI()
{    
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
        if (parityEnv && packedReadInfo.data.empty()) {
            resetPackedStorage(nReadsInput);
        }
        resetDebugStatusCounters();
        parityEnabled = parityEnv;
    }
#endif

    if (pSolo.readInfoYes[featureType] && !(pSolo.soloFlexMinimalMemory && pSolo.inlineHashMode)) {
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
        
        // Direct hash consumption: no materialization/legacy collapse
        collapseUMIall_fromHash();
        
        // Export readId/CB/UB/status TSV table (env var gated, after CB/UB finalized)
        writeReadIdTagTable();
        
        time(&rawTime);
        P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished collapsing UMIs (direct hash mode)" << endl;
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
            loader.load(*readFeatAll[ii], SoloReadInfoMode::Counting,
                        [&](const ReadInfoRecord &rec){ sink.onRecord(*this, rec); },
                        readBarSum->cbReadCountExact,
                        readFlagCounts,
                        nReadPerCBunique1, nReadPerCBmulti1);
            readFeatSum->addStats(*readFeatAll[ii]);
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

        // Run clique correction if enabled (before collapseUMIall)
        if (pSolo.umiCorrectionMode > 0) {
            runCliqueCorrection();
        }

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
        collapseUMIall();

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
        return;
    }
};
