#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "SoloReadInfoLoader.h"
#include "SoloReadInfoSink.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include <cstdio>

static const bool g_debugPrepare = (std::getenv("STAR_DEBUG_PREPARE") != nullptr);

void SoloFeature::prepareReadInfoOnly()
{
    // Minimal processing to populate readInfo without heavy matrix allocations
    // This is called in skipProcessing mode when tag table or CB/UB injection is needed
    
    time_t rawTime;
    
    if (!pSolo.readInfoYes[featureType]) {
        // readInfo not needed for this feature type
        return;
    }
    
    // Allocate and initialize readInfo (same as countCBgeneUMI init)
    resetPackedStorage(nReadsInput);
    
    // Set rguStride (needed for collapseUMIall)
    rguStride=2;
    if (pSolo.readIndexYes[featureType])
        rguStride=3; //to keep readI column
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) 
                     << " ... Allocated packed readInfo array for skipProcessing mode, nReadsInput = " 
                     << nReadsInput << endl;
    
    // Use loader + MinimalSink to populate readInfo directly without allocating rGeneUMI/CB arrays
    readFlagCounts.flagCounts.clear();
    readFlagCounts.flagCountsNoCB = {};

    SoloReadInfoLoader loader;
    
    // For clique correction: use CountingSink to collect UR histograms and build rGeneUMI
    // This allows clique correction to run and correct UB tags before writing to readInfo
    if (pSolo.umiCorrectionMode > 0) {
        if (g_debugPrepare) {
            fprintf(stderr, "[PREPARE] Starting clique correction path\n");
        }
        
        // Use CountingSink to collect histograms and build rGeneUMI structure
        readFlagCounts.flagCounts.reserve((pSolo.cbWLsize ? pSolo.cbWLsize : 1)*3/2);
        readFlagCounts.flagCountsNoCB = {};
        vector<uint32> nReadPerCBunique1(pSolo.cbWLsize), nReadPerCBmulti1(pSolo.cbWLsize);
        
        if (g_debugPrepare) {
            fprintf(stderr, "[PREPARE] Before CountingSink loop: P.runThreadN=%d\n", P.runThreadN);
        }
        
        CountingSink sink;
        for (int ii=0; ii<P.runThreadN; ii++) {
            if (g_debugPrepare && ii == 0) {
                fprintf(stderr, "[PREPARE] Loading thread %d\n", ii);
            }
            loader.load(*readFeatAll[ii], SoloReadInfoMode::Counting,
                        [&](const ReadInfoRecord &rec){ sink.onRecord(*this, rec); },
                        readBarSum->cbReadCountExact,
                        readFlagCounts,
                        nReadPerCBunique1, nReadPerCBmulti1);
            readFeatSum->addStats(*readFeatAll[ii]);
        }
        
        if (g_debugPrepare) {
            fprintf(stderr, "[PREPARE] Before sink.finalize\n");
        }
        
        readFlagCounts.countsAddNoCBarray(readFeatSum->readFlag.flagCountsNoCB);
        sink.finalize(*this);
        
        if (g_debugPrepare) {
            fprintf(stderr, "[PREPARE] After sink.finalize: nCB=%u rGeneUMI=%p rCBp=%p rCBn=%p\n", 
                    nCB, (void*)rGeneUMI, (void*)rCBp, (void*)rCBn);
        }
        
        // setRGUStride is called by CountingSink::finalize(), but ensure it's set
        setRGUStride(pSolo.readIndexYes[featureType] ? 3u : 2u);
        
        // Skip if no cells detected
        if (nCB == 0 || rGeneUMI == nullptr) {
            P.inOut->logMain << "No cells detected, skipping clique correction" << endl;
            return;
        }
        
        if (g_debugPrepare) {
            fprintf(stderr, "[PREPARE] Before runCliqueCorrection: nCB=%u featuresNumber=%u rguStride=%u\n",
                    nCB, featuresNumber, rguStride);
        }
        
        // Run clique correction if enabled (before applying corrections)
        runCliqueCorrection();
        
        if (g_debugPrepare) {
            fprintf(stderr, "[PREPARE] After runCliqueCorrection\n");
        }
        
        // Initialize per-CB arrays needed for collapseUMIall
        nReadPerCB.resize(nCB);
        nReadPerCBmax=0;
        for (uint32 iCB=0; iCB<nCB; iCB++) {
            nReadPerCB[iCB] = rCBn[iCB];
            if (nReadPerCB[iCB] > nReadPerCBmax) nReadPerCBmax = nReadPerCB[iCB];
        }
        
        // Validate rCBp / nReadPerCB layout before calling collapseUMIall
        if (g_debugPrepare) {
            fprintf(stderr, "[PREPARE] Validating rCBp layout before collapseUMIall\n");
        }
        
        // Calculate total reads for bounds checking
        uint64_t totalReads = 0;
        for (uint32 iCB=0; iCB<nCB; iCB++) {
            totalReads += nReadPerCB[iCB];
        }
        uint64_t totalWords = totalReads * rguStride;
        
        bool layoutValid = true;
        for (uint32 iCB=0; iCB<nCB; iCB++) {
            // Check rCBp[i+1] > rCBp[i] (except for last)
            if (iCB < nCB-1) {
                if (rCBp[iCB+1] <= rCBp[iCB]) {
                    fprintf(stderr, "[PREPARE] ERROR: rCBp[%u+1]=%p <= rCBp[%u]=%p\n",
                            iCB, (void*)rCBp[iCB+1], iCB, (void*)rCBp[iCB]);
                    layoutValid = false;
                    break;
                }
            }
            
            // Check stride matches nReadPerCB
            uint32 expectedStride = nReadPerCB[iCB] * rguStride;
            uint32 actualStride = (iCB < nCB-1) ? (rCBp[iCB+1] - rCBp[iCB]) : 
                                  (rGeneUMI + totalWords - rCBp[iCB]);
            if (actualStride != expectedStride) {
                fprintf(stderr, "[PREPARE] ERROR: CB %u stride mismatch: expected=%u actual=%u\n",
                        iCB, expectedStride, actualStride);
                layoutValid = false;
                break;
            }
            
            // Check rCBp[i] points inside rGeneUMI bounds
            if (rCBp[iCB] < rGeneUMI || rCBp[iCB] >= rGeneUMI + totalWords) {
                fprintf(stderr, "[PREPARE] ERROR: CB %u rCBp[%u]=%p out of bounds [%p, %p)\n",
                        iCB, iCB, (void*)rCBp[iCB], (void*)rGeneUMI, (void*)(rGeneUMI + totalWords));
                layoutValid = false;
                break;
            }
            
            // Dump first few entries for first 3 CBs
            if (iCB < 3 && rCBp[iCB] != nullptr) {
                fprintf(stderr, "[PREPARE] CB %u first 3 entries: [%u,%u,%u] [%u,%u,%u] [%u,%u,%u]\n",
                        iCB,
                        rCBp[iCB][0], rCBp[iCB][1], rCBp[iCB][2],
                        rCBp[iCB][3], rCBp[iCB][4], rCBp[iCB][5],
                        rCBp[iCB][6], rCBp[iCB][7], rCBp[iCB][8]);
            }
        }
        
        if (!layoutValid) {
            P.inOut->logMain << "ERROR: rCBp layout validation failed, aborting collapseUMIall" << endl;
            return;
        }
        
        if (g_debugPrepare) {
            fprintf(stderr, "[PREPARE] rCBp layout validation passed\n");
        }
        
        // Apply corrections via collapseUMIall in minimal mode (writes corrected readInfo)
        try {
            collapseUMIall(true); // minimalMode=true: writes readInfo but doesn't allocate matrices
        } catch (...) {
            P.inOut->logMain << "ERROR: Exception in collapseUMIall" << endl;
            throw;
        }
        
        if (g_debugPrepare) {
            fprintf(stderr, "[PREPARE] After collapseUMIall\n");
        }
        
        // Clean up rGeneUMI/rCBp after writing to readInfo
        if (rGeneUMI) { delete[] rGeneUMI; rGeneUMI=nullptr; }
        if (rCBp) { delete[] rCBp; rCBp=nullptr; }
        if (rCBn) { delete[] rCBn; rCBn=nullptr; }
    } else {
        // No clique correction: use MinimalSink for faster path
    MinimalSink sink;
    for (int ii = 0; ii < P.runThreadN; ii++) {
        loader.loadMinimal(*readFeatAll[ii],
                           [&](const ReadInfoRecord &rec){ sink.onRecord(*this, rec); },
                           readBarSum->cbReadCountExact);
        readFeatSum->addStats(*readFeatAll[ii]);
    }
    sink.finalize(*this);
    }

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime)
                     << " ... Finished populating readInfo via loader/minimal sink (skipProcessing mode)" << endl;
}

