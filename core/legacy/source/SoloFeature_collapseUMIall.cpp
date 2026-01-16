#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "serviceFuns.cpp"
#include <unordered_map>
#include "SoloCommon.h"
#include "UmiCodec.h"
#include "UMICorrector.h"
#include "ErrorWarning.h"
#include <cstdio>
#include <vector>
#include <algorithm>

static const bool g_debugCollapse = (std::getenv("STAR_DEBUG_COLLAPSE") != nullptr);

#ifdef DEBUG_CB_UB_PARITY
// Optional tracing of specific readIds via STAR_DEBUG_TRACE_READS=1,2,3
static std::unordered_set<uint32_t> buildTraceReadSet() {
    std::unordered_set<uint32_t> out;
    const char* env = std::getenv("STAR_DEBUG_TRACE_READS");
    if (!env || env[0]=='\0') return out;
    std::string s(env);
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (tok.empty()) continue;
        try {
            out.insert(static_cast<uint32_t>(std::stoul(tok)));
        } catch (...) {
            // ignore malformed entries
        }
    }
    return out;
}
static const std::unordered_set<uint32_t> g_traceReads = buildTraceReadSet();
#endif

inline int funCompareSolo1 (const void *a, const void *b);          //defined below
inline int funCompare_uint32_1_2_0 (const void *a, const void *b);

// Helper: validate packedReadInfo immediately after a write
#ifdef DEBUG_CB_UB_PARITY
static inline bool dbgCheckReadAfterWrite(SoloFeature *self, uint32_t readId, uint32_t cbExp, uint32_t umiExp, uint8_t statusExp, const char *site) {
    if (!self || !self->parityEnabled) return true;
    uint32_t cbGot = self->getPackedCB(readId);
    uint32_t umiGot = self->getPackedUMI(readId);
    uint8_t statusGot = self->getPackedStatus(readId);
    bool ok = (cbGot==cbExp && umiGot==umiExp && statusGot==statusExp);
    if (ok) {
        ++self->dbgRAWPass;
    } else {
        ++self->dbgRAWFail;
        if (self->dbgRAWFail <= 20) {
            fprintf(stderr, "[RAW-FAIL] site=%s read=%u exp(cb=%u umi=%u s=%u) got(cb=%u umi=%u s=%u)\n",
                    site, readId, cbExp, umiExp, statusExp, cbGot, umiGot, statusGot);
        }
    }
    return ok;
}
#endif

void SoloFeature::collapseUMIall(bool minimalMode)
{
    if (countMatStride == 0 || countMatStride > 1024) {
        // countMatStride may be unused in minimal/skip paths; enforce a sane non-zero stride
        countMatStride = 1;
    }
    // FORCED DEBUG: Always print to stderr to catch crashes early
    fprintf(stderr, "[COLLAPSE-FORCED] Entering collapseUMIall: nCB=%u nReadPerCBmax=%u featuresNumber=%u minimalMode=%d\n",
            nCB, nReadPerCBmax, featuresNumber, minimalMode ? 1 : 0);
    fflush(stderr);
    
    // Dump all member values before allocation to check for uninitialized values
    fprintf(stderr, "[COLLAPSE-FORCED] Member values before allocation:\n");
    fprintf(stderr, "  nCB=%u\n", nCB);
    fprintf(stderr, "  nReadPerCBmax=%u\n", nReadPerCBmax);
    fprintf(stderr, "  featuresNumber=%u\n", featuresNumber);
    fprintf(stderr, "  rguStride=%u\n", rguStride);
    fprintf(stderr, "  umiArrayStride=%u\n", umiArrayStride);
    fprintf(stderr, "  minimalMode=%d\n", minimalMode ? 1 : 0);
    fprintf(stderr, "  pSolo.multiMap.yes.multi=%d\n", pSolo.multiMap.yes.multi ? 1 : 0);
    fprintf(stderr, "  pSolo.umiFiltering.MultiGeneUMI=%d\n", pSolo.umiFiltering.MultiGeneUMI ? 1 : 0);
    fflush(stderr);
    
    // Calculate vector sizes and check for overflow
    uint32 umiArraySize = nReadPerCBmax * umiArrayStride;
    uint32 gIDSize = min(2*featuresNumber, nReadPerCBmax) + 1;
    uint32 gReadSSize = min(2*featuresNumber, nReadPerCBmax) + 1;
    
    fprintf(stderr, "[COLLAPSE-FORCED] Calculated vector sizes:\n");
    fprintf(stderr, "  umiArraySize=%u (bytes=%zu)\n", umiArraySize, (size_t)umiArraySize * sizeof(uint32));
    fprintf(stderr, "  gIDSize=%u (bytes=%zu)\n", gIDSize, (size_t)gIDSize * sizeof(uint32));
    fprintf(stderr, "  gReadSSize=%u (bytes=%zu)\n", gReadSSize, (size_t)gReadSSize * sizeof(uint32));
    fprintf(stderr, "  Total stack allocation: %zu bytes\n", 
            (size_t)(umiArraySize + gIDSize + gReadSSize) * sizeof(uint32));
    fflush(stderr);
    
    if (g_debugCollapse) {
        fprintf(stderr, "[COLLAPSE] Entering collapseUMIall: nCB=%u nReadPerCBmax=%u featuresNumber=%u minimalMode=%d\n",
                nCB, nReadPerCBmax, featuresNumber, minimalMode ? 1 : 0);
    }
    
    fprintf(stderr, "[COLLAPSE-FORCED] Allocating vectors: umiArray=%u gID=%u gReadS=%u\n",
            nReadPerCBmax*umiArrayStride, min(2*featuresNumber,nReadPerCBmax)+1, min(2*featuresNumber,nReadPerCBmax)+1);
    fflush(stderr);
    
    vector<uint32> umiArray(nReadPerCBmax*umiArrayStride);//temp array for collapsing UMI
    vector<uint32> gID(min(2*featuresNumber,nReadPerCBmax)+1); //gene IDs, 2* is needed because each gene can have unique and multi-mappers
    vector<uint32> gReadS(min(2*featuresNumber,nReadPerCBmax)+1); //start of gene reads TODO: allocate this array in the 2nd half of rGU

    fprintf(stderr, "[COLLAPSE-FORCED] Vectors allocated successfully\n");
    fflush(stderr);

    if (!minimalMode) {
        // Ensure per-CB bookkeeping vectors are sized for all cells before we start the main loop.
        if (countCellGeneUMIindex.size() < nCB + 1) {
            countCellGeneUMIindex.assign(nCB + 1, 0);
        }
        if (nGenePerCB.size() < nCB) {
            nGenePerCB.assign(nCB, 0);
        }
        if (nUMIperCB.size() < nCB) {
            nUMIperCB.assign(nCB, 0);
        }
        if (countCellGeneUMI.empty()) {
            // Seed count matrix with a conservative size so countCellGeneUMIindex arithmetic stays valid.
            // Use at least nCB entries of countMatStride.
            uint32_t stride = (countMatStride > 0) ? countMatStride : 1;
            size_t seed = std::max<size_t>(nCB * stride, (size_t)(nCB + 1) * stride);
            countCellGeneUMI.assign(seed, 0);
            fprintf(stderr, "[COLLAPSE-DUMP] initialized countCellGeneUMI to %zu entries (stride=%u)\n", seed, stride);
            fflush(stderr);
        }
    }

    if (g_debugCollapse) {
        fprintf(stderr, "[COLLAPSE] Allocated arrays: umiArray.size()=%zu gID.size()=%zu gReadS.size()=%zu\n",
                umiArray.size(), gID.size(), gReadS.size());
    }

    for (uint32 icb=0; icb<nCB; icb++) {//main collapse cycle
        // FORCED DEBUG: Print before each collapseUMIperCB call
        uint32 rN = nReadPerCB[icb];
        uint32 expectedStride = (icb < nCB-1) ? ((rCBp[icb+1] - rCBp[icb]) / rguStride) : rN;
        fprintf(stderr, "[COLLAPSE-FORCED] CB %u/%u: rCBp[%u]=%p rN=%u expectedStride=%u\n",
                icb, nCB, icb, (void*)rCBp[icb], rN, expectedStride);
        if (icb < 3 && rCBp[icb] != nullptr) {
            fprintf(stderr, "[COLLAPSE-FORCED]   First 3 words: %u %u %u\n",
                    rCBp[icb][0], rCBp[icb][1], rCBp[icb][2]);
        }
        fflush(stderr);
        
        if (g_debugCollapse && icb < 3) {
            fprintf(stderr, "[COLLAPSE] Processing CB %u/%u: rCBp[%u]=%p rN=%u\n",
                    icb, nCB, icb, (void*)rCBp[icb], nReadPerCB[icb]);
        }
        
        collapseUMIperCB(icb, umiArray, gID, gReadS, minimalMode);
        
        if (!minimalMode) {
            readFeatSum->stats.V[readFeatSum->stats.yesUMIs] += nUMIperCB[icb];
            if (nGenePerCB[icb]>0) //nGenePerCB contains only unique
                ++readFeatSum->stats.V[readFeatSum->stats.yesCellBarcodes];
            
            readFeatSum->stats.V[readFeatSum->stats.yesWLmatch] += nReadPerCBtotal[icb];        
            readFeatSum->stats.V[readFeatSum->stats.yessubWLmatch_UniqueFeature ] += nReadPerCBunique[icb];        
        }
    };
    
    // After all CB/UB values are finalized, finalize tag table if export is enabled
    if (!minimalMode && pSolo.writeTagTableEnabled) {
        finalizeTagTableFromReadInfo();
    }
    
    // Export readId/CB/UB/status TSV table (env var gated, after CB/UB finalized, before packedReadInfo cleared)
    if (!minimalMode) {
        writeReadIdTagTable();
    }
    
#ifdef DEBUG_CB_UB_PARITY
    // Write parity counters after all CBs processed
    if (parityEnabled) {
        writeParityCounters();
        if (mismatchLog.is_open()) {
            mismatchLog.flush();
        }
    }
#endif
};

void SoloFeature::collapseUMIperCB(uint32 iCB, vector<uint32> &umiArray, vector<uint32> &gID,  vector<uint32> &gReadS, bool minimalMode)
{
    // FORCED DEBUG: Always print at entry
    fprintf(stderr, "[COLLAPSE-FORCED] collapseUMIperCB entry: iCB=%u\n", iCB);
    fflush(stderr);

    if (g_debugCollapse && iCB < 3) {
        fprintf(stderr, "[COLLAPSE] collapseUMIperCB: iCB=%u rCBp[%u]=%p nReadPerCB[%u]=%u rguStride=%u\n",
                iCB, iCB, (void*)rCBp[iCB], iCB, nReadPerCB[iCB], rguStride);
    }

    uint32 *rGU=rCBp[iCB];
    uint32 rN=nReadPerCB[iCB]; //with multimappers, this is the number of all aligns, not reads
    
    fprintf(stderr, "[COLLAPSE-FORCED] Before qsort: rGU=%p rN=%u rguStride=%u\n", (void*)rGU, rN, rguStride);
    fflush(stderr);
    
    // Capture data before qsort for comparison
    const uint32 maxDumpTriplets = (rN < 10) ? rN : 10;
    uint32 beforeTriplets[30]; // 10 triplets max
    for (uint32 i = 0; i < maxDumpTriplets * rguStride && i < 30; i++) {
        beforeTriplets[i] = rGU[i];
    }
    fprintf(stderr, "[COLLAPSE-FORCED] Before qsort - first %u triplets:\n", maxDumpTriplets);
    for (uint32 i = 0; i < maxDumpTriplets; i++) {
        fprintf(stderr, "  [%u] g=%u u=%u r=%u\n", i, 
                rGU[i*rguStride + 0], rGU[i*rguStride + 1], 
                (rguStride == 3) ? rGU[i*rguStride + 2] : 0);
    }
    fflush(stderr);
    
    // Comment out qsort temporarily to test std::sort
    // qsort(rGU,rN,rguStride*sizeof(uint32),funCompareNumbers<uint32>); //sort by gene index
    
    // Test with std::sort using struct view
    struct Triplet {
        uint32 g, u, r;
        Triplet(uint32* base, uint32 stride) : g(base[0]), u(base[1]), r(stride == 3 ? base[2] : 0) {}
        void writeBack(uint32* base, uint32 stride) const {
            base[0] = g;
            base[1] = u;
            if (stride == 3) base[2] = r;
        }
    };
    
    fprintf(stderr, "[COLLAPSE-FORCED] Converting to Triplet vector for std::sort\n");
    fflush(stderr);
    
    std::vector<Triplet> triplets;
    triplets.reserve(rN);
    for (uint32 i = 0; i < rN; i++) {
        triplets.push_back(Triplet(rGU + i*rguStride, rguStride));
    }
    
    fprintf(stderr, "[COLLAPSE-FORCED] Sorting %zu triplets with std::sort\n", triplets.size());
    fflush(stderr);
    
    std::sort(triplets.begin(), triplets.end(), 
              [](const Triplet& a, const Triplet& b) { return a.g < b.g; });
    
    fprintf(stderr, "[COLLAPSE-FORCED] std::sort completed, writing back\n");
    fflush(stderr);
    
    for (uint32 i = 0; i < rN; i++) {
        triplets[i].writeBack(rGU + i*rguStride, rguStride);
    }
    
    fprintf(stderr, "[COLLAPSE-FORCED] std::sort writeback completed\n");
    fflush(stderr);
    
    fprintf(stderr, "[COLLAPSE-FORCED] After qsort, dumping same slice:\n");
    for (uint32 i = 0; i < maxDumpTriplets; i++) {
        fprintf(stderr, "  [%u] g=%u u=%u r=%u", i, 
                rGU[i*rguStride + 0], rGU[i*rguStride + 1], 
                (rguStride == 3) ? rGU[i*rguStride + 2] : 0);
        // Check if values changed unexpectedly
        bool changed = false;
        for (uint32 j = 0; j < rguStride; j++) {
            if (rGU[i*rguStride + j] != beforeTriplets[i*rguStride + j]) {
                changed = true;
                break;
            }
        }
        if (changed) {
            fprintf(stderr, " *** CHANGED ***");
            for (uint32 j = 0; j < rguStride; j++) {
                if (rGU[i*rguStride + j] != beforeTriplets[i*rguStride + j]) {
                    fprintf(stderr, " [%u: %u->%u]", j, beforeTriplets[i*rguStride + j], rGU[i*rguStride + j]);
                }
            }
        }
        fprintf(stderr, "\n");
    }
    fflush(stderr);
    
    fprintf(stderr, "[COLLAPSE-FORCED] After qsort, about to start loop\n");
    fflush(stderr);

    uint32 gid1 = -1;//current gID
    uint32 nGenes = 0, nGenesMult = 0; //number of genes
    
    fprintf(stderr, "[COLLAPSE-FORCED] Loop vars initialized\n");
    fflush(stderr);
    
    uint32 totalWords = rN * rguStride;
    fprintf(stderr, "[COLLAPSE-FORCED] totalWords=%u\n", totalWords);
    fflush(stderr);
    
    // Check pSolo.multiMap before using it
    bool multiMapEnabled = pSolo.multiMap.yes.multi;
    fprintf(stderr, "[COLLAPSE-FORCED] multiMapEnabled=%d\n", multiMapEnabled ? 1 : 0);
    fflush(stderr);
    
    for (uint32 iR=0; iR<totalWords; iR+=rguStride) {
        if (iR == 0) {
            fprintf(stderr, "[COLLAPSE-FORCED] First iteration: iR=0 rGU[0]=%u rGU[1]=%u rGU[2]=%u\n",
                    rGU[0], rGU[1], rGU[2]);
            fflush(stderr);
        }
        if (rGU[iR+rguG]!=gid1) {//record gene boundary
            // Bounds check before writing
            if (nGenes >= gReadS.size()) {
                fprintf(stderr, "[COLLAPSE-FORCED] ERROR in loop: nGenes=%u >= gReadS.size()=%zu at iR=%u\n",
                        nGenes, gReadS.size(), iR);
                fflush(stderr);
                return; // Early return
            }
            gReadS[nGenes]=iR;
            gid1=rGU[iR+rguG];
            gID[nGenes]=gid1;
            
            ++nGenes;            
            if (multiMapEnabled && (gid1 & geneMultMark) )
                ++nGenesMult;
        };
    };
    
    fprintf(stderr, "[COLLAPSE-FORCED] Gene boundary loop complete: nGenes=%u nGenesMult=%u\n", nGenes, nGenesMult);
    fflush(stderr);
    
    // Try to access gReadS more carefully - check if it's valid first
    fprintf(stderr, "[COLLAPSE-FORCED] About to access gReadS vector\n");
    fflush(stderr);
    
    // Store size in a local variable before any other operations
    const size_t gReadSSize = gReadS.size();
    fprintf(stderr, "[COLLAPSE-FORCED] gReadS.size()=%zu nGenes=%u\n", gReadSSize, nGenes);
    fflush(stderr);
    
    if (nGenes >= gReadSSize) {
        fprintf(stderr, "[COLLAPSE-FORCED] ERROR: nGenes=%u >= gReadS.size()=%zu\n", nGenes, gReadSSize);
        fflush(stderr);
        return; // Early return to avoid crash
    }
    
    fprintf(stderr, "[COLLAPSE-FORCED] About to set gReadS[%u]=%u\n", nGenes, rguStride*rN);
    fflush(stderr);
    
    gReadS[nGenes]=rguStride*rN;//so that gReadS[nGenes]-gReadS[nGenes-1] is the number of reads for nGenes, see below in qsort
    
    fprintf(stderr, "[COLLAPSE-FORCED] gReadS[%u] set successfully\n", nGenes);
    fflush(stderr);
    
    nGenes -= nGenesMult;//unique only gene
    
    bool multiGeneUMI = pSolo.umiFiltering.MultiGeneUMI;
#ifdef DEBUG_CB_UB_PARITY
    // For parity instrumentation, skip the MultiGeneUMI path to avoid noisy debug code crashes
    if (parityEnabled) {
        multiGeneUMI = false;
    }
#endif
    
    //unordered_map <uintUMI, unordered_set<uint32>> umiGeneMap;
    unordered_map <uintUMI, unordered_map<uint32,uint32>> umiGeneMapCount, umiGeneMapCount0;
                   //UMI                 //Gene //Count

    if (multiGeneUMI) {
        fprintf(stderr, "[COLLAPSE-FORCED] Entering MultiGeneUMI block, nGenes=%u\n", nGenes);
        fflush(stderr);
        
        if (nGenes >= gReadS.size()) {
            fprintf(stderr, "[COLLAPSE-FORCED] ERROR: nGenes=%u >= gReadS.size()=%zu\n", nGenes, gReadS.size());
            fflush(stderr);
            return;
        }
        
        fprintf(stderr, "[COLLAPSE-FORCED] About to loop: iR from 0 to gReadS[%u]=%u\n", nGenes, gReadS[nGenes]);
        fflush(stderr);
        
        for (uint32 iR=0; iR<gReadS[nGenes]; iR+=rguStride) {
            umiGeneMapCount[rGU[iR+1]][rGU[iR]]++; 
        };

        for (auto &iu : umiGeneMapCount) {//loop over all UMIs
            if (iu.second.size()==1)
                continue;
            uint32 maxu=0;
            for (const auto &ig : iu.second) {//loop over genes for a given UMI
                if (maxu<ig.second)
                    maxu=ig.second; //find gene with maximum count
            };
            if (maxu==1)
                maxu=2;//to kill UMIs with 1 read to one gene, 1 read to another gene
            for (auto &ig : iu.second) {
                if (maxu>ig.second)
                    ig.second=0; //kills Gene with read count *strictly* < maximum count
            };
        };
    };
    
    if (pSolo.umiFiltering.MultiGeneUMI_All) {
        for (uint32 iR=0; iR<gReadS[nGenes]; iR+=rguStride) {
            umiGeneMapCount[rGU[iR+1]][rGU[iR]]++; 
        };

        for (auto &iu : umiGeneMapCount) {//loop over all UMIs
            if (iu.second.size()>1) {
                for (auto &ig : iu.second)
                    ig.second=0; //kill all genes for this UMI
            };
        };
    };    
    
    vector<unordered_map <uintUMI,uintUMI>> umiCorrected(nGenes);

    if (!minimalMode) {
        if (countCellGeneUMI.size() < countCellGeneUMIindex[iCB] + nGenes*countMatStride)
            countCellGeneUMI.resize((countCellGeneUMI.size() + nGenes*countMatStride )*2);//allocated vector too small
        
        nGenePerCB[iCB]=0;
        nUMIperCB[iCB]=0;
        countCellGeneUMIindex[iCB+1]=countCellGeneUMIindex[iCB];
    }
    
    /////////////////////////////////////////////
    /////////// main cycle over genes with unique-gene-mappers
    for (uint32 iG=0; iG<nGenes; iG++) {//collapse UMIs for each gene
        uint32 *rGU1=rGU+gReadS[iG];
            
        uint32 nR0 = (gReadS[iG+1]-gReadS[iG])/rguStride; //total number of reads
        if (nR0==0)
            continue; //no reads - this should not happen?
            
        qsort(rGU1, nR0, rguStride*sizeof(uint32), 
              (int (*)(const void*, const void*))funCompareTypeShift<uint32,rguU>);
            
        //exact collapse
        uint32 iR1=-umiArrayStride; //number of distinct UMIs for this gene
        uint32 u1=-1;
        for (uint32 iR=rguU; iR<gReadS[iG+1]-gReadS[iG]; iR+=rguStride) {//count and collapse identical UMIs
            if (pSolo.umiFiltering.MultiGeneUMI && umiGeneMapCount[rGU1[iR]][gID[iG]]==0) {//multigene UMI is not recorded
                if ( pSolo.umiDedup.typeMain != UMIdedup::typeI::NoDedup ) //for NoDedup, the UMI filtering is not done
                    rGU1[iR] = (uintUMI) -1; //mark multigene UMI, so that UB tag will be set to -
                continue;
            };            
                
            if (rGU1[iR]!=u1) {
                iR1 += umiArrayStride;
                u1=rGU1[iR];
                umiArray[iR1]=u1;
                umiArray[iR1+1]=0;
            };
            umiArray[iR1+1]++;
            //if ( umiArray[iR1+1]>nRumiMax) nRumiMax=umiArray[iR1+1];
        };

        uint32 nU0 = (iR1+umiArrayStride)/umiArrayStride;//number of UMIs after simple exact collapse
       
        if (pSolo.umiFiltering.MultiGeneUMI_CR) {
            if (nU0==0)
                continue; //nothing to count
                
            for (uint64 iu=0; iu<nU0*umiArrayStride; iu+=umiArrayStride) {
                umiGeneMapCount0[umiArray[iu+0]][iG]+=umiArray[iu+1];//this sums read counts over UMIs that were collapsed
            };
                
            umiArrayCorrect_CR(nU0, umiArray.data(), true, false, umiCorrected[iG]);
                
            for (uint64 iu=0; iu<nU0*umiArrayStride; iu+=umiArrayStride) {//just fill the umiGeneMapCount - will calculate UMI counts later
                umiGeneMapCount[umiArray[iu+2]][iG]+=umiArray[iu+1];//this sums read counts over UMIs that were collapsed
            };
                
            continue; //done with MultiGeneUMI_CR, readInfo will be filled later
        };        
            
            
        if (!minimalMode) {
            if (pSolo.umiDedup.yes.NoDedup)
                countCellGeneUMI[countCellGeneUMIindex[iCB+1] + pSolo.umiDedup.countInd.NoDedup] = nR0;
        }

        if (nU0>0) {
            // Run dedup corrections but write matrices only if not minimal
            if (pSolo.umiDedup.yes.CR)
                umiArrayCorrect_CR(nU0, umiArray.data(), pSolo.umiDedup.typeMain==UMIdedup::typeI::CR, !minimalMode, umiCorrected[iG]);
            if (pSolo.umiDedup.yes.Directional)
                umiArrayCorrect_Directional(nU0, umiArray.data(), pSolo.umiDedup.typeMain==UMIdedup::typeI::Directional, !minimalMode, umiCorrected[iG], 0);
            if (pSolo.umiDedup.yes.Directional_UMItools)
                umiArrayCorrect_Directional(nU0, umiArray.data(), pSolo.umiDedup.typeMain==UMIdedup::typeI::Directional_UMItools, !minimalMode, umiCorrected[iG], -1);
            if (pSolo.umiDedup.yes.All)
                umiArrayCorrect_Graph(nU0, umiArray.data(), pSolo.umiDedup.typeMain==UMIdedup::typeI::All, !minimalMode, umiCorrected[iG]);
            if (!minimalMode) {
                if (pSolo.umiDedup.yes.Exact)
                    countCellGeneUMI[countCellGeneUMIindex[iCB+1] + pSolo.umiDedup.countInd.Exact] = nU0;
            }
        }
        
        if (!minimalMode) {
            uint32 totcount=0;
            for (uint32 ii=countCellGeneUMIindex[iCB+1]+1; ii<countCellGeneUMIindex[iCB+1]+countMatStride; ii++) {
                totcount += countCellGeneUMI[ii];
            };
            if (totcount>0) {
                countCellGeneUMI[countCellGeneUMIindex[iCB+1] + 0] = gID[iG];
                nGenePerCB[iCB]++;
                nUMIperCB[iCB] += countCellGeneUMI[countCellGeneUMIindex[iCB+1] + pSolo.umiDedup.countInd.main];
                countCellGeneUMIindex[iCB+1] = countCellGeneUMIindex[iCB+1] + countMatStride;
            };
        }
        
        {
            for (uint32 iR=0; iR<gReadS[iG+1]-gReadS[iG]; iR+=rguStride) {//cycle over reads
                uint64 iread1 = rGU1[iR+rguR];
                uint32 umi=rGU1[iR+rguU];
                if (umiCorrected[iG].count(umi)>0)
                    umi=umiCorrected[iG][umi]; //correct UMI
#ifdef DEBUG_CB_UB_PARITY
                if (parityEnabled) {
                    ++dbgWriteSiteUnique;
                }
#endif
                recordReadInfo((uint32_t)iread1, indCB[iCB], umi, 1);
#ifdef DEBUG_CB_UB_PARITY
                dbgCheckReadAfterWrite(this, (uint32_t)iread1, indCB[iCB], umi, 1, "unique");
#endif
            };
        }
    };

    if (pSolo.umiFiltering.MultiGeneUMI_CR) {
        
        vector<uint32> geneCounts(nGenes,0);
        
        vector<unordered_set<uintUMI>> geneUmiHash;
        geneUmiHash.resize(nGenes);
        
        for (const auto &iu: umiGeneMapCount) {//loop over UMIs for all genes
                       
            uint32 maxu=0, maxg=-1;
            for (const auto &ig : iu.second) {
                if (ig.second>maxu) {
                    maxu=ig.second;
                    maxg=ig.first;
                } else if (ig.second==maxu) {
                    maxg=-1;
                };
            };

            if ( maxg+1==0 )
                continue; //this umi is not counted for any gene, because two genes have the same read count for this UMI
            
            for (const auto &ig : umiGeneMapCount0[iu.first]) {//check that this umi/gene had also top count for uncorrected umis
                if (ig.second>umiGeneMapCount0[iu.first][maxg]) {
                    maxg=-1;
                    break;
                };
            };

            if ( maxg+1!=0 ) {//this UMI is counted
                geneCounts[maxg]++;
                geneUmiHash[maxg].insert(iu.first);
            };
        };

        if (!minimalMode) {
            for (uint32 ig=0; ig<nGenes; ig++) {
                if (geneCounts[ig] == 0)
                    continue; //no counts for this gene
                nGenePerCB[iCB]++;
                nUMIperCB[iCB] += geneCounts[ig];
                countCellGeneUMI[countCellGeneUMIindex[iCB+1] + 0] = gID[ig];
                countCellGeneUMI[countCellGeneUMIindex[iCB+1] + pSolo.umiDedup.countInd.CR] = geneCounts[ig];
                countCellGeneUMIindex[iCB+1] = countCellGeneUMIindex[iCB+1] + countMatStride;//iCB+1 accumulates the index
            };
        }
        
#ifdef DEBUG_CB_UB_PARITY
        // Snapshot legacy CB/UB before clique correction mutates them
        if (parityEnabled) {
            snapshotLegacyCBUB(iCB, rGU, gReadS, nGenes, gID, rguStride);
        }
#endif
        
        {
            // Apply clique corrections if enabled
            if (pSolo.umiCorrectionMode > 0 && !umiCorrections.empty()) {
                // Get CB16 and sample assignment for this CB
                if (iCB >= indCB.size()) {
                    P.inOut->logMain << "ERROR: iCB=" << iCB << " >= indCB.size()=" << indCB.size() << " in collapseUMIperCB" << endl;
                    ostringstream errOut;
                    errOut << "FATAL: indCB index out of bounds: iCB=" << iCB << " >= " << indCB.size();
                    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                }
                uint32_t wlIdx = indCB[iCB];
                if (wlIdx >= pSolo.cbWLstr.size()) {
                    P.inOut->logMain << "ERROR: wlIdx=" << wlIdx << " >= cbWLstr.size()=" << pSolo.cbWLstr.size() << " in collapseUMIperCB" << endl;
                    ostringstream errOut;
                    errOut << "FATAL: cbWLstr index out of bounds: wlIdx=" << wlIdx << " >= " << pSolo.cbWLstr.size();
                    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                }
                std::string cb16 = pSolo.cbWLstr[wlIdx].substr(0, 16);
                uint8_t sampleIdx = 0;
                auto assignIt = assignmentsMap.find(cb16);
                if (assignIt != assignmentsMap.end()) {
                    AssignmentInfo& info = assignIt->second;
                    if (info.status != "ambiguous" && info.status != "no_call") {
                        sampleIdx = info.sampleIdx;
                    }
                }
                
                for (uint32 iG=0; iG<nGenes; iG++) {
                    uint32 geneIdx = gID[iG];
                    // Build group key
                    uint64_t groupKey = (static_cast<uint64_t>(wlIdx) << 24) | 
                                      (static_cast<uint64_t>(sampleIdx) << 16) | 
                                      static_cast<uint64_t>(geneIdx);
                    
                    // Note: Ledger-based UMI correction path removed
                    // The inline flex path uses hash-based correction in buildHistogramsFromHash/applyCliqueCorrectionsToHash
                    (void)groupKey; // suppress unused warning
                }
            }
        }
        
            for (uint32 iG=0; iG<nGenes; iG++) {//cycle over genes
                uint32 *rGU1=rGU+gReadS[iG];            
                for (uint32 iR=0; iR<gReadS[iG+1]-gReadS[iG]; iR+=rguStride) {//cycle over reads
                    uint64 iread1 = rGU1[iR+rguR];
                    uint32 umi=rGU1[iR+rguU];
                    if (umiCorrected[iG].count(umi)>0)
                        umi=umiCorrected[iG][umi]; //correct UMI
                    if (geneUmiHash[iG].count(umi)>0) {
#ifdef DEBUG_CB_UB_PARITY
                        if (parityEnabled) {
                            ++dbgWriteSiteCliqueGood;
                        }
#endif
                        recordReadInfo((uint32_t)iread1, indCB[iCB], umi, 1);
#ifdef DEBUG_CB_UB_PARITY
                        dbgCheckReadAfterWrite(this, (uint32_t)iread1, indCB[iCB], umi, 1, "clique_good");
#endif
                    } else {
#ifdef DEBUG_CB_UB_PARITY
                        if (parityEnabled) {
                            ++dbgWriteSiteCliqueDrop;
                        }
#endif
                        recordReadInfo((uint32_t)iread1, indCB[iCB], (uint32_t)-1, 2);
#ifdef DEBUG_CB_UB_PARITY
                        dbgCheckReadAfterWrite(this, (uint32_t)iread1, indCB[iCB], (uint32_t)-1, 2, "clique_drop");
#endif
                    };
                };
            };
        }
        
#ifdef DEBUG_CB_UB_PARITY
        // Compare legacy vs inline CB/UB after all recordReadInfo writes for this CB
        if (parityEnabled) {
            compareCBUBParity(iCB, rGU, gReadS, nGenes, gID, rguStride);
        }
#endif
    
    //////////////////////////////////////////multi-gene reads to the end of function
    //////////////////////////////////////////
    if (!minimalMode && pSolo.multiMap.yes.multi)
        countMatMult.i[iCB+1] = countMatMult.i[iCB];
    
    if (nGenesMult>0) {//process multigene reads
        
        {
            for (uint32 iR=gReadS[nGenes]; iR<gReadS[nGenes+nGenesMult]; iR+=rguStride) {//cycle over multi-gene reads to record their CB and UMI, no corrections
                uint64 iread1 = rGU[iR+rguR];
#ifdef DEBUG_CB_UB_PARITY
                if (parityEnabled) {
                    ++dbgWriteSiteMulti;
                }
#endif
                recordReadInfo((uint32_t)iread1, indCB[iCB], rGU[iR+rguU], 1);
#ifdef DEBUG_CB_UB_PARITY
                dbgCheckReadAfterWrite(this, (uint32_t)iread1, indCB[iCB], rGU[iR+rguU], 1, "multi");
#endif
            };
        }
        
        if (!minimalMode) {
            std::vector<vector<uint32>> umiGenes;
            umiGenes.reserve(256);
            {//for each umi, count number of reads per gene. 
             //Output umiGenes: only genes with nReads = nReads-for-this-UMI will be kept for this UMI
                uint32 *rGUm = rGU + gReadS[nGenes];
                uint32 nRm=( gReadS[nGenes+nGenesMult] - gReadS[nGenes] ) / rguStride;
                
                //sort by UMI, then by read, then by gene
                // CRITICAL: funCompare_uint32_1_2_0 reads va[2] (readId), requires rguStride >= 3
                if (rguStride >= 3) {
                    int (*cmp)(const void*, const void*) = funCompare_uint32_1_2_0;
                    qsort(rGUm, nRm, rguStride*sizeof(uint32), cmp);
                } else {
                    // Fallback: sort by UMI only (va[1]) when stride == 2
                    // Use funCompareTypeShift to sort by UMI (offset 1)
                    qsort(rGUm, nRm, rguStride*sizeof(uint32), 
                          (int (*)(const void*, const void*))funCompareTypeShift<uint32,rguU>);
                }
                
                std::unordered_map<uint32, uint32> geneReadCount; //number of reads per gene
                uint32 nRumi=0;
                bool skipUMI=false;
                uintUMI umiPrev = (uintUMI)-1;
                uintRead readPrev = (uintRead)-1;
                for (uint32 iR=0; iR<nRm*rguStride; iR+=rguStride) {//for each umi, find intersection of genes from each read
                    uintUMI umi1 = rGUm[iR+1];
                    if (umi1!=umiPrev) {//starting new UMI
                        umiPrev = umi1;                
                        if (umiGeneMapCount.count(umi1)>0) {
                            skipUMI = true;//this UMI is skipped because it was among uniquely mapped
                        } else {
                            skipUMI = false;//new good umi
                            geneReadCount.clear();
                            nRumi=0;
                            readPrev = (uintRead)-1;
                        };
                    };
                    
                    if (skipUMI)
                        continue; //this UMI is skipped because it was among uniquely mapped
                    
                    uintRead read1 = rGUm[iR+2];
                    if (read1 != readPrev) {
                        ++nRumi;
                        readPrev = read1;
                    };
                    
                    uint32 g1 = rGUm[iR+0] ^ geneMultMark; //XOR to unset the geneMultMark bit
                    geneReadCount[g1]++;
                    
                    if (iR == nRm*rguStride-rguStride || umi1 != rGUm[iR+1+rguStride]) {//record this umi
                        uint32 ng=0;
                        for (const auto &gg: geneReadCount) {
                            if (gg.second == nRumi)
                                ++ng;
                        };
                        vector<uint32> vg;
                        vg.reserve(ng);//this and above is to construct vector of precise size, for efficiency?
                        for (const auto &gg: geneReadCount) {
                            if (gg.second == nRumi)
                                vg.push_back(gg.first);
                        };                
                        umiGenes.push_back(vg);
                    };
                };
            };
            
            std::map<uint32,uint32> genesM; //genes to quantify
            
            {//collect all genes, replace geneID with index in umiGenes
                uint32 ng = 0;
                for (auto &uu: umiGenes) {
                    for (auto &gg: uu) {
                        if (genesM.count(gg) == 0) {//new gene
                            genesM[gg]=ng;
                            ++ng;
                        };
                        gg = genesM[gg];
                    };
                };
            };
        
       
        vector<double> gEuniform(genesM.size(), 0);
        {//gEuniform=uniformly distribute multigene UMIs
            for (auto &ug: umiGenes) {
                for (auto &gg: ug) {
                    gEuniform[gg] += 1.0 / double(ug.size()); // 1/n_genes_umi
                };
            };
        };

            vector<vector<double>> gErescue(pSolo.umiDedup.yes.N);
            // Legacy multi-map "Rescue" mode is only safe when the full count matrices are present.
            // For Flex inline paths (umiCorrectionMode>0) or skip/minimal runs we gate it off.
            bool rescueEnabled = pSolo.multiMap.yes.Rescue && (pSolo.umiCorrectionMode == 0);
        if (!rescueEnabled && pSolo.multiMap.yes.Rescue) {
            static bool warnedRescue = false;
            if (!warnedRescue) {
                P.inOut->logMain << "NOTE: multi-map rescue mode is disabled for Flex/inline clique runs and when minimalMode is active.\n";
                warnedRescue = true;
            }
        }
        if (rescueEnabled) {
            
            for (uint32 indDedup=0; indDedup < pSolo.umiDedup.yes.N; indDedup++) {
                vector<double> gEu(genesM.size(), 0);
                {//collect unique gene counts
                    for (uint32 igm=countCellGeneUMIindex[iCB]; igm<countCellGeneUMIindex[iCB+1]; igm+=countMatStride) {
                        uint32 g1 = countCellGeneUMI[igm];
                        if (genesM.count(g1)>0)
                            gEu[genesM[g1]]=(double)countCellGeneUMI[igm+1+indDedup];
                    };
                };
                
                gErescue[indDedup].resize(genesM.size(), 0);
                {//gErescue=distribute UMI proportionally to gEuniform+gEu
                    for (auto &ug: umiGenes) {
                        double norm1 = 0.0;
                        for (auto &gg: ug)
                            norm1 += gEuniform[gg]+gEu[gg];
                        
                        if (norm1==0.0)
                            continue; //this should not happen since gEuniform is non-zero for all genes involved
                        norm1 = 1.0 / norm1;
                        
                        for (auto &gg: ug) {
                            gErescue[indDedup][gg] += (gEuniform[gg]+gEu[gg])*norm1;
                        };
                    };
                };
            };
        };
        
            vector<vector<double>> gEpropUnique(pSolo.umiDedup.yes.N);
            if (pSolo.multiMap.yes.PropUnique) {
            
            for (uint32 indDedup=0; indDedup < pSolo.umiDedup.yes.N; indDedup++) {
                vector<double> gEu(genesM.size(), 0);
                {//collect unique gene counts
                    for (uint32 igm=countCellGeneUMIindex[iCB]; igm<countCellGeneUMIindex[iCB+1]; igm+=countMatStride) {
                        uint32 g1 = countCellGeneUMI[igm];
                        if (genesM.count(g1)>0)
                            gEu[genesM[g1]]=(double)countCellGeneUMI[igm+1+indDedup];
                    };
                };
                
                gEpropUnique[indDedup].resize(genesM.size(), 0);
                {//gErescue=distribute UMI proportionally to gEuniform+gEu
                    for (auto &ug: umiGenes) {
                        double norm1 = 0.0;
                        for (auto &gg: ug)
                            norm1 += gEu[gg];
                        
                        if (norm1==0.0) {//this UMI has no genes with unique mappers - distribute it uniformly
                            for (auto &gg: ug)
                                gEpropUnique[indDedup][gg] += 1.0 / double(ug.size());
                        } else {//this UMI has genes with unique mappers - distribute it proportionally to unique mappers
                            norm1 = 1.0 / norm1;
                            for (auto &gg: ug)
                                gEpropUnique[indDedup][gg] += gEu[gg]*norm1;
                        };
                    };
                };
            };
        };
        
            vector<vector<double>> gEem(pSolo.umiDedup.yes.N);
            if (pSolo.multiMap.yes.EM) {
                
            for (uint32 indDedup=0; indDedup < pSolo.umiDedup.yes.N; indDedup++) {
                vector<double> gEu(genesM.size(), 0);
                {//collect unique gene counts
                    for (uint32 igm=countCellGeneUMIindex[iCB]; igm<countCellGeneUMIindex[iCB+1]; igm+=countMatStride) {
                        uint32 g1 = countCellGeneUMI[igm];
                        if (genesM.count(g1)>0)
                            gEu[genesM[g1]]=(double)countCellGeneUMI[igm+1+indDedup];
                    };
                };
                
                {//gEem = EM
                    
                    vector<double> gEM1 = gEuniform;
                    for (uint32 ii=0; ii<gEM1.size(); ii++)
                        gEM1[ii] += gEu[ii]; //start with sum of unique and uniform
                        
                    vector<double> gEM2(genesM.size(), 0);
                        
                    auto *gEM1p=&gEM1;
                    auto *gEM2p=&gEM2;
                                                
                    double maxAbsChange=1;
                    uint32 iterI=0;
                    while(true) {
                        ++iterI;
                        
                        auto &gEMold=*gEM1p; //for convenience - to use instead of pointer
                        auto &gEMnew=*gEM2p;
                        
                        std::copy(gEu.begin(), gEu.end(), gEMnew.begin());//gEMnew is initialized with unique counts
                        
                        for (auto &gg: gEMold) {//zero-out very small counts
                            if (gg<0.01) //hardcoded
                                gg=0;
                        };
                        
                        for (auto &ug: umiGenes) {
                            double norm1 = 0.0;
                            for (auto &gg: ug) //cycle over genes for this umi
                                norm1 += gEMold[gg];
                                
                            norm1 = 1.0 / norm1;
                                
                            for (auto &gg: ug) {
                                gEMnew[gg] += gEMold[gg]*norm1;
                            };
                        };
                        
                        maxAbsChange=0.0;
                        for (uint32 ii=0; ii<gEMnew.size(); ii++) {
                            double change1 = abs(gEMnew[ii]-gEMold[ii]);
                            if (change1 > maxAbsChange)
                                maxAbsChange = change1;
                        };
                        
                        if (maxAbsChange < 0.01 || iterI>100) {//hardcoded
                            gEem[indDedup] = gEMnew;
                            break;
                        };
                        
                        swap(gEM1p, gEM2p); //swap new and old for the next interation of EM
                    };
                        
                    for (uint32 ii=0; ii<gEM1.size(); ii++)
                        gEem[indDedup][ii] -= gEu[ii]; //gEem contains only multimapper counts
                };
                
            };
        };        
        
        #ifndef DEBUG_CB_UB_PARITY
            {//write to countMatMult (only in full mode)
            if (countMatMult.m.size() < countMatMult.i[iCB+1] + genesM.size()*countMatMult.s*pSolo.umiDedup.yes.N + 100) //+100 just in case
                countMatMult.m.resize((countMatMult.i[iCB+1] + genesM.size()*countMatMult.s*pSolo.umiDedup.yes.N + 100)*2);

            for (const auto &gm: genesM) {
                countMatMult.m[countMatMult.i[iCB+1] + 0] = gm.first;
                    
                for (uint32 indDedup=0; indDedup < pSolo.umiDedup.yes.N; indDedup++) {
                    uint32 ind1 = countMatMult.i[iCB+1] + indDedup;
                    
                    if (pSolo.multiMap.yes.Uniform)
                        countMatMult.m[ind1 + pSolo.multiMap.countInd.Uniform] = gEuniform[gm.second];
                        
                    if (pSolo.multiMap.yes.Rescue)
                        countMatMult.m[ind1 + pSolo.multiMap.countInd.Rescue] = gErescue[indDedup][gm.second];
                        
                    if (pSolo.multiMap.yes.PropUnique)
                        countMatMult.m[ind1 + pSolo.multiMap.countInd.PropUnique] = gEpropUnique[indDedup][gm.second];
                        
                    if (pSolo.multiMap.yes.EM)
                        countMatMult.m[ind1 + pSolo.multiMap.countInd.EM] = gEem[indDedup][gm.second];                    
                    
                    countMatMult.i[iCB+1] += countMatMult.s;
                };
            };
        } // write to countMatMult
        #endif
        } // !minimalMode block for multigene rescues
    };

};

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// sorting functions
inline int funCompareSolo1 (const void *a, const void *b) {
    uint32 *va= (uint32*) a;
    uint32 *vb= (uint32*) b;

    if (va[1]>vb[1]) {
        return 1;
    } else if (va[1]<vb[1]) {
        return -1;
    } else if (va[0]>vb[0]){
        return 1;
    } else if (va[0]<vb[0]){
        return -1;
    } else {
        return 0;
    };
};

inline int funCompare_uint32_1_2_0 (const void *a, const void *b) {
    uint32 *va= (uint32*) a;
    uint32 *vb= (uint32*) b;

    if (va[1]>vb[1]) {
        return 1;
    } else if (va[1]<vb[1]) {
        return -1;
    } else if (va[2]>vb[2]){
        return 1;
    } else if (va[2]<vb[2]){
        return -1;
    } else if (va[0]>vb[0]){
        return 1;
    } else if (va[0]<vb[0]){
        return -1;
    } else {
        return 0;
    };
};


////////////////////////////////////////////////////////////////////////////////////////////////
uint32 SoloFeature::umiArrayCorrect_CR(const uint32 nU0, uintUMI *umiArr, const bool readInfoRec, const bool nUMIyes, unordered_map <uintUMI,uintUMI> &umiCorr)
{
    int (*cmpSolo1)(const void*, const void*) = funCompareSolo1;
    qsort(umiArr, nU0, umiArrayStride*sizeof(uint32), cmpSolo1);
    
    for (uint64 iu=0; iu<nU0*umiArrayStride; iu+=umiArrayStride) {
        
        umiArr[iu+2] = umiArr[iu+0]; //stores corrected UMI for 1MM_CR and 1MM_Directional
        for (uint64 iuu=(nU0-1)*umiArrayStride; iuu>iu; iuu-=umiArrayStride) {

            uint32 uuXor = umiArr[iu+0] ^ umiArr[iuu+0];

            if ( (uuXor >> (__builtin_ctz(uuXor)/2)*2) <= 3 ) {//1MM                 
                umiArr[iu+2]=umiArr[iuu+0];//replace iu with iuu
                break;
            };
        };
    };
    
    if (readInfoRec) {//record corrections
        for (uint64 iu=0; iu<nU0*umiArrayStride; iu+=umiArrayStride) {
            if (umiArr[iu+0] != umiArr[iu+2])
                umiCorr[umiArr[iu+0]]=umiArr[iu+2];
        };
    };
    
    if (!nUMIyes) {
        return 0;
    } else {
        unordered_set<uintUMI> umiC;
        for (uint64 iu=0; iu<nU0*umiArrayStride; iu+=umiArrayStride) {
            umiC.insert(umiArr[iu+2]);
        };
       return umiC.size();
    };
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
uint32 SoloFeature::umiArrayCorrect_Directional(const uint32 nU0, uintUMI *umiArr, const bool readInfoRec, const bool nUMIyes, unordered_map <uintUMI,uintUMI> &umiCorr, const int32 dirCountAdd)
{
    int (*cmpReverseShift)(const void*, const void*) = funCompareNumbersReverseShift<uint32, 1>;
    qsort(umiArr, nU0, umiArrayStride*sizeof(uint32), cmpReverseShift);//TODO no need to sort by sequence here, only by count. 
    
    for (uint64 iu=0; iu<nU0*umiArrayStride; iu+=umiArrayStride)
        umiArr[iu+2] = umiArr[iu+0]; //initialized - it will store corrected UMI for 1MM_CR and 1MM_Directional

    uint32 nU1 = nU0;
    for (uint64 iu=umiArrayStride; iu<nU0*umiArrayStride; iu+=umiArrayStride) {
        
        for (uint64 iuu=0; iuu<iu; iuu+=umiArrayStride) {

            uint32 uuXor = umiArr[iu+0] ^ umiArr[iuu+0];

            if ( (uuXor >> (__builtin_ctz(uuXor)/2)*2) <= 3 && umiArr[iuu+1] >= (2*umiArr[iu+1]+dirCountAdd) ) {//1MM && directional condition
                umiArr[iu+2]=umiArr[iuu+2];//replace iuu with iu-corrected
                nU1--;
                break;
            };
        };
    };
    
    if (readInfoRec) {//record corrections
        for (uint64 iu=0; iu<nU0*umiArrayStride; iu+=umiArrayStride) {
            if (umiArr[iu+0] != umiArr[iu+2])
                umiCorr[umiArr[iu+0]]=umiArr[iu+2];
        };
    };
    
    if (!nUMIyes) {
        return 0;
    } else {
        unordered_set<uintUMI> umiC;
        for (uint64 iu=0; iu<nU0*umiArrayStride; iu+=umiArrayStride) {
            umiC.insert(umiArr[iu+2]);
        };
        if (umiC.size()!=nU1)
            cout << nU1 <<" "<< umiC.size()<<endl;
        return umiC.size();
    };
};

void SoloFeature::finalizeTagTableFromReadInfo() {
    // This method is called after collapseUMIall() has finalized all CB/UB values
    // No longer needs to call finalizeTag - the writeTagBinary method will derive
    // CB/UB/status directly from readInfo when needed
};

#ifdef DEBUG_CB_UB_PARITY
static inline void trimTrailingDash(std::string &s) {
    while (!s.empty() && s.back()=='-') {
        s.pop_back();
    }
}

#if 0
void SoloFeature::snapshotLegacyCBUB(uint32_t iCB, uint32_t *rGU, vector<uint32> &gReadS, uint32_t nGenes, vector<uint32> &gID, uint32_t rguStride) {
    // Clear previous snapshot for this CB
    legacyCBUB.clear();
    
    // Iterate through all reads for this CB, extract CB/UB before clique correction
    for (uint32 iG=0; iG<nGenes; iG++) {
        uint32 *rGU1=rGU+gReadS[iG];
        uint32 nR0 = (gReadS[iG+1]-gReadS[iG])/rguStride;
        
        for (uint32 iR=0; iR<nR0; iR++) {
            uint64 iread1 = rGU1[iR+rguStride*rguR];

            // Safety: skip if readId is out of range; defensive against malformed gReadS/count indices
            if (iread1 >= packedReadInfo.data.size()) {
                continue;
            }
            
            // Extract CB/UB from packed read info
            uint32_t cbIdx = getPackedCB((uint32_t)iread1);
            uint32_t umiPacked = getPackedUMI((uint32_t)iread1);
            uint8_t status = getPackedStatus((uint32_t)iread1);
            
            LegacyCBUB entry;
            entry.iCB = iCB;
            entry.iRead = (uint32_t)iread1;
            entry.geneIdx = gID[iG];
            
            // Decode CB
            if (status == 1 && cbIdx < pSolo.cbWLstr.size()) {
                entry.cb = pSolo.cbWLstr[cbIdx];
            } else {
                entry.cb = "-";
            }
            
            // Decode UMI
            if (status == 1) {
                entry.ub.reserve(pSolo.umiL);
                uint32_t tmp = umiPacked;
                for (int i = pSolo.umiL - 1; i >= 0; --i) {
                    uint8_t b = tmp & 0x3; tmp >>= 2;
                    char c = 'N';
                    switch (b) { case 0: c='A'; break; case 1: c='C'; break; case 2: c='G'; break; case 3: c='T'; }
                    entry.ub.push_back(c);
                }
                std::reverse(entry.ub.begin(), entry.ub.end());
                trimTrailingDash(entry.ub);
            } else {
                entry.ub = "-";
            }
            
            legacyCBUB.push_back(entry);
            parityReadsLegacy++;
        }
    }
}
#endif

// Replacement parity helpers using legacyCBUBByRead (keyed by readId)
void SoloFeature::snapshotLegacyCBUB(uint32_t iCB, uint32_t *rGU, vector<uint32> &gReadS, uint32_t nGenes, vector<uint32> &gID, uint32_t rguStride) {
    legacyCBUBByRead.clear();

    // Require readIndex in rGU for keying
    if (rguStride < 3) {
        return;
    }

    auto decodeUmi = [&](uint32_t packed) {
        std::string ub;
        ub.reserve(pSolo.umiL);
        uint32_t tmp = packed;
        for (int i = pSolo.umiL - 1; i >= 0; --i) {
            uint8_t b = tmp & 0x3; tmp >>= 2;
            char c = 'N';
            switch (b) { case 0: c='A'; break; case 1: c='C'; break; case 2: c='G'; break; case 3: c='T'; }
            ub.push_back(c);
        }
        std::reverse(ub.begin(), ub.end());
        trimTrailingDash(ub);
        return ub;
    };

    for (uint32 iG=0; iG<nGenes; iG++) {
        uint32 *rGU1=rGU+gReadS[iG];
        uint32 nR0 = (gReadS[iG+1]-gReadS[iG])/rguStride;

        for (uint32 iR=0; iR<nR0; iR++) {
            uint64 iread1 = rGU1[iR+rguStride*rguR];

            if (iread1 >= nReadsInput) {
                continue;
            }

            LegacyCBUB entry;
            entry.iCB = iCB;
            entry.iRead = (uint32_t)iread1;
            entry.geneIdx = gID[iG];

            if (iCB < indCB.size()) {
                uint32_t wlIdx = indCB[iCB];
                if (wlIdx < pSolo.cbWLstr.size()) {
                    entry.cb = pSolo.cbWLstr[wlIdx];
                } else {
                    entry.cb = "-";
                }
            } else {
                entry.cb = "-";
            }

            entry.ub = decodeUmi(rGU1[iR+rguU]);

            legacyCBUBByRead[entry.iRead] = entry;
            parityReadsLegacy++;

#ifdef DEBUG_CB_UB_PARITY
            if (!g_traceReads.empty() && g_traceReads.count(entry.iRead)) {
                fprintf(stderr, "[TRACE legacy] read=%u cb=%s ub=%s gene=%u iCB=%u\n",
                        entry.iRead, entry.cb.c_str(), entry.ub.c_str(), entry.geneIdx, entry.iCB);
            }
#endif
        }
    }
}

void SoloFeature::compareCBUBParity(uint32_t iCB, uint32_t *rGU, vector<uint32> &gReadS, uint32_t nGenes, vector<uint32> &gID, uint32_t rguStride) {
    // Get CB16 and sample assignment for clique size lookup
    uint32_t wlIdx = (iCB < indCB.size()) ? indCB[iCB] : 0;
    uint8_t sampleIdx = 0;
    if (wlIdx < pSolo.cbWLstr.size()) {
        std::string cb16 = pSolo.cbWLstr[wlIdx].substr(0, 16);
        auto assignIt = assignmentsMap.find(cb16);
        if (assignIt != assignmentsMap.end() && assignIt->second.status != "ambiguous" && assignIt->second.status != "no_call") {
            sampleIdx = assignIt->second.sampleIdx;
        }
    }

    auto decodeUmi = [&](uint32_t packed) {
        std::string ub;
        ub.reserve(pSolo.umiL);
        uint32_t tmp = packed;
        for (int i = pSolo.umiL - 1; i >= 0; --i) {
            uint8_t b = tmp & 0x3; tmp >>= 2;
            char c = 'N';
            switch (b) { case 0: c='A'; break; case 1: c='C'; break; case 2: c='G'; break; case 3: c='T'; }
            ub.push_back(c);
        }
        std::reverse(ub.begin(), ub.end());
        trimTrailingDash(ub);
        return ub;
    };

    for (uint32 iG=0; iG<nGenes; iG++) {
        uint32 *rGU1=rGU+gReadS[iG];
        uint32 nR0 = (gReadS[iG+1]-gReadS[iG])/rguStride;
        uint32 geneIdx = gID[iG];

        uint64_t groupKey = (static_cast<uint64_t>(wlIdx) << 24) |
                            (static_cast<uint64_t>(sampleIdx) << 16) |
                            static_cast<uint64_t>(geneIdx);
        uint32_t cliqueSize = 0;
        auto groupIt = umiCorrections.find(groupKey);
        if (groupIt != umiCorrections.end()) {
            cliqueSize = groupIt->second.size();
        }

        for (uint32 iR=0; iR<nR0; iR++) {
            uint64 iread1 = rGU1[iR+rguStride*rguR];

#ifdef DEBUG_CB_UB_PARITY
            if (!g_traceReads.empty() && g_traceReads.count((uint32_t)iread1)) {
                uint32_t dbgCb = getPackedCB((uint32_t)iread1);
                uint32_t dbgUmi = getPackedUMI((uint32_t)iread1);
                uint8_t dbgStatus = getPackedStatus((uint32_t)iread1);
                fprintf(stderr, "[TRACE packed_before_decode] read=%u cb=%u umi=%u status=%u wlIdx=%u iCB=%u gene=%u rguStride=%u\n",
                        (uint32_t)iread1, dbgCb, dbgUmi, dbgStatus, wlIdx, iCB, geneIdx, rguStride);
            }
#endif

            if (iread1 >= packedReadInfo.data.size()) {
                continue;
            }

            uint32_t cbIdx = getPackedCB((uint32_t)iread1);
            uint32_t umiPacked = getPackedUMI((uint32_t)iread1);
            uint8_t status = getPackedStatus((uint32_t)iread1);
#ifdef DEBUG_CB_UB_PARITY
            if (parityEnabled && status==0) {
                dbgReadStatus0++;
            }
#endif

            std::string inlineCB = "-";
            std::string inlineUB = "-";

            if (status == 1 && cbIdx < pSolo.cbWLstr.size()) {
                inlineCB = pSolo.cbWLstr[cbIdx];
            }
            if (status == 1) {
                inlineUB = decodeUmi(umiPacked);
            }

            parityReadsInline++;

            auto legacyIt = legacyCBUBByRead.find((uint32_t)iread1);
            if (legacyIt == legacyCBUBByRead.end()) {
                parityMismatches++;
                if (mismatchLog.is_open()) {
                    mismatchLog << iCB << "\t" << iread1 << "\t"
                                << "-" << "\t" << inlineCB << "\t"
                                << "-" << "\t" << inlineUB << "\t"
                                << "-1"
                                << "\t" << geneIdx << "\t"
                                << inlineUB << "\t"
                                << cliqueSize << "\t"
                                << (uint32_t)status << "\t" << cbIdx << "\t"
                                << "cb_added" << "\n";
                }
                continue;
            }

            LegacyCBUB &legacy = legacyIt->second;
            legacy.matched = true;

            bool cbMatch = (legacy.cb == inlineCB);
            bool ubMatch = (legacy.ub == inlineUB);
            bool geneMatch = (legacy.geneIdx == geneIdx);

#ifdef DEBUG_CB_UB_PARITY
            if (!g_traceReads.empty() && g_traceReads.count(legacy.iRead)) {
                fprintf(stderr, "[TRACE inline] read=%u cb_legacy=%s cb_inline=%s ub_legacy=%s ub_inline=%s gene_legacy=%u gene_inline=%u clique=%u reason=%s\n",
                        legacy.iRead, legacy.cb.c_str(), inlineCB.c_str(),
                        legacy.ub.c_str(), inlineUB.c_str(),
                        legacy.geneIdx, geneIdx, cliqueSize,
                        (cbMatch && ubMatch && geneMatch) ? "match" :
                        (!cbMatch ? "cb_changed" : (!ubMatch ? "clique_corrected" : "gene_changed")));
            }
#endif

            if (cbMatch && ubMatch && geneMatch) {
                parityMatches++;
            } else {
                parityMismatches++;
                std::string reason = "unknown";
                if (!cbMatch) {
                    reason = "cb_changed";
                } else if (!ubMatch) {
                    reason = "clique_corrected";
                } else if (!geneMatch) {
                    reason = "gene_changed";
                }
                std::string umiStr = (legacy.ub != "-") ? legacy.ub : inlineUB;
                trimTrailingDash(umiStr);

                if (mismatchLog.is_open()) {
                    mismatchLog << iCB << "\t" << legacy.iRead << "\t"
                                << legacy.cb << "\t" << inlineCB << "\t"
                                << legacy.ub << "\t" << inlineUB << "\t"
                                << legacy.geneIdx << "\t" << geneIdx << "\t"
                                << umiStr << "\t"
                                << cliqueSize << "\t"
                                << (uint32_t)status << "\t" << cbIdx << "\t"
                                << reason << "\n";
                }
            }
        }
    }

    // Any legacy reads that were never matched → dropped
    for (const auto &kv : legacyCBUBByRead) {
        const LegacyCBUB &legacy = kv.second;
        if (legacy.matched) continue;
        parityMismatches++;
        if (mismatchLog.is_open()) {
            mismatchLog << legacy.iCB << "\t" << legacy.iRead << "\t"
                        << legacy.cb << "\t" << "-" << "\t"
                        << legacy.ub << "\t" << "-" << "\t"
                        << legacy.geneIdx << "\t" << "-1" << "\t"
                        << legacy.ub << "\t"
                        << 0 << "\t" << (uint32_t)0 << "\t" << (uint32_t)0 << "\t" << "cb_dropped" << "\n";
        }
#ifdef DEBUG_CB_UB_PARITY
    if (!g_traceReads.empty() && g_traceReads.count(legacy.iRead)) {
        fprintf(stderr, "[TRACE drop] read=%u cb=%s ub=%s gene=%u iCB=%u\n",
                legacy.iRead, legacy.cb.c_str(), legacy.ub.c_str(), legacy.geneIdx, legacy.iCB);
    }
#endif
#ifdef DEBUG_CB_UB_PARITY
    // Reset matched flag so repeated calls don't carry state across CBs
    for (auto &kv : legacyCBUBByRead) {
        kv.second.matched = false;
    }
#endif

}

}

void SoloFeature::writeParityCounters() {
    if (parityEnabled) {
        P.inOut->logMain << "[CBUB-PARITY] reads_legacy=" << parityReadsLegacy
                         << " reads_inline=" << parityReadsInline
                         << " matches=" << parityMatches
                         << " mismatches=" << parityMismatches << endl;
    }
}
#endif
