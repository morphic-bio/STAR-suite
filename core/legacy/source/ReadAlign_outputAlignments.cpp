#include "ReadAlign.h"
#include "SlamQuant.h"
#include "SlamCompat.h"
#include "SampleDetector.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"
#include "SequenceFuns.h"
#include "UmiCodec.h"
#include "solo/CbCorrector.h"
#include "TranscriptQuantEC.h"
#include <atomic>
#include <mutex>
#include <cstring>
#include <cstdlib>
#include <unordered_map>
#include <limits>
#include <fstream>
#include <chrono>

// Forward declaration for Flex-only reject logging hook (defined in flex/SoloReadFeature_record_flex.cpp)
extern "C" void storeQnameMapping(uint64_t iRead, const char* qname);

namespace {
inline const char* sanitizeQname(char** readNameMates, uint32_t readNmates, uint32_t mateIdx) {
    if (mateIdx >= readNmates) {
        return nullptr;
    }
    const char* qptr = readNameMates[mateIdx];
    if (qptr && qptr[0] == '@') {
        ++qptr;
    }
    return qptr;
}

inline void writeNormalizedQname(std::ostream& out, const char* raw) {
    if (raw == nullptr) {
        return;
    }
    const char* start = raw;
    size_t len = std::strlen(raw);
    if (len > 0 && start[0] == '@') {
        ++start;
        --len;
    }
    size_t end = 0;
    while (end < len && start[end] != ' ' && start[end] != '\t' &&
           start[end] != '\n' && start[end] != '\r') {
        ++end;
    }
    len = end;
    if (len >= 2 && start[len - 2] == '/' && (start[len - 1] == '1' || start[len - 1] == '2')) {
        len -= 2;
    }
    if (len > 0) {
        out.write(start, static_cast<std::streamsize>(len));
    }
    out.put('\n');
}

static const bool g_debugAmbigResolve =
    (std::getenv("STAR_DEBUG_AMBIG_CB_RESOLVE") != nullptr);
static std::atomic<uint64_t> g_debugAmbigCreateCount{0};
static const uint64_t kMaxAmbigCreateLogs = 50;

static const bool g_dumpAmbigReadNames =
    (std::getenv("STAR_DUMP_AMBIG_READNAMES") != nullptr);
static std::once_flag g_ambigReadNamesInit;
static std::mutex g_ambigReadNamesMutex;
static std::ofstream g_ambigReadNamesStream;
static std::string g_ambigReadNamesPath;

static void initAmbigReadNamesStream(const Parameters &P) {
    const char* envPath = std::getenv("STAR_AMBIG_READNAMES_PATH");
    if (envPath != nullptr && envPath[0] != '\0') {
        g_ambigReadNamesPath = envPath;
    } else {
        g_ambigReadNamesPath = P.outFileNamePrefix + "ambig_cb_readnames.txt";
    }
    g_ambigReadNamesStream.open(g_ambigReadNamesPath.c_str(), std::ios::out | std::ios::trunc);
    g_ambigReadNamesStream.setf(std::ios::unitbuf);
}

static void logAmbigReadName(const Parameters &P, const char* rawName) {
    if (!g_dumpAmbigReadNames || rawName == nullptr) {
        return;
    }
    std::call_once(g_ambigReadNamesInit, [&]() { initAmbigReadNamesStream(P); });
    if (!g_ambigReadNamesStream.good()) {
        return;
    }
    std::lock_guard<std::mutex> lock(g_ambigReadNamesMutex);
    writeNormalizedQname(g_ambigReadNamesStream, rawName);
}

enum class CrRescueRegion : uint8_t { Intergenic = 0, Intronic = 1, Exonic = 2 };

struct CrRescueDecision {
    bool rescued = false;
    bool intronicFallback = false;
    uint64_t winnerAlignIndex = 0;
    uint64_t exonicCount = 0;
    uint64_t intronicCount = 0;
};

static int64_t findLastStartLE(const uint* starts, uint64_t n, uint64_t key) {
    if (starts == nullptr || n == 0) {
        return -1;
    }
    uint64_t lo = 0;
    uint64_t hi = n;
    while (lo < hi) {
        const uint64_t mid = lo + (hi - lo) / 2;
        if (starts[mid] <= key) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if (lo == 0) {
        return -1;
    }
    return static_cast<int64_t>(lo - 1);
}

static bool isSenseForTranscript(const Transcript& aln, uint8_t trStr, int32 strandType) {
    if (strandType == -1) {
        return true;
    }
    if (trStr == 0) {
        return false;
    }
    const uint32_t readStr = (strandType == 0) ? static_cast<uint32_t>(aln.Str) : (1U - static_cast<uint32_t>(aln.Str));
    return readStr == static_cast<uint32_t>(trStr - 1);
}

static bool isSenseForGene(const Transcript& aln, uint8_t geneStr, int32 strandType) {
    if (strandType == -1) {
        return true;
    }
    if (geneStr == 0) {
        return false;
    }
    const uint32_t readStr = (strandType == 0) ? static_cast<uint32_t>(aln.Str) : (1U - static_cast<uint32_t>(aln.Str));
    return readStr == static_cast<uint32_t>(geneStr - 1);
}

static bool alignmentHasSenseGeneOverlap(const Transcriptome& tr,
                                         const Transcript& aln,
                                         int32 strandType) {
    for (uint32_t ib = 0; ib < aln.nExons; ++ib) {
        if (aln.exons[ib][EX_L] == 0) {
            continue;
        }
        const uint64_t bStart = aln.exons[ib][EX_G];
        const uint64_t bEnd = bStart + aln.exons[ib][EX_L];
        const int64_t gi0 = findLastStartLE(tr.geneFull.s, tr.nGe, bEnd - 1);
        for (int64_t gi = gi0; gi >= 0; --gi) {
            if (tr.geneFull.eMax[gi] < bStart) {
                break;
            }
            if (tr.geneFull.e[gi] < bStart) {
                continue;
            }
            const uint64_t gStart = tr.geneFull.s[gi];
            const uint64_t gEnd = tr.geneFull.e[gi] + 1;
            const uint64_t ov = (std::min(bEnd, gEnd) > std::max(bStart, gStart))
                                    ? (std::min(bEnd, gEnd) - std::max(bStart, gStart))
                                    : 0;
            if (ov == 0) {
                continue;
            }
            if (!isSenseForGene(aln, tr.geneFull.str[gi], strandType)) {
                continue;
            }
            return true;
        }
    }
    return false;
}

static uint64_t alignmentMappedLength(const Transcript& aln) {
    if (aln.mappedLength > 0) {
        return static_cast<uint64_t>(aln.mappedLength);
    }
    uint64_t len = 0;
    for (uint32_t ib = 0; ib < aln.nExons; ++ib) {
        len += aln.exons[ib][EX_L];
    }
    return len;
}

static void alignmentBounds(const Transcript& aln, uint64_t& startOut, uint64_t& endOut) {
    startOut = std::numeric_limits<uint64_t>::max();
    endOut = 0;
    for (uint32_t ib = 0; ib < aln.nExons; ++ib) {
        if (aln.exons[ib][EX_L] == 0) {
            continue;
        }
        const uint64_t bStart = aln.exons[ib][EX_G];
        const uint64_t bEnd = bStart + aln.exons[ib][EX_L] - 1;
        startOut = std::min(startOut, bStart);
        endOut = std::max(endOut, bEnd);
    }
}

static uint64_t overlapBlocksWithTranscriptExons(const Transcript& aln,
                                                 uint64_t trStart,
                                                 uint16_t exN,
                                                 const uint32_t* exSE) {
    if (exSE == nullptr || exN == 0) {
        return 0;
    }
    uint64_t overlap = 0;
    uint32_t ib = 0;
    uint16_t ie = 0;
    while (ib < aln.nExons && ie < exN) {
        const uint64_t bStart = aln.exons[ib][EX_G];
        const uint64_t bEnd = bStart + aln.exons[ib][EX_L];
        const uint64_t eStart = trStart + exSE[2 * ie];
        const uint64_t eEnd = trStart + exSE[2 * ie + 1] + 1;

        if (bEnd <= eStart) {
            ++ib;
            continue;
        }
        if (eEnd <= bStart) {
            ++ie;
            continue;
        }

        overlap += std::min(bEnd, eEnd) - std::max(bStart, eStart);
        if (bEnd <= eEnd) {
            ++ib;
        }
        if (eEnd <= bEnd) {
            ++ie;
        }
    }
    return overlap;
}

static CrRescueRegion classifyAlignmentCrRescue(const Transcriptome& tr,
                                                const Transcript& aln,
                                                int32 strandType) {
    const uint64_t mappedLen = alignmentMappedLength(aln);
    if (mappedLen == 0) {
        return CrRescueRegion::Intergenic;
    }

    uint64_t aStart = 0, aEnd = 0;
    alignmentBounds(aln, aStart, aEnd);
    if (aStart == std::numeric_limits<uint64_t>::max()) {
        return CrRescueRegion::Intergenic;
    }

    bool hasSenseExonic = false;
    const int64_t trIdx0 = findLastStartLE(tr.trS, tr.nTr, aStart);
    for (int64_t trIdx = trIdx0; trIdx >= 0; --trIdx) {
        if (tr.trEmax[trIdx] < aEnd) {
            break;
        }
        if (tr.trE[trIdx] < aEnd) {
            continue;
        }
        if (!isSenseForTranscript(aln, tr.trStr[trIdx], strandType)) {
            continue;
        }
        const uint16_t exN = tr.trExN[trIdx];
        const uint32_t* exSE = tr.exSE + 2 * tr.trExI[trIdx];
        const uint64_t exonOverlap = overlapBlocksWithTranscriptExons(aln, tr.trS[trIdx], exN, exSE);
        if (2 * exonOverlap >= mappedLen) {
            hasSenseExonic = true;
            break;
        }
    }
    if (hasSenseExonic) {
        return CrRescueRegion::Exonic;
    }

    std::unordered_map<uint32_t, uint64_t> senseGeneOverlap;
    senseGeneOverlap.reserve(8);
    for (uint32_t ib = 0; ib < aln.nExons; ++ib) {
        if (aln.exons[ib][EX_L] == 0) {
            continue;
        }
        const uint64_t bStart = aln.exons[ib][EX_G];
        const uint64_t bEnd = bStart + aln.exons[ib][EX_L];
        const int64_t gi0 = findLastStartLE(tr.geneFull.s, tr.nGe, bEnd - 1);
        for (int64_t gi = gi0; gi >= 0; --gi) {
            if (tr.geneFull.eMax[gi] < bStart) {
                break;
            }
            if (tr.geneFull.e[gi] < bStart) {
                continue;
            }
            const uint64_t gStart = tr.geneFull.s[gi];
            const uint64_t gEnd = tr.geneFull.e[gi] + 1;
            const uint64_t ov = (std::min(bEnd, gEnd) > std::max(bStart, gStart))
                                    ? (std::min(bEnd, gEnd) - std::max(bStart, gStart))
                                    : 0;
            if (ov == 0) {
                continue;
            }
            if (!isSenseForGene(aln, tr.geneFull.str[gi], strandType)) {
                continue;
            }
            senseGeneOverlap[tr.geneFull.g[gi]] += ov;
        }
    }

    uint64_t maxSenseGeneOverlap = 0;
    for (const auto& kv : senseGeneOverlap) {
        maxSenseGeneOverlap = std::max(maxSenseGeneOverlap, kv.second);
    }
    if (2 * maxSenseGeneOverlap > mappedLen) {
        return CrRescueRegion::Intronic;
    }
    return CrRescueRegion::Intergenic;
}

static CrRescueDecision evaluateCrRescueDecision(const std::vector<CrRescueRegion>& states,
                                                 bool allowIntronicFallback) {
    CrRescueDecision d;
    for (uint64_t ia = 0; ia < states.size(); ++ia) {
        if (states[ia] == CrRescueRegion::Exonic) {
            ++d.exonicCount;
            d.winnerAlignIndex = ia;
        } else if (states[ia] == CrRescueRegion::Intronic) {
            ++d.intronicCount;
            d.winnerAlignIndex = ia;
        }
    }
    if (d.exonicCount == 1) {
        d.rescued = true;
        d.intronicFallback = false;
        return d;
    }
    if (allowIntronicFallback && d.exonicCount == 0 && d.intronicCount == 1) {
        d.rescued = true;
        d.intronicFallback = true;
        return d;
    }
    d.rescued = false;
    d.intronicFallback = false;
    d.winnerAlignIndex = 0;
    return d;
}
}

void ReadAlign::detectSampleFromRawR2() {
    detectedSampleByte_ = 0xFFu;
    if (!(sampleDetReady_ && readNmates > 0 && readLengthOriginal[0] >= 8)) {
        return;
    }

    const uint32_t seqLen = readLengthOriginal[0];
    const uint32_t packedLen = (seqLen + 1) / 2;
    std::vector<uint8_t> packedSeq(packedLen);
    nuclPackBAM(Read0[0], reinterpret_cast<char*>(packedSeq.data()), seqLen);

    const uint32_t detectedIdx = sampleDet_->detectSampleIndex(packedSeq.data(), seqLen, false);
    if (detectedIdx > 0) {
        detectedSampleByte_ = static_cast<uint8_t>(detectedIdx & 0x1Fu);
    }
}

void ReadAlign::outputAlignments() {
  
    outBAMbytes=0;
    readAnnot.reset();
    
    if (mapGen.pGe.gType==101) {//temporary
        ReadAlign::spliceGraphWriteSAM();
        return;
    };    

    ReadAlign::outFilterBySJout();//sets outFilterBySJoutPass=false if read is held for the 2nd stage of outFilterBySJout
    
    if (outFilterBySJoutPass) {//otherwise align is held for the 2nd stage of outFilterBySJout
        ////////////////////////////////////
        if (unmapType<0) {//passed mappedFilter. Unmapped reads can have nTr>0

            crMultiMapRescued_ = false;
            crMultiMapRescuedIntronic_ = false;
            if (P.pSolo.crMultimapRescue && nTr > 1) {
                statsRA.crRescueTotal++;
                const uint64_t nTrBefore = nTr;

                uint64_t geneSenseCount = 0;
                uint64_t geneSenseWinner = 0;
                for (uint64_t ia = 0; ia < nTr; ++ia) {
                    if (trMult[ia] == nullptr) {
                        continue;
                    }
                    if (alignmentHasSenseGeneOverlap(*chunkTr, *trMult[ia], P.pSolo.strand)) {
                        ++geneSenseCount;
                        geneSenseWinner = ia;
                    }
                }

                CrRescueDecision decision;
                if (geneSenseCount == 1) {
                    CrRescueRegion winnerRegion = classifyAlignmentCrRescue(
                        *chunkTr, *trMult[geneSenseWinner], P.pSolo.strand);
                    if (winnerRegion == CrRescueRegion::Exonic) {
                        decision.rescued = true;
                        decision.intronicFallback = false;
                        decision.winnerAlignIndex = geneSenseWinner;
                        statsRA.crRescueGeneVsNonGene++;
                    } else if (winnerRegion == CrRescueRegion::Intronic
                               && P.pSolo.crMultimapRescueIntronic) {
                        decision.rescued = true;
                        decision.intronicFallback = true;
                        decision.winnerAlignIndex = geneSenseWinner;
                        statsRA.crRescueGeneVsNonGene++;
                    } else if (winnerRegion == CrRescueRegion::Intergenic) {
                        statsRA.crRescueFastPathRejected50pct++;
                    } else {
                        statsRA.crRescueFastPathIntronicFallbackOff++;
                    }
                } else {
                    std::vector<CrRescueRegion> states;
                    states.reserve(static_cast<size_t>(nTr));
                    for (uint64_t ia = 0; ia < nTr; ++ia) {
                        if (trMult[ia] == nullptr) {
                            states.push_back(CrRescueRegion::Intergenic);
                            continue;
                        }
                        states.push_back(classifyAlignmentCrRescue(*chunkTr, *trMult[ia], P.pSolo.strand));
                    }
                    decision = evaluateCrRescueDecision(states, P.pSolo.crMultimapRescueIntronic);
                    if (decision.rescued) {
                        if (decision.intronicFallback) {
                            statsRA.crRescueIntronicFallback++;
                        } else {
                            statsRA.crRescueExonicWinner++;
                        }
                    } else {
                        if (decision.exonicCount > 1) {
                            statsRA.crRescueMultiExonicNoRescue++;
                        } else if (decision.intronicCount > 1) {
                            statsRA.crRescueMultiIntronicNoRescue++;
                        } else if (decision.exonicCount == 0 && decision.intronicCount == 1
                                   && !P.pSolo.crMultimapRescueIntronic) {
                            statsRA.crRescueIntronicFallbackOffNoRescue++;
                        } else {
                            statsRA.crRescueAllIntergenicNoRescue++;
                        }
                    }
                }

                if (decision.rescued && decision.winnerAlignIndex < nTrBefore
                    && trMult[decision.winnerAlignIndex] != nullptr) {
                    Transcript* winner = trMult[decision.winnerAlignIndex];
                    for (uint64_t ia = 0; ia < nTrBefore; ++ia) {
                        if (trMult[ia] != nullptr) {
                            trMult[ia]->primaryFlag = (ia == decision.winnerAlignIndex);
                        }
                    }
                    trMult[0] = winner;
                    nTr = 1;
                    trBest = winner;
                    crMultiMapRescued_ = true;
                    crMultiMapRescuedIntronic_ = decision.intronicFallback;

                    if (alignsGenOut.yes && alignsGenOut.alN == nTrBefore
                        && decision.winnerAlignIndex < alignsGenOut.alN
                        && alignsGenOut.alMult[decision.winnerAlignIndex] != nullptr) {
                        Transcript* winnerOut = alignsGenOut.alMult[decision.winnerAlignIndex];
                        for (uint64_t ia = 0; ia < alignsGenOut.alN; ++ia) {
                            if (alignsGenOut.alMult[ia] != nullptr) {
                                alignsGenOut.alMult[ia]->primaryFlag = (ia == decision.winnerAlignIndex);
                            }
                        }
                        alignsGenOut.alMult[0] = winnerOut;
                        alignsGenOut.alN = 1;
                        alignsGenOut.alBest = winnerOut;
                    }
                }
            }

            auto nTr1 = nTr;
            auto trOut1 = trMult[0];

            if (P.pGe.transform.outYes) {
                nTr1 = alignsGenOut.alN;
                trOut1 = alignsGenOut.alMult[0];
            };

            if (nTr1>1) {//multimappers
                statsRA.mappedReadsM++;
                unmapType = -2; //not sure if this used
            } else if (nTr1==1) {//unique mappers
                statsRA.mappedReadsU++;
                statsRA.transcriptStats(*trOut1, Lread);
            };

            if (P.pGe.transform.outSAM && (!P.twoPass.yes || P.twoPass.pass2) ) {//transform genome only on 2nd pass
                ReadAlign::recordSJ(alignsGenOut.alN, alignsGenOut.alMult, chunkOutSJ);
            } else {
                ReadAlign::recordSJ(nTr, trMult, chunkOutSJ); //this will set mateMapped
            };            
            
            ReadAlign::alignedAnnotation();
            
            if (P.quant.slam.yes && slamQuant != nullptr) {
                auto &annFeat = readAnnot.annotFeatures[SoloFeatureTypes::Gene];
                auto &annFeatFull = readAnnot.annotFeatures[SoloFeatureTypes::GeneFull_ExonOverIntron];
                size_t nAlign = annFeat.fAlign.size();
                if (nTr > 0) {
                    size_t nTrSize = static_cast<size_t>(nTr);
                    // Ensure we process exactly nTr alignments
                    // fAlign should be resized to nTr by classifyAlign, but be safe
                    if (nAlign != nTrSize) {
                        if (nAlign < nTrSize) {
                            annFeat.fAlign.resize(nTrSize);
                        } else {
                            nAlign = nTrSize;
                        }
                    }
                    
                    // Track nTr distribution
                    slamQuant->diagnostics().nTrDistribution[nTr]++;
                    
                    // Compute read-level gene count for weighting (union across all alignments)
                    size_t readLevelGeneCount = annFeat.fSet.size();
                    
                    // Count alignments with gene assignments
                    size_t nAlignWithGene = 0;
                    double sumWeight = 0.0;
                    bool trackIntronic = annFeat.fSet.empty() &&
                        annFeatFull.ovType == ReadAnnotFeature::overlapTypes::intronic;
                    
                    // Iterate over all nTr alignments
                    for (size_t ia = 0; ia < nTrSize; ++ia) {
                        // Weight is applied only to alignments that have gene assignments
                        double weight = 1.0;
                        if (P.quant.slam.weightMode == 0) {
                            weight = 1.0 / static_cast<double>(nTr);
                        }
                        
                        // Apply compat overlap weighting if enabled (per alignment, using read-level count)
                        if (slamCompat && slamCompat->cfg().overlapWeight && readLevelGeneCount > 1) {
                            weight = slamCompat->compatOverlapWeight(weight, readLevelGeneCount);
                            slamQuant->diagnostics().compatAlignsOverlapWeightApplied++;
                        }
                        if (ia < annFeat.fAlign.size()) {
                            const auto &exonicGenes = annFeat.fAlign[ia];
                            // Track gene set size distribution
                            slamQuant->diagnostics().geneSetSizeDistribution[exonicGenes.size()]++;
                            if (!exonicGenes.empty()) {
                                nAlignWithGene++;
                                sumWeight += weight;
                                // Apply same weight to all genes in the set (current approach)
                                slamCollect(*trMult[ia], exonicGenes, weight, false);
                            } else if (slamCompat && slamCompat->cfg().intronic) {
                                // No exonic assignment - try compat intronic classification
                                // Candidate genes from GeneFull_ExonOverIntron (gene body overlap)
                                if (ia < annFeatFull.fAlign.size()) {
                                    const auto &candidateGenes = annFeatFull.fAlign[ia];
                                    if (!candidateGenes.empty()) {
                                        std::set<uint32_t> intronicGenes;
                                        if (slamCompat->compatIsIntronic(*trMult[ia], candidateGenes, intronicGenes)) {
                                            slamCollect(*trMult[ia], intronicGenes, weight, true);
                                            slamQuant->diagnostics().compatAlignsReclassifiedIntronic++;
                                        }
                                    }
                                }
                            } else if (trackIntronic && ia < annFeatFull.fAlign.size()) {
                                // Fallback to existing intronic logic (non-compat mode)
                                const auto &intronicGenes = annFeatFull.fAlign[ia];
                                if (!intronicGenes.empty()) {
                                    slamCollect(*trMult[ia], intronicGenes, weight, true);
                                }
                            }
                        }
                    }
                    
                    // Track lenient acceptance count (per-alignment, not per-transcript)
                    for (size_t la = 0; la < annFeat.lenientAcceptByAlign.size(); ++la) {
                        if (annFeat.lenientAcceptByAlign[la]) {
                            slamQuant->diagnostics().compatAlignsLenientAccepted++;
                        }
                    }
                    
                    // Track weight denominator diagnostics
                    slamQuant->diagnostics().nAlignWithGeneDistribution[nAlignWithGene]++;
                    // Bucket sumWeight for distribution (0.0-0.1, 0.1-0.2, ..., 0.9-1.0, 1.0+)
                    size_t weightBucket = static_cast<size_t>(sumWeight * 10.0);
                    if (weightBucket > 10) weightBucket = 10;
                    slamQuant->diagnostics().sumWeightDistribution[weightBucket]++;
                    
                    if (nAlignWithGene == 0) {
                        slamQuant->diagnostics().readsNAlignWithGeneZero++;
                    }
                    if (sumWeight < 0.999) {  // Allow small floating point tolerance
                        slamQuant->diagnostics().readsSumWeightLessThanOne++;
                    }
                } else if (nTr > 0 && nAlign == 0) {
                    slamQuant->diagnostics().readsZeroGenes++;
                    slamQuant->diagnostics().readsNAlignWithGeneZero++;
                }
            }
        };
        
        // Y-chromosome routing decision: follow emitted-alignment semantics.
        // For mapped reads/pairs, route to Y if any emitted alignment touches Y.
        // For non-BAM failure classes (unmapType 0..3), keep reads in noY.
        // This matches the external stepwise path, which derives Y read names from
        // the emitted Y BAM rather than from filtered-out candidate alignments.
        hasYAlignment_ = false;
        if (P.emitNoYBAMyes || P.emitYReadNamesyes || P.emitYNoYFastqyes) {
            // Use transformed genome alignments if available, otherwise use original.
            if (unmapType < 0) {
                uint64 nTrCheck = nTr;
                Transcript **trCheck = trMult;
                if (P.pGe.transform.outSAM && (!P.twoPass.yes || P.twoPass.pass2)) {
                    nTrCheck = alignsGenOut.alN;
                    trCheck = alignsGenOut.alMult;
                }

                // Check all mapped alignments that can contribute to emitted BAM output.
                for (uint iTr = 0; iTr < nTrCheck && !hasYAlignment_; iTr++) {
                    if (trCheck[iTr] == nullptr) continue;
                    if (mapGen.yTids.count(trCheck[iTr]->Chr)) {
                        hasYAlignment_ = true;
                        break;
                    }
                }
            }
        }

        if (P.emitYReadNamesyes && hasYAlignment_) {
            writeNormalizedQname(chunkOutYReadNames, readNameMates[0]);
        }
        
        // Write Y/noY FASTQ output
        if (P.emitYNoYFastqyes) {
            const uint32 yFastqEmitCount = P.yFastqEmitReadCount();
            // Determine if we're in Flex/solo mode (route per read) vs bulk PE (route both mates together)
            bool isSingleCellOrFlex = (P.pSolo.type != 0) || P.pSolo.flexMode;
            bool routePerRead = isSingleCellOrFlex;  // Flex: route each read independently
            
            if (routePerRead) {
                // Flex/solo mode: route each read independently based on its own alignments
                // Note: hasYAlignment_ reflects current read's alignments in Flex mode
                // When scRNA/Flex has a separate barcode read, emit it to the same side
                // as the aligned cDNA read by iterating over all emitted input ends.
                for (uint im = 0; im < yFastqEmitCount; im++) {
                    if (hasYAlignment_) {
                        writeFastxRecord(im, true);  // Write to Y FASTQ
                    } else {
                        writeFastxRecord(im, false); // Write to noY FASTQ
                    }
                }
            } else {
                // Bulk PE mode: route both mates together based on hasYAlignment_ (which reflects both mates)
                // If either mate has Y alignment, both mates go to Y FASTQ
                for (uint im = 0; im < yFastqEmitCount; im++) {
                    if (hasYAlignment_) {
                        writeFastxRecord(im, true);  // Write both mates to Y FASTQ
                    } else {
                        writeFastxRecord(im, false); // Write both mates to noY FASTQ
                    }
                }
            }
        }

        //the operations below are both for mapped and unmapped reads
        soloRead->readBar->getCBandUMI(Read0, Qual0, readLengthOriginal, readNameExtra[0], readFilesIndex, readName);
        
        // Extract CB/UMI/Sample metadata from Solo structures (upstream detection)
        // This avoids re-parsing BAM tags and ensures we use the original FASTQ data
        extractedCbIdxPlus1_ = 0;
        extractedUmi24_ = 0;
        
        // Store CB sequence for Phase 2 resolution (will be looked up in BAMoutput if needed)
        // Note: For ambiguous CBs, the sequence is stored in pendingAmbiguous_ and will be
        // looked up via resolvedCbByKey_ when writing keys
        
        // Use CbCorrector for inline CB correction if available
        SoloReadBarcode *readBar = soloRead->readBar;
        static int debugCount = 0;
        const int MAX_DEBUG_LOGS = 10;
        
        if (P.pSolo.cbCorrector && !readBar->cbSeq.empty()) {
            CbMatch match = P.pSolo.cbCorrector->correct(readBar->cbSeq);
            // Extract UMI first (needed for ambiguous accumulation)
            uint32_t umi24 = 0;
            if (readBar->urValid && readBar->urPacked != UINT32_MAX) {
                umi24 = readBar->urPacked & 0xFFFFFFu;
            } else if (!readBar->umiSeq.empty() && readBar->umiSeq.length() == 12) {
                uint32_t encoded = encodeUMI12(readBar->umiSeq);
                if (encoded != UINT32_MAX) {
                    umi24 = encoded & 0xFFFFFFu;
                }
            }
            
            static const bool disableAmbigResolve =
                (std::getenv("STAR_DISABLE_AMBIG_CB_RESOLVE") != nullptr);

            if (disableAmbigResolve && match.ambiguous) {
                extractedCbIdxPlus1_ = 0;
                extractedCbSeq_ = readBar->cbSeq;
            } else if (match.ambiguous && !match.ambiguousIdx.empty() && umi24 != 0) {
                // Ambiguous CB: accumulate UMI counts for Bayesian resolution
                // This matches process_features: capture raw CB sequence (may contain Ns),
                // quality scores, candidate whitelist indices, and UMI counts
                AmbigKey key = hashCbSeq(readBar->cbSeq);
                auto &entry = pendingAmbiguous_[key];
                bool entryWasEmpty = entry.candidateIdx.empty();
                
                if (entry.candidateIdx.empty()) {
                    // First time seeing this ambiguous CB: initialize candidates
                    // Store raw observed CB sequence (e.g., "ACGNT...") and quality scores
                    entry.candidateIdx = match.ambiguousIdx; // 1-based whitelist indices
                    entry.cbSeq = readBar->cbSeq; // Raw observed sequence (may contain Ns)
                    entry.cbQual = readBar->cbQual; // Phred quality scores (same length as cbSeq)
                    entry.umiCounts.reserve(32);
                    
                    // Validate quality scores match sequence length
                    if (entry.cbQual.length() != entry.cbSeq.length()) {
                        // Pad with default quality if needed (shouldn't happen, but be defensive)
                        if (entry.cbQual.length() < entry.cbSeq.length()) {
                            entry.cbQual.append(entry.cbSeq.length() - entry.cbQual.length(), 'H'); // Q39 default
                        } else {
                            entry.cbQual = entry.cbQual.substr(0, entry.cbSeq.length());
                        }
                    }
                }
                
                // Accumulate UMI count (24-bit packed UMI -> count)
                entry.umiCounts[umi24]++;

                // Debug: dump read name for ambiguous CBs so we can reconstruct a minimal FASTQ set.
                if (g_dumpAmbigReadNames) {
                    const char* rawName = (readNameMates[0] != nullptr) ? readNameMates[0] : readName;
                    logAmbigReadName(P, rawName);
                }

                if (g_debugAmbigResolve) {
                    bool badSeq = entry.cbSeq.empty() || entry.cbQual.empty() ||
                                  (entry.cbSeq.size() != entry.cbQual.size());
                    uint64_t logNo = g_debugAmbigCreateCount.fetch_add(1);
                    if (logNo < kMaxAmbigCreateLogs || badSeq) {
                        if (P.inOut && P.inOut->logMain.good()) {
                            P.inOut->logMain << "[AMBIG-CB-DEBUG] accumulate entry=" << logNo
                                             << " key=0x" << std::hex << key << std::dec
                                             << " first=" << (entryWasEmpty ? 1 : 0)
                                             << " cbSeqLen=" << entry.cbSeq.size()
                                             << " cbQualLen=" << entry.cbQual.size()
                                             << " cbSeq=" << entry.cbSeq
                                             << " cbQual=" << entry.cbQual
                                             << " umi24=0x" << std::hex << umi24 << std::dec
                                             << " candidates=" << entry.candidateIdx.size()
                                             << " umiCounts=" << entry.umiCounts.size()
                                             << " badSeq=" << badSeq
                                             << endl;
                        }
                    }
                }
                
                // Leave as 0 for now; will be resolved after mapping completes
                // Store CB sequence for lookup when writing keys
                extractedCbIdxPlus1_ = 0;
                extractedCbSeq_ = readBar->cbSeq;
            } else if (match.whitelistIdx != 0) {
                // Immediate resolution (non-ambiguous)
                extractedCbIdxPlus1_ = match.whitelistIdx;
                extractedCbSeq_.clear(); // Not needed for resolved CBs
                cbResolutionStats_.resolvedImmediate++;
            } else {
                // No match
                extractedCbIdxPlus1_ = 0;
                extractedCbSeq_.clear();
                cbResolutionStats_.noMatch++;
            }
            
            // Debug logging for first few reads and specific failing CBs
            if (debugCount++ < MAX_DEBUG_LOGS || readBar->cbSeq == "CNGTATTTCGGGCAGT") {
                P.inOut->logMain << "CbCorrector::correct: cbSeq=" << readBar->cbSeq
                                 << ", whitelistIdx=" << match.whitelistIdx
                                 << ", hammingDist=" << (int)match.hammingDist
                                 << ", ambiguous=" << match.ambiguous
                                 << ", extractedCbIdxPlus1_=" << extractedCbIdxPlus1_ << endl;
            }
        } else {
            // Debug logging for why CbCorrector wasn't used
            if (debugCount++ < MAX_DEBUG_LOGS) {
                P.inOut->logMain << "CbCorrector NOT used: cbCorrector=" << (P.pSolo.cbCorrector ? "non-null" : "null")
                                 << ", cbSeq.empty()=" << readBar->cbSeq.empty()
                                 << ", cbSeq=" << (readBar->cbSeq.empty() ? "(empty)" : readBar->cbSeq) << endl;
            }
            // Fallback to original SoloReadBarcode matching logic
            // cbMatch: 0=exact match, 1=one match with 1MM, -1=no match, -3=multiple matches (not allowed)
            if ((readBar->cbMatch == 0 || readBar->cbMatch == 1) && !readBar->cbMatchInd.empty()) {
                // Exact match or single mismatch match: use first match index (0-based) + 1
                extractedCbIdxPlus1_ = static_cast<uint32_t>(readBar->cbMatchInd[0] + 1);
            }
        }
        
        // Extract UMI (24-bit packed) from SoloReadBarcode (if not already extracted above)
        // Also determine UMI validity based on extraction success and umiCheck status
        extractedUmiValid_ = false;  // Default: invalid
        if (readBar) {  // Guard against null readBar
            if (extractedUmi24_ == 0) {
                if (readBar->urValid && readBar->urPacked != UINT32_MAX) {
                    extractedUmi24_ = readBar->urPacked & 0xFFFFFFu;
                    extractedUmiValid_ = true;  // urValid indicates successful extraction
                } else if (!readBar->umiSeq.empty() && readBar->umiSeq.length() == 12) {
                    uint32_t encoded = encodeUMI12(readBar->umiSeq);
                    if (encoded != UINT32_MAX) {
                        extractedUmi24_ = encoded & 0xFFFFFFu;
                        // UMI is valid if umiCheck >= 0 (umiCheck < 0 means N in UMI or other failure)
                        extractedUmiValid_ = (readBar->umiCheck >= 0);
                    }
                }
            } else {
                // UMI was already extracted (e.g., from previous code path)
                // Consider it valid if umiCheck >= 0
                extractedUmiValid_ = (readBar->umiCheck >= 0);
            }
        }
        
        // Debug: log first few extractions
        static int extractCount = 0;
        if (extractCount++ < 5) {
            P.inOut->logMain << "ReadAlign::outputAlignments: extractedCbIdxPlus1_=" << extractedCbIdxPlus1_
                             << ", extractedUmi24_=0x" << std::hex << extractedUmi24_ << std::dec
                             << ", cbMatch=" << readBar->cbMatch
                             << ", umiSeq.length()=" << readBar->umiSeq.length() << std::endl;
        }
        
        const auto sampleDetectStart = std::chrono::steady_clock::now();
        detectSampleFromRawR2();
        const auto sampleDetectEnd = std::chrono::steady_clock::now();
        statsRA.sampleDetectOutputCalls++;
        statsRA.sampleDetectOutputNs += static_cast<uint64>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(sampleDetectEnd - sampleDetectStart).count());

        //transcripts: need to be run after CB/UMI are obtained to output CR/UR tags
        uint nAlignT = 0;
        if ((P.quant.trSAM.yes || P.quant.transcriptVB.yes) && unmapType<0 && chunkTr != nullptr && alignTrAll != nullptr) {//Aligned.toTranscriptome output, only for mapped
            if (P.pGe.transform.outQuant) {
                nAlignT = quantTranscriptome(chunkTr, alignsGenOut.alN,  alignsGenOut.alMult,  alignTrAll);
            } else {
                nAlignT = quantTranscriptome(chunkTr, nTr, trMult,  alignTrAll);
            };
        };
        
        // Build ECs for transcript quantification (TranscriptVB mode)
        if (P.quant.transcriptVB.yes && unmapType<0 && quantEC != nullptr && alignTrAll != nullptr && nAlignT > 0) {
            const char* qname = sanitizeQname(readNameMates, readNmates, 0);
            
            // Pack read sequences for error model (4-bit BAM encoding)
            // Pack once per read, not per alignment
            std::vector<uint8_t> packed_read1, packed_read1_rc, packed_read2, packed_read2_rc;
            std::vector<char> read1_rc, read2_rc;
            uint32_t read1_len = 0, read2_len = 0;
            
            if (readNmates > 0 && readLengthOriginal[0] > 0) {
                read1_len = readLengthOriginal[0];
                uint32_t packed_len1 = (read1_len + 1) / 2;
                packed_read1.resize(packed_len1);
                nuclPackBAM(Read0[0], reinterpret_cast<char*>(packed_read1.data()), read1_len);
                
                read1_rc.resize(read1_len);
                revComplementNucleotides(Read0[0], read1_rc.data(), read1_len);
                packed_read1_rc.resize(packed_len1);
                nuclPackBAM(read1_rc.data(), reinterpret_cast<char*>(packed_read1_rc.data()), read1_len);
            }
            
            if (readNmates > 1 && readLengthOriginal[1] > 0) {
                read2_len = readLengthOriginal[1];
                uint32_t packed_len2 = (read2_len + 1) / 2;
                packed_read2.resize(packed_len2);
                nuclPackBAM(Read0[1], reinterpret_cast<char*>(packed_read2.data()), read2_len);
                
                read2_rc.resize(read2_len);
                revComplementNucleotides(Read0[1], read2_rc.data(), read2_len);
                packed_read2_rc.resize(packed_len2);
                nuclPackBAM(read2_rc.data(), reinterpret_cast<char*>(packed_read2_rc.data()), read2_len);
            }
            
            quantEC->addReadAlignments(alignTrAll, nAlignT, {}, qname,
                                      packed_read1.empty() ? nullptr : packed_read1.data(),
                                      packed_read1_rc.empty() ? nullptr : packed_read1_rc.data(), read1_len,
                                      packed_read2.empty() ? nullptr : packed_read2.data(),
                                      packed_read2_rc.empty() ? nullptr : packed_read2_rc.data(), read2_len);
            
            // Compute uniformWeight for GC bias collection
            double uniformWeight = (nAlignT > 0) ? 1.0 / nAlignT : 1.0;
            
            // Collect GC observations from properly-paired transcriptomic alignments
            if (P.quant.transcriptVB.gcBias && P.readNmates == 2) {
                for (uint i = 0; i < nAlignT; i++) {
                    Transcript &aT = alignTrAll[i];
                    
                    // Check if this is a properly-paired alignment (both mates present)
                    // For PE, first exon has EX_iFrag=0, last exon has EX_iFrag=1
                    if (aT.nExons < 2) continue;
                    uint firstFrag = aT.exons[0][EX_iFrag];
                    uint lastFrag = aT.exons[aT.nExons-1][EX_iFrag];
                    if (firstFrag == lastFrag) continue;  // Both mates from same fragment = not properly paired
                    
                    // Get fragment boundaries in transcriptomic coordinates
                    // Find the leftmost and rightmost positions
                    uint64 fragStart = aT.exons[0][EX_G];
                    uint64 fragEnd = fragStart;
                    for (uint32 iex = 0; iex < aT.nExons; iex++) {
                        uint64 exStart = aT.exons[iex][EX_G];
                        uint64 exEnd = exStart + aT.exons[iex][EX_L];
                        if (exStart < fragStart) fragStart = exStart;
                        if (exEnd > fragEnd) fragEnd = exEnd;
                    }
                    
                    uint64 fragLen = fragEnd - fragStart;
                    // Only truncate very long fragments (Salmon doesn't drop small lengths)
                    // Keep fragments >= 1 and <= MAX_FRAG_LEN (2000)
                    if (fragLen == 0 || fragLen > FLDAccumulator::MAX_FRAG_LEN) continue;
                    
                    // Get transcript ID and compute GC from transcript sequence
                    uint32_t trId = aT.Chr;
                    if (trId >= chunkTr->nTr) continue;
                    
                    // Get genomic coordinates for this transcript
                    uint64 trStart = chunkTr->trS[trId];
                    uint64 trEnd = chunkTr->trE[trId];
                    
                    // Convert transcriptomic fragment coords to genomic
                    // For single-exon transcripts, this is straightforward
                    // For multi-exon, we'd need to map through exons (simplified here)
                    uint64 genomicStart = trStart + fragStart;
                    uint64 genomicEnd = trStart + fragEnd;
                    
                    // Safety check
                    if (genomicEnd > trEnd + 1) continue;
                    if (genomicEnd > mapGen.nGenome) continue;
                    
                    // Count GC bases in the fragment
                    uint64 gcCount = 0;
                    for (uint64 pos = genomicStart; pos < genomicEnd && pos < mapGen.nGenome; pos++) {
                        char base = mapGen.G[pos];
                        // In STAR encoding: 0=A, 1=C, 2=G, 3=T
                        if (base == 1 || base == 2) {  // C or G
                            gcCount++;
                        }
                    }
                    
                    // Compute GC percentage
                    int32_t gcPct = static_cast<int32_t>((gcCount * 100) / fragLen);
                    if (gcPct > 100) gcPct = 100;
                    
                    // Add GC observation with uniform weight
                    quantEC->addGCObservation(static_cast<int32_t>(fragStart), 
                                              static_cast<int32_t>(fragEnd), 
                                              gcPct, 
                                              uniformWeight);
                    
                    // FLD updates are now done in TranscriptQuantEC::addReadAlignments
                    // using Salmon-style stochastic acceptance based on exp(aln->logProb)
                    // (Removed uniform FLD updates for Salmon parity)
                }
            }
            // Note: FLD collection without GC bias is also done in TranscriptQuantEC
        }
        
        // Set detected sample token in SoloReadBarcode for tag extraction in inline hash capture
        if ((P.pSolo.flexMode || P.pSolo.inlineHashMode) && soloRead && soloRead->readBar) {
            soloRead->readBar->detectedSampleToken = detectedSampleByte_;
        }

        // Populate optional MAPQ/CIGAR/score on transcripts for downstream consumers (Flex resolver)
        // Gate Flex-only side effects so non-Flex STARsolo behaves like upstream
        if ((P.pSolo.flexMode || P.pSolo.inlineHashMode) && unmapType < 0 && nTr > 0) {
            int mapqUnique = P.outSAMmapqUnique;
            for (uint64 i = 0; i < nTr; i++) {
                if (trMult[i] == nullptr) continue;
                // Heuristic MAPQ: unique vs multimapper
                trMult[i]->mapq = (nTr == 1) ? mapqUnique : 1;
                // Alignment score / mismatch analogues for downstream ranking (AS/NM-like)
                trMult[i]->asScore = static_cast<int>(trMult[i]->maxScore);
                trMult[i]->nm = static_cast<int>(trMult[i]->nMM);
                trMult[i]->cigarString = trMult[i]->generateCigarP();
                // store chrName for probe chromosome-derived gene IDs
                trMult[i]->chrName = genOut.chrName[trMult[i]->Chr];
            }
        }

        // Store qname mapping for reject logging if enabled (Flex-only)
        // Forward declaration - function defined in flex/SoloReadFeature_record_flex.cpp
        if (P.pSolo.flexMode || P.pSolo.inlineHashMode) {
            if (readName) {
                storeQnameMapping(iReadAll, readName);
            }
        }
        
        if (!hashCacheSynthProbe_) {
            soloRead->record((unmapType<0 ? nTr : 0), trMult, iReadAll, readAnnot); //need to supply nTr=0 for unmapped reads
        }

        if (P.pGe.transform.outSAM) {
            ReadAlign::writeSAM(alignsGenOut.alN, alignsGenOut.alMult, alignsGenOut.alBest);
        } else {
            ReadAlign::writeSAM(nTr, trMult, trBest); //this will set mateMapped
        };
    };    

    if (unmapType>=0) {//unmapped reads
        statsRA.unmappedAll++; //include unmapType==4, i.e. one-mate alignments of PE reads - which may have been set in writeSAM above
        ReadAlign::outReadsUnmapped(); //uses mateMapped that was set in writeSAM above
    };    
};


///////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::recordSJ(uint64 nTrO, Transcript **trO, OutSJ *cSJ)
{//junction output for mapped reads (i.e. passed BySJout filtering)
    if (!P.outSJ.yes)
        return; //no SJ output
    
    if ( P.outSJfilterReads=="All" || nTrO==1 ) {
        uint64 sjReadStartN=cSJ->N;
        for (uint64 iTr=0; iTr<nTrO; iTr++) {//write all transcripts junctions
            outputTranscriptSJ (*(trO[iTr]), nTrO, cSJ, sjReadStartN);
        };
    };
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::outFilterBySJout()
{//filtering by SJout
    outFilterBySJoutPass=true;//only false if the alignment is held for outFilterBySJoutStage. True even if unmapped
    
    if (unmapType>0 || P.outFilterBySJoutStage!=1)
        return; //unmapped, or 2nd stage
   
    for (uint iTr=0;iTr<nTr;iTr++) {//check transcript for unannotated junctions
        for (uint iex=0;iex<trMult[iTr]->nExons-1;iex++) {//check all junctions
            if (trMult[iTr]->canonSJ[iex]>=0 && trMult[iTr]->sjAnnot[iex]==0) {
                outFilterBySJoutPass=false;
                break;
            };
        };
        if (!outFilterBySJoutPass) 
            break;
    };
    
    if (!outFilterBySJoutPass) {//this read is held for further filtering BySJout, record fastq
        unmapType=-3; //the read is not conisdered mapped
        statsRA.readN--;
        statsRA.readBases -= readLength[0]+readLength[1];

        for (uint im=0;im<P.readNends;im++) {
            chunkOutFilterBySJoutFiles[im] << readNameMates[im] <<" "<< iReadAll <<" "<< readFilter <<" "<< readFilesIndex;
            if (!readNameExtra[im].empty())
                chunkOutFilterBySJoutFiles[im]<<" "<< readNameExtra[im];
            chunkOutFilterBySJoutFiles[im] <<"\n";
            chunkOutFilterBySJoutFiles[im] << Read0[im] <<"\n";
            if (readFileType==2) {//fastq
                chunkOutFilterBySJoutFiles[im] << "+\n";
                chunkOutFilterBySJoutFiles[im] << Qual0[im] <<"\n";
            };
        };
    };
    
    //SJ output for all reads, including those not passed bySJout filtering. This only needs to be at the 1st stage of BySJout filtering
    ReadAlign::recordSJ(nTr, trMult, chunkOutSJ1); //this will set mateMapped
         
};

////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::writeSAM(uint64 nTrOutSAM, Transcript **trOutSAM, Transcript *trBestSAM)
{
    outBAMbytes=0;
    mateMapped[0] = mateMapped[1] = false; //mateMapped = are mates present in any of the transcripts?

    if (unmapType < 0 && outFilterBySJoutPass) {//write to SAM/BAM
        
        //////////////////////////////////////////////////////////////////////////////////
        /////////////outSAMfilter
        if (P.outSAMfilter.yes) {
            if (P.outSAMfilter.KeepOnlyAddedReferences) {
                for (uint itr=0;itr<nTrOutSAM;itr++) {//check if transcripts map to chr other than added references
                    if (trOutSAM[itr]->Chr<mapGen.genomeInsertChrIndFirst) {
                        return;//no SAM output
                    };
                };
            } else if (P.outSAMfilter.KeepAllAddedReferences) {
                uint64 nTrOutSAM1=0;
                for (uint itr=0;itr<nTrOutSAM;itr++) {//check if transcripts map to chr other than added references
                    if (trOutSAM[itr]->Chr>=mapGen.genomeInsertChrIndFirst) {
                        trOutSAM[nTrOutSAM1]=trOutSAM[itr];
                        trOutSAM[nTrOutSAM1]->primaryFlag=false;
                        ++nTrOutSAM1;
                    };
                };
                if (nTrOutSAM1==0) {
                   return;//no SAM output
                } else {
                    trOutSAM[0]->primaryFlag=true;
                };
                nTrOutSAM = nTrOutSAM1;
            };
        };
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////// write SAM/BAM 
        auto nTrOutWrite=min(P.outSAMmultNmax,nTrOutSAM); //number of aligns to write to SAM/BAM files            
        
        for (uint iTr=0;iTr<nTrOutWrite;iTr++) {//write transcripts
            //mateMapped1 = true if a mate is present in this transcript
            bool mateMapped1[2]={false,false};
            mateMapped1[trOutSAM[iTr]->exons[0][EX_iFrag]]=true;
            mateMapped1[trOutSAM[iTr]->exons[trOutSAM[iTr]->nExons-1][EX_iFrag]]=true;

            if (P.outSAMbool) {//SAM output
                outBAMbytes+=outputTranscriptSAM(*(trOutSAM[iTr]), nTrOutSAM, iTr, (uint) -1, (uint) -1, 0, -1, NULL, outSAMstream);
                if (P.outSAMunmapped.keepPairs && P.readNmates>1 && ( !mateMapped1[0] || !mateMapped1[1] ) ) {//keep pairs && paired reads && one of the mates not mapped in this transcript //not readNends: this is alignment
                    outBAMbytes+= outputTranscriptSAM(*(trOutSAM[iTr]), 0, 0, (uint) -1, (uint) -1, 0, 4, mateMapped1, outSAMstream);
                };
            };

            if (P.outBAMunsorted || P.outBAMcoord) {//BAM output
                alignBAM(*(trOutSAM[iTr]), nTrOutSAM, iTr, mapGen.chrStart[trOutSAM[iTr]->Chr], (uint) -1, (uint) -1, 0, -1, NULL, P.outSAMattrOrder,outBAMoneAlign, outBAMoneAlignNbytes);

                if (P.outBAMunsorted && outBAMunsorted != NULL) {//unsorted mode
                    for (uint imate=0; imate<P.readNmates; imate++) {//output each mate //not readNends: this is alignment
                        // Pass CB/UB for both mates so both get tags injected (consistent with sorted BAM behavior)
                        outBAMunsorted->unsortedOneAlign(
                            outBAMoneAlign[imate],
                            outBAMoneAlignNbytes[imate],
                            (imate>0 || iTr>0) ? 0 : (outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1])*2*nTrOutWrite,
                            iReadAll,
                            (imate==0 ? detectedSampleByte_ : 0xFFu),
                            extractedCbIdxPlus1_,  // Same CB for both mates
                            extractedUmi24_,       // Same UMI for both mates
                            extractedUmiValid_,    // UMI validity flag
                            (extractedCbIdxPlus1_ == 0 ? extractedCbSeq_ : std::string()),
                            hasYAlignment_);
                    };
                    if (P.outSAMunmapped.keepPairs && P.readNmates>1 && ( !mateMapped1[0] || !mateMapped1[1] ) ) {//keep pairs && paired reads && one of the mates not mapped in this transcript //not readNends: this is alignment
                        alignBAM(*trOutSAM[iTr], 0, 0, mapGen.chrStart[trOutSAM[iTr]->Chr], (uint) -1, (uint) -1, 0, 4, mateMapped1, P.outSAMattrOrder, outBAMoneAlign, outBAMoneAlignNbytes);
                        for (uint imate=0; imate<P.readNmates; imate++) {//output each mate //not readNends: this is alignment
                            // Pass CB/UB for unmapped pairs too (consistent tagging)
                            outBAMunsorted->unsortedOneAlign(
                                outBAMoneAlign[imate],
                                outBAMoneAlignNbytes[imate],
                                (imate>0 || iTr>0) ? 0 : (outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1])*2*nTrOutWrite,
                                iReadAll,
                                (imate==0 ? detectedSampleByte_ : 0xFFu),
                                extractedCbIdxPlus1_,  // Same CB for both mates
                                extractedUmi24_,       // Same UMI for both mates
                                extractedUmiValid_,    // UMI validity flag
                                (extractedCbIdxPlus1_ == 0 ? extractedCbSeq_ : std::string()),
                                hasYAlignment_);
                        };
                    };
                };

                if (P.outBAMcoord) {//coordinate sorted
                    for (uint imate=0; imate<P.readNmates; imate++) {//output each mate //not readNends: this is alignment
                        outBAMcoord->coordOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], (iReadAll<<32) | (iTr<<8) | trOutSAM[iTr]->exons[0][EX_iFrag], hasYAlignment_);
                    };
                };
            };
        };

        /////////////////////////////////////////////////////////////////////////////////////////////
        //////// write unmapped ends
        //TODO it's better to check all transcripts in the loop above for presence of both mates
        mateMapped[trBestSAM->exons[0][EX_iFrag]] = true;
        mateMapped[trBestSAM->exons[trBestSAM->nExons-1][EX_iFrag]] = true;

        if (P.readNmates>1 && !(mateMapped[0] && mateMapped[1]) ) {//not readNends: this is alignment
            unmapType=4;
        };

        if (unmapType==4 && P.outSAMunmapped.within) {//output unmapped ends for single-end alignments of PE reads
            if (P.outSAMbool && !P.outSAMunmapped.keepPairs ) {
                outBAMbytes+= outputTranscriptSAM(*trBestSAM, 0, 0, (uint) -1, (uint) -1, 0, unmapType, mateMapped, outSAMstream);
            };

            if ( P.outBAMcoord || (P.outBAMunsorted && !P.outSAMunmapped.keepPairs) ) {//BAM output
                alignBAM(*trBestSAM, 0, 0, mapGen.chrStart[trBestSAM->Chr], (uint) -1, (uint) -1, 0, unmapType, mateMapped, P.outSAMattrOrder, outBAMoneAlign, outBAMoneAlignNbytes);
                for (uint imate=0; imate<P.readNmates; imate++) {//alignBAM output is empty for mapped mate, but still need to scan through it //not readNends: this is alignment
                    if (P.outBAMunsorted && outBAMunsorted != NULL && !P.outSAMunmapped.keepPairs) {
                        // Pass CB/UB for both mates (consistent tagging)
                        outBAMunsorted->unsortedOneAlign(outBAMoneAlign[imate],
                                                          outBAMoneAlignNbytes[imate],
                                                          imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1],
                                                          iReadAll,
                                                          (imate==0 ? detectedSampleByte_ : 0xFFu),
                                                          extractedCbIdxPlus1_,  // Same CB for both mates
                                                          extractedUmi24_,       // Same UMI for both mates
                                                          extractedUmiValid_,    // UMI validity flag
                                                          (extractedCbIdxPlus1_ == 0 ? extractedCbSeq_ : std::string()),
                                                          hasYAlignment_);
                    };
                    if (P.outBAMcoord) {//KeepPairs option does not affect for sorted BAM since we do not want multiple entries for the same unmapped read
                        outBAMcoord->coordOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], iReadAll<<32, hasYAlignment_);
                    };
                };
            };
        };  
        
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////write completely unmapped
    } else if (unmapType>=0 && P.outSAMunmapped.within) {//output unmapped within && unmapped read && both mates unmapped
        if (P.outBAMcoord || P.outBAMunsorted || P.quant.trSAM.bamYes) {//BAM output
            alignBAM(*trBestSAM, 0, 0, mapGen.chrStart[trBestSAM->Chr], (uint) -1, (uint) -1, 0, unmapType, mateMapped, P.outSAMattrOrder, outBAMoneAlign, outBAMoneAlignNbytes);
            for (uint imate=0; imate<P.readNmates; imate++) {//output each mate //not readNends: this is alignment
                if (P.outBAMunsorted && outBAMunsorted != NULL) {
                    // Pass CB/UB for both mates (consistent tagging)
                    outBAMunsorted->unsortedOneAlign(outBAMoneAlign[imate],
                                                      outBAMoneAlignNbytes[imate],
                                                      imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1],
                                                      iReadAll,
                                                      (imate==0 ? detectedSampleByte_ : 0xFFu),
                                                      extractedCbIdxPlus1_,  // Same CB for both mates
                                                      extractedUmi24_,       // Same UMI for both mates
                                                      extractedUmiValid_,    // UMI validity flag
                                                      (extractedCbIdxPlus1_ == 0 ? extractedCbSeq_ : std::string()),
                                                      hasYAlignment_);
                };
                if (P.quant.trSAM.bamYes) {
                    // Pass CB/UB for both mates (consistent tagging)
                    outBAMquant->unsortedOneAlign(outBAMoneAlign[imate],
                                                  outBAMoneAlignNbytes[imate],
                                                  imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1],
                                                  iReadAll,
                                                  (imate==0 ? detectedSampleByte_ : 0xFFu),
                                                  extractedCbIdxPlus1_,  // Same CB for both mates
                                                  extractedUmi24_,       // Same UMI for both mates
                                                  extractedUmiValid_,    // UMI validity flag
                                                  (extractedCbIdxPlus1_ == 0 ? extractedCbSeq_ : std::string()));
                };
                if (P.outBAMcoord) {
                    outBAMcoord->coordOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], iReadAll<<32, hasYAlignment_);
                };
            };
        };

        if (P.outSAMbool) {//output SAM
            outBAMbytes+= outputTranscriptSAM(*trBestSAM, 0, 0, (uint) -1, (uint) -1, 0, unmapType, mateMapped, outSAMstream);
        };
    };       
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::writeFastxRecord(uint imate, bool isY)
{
    // Write FASTQ/FASTA record to Y or noY stream, preserving input header tokens
    const char* rawName = readNameMates[imate];
    if (rawName == nullptr || rawName[0] == 0) {
        return;
    }
    char headerPrefix = (readFileType == 2) ? '@' : '>';
    const char* name = rawName;
    if (name[0] == '@' || name[0] == '>') {
        ++name;
    }
    const string& extra = readNameExtra[imate];

    if (P.emitYNoYFastqCompression == "gz") {
        gzFile stream = isY ? chunkOutYFastqGz[imate] : chunkOutNoYFastqGz[imate];
        if (stream == nullptr) return;
        
        // Write header
        if (!extra.empty()) {
            gzprintf(stream, "%c%s %s\n", headerPrefix, name, extra.c_str());
        } else {
            gzprintf(stream, "%c%s\n", headerPrefix, name);
        }
        gzprintf(stream, "%s\n", Read0[imate]);
        if (readFileType == 2) { // fastq
            gzprintf(stream, "+\n");
            gzprintf(stream, "%s\n", Qual0[imate]);
        }
    } else {
        fstream &stream = isY ? chunkOutYFastqStream[imate] : chunkOutNoYFastqStream[imate];
        // Check if stream is open (not in bad state from gzip mode)
        if (stream.rdstate() & ios::badbit) return;
        
        // Write header
        if (!extra.empty()) {
            stream << headerPrefix << name << " " << extra << "\n";
        } else {
            stream << headerPrefix << name << "\n";
        }
        stream << Read0[imate] << "\n";
        if (readFileType == 2) { // fastq
            stream << "+\n";
            stream << Qual0[imate] << "\n";
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::outReadsUnmapped()
{
    if (P.outReadsUnmapped=="Fastx" ) {//output to fasta/q files. Include unmapType==4, i.e. one-mate alignments of PE reads
       for (uint im=0;im<P.readNends;im++) {
           chunkOutUnmappedReadsStream[im] << readNameMates[im]  <<" "<<im<<":"<< readFilter <<": "<< readNameExtra[im];
           if (P.readNmates>1) //not readNends: this is alignment
               chunkOutUnmappedReadsStream[im] <<" "<< int(mateMapped[0]) <<  int(mateMapped[1]);
           chunkOutUnmappedReadsStream[im] <<"\n";
           chunkOutUnmappedReadsStream[im] << Read0[im] <<"\n";
            if (readFileType==2) {//fastq
                chunkOutUnmappedReadsStream[im] << "+\n";
                chunkOutUnmappedReadsStream[im] << Qual0[im] <<"\n";
            };
       };
    };
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::spliceGraphWriteSAM()
{//temporary: SAM output for SpliceGraph
    outBAMbytes=0;
    uint64 nTrOutSAM = nTr;
    if (mapGen.genomeOut.convYes) {//convert to new genome
        nTrOutSAM=0;
        for (uint32 iTr=0; iTr<nTrOutSAM; iTr++) {//convert output transcripts into new genome
            *alignsGenOut.alMult[nTrOutSAM] = *trMult[iTr];//copy information before conversion
            if (trMult[iTr]->convertGenomeCigar(*mapGen.genomeOut.g, *alignsGenOut.alMult[nTrOutSAM])) {
                ++nTrOutSAM;
                trMult[nTrOutSAM-1] = alignsGenOut.alMult[nTrOutSAM-1]; //point to new transcsript
            };
        };
        nTrOutSAM=nTrOutSAM;
    };

    for (uint iTr=0; iTr<nTrOutSAM; iTr++) {//write all transcripts            
        outBAMbytes += outputSpliceGraphSAM(*(trMult[iTr]), nTrOutSAM, iTr, outSAMstream);
    };
};

void ReadAlign::alignedAnnotation()
{
    //TODO maybe initialize readAnnot to all empty?
    //genes
    if ( P.quant.geCount.yes ) {
        if (P.pGe.transform.outQuant) {
            chunkTr->geneCountsAddAlign(alignsGenOut.alN, alignsGenOut.alMult, readAnnot.geneExonOverlap);
        } else {
            chunkTr->geneCountsAddAlign(nTr, trMult, readAnnot.geneExonOverlap);
        };        
    };
    //solo-GeneFull
    if ( P.quant.geneFull.yes ) {
        chunkTr->geneFullAlignOverlap(nTr, trMult, P.pSolo.strand, readAnnot.annotFeatures[SoloFeatureTypes::GeneFull]);
    };   
    //solo-Gene (also needed for CR-compat exonic-over-intronic filter)
    bool geneClassified = false;
    if ( P.quant.gene.yes || (P.pSolo.crMultimapRescue && P.quant.geneFull.yes) ) {
        chunkTr->classifyAlign(trMult, nTr, readAnnot, slamCompat);
        geneClassified = true;
    };
    // SLAM quantification needs gene annotations for alignments (if not already done above)
    if ( P.quant.slam.yes && !geneClassified ) {
        chunkTr->classifyAlign(trMult, nTr, readAnnot, slamCompat);
        geneClassified = true;
    }
    // CR-compat: exonic-over-intronic filter for GeneFull
    // When a GeneFull alignment overlaps multiple gene bodies but only some have
    // exonic (transcript-concordant) evidence, keep only the exonic genes.
    // This resolves nested small-gene-inside-large-gene-intron cases (sentinel pattern).
    if ( P.pSolo.crMultimapRescue && P.quant.geneFull.yes && geneClassified ) {
        ReadAnnotFeature &gfAnn = readAnnot.annotFeatures[SoloFeatureTypes::GeneFull];
        const ReadAnnotFeature &gAnn = readAnnot.annotFeatures[SoloFeatureTypes::Gene];

        if ( gfAnn.fSet.size() > 1 && gfAnn.fAlign.size() == gAnn.fAlign.size() ) {
            bool anyFiltered = false;
            bool anyMultiGeneAlign = false;
            for (uint32 iA = 0; iA < (uint32)gfAnn.fAlign.size(); iA++) {
                if (gfAnn.fAlign[iA].size() <= 1)
                    continue;
                anyMultiGeneAlign = true;

                set<uint32> exonicInAlign;
                for (uint32 g : gfAnn.fAlign[iA]) {
                    if (gAnn.fAlign[iA].count(g))
                        exonicInAlign.insert(g);
                }

                if (!exonicInAlign.empty() && exonicInAlign.size() < gfAnn.fAlign[iA].size()) {
                    gfAnn.fAlign[iA] = exonicInAlign;
                    anyFiltered = true;
                    statsRA.crGeneFullExonicOverIntronicFiltered++;
                }
            }

            if (!anyMultiGeneAlign)
                statsRA.crGeneFullCrossAlignMultiGene++;

            if (anyFiltered) {
                gfAnn.fSet.clear();
                for (uint32 iA = 0; iA < (uint32)gfAnn.fAlign.size(); iA++)
                    for (uint32 g : gfAnn.fAlign[iA])
                        gfAnn.fSet.insert(g);

                if (gfAnn.fSet.size() == 1)
                    statsRA.crGeneFullResolvedToUniqueAfterFilter++;
                else if (gfAnn.fSet.size() > 1)
                    statsRA.crGeneFullStillMultiExonic++;
            }
        }
    }

    //solo-GeneFull_ExonOverIntron
    if ( P.quant.geneFull_ExonOverIntron.yes ) {
        chunkTr->geneFullAlignOverlap_ExonOverIntron(nTr, trMult, P.pSolo.strand, readAnnot.annotFeatures[SoloFeatureTypes::GeneFull_ExonOverIntron], readAnnot.annotFeatures[SoloFeatureTypes::Gene]);
    };
    // SLAM needs intronic/exonic classification for mismatch reporting
    if ( P.quant.slam.yes && !P.quant.geneFull_ExonOverIntron.yes ) {
        chunkTr->geneFullAlignOverlap_ExonOverIntron(nTr, trMult, P.pSolo.strand, readAnnot.annotFeatures[SoloFeatureTypes::GeneFull_ExonOverIntron], readAnnot.annotFeatures[SoloFeatureTypes::Gene]);
    };
    //solo-GeneFull_Ex50pAS
    if ( P.quant.geneFull_Ex50pAS.yes ) {
        chunkTr->alignExonOverlap(nTr, trMult, P.pSolo.strand, readAnnot.annotFeatures[SoloFeatureTypes::GeneFull_Ex50pAS]);
    };    
};
