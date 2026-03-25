#include "SoloReadFeature.h"
#include "SoloFeature.h"
#include "Transcriptome.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"
#include "ReadAnnotations.h"
#include "SoloReadBarcode.h"
#include "hash_shims_cpp_compat.h"
#include "ReadAlign.h"
#include "SoloBinarySpool.h"
#include "SoloReadFeature_record_shared.h"
#include "ErrorWarning.h"
#include <unordered_set>
#include <sstream>
#include <string>
#include <cstdlib>
#include <fstream>

namespace {
std::vector<uint8_t> bridgeCandidateQuals(const SoloReadBarcode &soloBar);

bool bridgeNonFlexAmbigSidecar(const SoloReadFeature *soloReadFeat, const SoloReadBarcode &soloBar)
{
    return soloReadFeat != nullptr
        && std::getenv("STAR_SOLO_NONFLEX_HASH_BRIDGE") != nullptr
        && soloBar.pSolo.inlineHashMode
        && !soloBar.pSolo.flexMode
        && soloBar.cbMatch > 1
        && !soloBar.cbMatchInd.empty();
}

int64_t bridgeCbQualTotalScoreLocal(const std::string &qual, char qsBase, uint32_t qsMax)
{
    int64_t s = 0;
    for (unsigned char uc : qual) {
        int v = static_cast<int>(uc) - static_cast<int>(qsBase);
        if (v < 0) {
            v = 0;
        }
        if (static_cast<uint32_t>(v) > qsMax) {
            v = static_cast<int>(qsMax);
        }
        s += v;
    }
    return s;
}

void ensureBridgeAmbigEntryBase(SoloReadFeature *soloReadFeat, const SoloReadBarcode &soloBar)
{
    ReadAlign::AmbigKey ambigKey = ReadAlign::hashCbSeq(soloBar.cbSeq);
    auto &entry = soloReadFeat->pendingAmbiguous_[ambigKey];
    if (!entry.candidateIdx.empty()) {
        return;
    }
    entry.candidateIdx.reserve(soloBar.cbMatchInd.size());
    for (auto idx : soloBar.cbMatchInd) {
        entry.candidateIdx.push_back(static_cast<uint32_t>(idx + 1));
    }
    entry.cbSeq = soloBar.cbSeq;
    entry.cbQual = soloBar.cbQual;
    if (entry.cbQual.length() != entry.cbSeq.length()) {
        if (entry.cbQual.length() < entry.cbSeq.length()) {
            entry.cbQual.append(entry.cbSeq.length() - entry.cbQual.length(), 'H');
        } else {
            entry.cbQual = entry.cbQual.substr(0, entry.cbSeq.length());
        }
    }
    const uint32_t kCand = static_cast<uint32_t>(entry.candidateIdx.size());
    const char qsBase = soloBar.pSolo.QSbase;
    const std::vector<uint8_t> quals = bridgeCandidateQuals(soloBar);
    entry.bridgeAmbigPinCandQuals_.resize(kCand);
    for (uint32_t ii = 0; ii < kCand; ++ii) {
        entry.bridgeAmbigPinCandQuals_[ii] = ii < quals.size() ? quals[ii] : static_cast<uint8_t>(qsBase);
    }
}

void maybeUpgradeBridgeAmbigRepresentativeQual(SoloReadFeature::ExtendedAmbiguousEntry &entry,
                                               const SoloReadBarcode &soloBar,
                                               char qsBase,
                                               uint32_t qsMax)
{
    std::string newQual = soloBar.cbQual;
    if (newQual.length() != soloBar.cbSeq.length()) {
        if (newQual.length() < soloBar.cbSeq.length()) {
            newQual.append(soloBar.cbSeq.length() - newQual.length(), 'H');
        } else {
            newQual = newQual.substr(0, soloBar.cbSeq.length());
        }
    }
    const int64_t oldScore = bridgeCbQualTotalScoreLocal(entry.cbQual, qsBase, qsMax);
    const int64_t newScore = bridgeCbQualTotalScoreLocal(newQual, qsBase, qsMax);
    if (newScore > oldScore || (newScore == oldScore && newQual < entry.cbQual)) {
        entry.cbQual = std::move(newQual);
        entry.cbSeq = soloBar.cbSeq;
        const uint32_t kCand = static_cast<uint32_t>(entry.candidateIdx.size());
        const std::vector<uint8_t> quals = bridgeCandidateQuals(soloBar);
        entry.bridgeAmbigPinCandQuals_.resize(kCand);
        for (uint32_t ii = 0; ii < kCand; ++ii) {
            entry.bridgeAmbigPinCandQuals_[ii] = ii < quals.size() ? quals[ii] : static_cast<uint8_t>(qsBase);
        }
    }
}

void maybeUpgradeBridgeAmbigOrphanQual(SoloReadFeature::BridgeAmbigReadInfoOrphanEntry &orph,
                                       const SoloReadBarcode &soloBar,
                                       char qsBase,
                                       uint32_t qsMax)
{
    std::string newQual = soloBar.cbQual;
    if (newQual.length() != soloBar.cbSeq.length()) {
        if (newQual.length() < soloBar.cbSeq.length()) {
            newQual.append(soloBar.cbSeq.length() - newQual.length(), 'H');
        } else {
            newQual = newQual.substr(0, soloBar.cbSeq.length());
        }
    }
    const int64_t oldScore = bridgeCbQualTotalScoreLocal(orph.cbQual, qsBase, qsMax);
    const int64_t newScore = bridgeCbQualTotalScoreLocal(newQual, qsBase, qsMax);
    if (newScore > oldScore || (newScore == oldScore && newQual < orph.cbQual)) {
        orph.cbQual = std::move(newQual);
        orph.cbSeq = soloBar.cbSeq;
        const uint32_t kCand = static_cast<uint32_t>(orph.candidateIdx.size());
        const std::vector<uint8_t> quals = bridgeCandidateQuals(soloBar);
        orph.pinCandQuals_.resize(kCand);
        for (uint32_t ii = 0; ii < kCand; ++ii) {
            orph.pinCandQuals_[ii] = ii < quals.size() ? quals[ii] : static_cast<uint8_t>(qsBase);
        }
    }
}

void accumulateBridgeAmbigGeneAggregates(SoloReadFeature *soloReadFeat,
                                         const SoloReadBarcode &soloBar,
                                         const SoloReadFlagClass &readFlag,
                                         bool multiFeature)
{
    if (!bridgeNonFlexAmbigSidecar(soloReadFeat, soloBar)) {
        return;
    }
    const ReadAlign::AmbigKey ambigKey = ReadAlign::hashCbSeq(soloBar.cbSeq);
    auto it = soloReadFeat->pendingAmbiguous_.find(ambigKey);
    if (it == soloReadFeat->pendingAmbiguous_.end() || it->second.candidateIdx.empty()) {
        return;
    }
    auto &entry = it->second;
    if (multiFeature) {
        entry.bridgeAmbigGeneFeatM_++;
        if (!entry.bridgeAmbigGeneHaveSampleM_) {
            entry.bridgeAmbigGeneSampleFlagM_ = readFlag.flag;
            entry.bridgeAmbigGeneHaveSampleM_ = true;
        }
    } else {
        entry.bridgeAmbigGeneFeatU_++;
        if (!entry.bridgeAmbigGeneHaveSampleU_) {
            entry.bridgeAmbigGeneSampleFlagU_ = readFlag.flag;
            entry.bridgeAmbigGeneHaveSampleU_ = true;
        }
    }
}

void accumulateBridgeAmbigReadInfoAggregates(SoloReadFeature *soloReadFeat,
                                             const SoloReadBarcode &soloBar,
                                             const SoloReadFlagClass &readFlag)
{
    if (!bridgeNonFlexAmbigSidecar(soloReadFeat, soloBar)) {
        return;
    }
    ReadAlign::AmbigKey ambigKey = ReadAlign::hashCbSeq(soloBar.cbSeq);
    auto pit = soloReadFeat->pendingAmbiguous_.find(ambigKey);
    if (pit != soloReadFeat->pendingAmbiguous_.end() && !pit->second.candidateIdx.empty()) {
        auto &entry = pit->second;
        maybeUpgradeBridgeAmbigRepresentativeQual(entry, soloBar, soloBar.pSolo.QSbase, soloBar.pSolo.QSmax);
        entry.bridgeAmbigReadInfoN_++;
        if (!entry.bridgeAmbigReadInfoHaveSample_) {
            entry.bridgeAmbigReadInfoSampleFlag_ = readFlag.flag;
            entry.bridgeAmbigReadInfoHaveSample_ = true;
        }
        return;
    }

    auto &orph = soloReadFeat->bridgeAmbigReadInfoOrphan_[ambigKey];
    if (orph.candidateIdx.empty()) {
        orph.candidateIdx.reserve(soloBar.cbMatchInd.size());
        for (auto idx : soloBar.cbMatchInd) {
            orph.candidateIdx.push_back(static_cast<uint32_t>(idx + 1));
        }
        orph.cbSeq = soloBar.cbSeq;
        orph.cbQual = soloBar.cbQual;
        if (orph.cbQual.length() != orph.cbSeq.length()) {
            if (orph.cbQual.length() < orph.cbSeq.length()) {
                orph.cbQual.append(orph.cbSeq.length() - orph.cbQual.length(), 'H');
            } else {
                orph.cbQual = orph.cbQual.substr(0, orph.cbSeq.length());
            }
        }
        const uint32_t kCand = static_cast<uint32_t>(orph.candidateIdx.size());
        const std::vector<uint8_t> quals = bridgeCandidateQuals(soloBar);
        orph.pinCandQuals_.resize(kCand);
        for (uint32_t ii = 0; ii < kCand; ++ii) {
            orph.pinCandQuals_[ii] = ii < quals.size() ? quals[ii] : static_cast<uint8_t>(soloBar.pSolo.QSbase);
        }
    } else {
        maybeUpgradeBridgeAmbigOrphanQual(orph, soloBar, soloBar.pSolo.QSbase, soloBar.pSolo.QSmax);
    }
    orph.readInfoN_++;
    if (!orph.haveSample_) {
        orph.sampleFlag_ = readFlag.flag;
        orph.haveSample_ = true;
    }
}

void accumulateBridgeImmediateReadFlags(SoloReadFeature *soloReadFeat,
                                        uint32_t cbIdx,
                                        SoloReadFlagClass::typeFlag baseFlag,
                                        bool featGood,
                                        bool multiFeature,
                                        int32_t cbMatch,
                                        bool readStatsEnabled)
{
    if (soloReadFeat == nullptr) {
        return;
    }

    auto &packedCounts = soloReadFeat->bridgeImmediateReadCounts_[cbIdx];
    if (featGood) {
        if (multiFeature) {
            packedCounts += (uint64_t{1} << 32);
        } else {
            packedCounts += 1;
        }
    }

    if (!readStatsEnabled) {
        return;
    }

    SoloReadFlagClass::typeFlag finalFlag = baseFlag;
    finalFlag |= (SoloReadFlagClass::typeFlag{1} << SoloReadFlagClass::cbMatch);
    if (cbMatch == 0) {
        finalFlag |= (SoloReadFlagClass::typeFlag{1} << SoloReadFlagClass::cbPerfect);
    } else if (cbMatch == 1) {
        finalFlag |= (SoloReadFlagClass::typeFlag{1} << SoloReadFlagClass::cbMMunique);
    }

    if (featGood) {
        if (multiFeature) {
            finalFlag |= (SoloReadFlagClass::typeFlag{1} << SoloReadFlagClass::countedM);
        } else {
            finalFlag |= (SoloReadFlagClass::typeFlag{1} << SoloReadFlagClass::countedU);
        }
    }

    auto &arr = soloReadFeat->readFlag.flagCounts[cbIdx];
    for (uint32_t bit = 0; bit < SoloReadFlagClass::nBits; ++bit) {
        arr[bit] += static_cast<uint64_t>((finalFlag >> bit) & SoloReadFlagClass::typeFlag{1});
    }
}

const std::unordered_set<std::string>& bridgeDebugBarcodeSet()
{
    static std::unordered_set<std::string> barcodes;
    static bool loaded = false;
    if (loaded) {
        return barcodes;
    }
    loaded = true;

    const char *path = std::getenv("STAR_SOLO_DEBUG_BARCODE_FILE");
    if (path == nullptr || path[0] == '\0') {
        return barcodes;
    }

    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty()) {
            barcodes.insert(line);
        }
    }
    return barcodes;
}

bool shouldTraceBridgeBarcode(const ParametersSolo &pSolo, uint32_t wlIdx)
{
    if (wlIdx >= pSolo.cbWLstr.size()) {
        return false;
    }
    const auto &debugSet = bridgeDebugBarcodeSet();
    if (debugSet.empty()) {
        return false;
    }
    return debugSet.count(pSolo.cbWLstr[wlIdx]) != 0;
}

void insertInlineHashEntry(SoloReadFeature *soloReadFeat, const SoloReadBarcode &soloBar, uint32_t geneIdx)
{
    if (soloReadFeat == nullptr
        || soloReadFeat->inlineHash_ == nullptr
        || soloBar.cbMatch < 0
        || soloBar.cbMatch > 1
        || soloBar.cbMatchInd.empty()) {
        return;
    }

    const uint32_t wlIdx = soloBar.cbMatchInd[0];
    uint32_t umi24 = soloBar.umiB & 0xFFFFFF;
    uint64_t key = 0;
    if (std::getenv("STAR_SOLO_NONFLEX_HASH_BRIDGE") != nullptr) {
        if (geneIdx > 0xFFFFu) {
            std::fprintf(stderr,
                         "FATAL ERROR: non-Flex Solo direct bridge: gene index %u exceeds 16-bit tuple field\n",
                         geneIdx);
            std::exit(1);
        }
        if (wlIdx >= 0x1000000u) {
            std::fprintf(stderr,
                         "FATAL ERROR: non-Flex Solo direct bridge: whitelist CB index %u exceeds 24-bit key field\n",
                         wlIdx);
            std::exit(1);
        }
        key = packBridgeWlUmiGeneKey(wlIdx, umi24, static_cast<uint16_t>(geneIdx));
    } else {
        key = packCgAggKey(wlIdx, umi24, static_cast<uint16_t>(geneIdx), 0);
    }
    int absent;
    khiter_t iter = kh_put(cg_agg, soloReadFeat->inlineHash_, key, &absent);
    if (absent) {
        kh_val(soloReadFeat->inlineHash_, iter) = 1;
    } else {
        kh_val(soloReadFeat->inlineHash_, iter)++;
    }

    if (shouldTraceBridgeBarcode(soloBar.pSolo, wlIdx)) {
        std::fprintf(stderr,
                     "[BRIDGE-CB-TRACE] mode=direct wlIdx=%u cb=%s cbMatch=%d gene=%u umi=%u\n",
                     wlIdx,
                     soloBar.pSolo.cbWLstr[wlIdx].c_str(),
                     soloBar.cbMatch,
                     geneIdx,
                     umi24);
    }
}

std::vector<uint8_t> bridgeCandidateQuals(const SoloReadBarcode &soloBar)
{
    std::vector<uint8_t> quals;
    if (soloBar.cbMatch <= 1) {
        return quals;
    }

    quals.reserve(soloBar.cbMatchInd.size());
    std::istringstream qualsIn(soloBar.cbMatchString);
    uint32_t cbIdx;
    char qualChar;
    while (qualsIn >> cbIdx >> qualChar) {
        quals.push_back(static_cast<uint8_t>(qualChar));
    }
    return quals;
}

void accumulateBridgeAmbiguousCB(SoloReadFeature *soloReadFeat,
                                 const SoloReadBarcode &soloBar,
                                 uint32_t geneIdx)
{
    bool isAmbiguous = (soloBar.cbMatchInd.size() > 1) || (soloBar.cbMatch > 1);
    if (!soloReadFeat || !soloReadFeat->inlineHash_ || !isAmbiguous || soloBar.cbMatchInd.empty()) {
        return;
    }

    if (geneIdx > 0xFFFFu) {
        std::fprintf(stderr,
                     "FATAL ERROR: non-Flex Solo direct bridge: gene index %u exceeds 16-bit ambiguous observation field\n",
                     geneIdx);
        std::exit(1);
    }

    ensureBridgeAmbigEntryBase(soloReadFeat, soloBar);
    ReadAlign::AmbigKey ambigKey = ReadAlign::hashCbSeq(soloBar.cbSeq);
    auto &entry = soloReadFeat->pendingAmbiguous_[ambigKey];
    maybeUpgradeBridgeAmbigRepresentativeQual(entry, soloBar, soloBar.pSolo.QSbase, soloBar.pSolo.QSmax);

    auto oit = soloReadFeat->bridgeAmbigReadInfoOrphan_.find(ambigKey);
    if (oit != soloReadFeat->bridgeAmbigReadInfoOrphan_.end()) {
        entry.bridgeAmbigReadInfoN_ += oit->second.readInfoN_;
        if (!entry.bridgeAmbigReadInfoHaveSample_ && oit->second.haveSample_) {
            entry.bridgeAmbigReadInfoSampleFlag_ = oit->second.sampleFlag_;
            entry.bridgeAmbigReadInfoHaveSample_ = true;
        }
        soloReadFeat->bridgeAmbigReadInfoOrphan_.erase(oit);
    }

    const uint32_t umi24 = soloBar.umiB & 0xFFFFFF;
    const uint64_t ugKey = (static_cast<uint64_t>(umi24 & 0xFFFFFFu) << 16)
        | static_cast<uint64_t>(static_cast<uint16_t>(geneIdx));
    entry.bridgeAmbigUmiGene_[ugKey] += 1u;
}
}

#ifdef DEBUG_CB_UB_PARITY
// Optional tracing of specific readIds via STAR_DEBUG_TRACE_READS=1,2,3
static std::unordered_set<uint32_t> buildTraceReadSetWriter() {
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
static const std::unordered_set<uint32_t> g_traceReadsWriter = buildTraceReadSetWriter();
#endif

void record_base(SoloReadFeature *soloReadFeat, SoloReadBarcode &soloBar, uint nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot)
{
    if (soloReadFeat->pSolo.type==0)
        return;

    if (soloReadFeat->pSolo.readStatsYes[soloReadFeat->featureType]) {//readFlag

        if (nTr==1) {
            soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.genomeU);
        } else if (nTr>1) {
            soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.genomeM);
        };

        for (uint64 itr=0; itr<nTr; itr++) {
            if (soloReadFeat->P.pGe.chrSet.mito.count(alignOut[itr]->Chr) == 1) {
                soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.mito);
            };
        };

        switch (readAnnot.annotFeatures[soloReadFeat->featureType].ovType) {
            case ReadAnnotFeature::overlapTypes::exonic : 
            case ReadAnnotFeature::overlapTypes::exonic50p :
                soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.exonic);
                break;
            case ReadAnnotFeature::overlapTypes::intronic : 
                soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.intronic);
                break;
            case ReadAnnotFeature::overlapTypes::exonicAS : 
            case ReadAnnotFeature::overlapTypes::exonic50pAS :
                soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.exonicAS);            
                break;
            case ReadAnnotFeature::overlapTypes::intronicAS :
                soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.intronicAS);            
        };

        if (soloBar.cbMatch<0 && soloReadFeat->pSolo.cbWLyes) {//no CB match in the WL
            if (readAnnot.annotFeatures[soloReadFeat->featureType].fSet.size()==1) {
                soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.featureU);
            } else if (readAnnot.annotFeatures[soloReadFeat->featureType].fSet.size()>1){
                soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.featureM);
            };
            soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.cbMatch);//this will counts reads with no CB match
            soloReadFeat->readFlag.countsAddNoCB();
        };
    };

    if (soloBar.cbMatch<0) {
        return;
    }

       
    ReadSoloFeatures reFe;
    reFe.alignOut=alignOut;
    reFe.indAnnotTr = 0;    

    uint32 nFeat=0; //number of features in this read (could be >1 for SJs)
    if (nTr==0) {//unmapped
        soloReadFeat->stats.V[soloReadFeat->stats.noUnmapped]++;
        
    } else {
        switch (soloReadFeat->featureType) {
            case SoloFeatureTypes::Gene :
            case SoloFeatureTypes::GeneFull :
            case SoloFeatureTypes::GeneFull_Ex50pAS :
            case SoloFeatureTypes::GeneFull_ExonOverIntron :
                {
                    auto *readGe = &readAnnot.annotFeatures[soloReadFeat->featureType].fSet;

                    if (soloBar.pSolo.type==soloBar.pSolo.SoloTypes::SmartSeq) {
                        for (int32 itr=nTr-1; itr>=0; itr--) {
                            if (readAnnot.annotFeatures[soloReadFeat->featureType].fAlign[itr].size() > 0) {
                                reFe.indAnnotTr = itr;
                                break;//indAnnotTr is the first genic align, and is used for Smart-seq deduplication, below in outputReadCB()
                                      //TODO this is ambiguous for reads that map to the same gene multiple times. *Last* align is chosen to match the previous version.
                            };
                        };
                    };
                        
                    if (readGe->size()==0) {//check genes
                        soloReadFeat->stats.V[soloReadFeat->stats.noNoFeature]++;//no gene
                    } else if (readGe->size()>1) {
                        soloReadFeat->stats.V[soloReadFeat->stats.MultiFeature]++;//multi-gene reads
                        soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.featureM);
                        if (nTr>1)
                            soloReadFeat->stats.V[soloReadFeat->stats.subMultiFeatureMultiGenomic]++;//multigene caused by multimapper
#ifdef DEBUG_CB_UB_PARITY
                        if (!g_traceReadsWriter.empty() && iRead!=(uint64)-1 && g_traceReadsWriter.count((uint32_t)iRead)) {
                            std::string genes;
                            genes.reserve(readGe->size()*8);
                            bool first=true;
                            for (auto gVal : *readGe) {
                                if (!first) genes.push_back(',');
                                first=false;
                                genes += std::to_string(gVal);
                            }
                            fprintf(stderr, "[TRACE annotate] read=%llu genes=%s nTr=%u cbMatch=%d cbMatchString=%s\n",
                                    (unsigned long long)iRead, genes.c_str(), (unsigned)nTr,
                                    soloBar.cbMatch, soloBar.cbMatchString.c_str());
                        }
#endif
                            
                        if (soloReadFeat->pSolo.multiMap.yes.multi) {//output multimappers
                            reFe.geneMult.resize(readGe->size());
                            uint32 ii=0;
                            for (auto &g: *readGe) {//set high bit to mark multimappers
                                reFe.geneMult[ii] = g | geneMultMark;
                                ++ii;
                            };
                                
                            nFeat = outputReadCB_base(soloReadFeat->streamReads, iRead, soloReadFeat->featureType, soloBar, reFe, readAnnot, soloReadFeat->readFlag, soloReadFeat);
                        };
                    } else {//unique-gene reads
                        reFe.gene = *readGe->begin();
                        soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.featureU);
                        nFeat = outputReadCB_base(soloReadFeat->streamReads, (soloReadFeat->readIndexYes ? iRead : (uint64)-1), soloReadFeat->featureType, soloBar, reFe, readAnnot, soloReadFeat->readFlag, soloReadFeat);
                    };
                };
                break;
        
            case SoloFeatureTypes::SJ : 
                if (nTr>1) {//reject all multimapping reads
                    soloReadFeat->stats.V[soloReadFeat->stats.subMultiFeatureMultiGenomic]++;
                    soloReadFeat->stats.V[soloReadFeat->stats.MultiFeature]++;
                } else {//one gene or no gene
                    alignOut[0]->extractSpliceJunctions(reFe.sj, reFe.sjAnnot);
                    if ( reFe.sj.empty() ) {
                        soloReadFeat->stats.V[soloReadFeat->stats.noNoFeature]++;
                    } else {//good junction
                        soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.featureU);
                        nFeat = outputReadCB_base(soloReadFeat->streamReads, (soloReadFeat->readIndexYes ? iRead : (uint64)-1), soloReadFeat->featureType, soloBar, reFe, readAnnot, soloReadFeat->readFlag, soloReadFeat);
                    };
                };                  
                break;
        
            case SoloFeatureTypes::Transcript3p :
                if (readAnnot.transcriptConcordant.size()==0 || soloBar.cbMatch>1) {//do not record ambiguous CB  
                    soloReadFeat->stats.V[soloReadFeat->stats.noNoFeature]++;
                } else {
                    nFeat = outputReadCB_base(soloReadFeat->streamReads, iRead, soloReadFeat->featureType, soloBar, reFe, readAnnot, soloReadFeat->readFlag, soloReadFeat);
                };
                if (readAnnot.transcriptConcordant.size()==1 && readAnnot.transcriptConcordant[0][1] < soloReadFeat->transcriptDistCount.size()) {
                    //read maps to one transcript - use for distTTS distribution function
                    soloReadFeat->transcriptDistCount[readAnnot.transcriptConcordant[0][1]]++;
                };
                
                break;

            case SoloFeatureTypes::Velocyto :
                //different record: iRead, nTr, tr1, type1, tr2, type2 ...
                if (readAnnot.trVelocytoType.size()>0) {//otherwise, no gene
                    
                    sort(readAnnot.trVelocytoType.begin(), readAnnot.trVelocytoType.end(),
                         [](const trTypeStruct &t1, const trTypeStruct &t2) {return t1.tr < t2.tr;});

                    if (soloReadFeat->streamReads) {
                        *soloReadFeat->streamReads << iRead <<' '<< readAnnot.trVelocytoType.size();
                        for (auto &tt: readAnnot.trVelocytoType)
                             *soloReadFeat->streamReads <<' '<< tt.tr <<' '<< (uint32) tt.type;
                        *soloReadFeat->streamReads <<'\n';
                    }
                    nFeat=1;
                } else {
                    soloReadFeat->stats.V[soloReadFeat->stats.noNoFeature]++;
                };
                break; //no need to go with downstream processing                
                
        };//switch (featureType)
    };//if (nTr==0)
    
    if ( nFeat==0 && (soloReadFeat->readInfoYes | soloReadFeat->pSolo.readStatsYes[soloReadFeat->featureType]) ) {//no feature, but readInfo requested
        outputReadCB_base(soloReadFeat->streamReads, iRead, (uint32)-1, soloBar, reFe, readAnnot, soloReadFeat->readFlag, soloReadFeat);
    };
    
    if (nFeat==0)
        return; //no need to record the number of reads per CB
    
    if (soloReadFeat->pSolo.cbWLyes) {//WL
        for (auto &cbi : soloBar.cbMatchInd)
            soloReadFeat->cbReadCount[cbi] += nFeat;
    } else {//no WL
        soloReadFeat->cbReadCountMap[soloBar.cbMatchInd[0]] += nFeat;
    };
    
    return;
};

uint32 outputReadCB_base(fstream *streamOut, const uint64 iRead, const int32 featureType, SoloReadBarcode &soloBar, 
                         const ReadSoloFeatures &reFe, const ReadAnnotations &readAnnot, const SoloReadFlagClass &readFlag,
                         SoloReadFeature *soloReadFeat)
{   
    /*format of the temp output file
     * UMI [iRead] type feature* cbMatchString
     *             0=exact match, 1=one non-exact match, 2=multipe non-exact matches
     *                   gene or sj[0] sj[1]
     *                         CB or nCB {CB Qual, ...}
     */
    
    if (soloBar.pSolo.type==soloBar.pSolo.SoloTypes::SmartSeq && featureType!=-1) {//need to calculate "UMI" from align start/end
        soloBar.umiB=reFe.alignOut[reFe.indAnnotTr]->chrStartLengthExtended();
    };

    const bool binarySpool = soloReadFeat != nullptr && soloReadFeat->binarySpool;
    std::vector<std::pair<uint32_t, uint8_t>> ambiguousCbCandidates;
    if (binarySpool && soloBar.cbMatch > 1) {
        ambiguousCbCandidates.reserve(static_cast<size_t>(soloBar.cbMatch));
        std::istringstream cbMatchStream(soloBar.cbMatchString);
        uint32_t cbIdx = 0;
        char qualChar = 0;
        while (cbMatchStream >> cbIdx >> qualChar) {
            ambiguousCbCandidates.emplace_back(cbIdx, static_cast<uint8_t>(qualChar));
        }
        if (static_cast<int32>(ambiguousCbCandidates.size()) != soloBar.cbMatch) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: could not encode ambiguous CB payload into experimental binary Solo spool\n"
                   << "cbMatch=" << soloBar.cbMatch << " parsedCandidates=" << ambiguousCbCandidates.size() << "\n";
            exitWithError(errOut.str(), std::cerr, soloReadFeat->P.inOut->logMain, EXIT_CODE_PARAMETER, soloReadFeat->P);
        }
    }

    auto writeCbPayload = [&]() {
        const bool binarySpoolInMemoryNow = binarySpool && soloReadFeat->binarySpoolInMemory;
        fstream *activeStream = binarySpool ? soloReadFeat->streamReads : streamOut;
        if (activeStream == nullptr && !binarySpoolInMemoryNow) {
            return;
        }
        if (binarySpool) {
            if (soloBar.cbMatch <= 1) {
                uint64_t cbValue = 0;
                if (!soloBar.cbMatchInd.empty()) {
                    cbValue = soloBar.cbMatchInd[0];
                } else if (!soloBar.cbMatchString.empty()) {
                    cbValue = static_cast<uint64_t>(std::stoull(soloBar.cbMatchString));
                } else {
                    ostringstream errOut;
                    errOut << "EXITING because of fatal ERROR: missing resolved CB payload for experimental binary Solo spool\n";
                    exitWithError(errOut.str(), std::cerr, soloReadFeat->P.inOut->logMain, EXIT_CODE_PARAMETER, soloReadFeat->P);
                }
                if (binarySpoolInMemoryNow) {
                    SoloBinarySpool::writeResolvedCb(soloReadFeat->binarySpoolBuffer, cbValue);
                } else {
                    SoloBinarySpool::writeResolvedCb(*activeStream, cbValue);
                }
            } else {
                for (const auto &candidate : ambiguousCbCandidates) {
                    if (binarySpoolInMemoryNow) {
                        SoloBinarySpool::writeAmbiguousCandidate(soloReadFeat->binarySpoolBuffer, candidate.first, candidate.second);
                    } else {
                        SoloBinarySpool::writeAmbiguousCandidate(*activeStream, candidate.first, candidate.second);
                    }
                }
            }
        } else {
            *activeStream << soloBar.cbMatchString;
        }
    };
    auto binaryHeaderBytes = [&](bool writeReadIndex) -> size_t {
        return sizeof(uint64_t) + sizeof(uint32_t) + sizeof(int32_t) +
               (writeReadIndex ? (sizeof(uint64_t) + sizeof(uint32_t)) : 0);
    };
    auto binaryPayloadBytes = [&]() -> size_t {
        if (soloBar.cbMatch <= 1) {
            return sizeof(uint64_t);
        }
        return static_cast<size_t>(soloBar.cbMatch) * (sizeof(uint32_t) + sizeof(uint8_t));
    };
    auto maybeSpillBeforeBinaryWrite = [&](bool writeReadIndex) {
        if (binarySpool && soloReadFeat->binarySpoolInMemory) {
            soloReadFeat->maybeSpillBinarySpool(binaryHeaderBytes(writeReadIndex) + binaryPayloadBytes());
        }
    };
    
    uint64 nout=1;
    
    switch (featureType) {
        case -1 : {
            // no feature, output for readInfo
            bool nonFlexHashBridge = soloReadFeat != nullptr
                && soloReadFeat->pSolo.inlineHashMode
                && !soloReadFeat->pSolo.flexMode
                && std::getenv("STAR_SOLO_NONFLEX_HASH_BRIDGE") != nullptr;
            if (nonFlexHashBridge) {
                if (bridgeNonFlexAmbigSidecar(soloReadFeat, soloBar)) {
                    accumulateBridgeAmbigReadInfoAggregates(soloReadFeat, soloBar, readFlag);
                } else if (soloBar.cbMatch >= 0 && soloBar.cbMatch <= 1 && !soloBar.cbMatchInd.empty()) {
                    accumulateBridgeImmediateReadFlags(soloReadFeat,
                                                       static_cast<uint32_t>(soloBar.cbMatchInd[0]),
                                                       readFlag.flag,
                                                       false,
                                                       false,
                                                       soloBar.cbMatch,
                                                       soloReadFeat->pSolo.readStatsYes[soloReadFeat->featureType]);
                }
                break;
            }
            if (streamOut || (binarySpool && soloReadFeat->binarySpoolInMemory)) {
                if (binarySpool) {
                    maybeSpillBeforeBinaryWrite(true);
                    fstream *activeStream = soloReadFeat->binarySpoolInMemory ? nullptr : soloReadFeat->streamReads;
                    if (soloReadFeat->binarySpoolInMemory) {
                        SoloBinarySpool::writeRecordHeader(soloReadFeat->binarySpoolBuffer, true, soloBar.umiB, iRead, readFlag.flag, static_cast<uint32_t>(-1), soloBar.cbMatch);
                    } else {
                        SoloBinarySpool::writeRecordHeader(*activeStream, true, soloBar.umiB, iRead, readFlag.flag, static_cast<uint32_t>(-1), soloBar.cbMatch);
                    }
                    writeCbPayload();
                } else {
                    *streamOut << soloBar.umiB <<' '<< iRead <<' '<< readFlag.flag <<' '<< -1 <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';
                }
            }
            break;
        }
            
        case SoloFeatureTypes::Gene :
        case SoloFeatureTypes::GeneFull :
        case SoloFeatureTypes::GeneFull_Ex50pAS :
        case SoloFeatureTypes::GeneFull_ExonOverIntron : {
            bool nonFlexHashBridge = soloReadFeat != nullptr
                && soloReadFeat->pSolo.inlineHashMode
                && !soloReadFeat->pSolo.flexMode
                && std::getenv("STAR_SOLO_NONFLEX_HASH_BRIDGE") != nullptr;
            if (nonFlexHashBridge) {
                const bool multiFeature = !reFe.geneMult.empty();
                if (soloBar.cbMatch >= 0 && soloBar.cbMatch <= 1 && !soloBar.cbMatchInd.empty()) {
                    accumulateBridgeImmediateReadFlags(soloReadFeat,
                                                       static_cast<uint32_t>(soloBar.cbMatchInd[0]),
                                                       readFlag.flag,
                                                       true,
                                                       multiFeature,
                                                       soloBar.cbMatch,
                                                       soloReadFeat->pSolo.readStatsYes[soloReadFeat->featureType]);
                    if (soloBar.cbMatch == 0) {
                        soloReadFeat->stats.V[soloReadFeat->stats.yessubWLmatchExact]++;
                    }
                }
                const bool isAmbiguous = (soloBar.cbMatchInd.size() > 1) || (soloBar.cbMatch > 1);
                if (!reFe.geneMult.empty()) {
                    for (uint32_t geneIdx : reFe.geneMult) {
                        uint32_t resolvedGeneIdx = static_cast<uint32_t>(geneIdx ^ geneMultMark);
                        if (isAmbiguous) {
                            accumulateBridgeAmbiguousCB(soloReadFeat, soloBar, resolvedGeneIdx);
                        } else {
                            insertInlineHashEntry(soloReadFeat, soloBar, resolvedGeneIdx);
                        }
                    }
                    nout = reFe.geneMult.size();
                } else {
                    if (isAmbiguous) {
                        accumulateBridgeAmbiguousCB(soloReadFeat, soloBar, reFe.gene);
                    } else {
                        insertInlineHashEntry(soloReadFeat, soloBar, reFe.gene);
                    }
                    nout = 1;
                }
                if (isAmbiguous && bridgeNonFlexAmbigSidecar(soloReadFeat, soloBar)) {
                    accumulateBridgeAmbigGeneAggregates(soloReadFeat, soloBar, readFlag, multiFeature);
                }
                break;
            }

            if (streamOut || (binarySpool && soloReadFeat->binarySpoolInMemory)) {
                if (!reFe.geneMult.empty()) {
                    for (uint32_t geneIdx : reFe.geneMult) {
                        if (binarySpool) {
                            maybeSpillBeforeBinaryWrite(iRead != (uint64)-1);
                            fstream *activeStream = soloReadFeat->binarySpoolInMemory ? nullptr : soloReadFeat->streamReads;
                            if (soloReadFeat->binarySpoolInMemory) {
                                SoloBinarySpool::writeRecordHeader(soloReadFeat->binarySpoolBuffer, iRead != (uint64)-1, soloBar.umiB, iRead, readFlag.flag, geneIdx, soloBar.cbMatch);
                            } else {
                                SoloBinarySpool::writeRecordHeader(*activeStream, iRead != (uint64)-1, soloBar.umiB, iRead, readFlag.flag, geneIdx, soloBar.cbMatch);
                            }
                            writeCbPayload();
                        } else {
                            *streamOut << soloBar.umiB <<' ';
                            if (iRead != (uint64)-1) {
                                *streamOut << iRead <<' '<< readFlag.flag <<' ';
                            }
                            *streamOut << geneIdx <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';
                        }
                    }
                    nout = reFe.geneMult.size();
                } else {
                    if (binarySpool) {
                        maybeSpillBeforeBinaryWrite(iRead != (uint64)-1);
                        fstream *activeStream = soloReadFeat->binarySpoolInMemory ? nullptr : soloReadFeat->streamReads;
                        if (soloReadFeat->binarySpoolInMemory) {
                            SoloBinarySpool::writeRecordHeader(soloReadFeat->binarySpoolBuffer, iRead != (uint64)-1, soloBar.umiB, iRead, readFlag.flag, reFe.gene, soloBar.cbMatch);
                        } else {
                            SoloBinarySpool::writeRecordHeader(*activeStream, iRead != (uint64)-1, soloBar.umiB, iRead, readFlag.flag, reFe.gene, soloBar.cbMatch);
                        }
                        writeCbPayload();
                    } else {
                        *streamOut << soloBar.umiB <<' ';
                        if (iRead != (uint64)-1) {
                            *streamOut << iRead <<' '<< readFlag.flag <<' ';
                        }
                        *streamOut << reFe.gene <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';
                    }
                    nout = 1;
                }
            }
            break;
        }

        case SoloFeatureTypes::SJ : {
            //sj - two numbers, multiple sjs per read
            if (streamOut) {
                for (auto &sj : reFe.sj) {
                    *streamOut << soloBar.umiB <<' ';//UMI
                    if ( iRead != (uint64)-1 )
                        *streamOut << iRead <<' '<< readFlag.flag <<' ';//iRead            
                    *streamOut << sj[0] <<' '<< sj[1] <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n' << flush;
                }
            }
            nout=reFe.sj.size();
            break;
        }

        case SoloFeatureTypes::Transcript3p : {
            //transcript,distToTTS structure
            if (streamOut) {
                *streamOut << soloBar.cbMatchString <<' ';            
                *streamOut << soloBar.umiB <<' ';
                *streamOut << readAnnot.transcriptConcordant.size();
                for (auto &tt: readAnnot.transcriptConcordant) {
                    *streamOut <<' '<< tt[0] <<' '<< tt[1];
                };
                if ( iRead != (uint64)-1 )
                    *streamOut  <<' '<< iRead;//iRead
                *streamOut  <<'\n';
            }
            nout=1;

            break;
        }
        default:
            break;
    }; //switch (featureType)
    
    return nout;
};
