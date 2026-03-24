#include "SoloReadFeature.h"
#include "streamFuns.h"
#include "SoloFeatureTypes.h"
#include <cstdlib>
#include <cstdio>

SoloReadFeature::SoloReadFeature(int32 feTy, Parameters &Pin, int iChunk)
             : featureType(feTy), P(Pin), pSolo(P.pSolo), streamReads(nullptr), inlineHash_(nullptr), readIdTracker_(nullptr)
{
    if (pSolo.type==0)
        return;
//     if (pSolo.type==pSolo.SoloTypes::CB_samTagOut)
//         return;
    
    readInfoYes = pSolo.readInfoYes[featureType];
    readIndexYes = pSolo.readIndexYes[featureType];
    
    if (pSolo.cbWLyes) {
        cbReadCount.resize(pSolo.cbWLsize,0);
    };

    if (pSolo.inlineHashMode) {
        // Initialize inline hash instead of opening temp stream file
        inlineHash_ = kh_init(cg_agg);
        streamReads = nullptr; // Do NOT open stream file in inline hash mode
        
        // Initialize parallel readId tracker for sorted BAM CB/UB tag injection
        if (pSolo.trackReadIdsForTags) {
            readIdTracker_ = kh_init(readid_cbumi);
        }
    } else if (iChunk>=0) {
        //open with flagDelete=false, i.e. try to keep file if it exists
        streamReads = &fstrOpen(P.outFileTmp+"/solo"+SoloFeatureTypes::Names[featureType]+'_'+std::to_string(iChunk), ERROR_OUT, P, false);
    };
    
    if (featureType==SoloFeatureTypes::Transcript3p)
        transcriptDistCount.resize(10000,0);
};

SoloReadFeature::~SoloReadFeature() {
    if (inlineHash_) {
        kh_destroy(cg_agg, inlineHash_);
        inlineHash_ = nullptr;
    }
    if (readIdTracker_) {
        kh_destroy(readid_cbumi, readIdTracker_);
        readIdTracker_ = nullptr;
    }
}

void SoloReadFeature::addCounts(const SoloReadFeature &rfIn)
{
    if (pSolo.cbWLyes) {//WL
        for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
            cbReadCount[ii] += rfIn.cbReadCount[ii];
        };
    } else {
        for (auto ii=rfIn.cbReadCountMap.cbegin(); ii!=rfIn.cbReadCountMap.cend(); ++ii) {
            cbReadCountMap[ii->first] += ii->second;
        };
    };
    
    if (transcriptDistCount.size()>0) {
        for (uint32 ii=0; ii<transcriptDistCount.size(); ii++)
            transcriptDistCount[ii] += rfIn.transcriptDistCount[ii];
    };
};

void SoloReadFeature::addStats(const SoloReadFeature &rfIn)
{
    for (uint32 ii=0; ii<stats.nStats; ii++)
        stats.V[ii] += rfIn.stats.V[ii];

    for (const auto &kv : rfIn.readFlag.flagCounts) {
        auto &dst = readFlag.flagCounts[kv.first];
        for (uint32 ii = 0; ii < readFlag.nBits; ++ii) {
            dst[ii] += kv.second[ii];
        }
    }
    for (uint32 ii=0; ii<readFlag.nBits; ii++)
        readFlag.flagCountsNoCB[ii] += rfIn.readFlag.flagCountsNoCB[ii];

    for (const auto &kv : rfIn.bridgeImmediateReadCounts_) {
        bridgeImmediateReadCounts_[kv.first] += kv.second;
    }
};

void SoloReadFeature::statsOut(ofstream &streamOut)
{
    //streamOut << setw(50) << "CELL BARCODES IN READS:\n"
    for (uint32 ii=0; ii<stats.nStats; ii++) {
        streamOut << setw(50) << stats.names[ii] << setw(15) << stats.V[ii] << '\n';
    };
    streamOut.flush();
};

uint32_t SoloReadFeature::getOrCreateBridgeCompactCb(uint32_t wlIdx)
{
    auto it = bridgeCbCompactByWl_.find(wlIdx);
    if (it != bridgeCbCompactByWl_.end()) {
        return it->second;
    }

    uint32_t compactIdx = bridgeCbWlByCompact_.size();
    bridgeCbCompactByWl_[wlIdx] = compactIdx;
    bridgeCbWlByCompact_.push_back(wlIdx);
    return compactIdx;
}

uint32_t SoloReadFeature::bridgeCompactToWl(uint32_t compactIdx) const
{
    if (compactIdx >= bridgeCbWlByCompact_.size()) {
        return static_cast<uint32_t>(-1);
    }
    return bridgeCbWlByCompact_[compactIdx];
}

uint16_t SoloReadFeature::getOrCreateBridgeCompactGene(uint32_t geneIdx)
{
    auto it = bridgeGeneCompactByFull_.find(geneIdx);
    if (it != bridgeGeneCompactByFull_.end()) {
        return it->second;
    }

    size_t compactIdx = bridgeGeneFullByCompact_.size();
    if (compactIdx > 0xFFFFu) {
        std::fprintf(stderr,
                     "FATAL ERROR: non-Flex Solo hash bridge observed more than 65536 distinct genes; packed gene field overflow\n");
        std::exit(1);
    }

    uint16_t compactGene = static_cast<uint16_t>(compactIdx);
    bridgeGeneCompactByFull_[geneIdx] = compactGene;
    bridgeGeneFullByCompact_.push_back(geneIdx);
    return compactGene;
}

uint32_t SoloReadFeature::bridgeCompactToGene(uint16_t compactIdx) const
{
    if (compactIdx >= bridgeGeneFullByCompact_.size()) {
        return static_cast<uint32_t>(-1);
    }
    return bridgeGeneFullByCompact_[compactIdx];
}

void SoloReadFeature::mergeInlineHash(SoloReadFeature &other)
{
    if (!inlineHash_ || !other.inlineHash_) {
        // Still merge non-hash sidecars below.
    } else {
        bool nonFlexHashBridge = pSolo.inlineHashMode
            && !pSolo.flexMode
            && std::getenv("STAR_SOLO_NONFLEX_HASH_BRIDGE") != nullptr;

        // Merge hash tables: iterate over source hash, add counts
        for (khiter_t iter = kh_begin(other.inlineHash_); iter != kh_end(other.inlineHash_); ++iter) {
            if (!kh_exist(other.inlineHash_, iter)) continue;

            uint64_t key = kh_key(other.inlineHash_, iter);
            uint32_t count = kh_val(other.inlineHash_, iter);

            if (nonFlexHashBridge) {
                uint32_t otherCompactCb = 0, umi24 = 0;
                uint16_t geneIdx = 0;
                unpackBridgeCgAggKey(key, &otherCompactCb, &umi24, &geneIdx);

                uint32_t wlIdx = other.bridgeCompactToWl(otherCompactCb);
                if (wlIdx == static_cast<uint32_t>(-1)) {
                    continue;
                }

                uint32_t compactCb = getOrCreateBridgeCompactCb(wlIdx);
                uint32_t fullGeneIdx = other.bridgeCompactToGene(geneIdx);
                if (fullGeneIdx == static_cast<uint32_t>(-1)) {
                    continue;
                }
                uint16_t compactGene = getOrCreateBridgeCompactGene(fullGeneIdx);
                key = packBridgeCgAggKey(compactCb, umi24, compactGene);
            }

            int absent;
            khiter_t dest_iter = kh_put(cg_agg, inlineHash_, key, &absent);
            if (absent) {
                kh_val(inlineHash_, dest_iter) = count;
            } else {
                kh_val(inlineHash_, dest_iter) += count;
            }
        }
    }
    
    // Merge readIdTracker_ if both have it
    // Note: For readId tracking, we keep ALL entries (no collision resolution needed)
    // Each readId should only appear in one thread's tracker
    if (readIdTracker_ && other.readIdTracker_) {
        for (khiter_t iter = kh_begin(other.readIdTracker_); iter != kh_end(other.readIdTracker_); ++iter) {
            if (!kh_exist(other.readIdTracker_, iter)) continue;
            
            uint32_t readId = kh_key(other.readIdTracker_, iter);
            uint64_t val = kh_val(other.readIdTracker_, iter);
            
            int absent;
            khiter_t dest_iter = kh_put(readid_cbumi, readIdTracker_, readId, &absent);
            // Should always be absent (each readId processed by one thread)
            kh_val(readIdTracker_, dest_iter) = val;
        }
    }
    
    mergePendingAmbiguous(other);
    mergeDeferredBridgeAccounting(other);
}

void SoloReadFeature::mergePendingAmbiguous(const SoloReadFeature &other)
{
    // Merge ambiguous CB structs: combine UMI counts and observations on key collision
    for (const auto &kv : other.pendingAmbiguous_) {
        ReadAlign::AmbigKey key = kv.first;
        const ExtendedAmbiguousEntry &otherEntry = kv.second;
        
        auto &entry = pendingAmbiguous_[key];
        if (entry.candidateIdx.empty()) {
            // First time seeing this ambiguous CB: copy entire entry
            entry.candidateIdx = otherEntry.candidateIdx;
            entry.cbSeq = otherEntry.cbSeq;
            entry.cbQual = otherEntry.cbQual;
            entry.umiCounts = otherEntry.umiCounts;
            entry.observations = otherEntry.observations;
        } else {
            // Merge UMI counts
            for (const auto &umiCount : otherEntry.umiCounts) {
                entry.umiCounts[umiCount.first] += umiCount.second;
            }
            // Merge observations (gene/tag/umi combinations)
            entry.observations.insert(entry.observations.end(), 
                                     otherEntry.observations.begin(), 
                                     otherEntry.observations.end());
        }
    }
}

void SoloReadFeature::mergeDeferredBridgeAccounting(const SoloReadFeature &other)
{
    if (!other.bridgeDeferredAccounting_.empty()) {
        const uint32_t candidateOffsetBase = static_cast<uint32_t>(bridgeDeferredCandidates_.size());
        bridgeDeferredCandidates_.insert(bridgeDeferredCandidates_.end(),
                                         other.bridgeDeferredCandidates_.begin(),
                                         other.bridgeDeferredCandidates_.end());
        bridgeDeferredAccounting_.reserve(bridgeDeferredAccounting_.size() + other.bridgeDeferredAccounting_.size());
        for (const auto &rec : other.bridgeDeferredAccounting_) {
            BridgeDeferredReadAccounting mergedRec = rec;
            mergedRec.candidateOffset += candidateOffsetBase;
            bridgeDeferredAccounting_.push_back(mergedRec);
        }
    }
}
