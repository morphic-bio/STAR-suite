#include "SoloReadFeature.h"
#include "ErrorWarning.h"
#include "SoloBinarySpool.h"
#include "streamFuns.h"
#include "SoloFeatureTypes.h"
#include "solo/CbBayesianResolver.h"
#include <cstdlib>
#include <cstdio>

SoloReadFeature::SoloReadFeature(int32 feTy, Parameters &Pin, int iChunk)
             : featureType(feTy), P(Pin), pSolo(P.pSolo), binarySpool(false), binarySpoolInMemory(false), binarySpoolMemoryLimitBytes(0), streamReads(nullptr), binarySpoolReadPos(0), inlineHash_(nullptr), readIdTracker_(nullptr)
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

    const bool wantBinarySpool = SoloBinarySpool::envEnabled();
    const bool wantBinarySpoolInMemory = SoloBinarySpool::envInMemoryEnabled();
    const uint64_t binarySpoolMemoryLimitBytesEnv = SoloBinarySpool::envInMemoryLimitBytes();
    if (wantBinarySpool && P.runRestart.type == 1) {
        ostringstream errOut;
        errOut << "EXITING because STAR_SOLO_BINARY_SPOOL is not compatible with restart mode in this experiment.\n"
               << "SOLUTION: rerun without restart files or disable STAR_SOLO_BINARY_SPOOL.\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    if (pSolo.inlineHashMode) {
        // Initialize inline hash instead of opening temp stream file
        inlineHash_ = kh_init(cg_agg);
        streamReads = nullptr; // Do NOT open stream file in inline hash mode
        
        // Initialize parallel readId tracker for sorted BAM CB/UB tag injection
        if (pSolo.trackReadIdsForTags) {
            readIdTracker_ = kh_init(readid_cbumi);
        }
    } else if (iChunk>=0) {
        binarySpool = wantBinarySpool && SoloBinarySpool::supportsFeature(featureType);
        binarySpoolInMemory = binarySpool && wantBinarySpoolInMemory;
        binarySpoolFileName = P.outFileTmp+"/solo"+SoloFeatureTypes::Names[featureType]+'_'+std::to_string(iChunk);
        binarySpoolMemoryLimitBytes = binarySpoolInMemory ? binarySpoolMemoryLimitBytesEnv : 0;
        if (binarySpoolInMemory) {
            streamReads = nullptr;
            SoloBinarySpool::writeFileHeader(binarySpoolBuffer, featureType, readIndexYes);
            if (iChunk == 0) {
                P.inOut->logMain << "Using experimental in-memory binary Solo temp spool for feature "
                                 << SoloFeatureTypes::Names[featureType] << endl;
                if (binarySpoolMemoryLimitBytes > 0) {
                    P.inOut->logMain << "Experimental in-memory Solo spool limit: "
                                     << (binarySpoolMemoryLimitBytes / (1024ull * 1024ull))
                                     << " MiB per thread; spilling overflow to disk" << endl;
                }
            }
        } else if (binarySpool) {
            streamReads = &fstrOpenBinary(binarySpoolFileName, ERROR_OUT, P, true);
            SoloBinarySpool::writeFileHeader(*streamReads, featureType, readIndexYes);
            if (iChunk == 0) {
                P.inOut->logMain << "Using experimental binary Solo temp spool for feature "
                                 << SoloFeatureTypes::Names[featureType] << endl;
            }
        } else {
            //open with flagDelete=false, i.e. try to keep file if it exists
            streamReads = &fstrOpen(binarySpoolFileName, ERROR_OUT, P, false);
            if (wantBinarySpool && iChunk == 0 && !SoloBinarySpool::supportsFeature(featureType)) {
                P.inOut->logMain << "WARNING: STAR_SOLO_BINARY_SPOOL is enabled, but feature "
                                 << SoloFeatureTypes::Names[featureType]
                                 << " is not supported by the experimental binary spool; using legacy text temp stream"
                                 << endl;
            }
        }
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

void SoloReadFeature::bridgeDirectTupleAdd(uint64_t tupleKey, uint32_t geneFull, uint32_t umi24, uint32_t delta)
{
    if (inlineHash_ == nullptr || delta == 0) {
        return;
    }
    if (geneFull >= (1u << 18)) {
        std::fprintf(stderr,
                     "FATAL ERROR: non-Flex bridge direct path: gene index %u exceeds 18-bit collapse key field\n",
                     geneFull);
        std::exit(1);
    }
    int absent = 0;
    khiter_t it = kh_put(cg_agg, inlineHash_, tupleKey, &absent);
    if (absent) {
        const bool over = delta > kBridgePackedSlotCountMax;
        const uint32_t cnt = over ? kBridgePackedSlotCountMax : delta;
        if (over) {
            bridgeSlotOverflowEvents_++;
        }
        const uint32_t sid = static_cast<uint32_t>(bridgePackedSlots_.size());
        bridgePackedSlots_.push_back(packBridgePackedSlot(umi24, geneFull, cnt, over));
        kh_val(inlineHash_, it) = sid;
        return;
    }
    const uint32_t sid = kh_val(inlineHash_, it);
    if (sid >= bridgePackedSlots_.size()) {
        return;
    }
    bridgePackedSlotAddCount(&bridgePackedSlots_[sid], delta, &bridgeSlotOverflowEvents_);
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
                uint32_t otherSid = count;
                if (otherSid >= other.bridgePackedSlots_.size()) {
                    continue;
                }
                const uint64_t srcPack = other.bridgePackedSlots_[otherSid];

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
                const uint64_t newKey = packBridgeCgAggKey(compactCb, umi24, compactGene);

                int absent = 0;
                khiter_t dest_iter = kh_put(cg_agg, inlineHash_, newKey, &absent);
                if (absent) {
                    const uint32_t sid = static_cast<uint32_t>(bridgePackedSlots_.size());
                    bridgePackedSlots_.push_back(srcPack);
                    kh_val(inlineHash_, dest_iter) = sid;
                } else {
                    const uint32_t sid = kh_val(inlineHash_, dest_iter);
                    if (sid < bridgePackedSlots_.size()) {
                        bridgePackedSlots_[sid] =
                            bridgePackedSlotMerge(bridgePackedSlots_[sid], srcPack, &bridgeSlotOverflowEvents_);
                    }
                }
                continue;
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
            entry.cbLogLikMatch = otherEntry.cbLogLikMatch;
            entry.cbLogLikMismatch = otherEntry.cbLogLikMismatch;
            entry.cbEvidenceReads = otherEntry.cbEvidenceReads;
            entry.umiCounts = otherEntry.umiCounts;
            entry.observations = otherEntry.observations;
        } else {
            cb_bayesian::mergeCbQualityEvidence(otherEntry.cbLogLikMatch,
                                                otherEntry.cbLogLikMismatch,
                                                otherEntry.cbEvidenceReads,
                                                entry.cbLogLikMatch,
                                                entry.cbLogLikMismatch,
                                                entry.cbEvidenceReads);
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

void SoloReadFeature::maybeSpillBinarySpool(size_t extraBytes)
{
    if (!binarySpoolInMemory || binarySpoolMemoryLimitBytes == 0) {
        return;
    }
    if (binarySpoolBuffer.totalBytes + extraBytes <= binarySpoolMemoryLimitBytes) {
        return;
    }

    streamReads = &fstrOpenBinary(binarySpoolFileName, ERROR_OUT, P, true);
    SoloBinarySpool::flushToStream(binarySpoolBuffer, *streamReads);
    streamReads->flush();
    binarySpoolBuffer.clear();
    binarySpoolReadPos = 0;
    binarySpoolInMemory = false;

    P.inOut->logMain << "Spilling experimental in-memory binary Solo temp spool to disk for feature "
                     << SoloFeatureTypes::Names[featureType]
                     << " after reaching "
                     << (binarySpoolMemoryLimitBytes / (1024ull * 1024ull))
                     << " MiB per-thread limit" << endl;
}
