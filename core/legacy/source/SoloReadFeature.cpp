#include "SoloReadFeature.h"
#include "ErrorWarning.h"
#include "SoloBinarySpool.h"
#include "streamFuns.h"
#include "SoloFeatureTypes.h"
#include "solo/CbBayesianResolver.h"
#include <algorithm>
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

    const bool nonFlexHashBridge =
        pSolo.inlineHashMode && !pSolo.flexMode
        && std::getenv("STAR_SOLO_NONFLEX_HASH_BRIDGE") != nullptr;
    const bool keepLegacyVelocytoStream =
        nonFlexHashBridge && featureType == SoloFeatureTypes::Velocyto;
    const bool useInlineHashStorage =
        pSolo.inlineHashMode && (pSolo.flexMode || nonFlexHashBridge) && !keepLegacyVelocytoStream;

    if (pSolo.bucketStoreEnabled && pSolo.flexMode && pSolo.inlineHashMode
        && featureType == SoloFeatureTypes::Gene) {
        bucketStore_ = pSolo.cbBucketStore;
        bucketSegments_.resize(pSolo.bucketCount);
        bucketWorkerIndex_ = iChunk >= 0
            ? static_cast<uint32_t>(iChunk)
            : static_cast<uint32_t>(-iChunk);
    }

    if (useInlineHashStorage) {
        // Bucket mode replaces the per-thread fused Gene hash. Other features
        // and the off path keep the established hash storage unchanged.
        if (!bucketStorageEnabled())
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

void SoloReadFeature::appendInlineObservation(uint64_t key, uint32_t value)
{
    if (bucketStore_) {
        const uint32_t cb = static_cast<uint32_t>((key >> 44) & 0xFFFFFu);
        uint32_t bucket = 0;
        try {
            bucket = bucketStore_->bucket_for_cb(cb);
        } catch (const std::exception &error) {
            ostringstream errOut;
            errOut << "EXITING because a fused Flex record could not be bucketed: "
                   << error.what() << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain,
                          EXIT_CODE_INCONSISTENT_DATA, P);
        }
        std::vector<star::solo::PackedCbRecord> &segment = bucketSegments_[bucket];
        star::solo::PackedCbRecord record;
        record.key = key;
        record.value = value;
        segment.push_back(record);
        // Large enough to amortize spill pwrite overhead while keeping bounded
        // producer staging (4096 records is 64 KiB in the in-memory vector).
        if (segment.size() >= 4096) {
            std::string error;
            std::vector<star::solo::PackedCbRecord> sealed;
            sealed.swap(segment);
            if (!bucketStore_->append_segment(bucketWorkerIndex_, bucket,
                                              std::move(sealed), &error)) {
                exitWithError("EXITING because a CB bucket append failed: " + error + "\n",
                              std::cerr, P.inOut->logMain,
                              EXIT_CODE_INCONSISTENT_DATA, P);
            }
        }
        return;
    }

    if (!inlineHash_)
        return;
    int absent;
    const khiter_t iter = kh_put(cg_agg, inlineHash_, key, &absent);
    if (absent) {
        kh_val(inlineHash_, iter) = value;
    } else {
        kh_val(inlineHash_, iter) = pSolo.flexMode
            ? flexGdnaMergeValue(kh_val(inlineHash_, iter), value)
            : kh_val(inlineHash_, iter) + value;
    }
}

void SoloReadFeature::flushBucketSegments()
{
    if (!bucketStore_)
        return;
    for (uint32_t bucket = 0; bucket < bucketSegments_.size(); ++bucket) {
        if (bucketSegments_[bucket].empty())
            continue;
        std::vector<star::solo::PackedCbRecord> sealed;
        sealed.swap(bucketSegments_[bucket]);
        std::string error;
        if (!bucketStore_->append_segment(bucketWorkerIndex_, bucket,
                                          std::move(sealed), &error)) {
            exitWithError("EXITING because a CB bucket append failed: " + error + "\n",
                          std::cerr, P.inOut->logMain,
                          EXIT_CODE_INCONSISTENT_DATA, P);
        }
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

void SoloReadFeature::mergeInlineHash(SoloReadFeature &other)
{
    other.flushBucketSegments();
    if (!inlineHash_ || !other.inlineHash_) {
        // Still merge non-hash sidecars below.
    } else {
        // Merge hash tables: iterate over source hash, add counts
        for (khiter_t iter = kh_begin(other.inlineHash_); iter != kh_end(other.inlineHash_); ++iter) {
            if (!kh_exist(other.inlineHash_, iter)) continue;

            uint64_t key = kh_key(other.inlineHash_, iter);
            uint32_t count = kh_val(other.inlineHash_, iter);

            int absent;
            khiter_t dest_iter = kh_put(cg_agg, inlineHash_, key, &absent);
            if (absent) {
                kh_val(inlineHash_, dest_iter) = count;
            } else {
                kh_val(inlineHash_, dest_iter) = pSolo.flexMode
                    ? flexGdnaMergeValue(kh_val(inlineHash_, dest_iter), count)
                    : kh_val(inlineHash_, dest_iter) + count;
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
    
    mergePendingAmbiguous(
        other, bucketStorageEnabled() && other.bucketStorageEnabled());
}

namespace {
int64_t bridgeCbQualTotalScore(const std::string &qual, char qsBase, uint32_t qsMax)
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

template <typename T>
bool bridgeVectorLess(const std::vector<T> &lhs, const std::vector<T> &rhs)
{
    return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

bool bridgeRepresentativeBetter(const std::string &newSeq,
                                const std::string &newQual,
                                const std::vector<uint8_t> &newPinCandQuals,
                                const std::vector<uint32_t> &newCandidateIdx,
                                const std::string &oldSeq,
                                const std::string &oldQual,
                                const std::vector<uint8_t> &oldPinCandQuals,
                                const std::vector<uint32_t> &oldCandidateIdx,
                                char qsBase,
                                uint32_t qsMax)
{
    const int64_t oldScore = bridgeCbQualTotalScore(oldQual, qsBase, qsMax);
    const int64_t newScore = bridgeCbQualTotalScore(newQual, qsBase, qsMax);
    if (newScore != oldScore) {
        return newScore > oldScore;
    }
    if (newQual != oldQual) {
        return newQual < oldQual;
    }
    if (newSeq != oldSeq) {
        return newSeq < oldSeq;
    }
    if (newPinCandQuals.empty() != oldPinCandQuals.empty()) {
        return !newPinCandQuals.empty();
    }
    if (newPinCandQuals != oldPinCandQuals) {
        return bridgeVectorLess(newPinCandQuals, oldPinCandQuals);
    }
    return bridgeVectorLess(newCandidateIdx, oldCandidateIdx);
}

bool bridgeSampleFlagBetter(bool oldHave,
                            SoloReadFlagClass::typeFlag oldFlag,
                            bool newHave,
                            SoloReadFlagClass::typeFlag newFlag)
{
    return newHave && (!oldHave || newFlag < oldFlag);
}

void bridgeMergePinCandQualsMax(std::vector<uint8_t> &dst,
                                const std::vector<uint8_t> &src,
                                char qsBase)
{
    if (src.empty()) {
        return;
    }
    if (dst.size() < src.size()) {
        dst.resize(src.size(), static_cast<uint8_t>(qsBase));
    }
    for (size_t ii = 0; ii < src.size(); ++ii) {
        if (src[ii] > dst[ii]) {
            dst[ii] = src[ii];
        }
    }
}
} // namespace

void SoloReadFeature::mergePendingAmbiguous(SoloReadFeature &other,
                                            bool takeOwnership)
{
    Parameters &P = other.P;
    const char qsBase = P.pSolo.QSbase;
    const uint32_t qsMax = P.pSolo.QSmax;

    for (auto &kv : other.pendingAmbiguous_) {
        ReadAlign::AmbigKey key = kv.first;
        ExtendedAmbiguousEntry &otherEntry = kv.second;

        auto found = pendingAmbiguous_.find(key);
        if (found == pendingAmbiguous_.end()) {
            if (takeOwnership) {
                pendingAmbiguous_.emplace(key, std::move(otherEntry));
            } else {
                ExtendedAmbiguousEntry &entry = pendingAmbiguous_[key];
                entry.candidateIdx = otherEntry.candidateIdx;
                entry.cbSeq = otherEntry.cbSeq;
                entry.cbQual = otherEntry.cbQual;
                entry.umiCounts = otherEntry.umiCounts;
                entry.observations = otherEntry.observations;
                entry.bridgeAmbigUmiGene_ = otherEntry.bridgeAmbigUmiGene_;
                entry.cbLogLikMatch = otherEntry.cbLogLikMatch;
                entry.cbLogLikMismatch = otherEntry.cbLogLikMismatch;
                entry.cbEvidenceReads = otherEntry.cbEvidenceReads;
                entry.bridgeAmbigGeneFeatU_ = otherEntry.bridgeAmbigGeneFeatU_;
                entry.bridgeAmbigGeneFeatM_ = otherEntry.bridgeAmbigGeneFeatM_;
                entry.bridgeAmbigGeneHaveSampleU_ = otherEntry.bridgeAmbigGeneHaveSampleU_;
                entry.bridgeAmbigGeneHaveSampleM_ = otherEntry.bridgeAmbigGeneHaveSampleM_;
                entry.bridgeAmbigGeneSampleFlagU_ = otherEntry.bridgeAmbigGeneSampleFlagU_;
                entry.bridgeAmbigGeneSampleFlagM_ = otherEntry.bridgeAmbigGeneSampleFlagM_;
                entry.bridgeAmbigReadInfoN_ = otherEntry.bridgeAmbigReadInfoN_;
                entry.bridgeAmbigReadInfoHaveSample_ = otherEntry.bridgeAmbigReadInfoHaveSample_;
                entry.bridgeAmbigReadInfoSampleFlag_ = otherEntry.bridgeAmbigReadInfoSampleFlag_;
                entry.bridgeAmbigPinCandQuals_ = otherEntry.bridgeAmbigPinCandQuals_;
            }
        } else {
            ExtendedAmbiguousEntry &entry = found->second;
            for (const auto &umiCount : otherEntry.umiCounts) {
                entry.umiCounts[umiCount.first] += umiCount.second;
            }
            entry.observations.insert(entry.observations.end(),
                                      otherEntry.observations.begin(),
                                      otherEntry.observations.end());
            for (const auto &ug : otherEntry.bridgeAmbigUmiGene_) {
                entry.bridgeAmbigUmiGene_[ug.first] += ug.second;
            }
            cb_bayesian::mergeCbQualityEvidence(otherEntry.cbLogLikMatch,
                                                otherEntry.cbLogLikMismatch,
                                                otherEntry.cbEvidenceReads,
                                                entry.cbLogLikMatch,
                                                entry.cbLogLikMismatch,
                                                entry.cbEvidenceReads);
            bridgeMergePinCandQualsMax(entry.bridgeAmbigPinCandQuals_,
                                       otherEntry.bridgeAmbigPinCandQuals_,
                                       qsBase);
            entry.bridgeAmbigGeneFeatU_ += otherEntry.bridgeAmbigGeneFeatU_;
            entry.bridgeAmbigGeneFeatM_ += otherEntry.bridgeAmbigGeneFeatM_;
            if (bridgeSampleFlagBetter(entry.bridgeAmbigGeneHaveSampleU_,
                                       entry.bridgeAmbigGeneSampleFlagU_,
                                       otherEntry.bridgeAmbigGeneHaveSampleU_,
                                       otherEntry.bridgeAmbigGeneSampleFlagU_)) {
                entry.bridgeAmbigGeneSampleFlagU_ = otherEntry.bridgeAmbigGeneSampleFlagU_;
                entry.bridgeAmbigGeneHaveSampleU_ = true;
            }
            if (bridgeSampleFlagBetter(entry.bridgeAmbigGeneHaveSampleM_,
                                       entry.bridgeAmbigGeneSampleFlagM_,
                                       otherEntry.bridgeAmbigGeneHaveSampleM_,
                                       otherEntry.bridgeAmbigGeneSampleFlagM_)) {
                entry.bridgeAmbigGeneSampleFlagM_ = otherEntry.bridgeAmbigGeneSampleFlagM_;
                entry.bridgeAmbigGeneHaveSampleM_ = true;
            }
            entry.bridgeAmbigReadInfoN_ += otherEntry.bridgeAmbigReadInfoN_;
            if (bridgeSampleFlagBetter(entry.bridgeAmbigReadInfoHaveSample_,
                                       entry.bridgeAmbigReadInfoSampleFlag_,
                                       otherEntry.bridgeAmbigReadInfoHaveSample_,
                                       otherEntry.bridgeAmbigReadInfoSampleFlag_)) {
                entry.bridgeAmbigReadInfoSampleFlag_ = otherEntry.bridgeAmbigReadInfoSampleFlag_;
                entry.bridgeAmbigReadInfoHaveSample_ = true;
            }

            if (bridgeRepresentativeBetter(otherEntry.cbSeq,
                                           otherEntry.cbQual,
                                           otherEntry.bridgeAmbigPinCandQuals_,
                                           otherEntry.candidateIdx,
                                           entry.cbSeq,
                                           entry.cbQual,
                                           entry.bridgeAmbigPinCandQuals_,
                                           entry.candidateIdx,
                                           qsBase,
                                           qsMax)) {
                entry.cbQual = otherEntry.cbQual;
                entry.cbSeq = otherEntry.cbSeq;
            }
        }
    }

    for (auto &okv : other.bridgeAmbigReadInfoOrphan_) {
        ReadAlign::AmbigKey okey = okv.first;
        SoloReadFeature::BridgeAmbigReadInfoOrphanEntry &oe = okv.second;

        auto pit = pendingAmbiguous_.find(okey);
        if (pit != pendingAmbiguous_.end() && !pit->second.candidateIdx.empty()) {
            auto &pent = pit->second;
            pent.bridgeAmbigReadInfoN_ += oe.readInfoN_;
            if (bridgeSampleFlagBetter(pent.bridgeAmbigReadInfoHaveSample_,
                                       pent.bridgeAmbigReadInfoSampleFlag_,
                                       oe.haveSample_,
                                       oe.sampleFlag_)) {
                pent.bridgeAmbigReadInfoSampleFlag_ = oe.sampleFlag_;
                pent.bridgeAmbigReadInfoHaveSample_ = true;
            }
            bridgeMergePinCandQualsMax(pent.bridgeAmbigPinCandQuals_, oe.pinCandQuals_, qsBase);
            if (bridgeRepresentativeBetter(oe.cbSeq,
                                           oe.cbQual,
                                           oe.pinCandQuals_,
                                           oe.candidateIdx,
                                           pent.cbSeq,
                                           pent.cbQual,
                                           pent.bridgeAmbigPinCandQuals_,
                                           pent.candidateIdx,
                                           qsBase,
                                           qsMax)) {
                pent.cbQual = oe.cbQual;
                pent.cbSeq = oe.cbSeq;
            }
            continue;
        }

        auto found = bridgeAmbigReadInfoOrphan_.find(okey);
        if (found == bridgeAmbigReadInfoOrphan_.end()) {
            if (takeOwnership) {
                bridgeAmbigReadInfoOrphan_.emplace(okey, std::move(oe));
            } else {
                bridgeAmbigReadInfoOrphan_[okey] = oe;
            }
        } else {
            SoloReadFeature::BridgeAmbigReadInfoOrphanEntry &dest = found->second;
            dest.readInfoN_ += oe.readInfoN_;
            if (bridgeSampleFlagBetter(dest.haveSample_,
                                       dest.sampleFlag_,
                                       oe.haveSample_,
                                       oe.sampleFlag_)) {
                dest.sampleFlag_ = oe.sampleFlag_;
                dest.haveSample_ = true;
            }
            if (bridgeRepresentativeBetter(oe.cbSeq,
                                           oe.cbQual,
                                           oe.pinCandQuals_,
                                           oe.candidateIdx,
                                           dest.cbSeq,
                                           dest.cbQual,
                                           dest.pinCandQuals_,
                                           dest.candidateIdx,
                                           qsBase,
                                           qsMax)) {
                dest.cbQual = oe.cbQual;
                dest.cbSeq = oe.cbSeq;
                dest.candidateIdx = oe.candidateIdx;
            }
            bridgeMergePinCandQualsMax(dest.pinCandQuals_, oe.pinCandQuals_, qsBase);
        }
    }
}

namespace {
void bridgeBulkAddFlagCounts(SoloReadFlagClass &rf,
                             uint32_t cb,
                             SoloReadFlagClass::typeFlag f,
                             uint64_t n)
{
    if (n == 0) {
        return;
    }
    auto cbInserted = rf.flagCounts.insert({static_cast<uintCB>(cb), {}});
    for (uint32_t ibit = 0; ibit < SoloReadFlagClass::nBits; ++ibit) {
        (*cbInserted.first).second[ibit]
            += n * static_cast<uint64_t>((f >> ibit) & SoloReadFlagClass::typeFlag{1});
    }
}

void bridgeBulkAddFlagCountsNoCB(SoloReadFlagClass &rf, SoloReadFlagClass::typeFlag f, uint64_t n)
{
    if (n == 0) {
        return;
    }
    for (uint32_t ibit = 0; ibit < SoloReadFlagClass::nBits; ++ibit) {
        rf.flagCountsNoCB[ibit] += n * static_cast<uint64_t>((f >> ibit) & SoloReadFlagClass::typeFlag{1});
    }
}

SoloReadFlagClass::typeFlag bridgeGeneFinalReadFlag(const SoloReadFlagClass &proto,
                                                    SoloReadFlagClass::typeFlag base,
                                                    bool pinOk)
{
    SoloReadFlagClass::typeFlag f = base;
    if (pinOk) {
        if ((f >> proto.featureU) & SoloReadFlagClass::typeFlag{1}) {
            f |= (SoloReadFlagClass::typeFlag{1} << proto.countedU);
        }
        if ((f >> proto.featureM) & SoloReadFlagClass::typeFlag{1}) {
            f |= (SoloReadFlagClass::typeFlag{1} << proto.countedM);
        }
    }
    f |= (SoloReadFlagClass::typeFlag{1} << proto.cbMatch);
    if (pinOk) {
        f |= (SoloReadFlagClass::typeFlag{1} << proto.cbMMmultiple);
    }
    return f;
}

SoloReadFlagClass::typeFlag bridgeReadInfoFinalReadFlag(const SoloReadFlagClass &proto,
                                                        SoloReadFlagClass::typeFlag base,
                                                        bool pinOk)
{
    SoloReadFlagClass::typeFlag f = base;
    f |= (SoloReadFlagClass::typeFlag{1} << proto.cbMatch);
    if (pinOk) {
        f |= (SoloReadFlagClass::typeFlag{1} << proto.cbMMmultiple);
    }
    return f;
}
} // namespace

void SoloReadFeature::applyBridgeAmbiguousAggregatedReadAccounting(Parameters &P,
                                                                   int32 featureType,
                                                                   const ExtendedAmbiguousEntry &entry,
                                                                   bool bayesResolved,
                                                                   uint32_t resolvedCbIdx0)
{
    const size_t wlN = P.pSolo.cbWLsize;
    if (bridgePinNreadUnique_.size() < wlN) {
        bridgePinNreadUnique_.assign(wlN, 0);
        bridgePinNreadMulti_.assign(wlN, 0);
    }

    const bool readStatsEnabled = P.pSolo.readStatsYes[featureType];
    const uint32_t geneFeatU = entry.bridgeAmbigGeneFeatU_;
    const uint32_t geneFeatM = entry.bridgeAmbigGeneFeatM_;
    const uint32_t readInfoN = entry.bridgeAmbigReadInfoN_;

    const uint32_t accountCb = bayesResolved ? resolvedCbIdx0 : 0;
    const bool pinOk = bayesResolved && accountCb < bridgePinNreadUnique_.size();

    if (!pinOk) {
        const uint64_t nGeneAmbig = static_cast<uint64_t>(geneFeatU) + static_cast<uint64_t>(geneFeatM);
        if (nGeneAmbig > 0) {
            stats.V[stats.noTooManyWLmatches] += nGeneAmbig;
        }
    }

    if (pinOk && accountCb < bridgePinNreadUnique_.size()) {
        bridgePinNreadUnique_[accountCb] += geneFeatU;
        bridgePinNreadMulti_[accountCb] += geneFeatM;
    }

    if (!readStatsEnabled) {
        return;
    }

    const SoloReadFlagClass &proto = readFlag;
    const bool cbStatsOk = pinOk && accountCb < bridgePinNreadUnique_.size();

    if (geneFeatU > 0 && entry.bridgeAmbigGeneHaveSampleU_) {
        const SoloReadFlagClass::typeFlag fU = bridgeGeneFinalReadFlag(proto, entry.bridgeAmbigGeneSampleFlagU_, pinOk);
        if (cbStatsOk) {
            bridgeBulkAddFlagCounts(readFlag, accountCb, fU, geneFeatU);
        } else {
            bridgeBulkAddFlagCountsNoCB(readFlag, fU, geneFeatU);
        }
    }
    if (geneFeatM > 0 && entry.bridgeAmbigGeneHaveSampleM_) {
        const SoloReadFlagClass::typeFlag fM = bridgeGeneFinalReadFlag(proto, entry.bridgeAmbigGeneSampleFlagM_, pinOk);
        if (cbStatsOk) {
            bridgeBulkAddFlagCounts(readFlag, accountCb, fM, geneFeatM);
        } else {
            bridgeBulkAddFlagCountsNoCB(readFlag, fM, geneFeatM);
        }
    }

    if (readInfoN > 0 && entry.bridgeAmbigReadInfoHaveSample_) {
        const SoloReadFlagClass::typeFlag fR = bridgeReadInfoFinalReadFlag(proto, entry.bridgeAmbigReadInfoSampleFlag_, pinOk);
        if (cbStatsOk) {
            bridgeBulkAddFlagCounts(readFlag, accountCb, fR, readInfoN);
        } else {
            bridgeBulkAddFlagCountsNoCB(readFlag, fR, readInfoN);
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

// ---------------------------------------------------------------------------
// Shared sharded ambiguous-barcode accumulation.
//
// See the comment on SoloReadFeature::AmbigShard. One striped structure that
// every fused thread writes into, so an ambiguous barcode seen by many threads
// is accumulated in one place instead of being reconciled field by field
// afterwards. Enabled only for the fused Flex path that has the expensive fold.
// ---------------------------------------------------------------------------
namespace {
SoloReadFeature::AmbigShard g_ambigShards[SoloReadFeature::kAmbigShardCount];
bool g_ambigSharedActive = false;

inline size_t ambigShardIndex(ReadAlign::AmbigKey key)
{
    // The key is already a hash of the CB sequence; mix once more so the low
    // bits used for striping are well distributed.
    uint64_t h = key;
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    return static_cast<size_t>(h & (SoloReadFeature::kAmbigShardCount - 1));
}
} // namespace

bool SoloReadFeature::sharedAmbigActive()
{
    return g_ambigSharedActive;
}

void SoloReadFeature::sharedAmbigEnable(bool on)
{
    g_ambigSharedActive = on;
}

SoloReadFeature::AmbigEntryRef SoloReadFeature::sharedAmbigEntry(ReadAlign::AmbigKey key)
{
    AmbigShard &shard = g_ambigShards[ambigShardIndex(key)];
    AmbigEntryRef ref;
    ref.lock = std::unique_lock<std::mutex>(shard.mutex);
    ref.entry = &shard.map[key];
    return ref;
}

size_t SoloReadFeature::sharedAmbigSize()
{
    size_t total = 0;
    for (size_t ii = 0; ii < kAmbigShardCount; ++ii) {
        std::lock_guard<std::mutex> lock(g_ambigShards[ii].mutex);
        total += g_ambigShards[ii].map.size();
    }
    return total;
}

void SoloReadFeature::sharedAmbigDrainInto(
    std::unordered_map<ReadAlign::AmbigKey, ExtendedAmbiguousEntry> &dest)
{
    // Stripes are disjoint by key, so this is a move of already-accumulated
    // entries: no field-wise reconciliation, which is the whole point.
    size_t total = 0;
    for (size_t ii = 0; ii < kAmbigShardCount; ++ii)
        total += g_ambigShards[ii].map.size();
    if (total == 0)
        return;
    dest.reserve(dest.size() + total);
    for (size_t ii = 0; ii < kAmbigShardCount; ++ii) {
        auto &src = g_ambigShards[ii].map;
        for (auto &kv : src) {
            const auto inserted = dest.emplace(kv.first, std::move(kv.second));
            if (!inserted.second) {
                // Unreachable by construction: when the shared store is active
                // it is the only writer, so the per-thread maps are empty and
                // `dest` cannot already hold this key. Reconciling here would
                // need the field-wise merge, and silently keeping one side
                // would corrupt counts, so fail rather than guess.
                std::cerr << "EXITING because a shared ambiguous-barcode key was "
                             "already present during drain (key "
                          << kv.first << "); per-thread and shared accumulation "
                             "must not both be active\n";
                exit(EXIT_CODE_INCONSISTENT_DATA);
            }
        }
        decltype(g_ambigShards[ii].map)().swap(src);
    }
}
