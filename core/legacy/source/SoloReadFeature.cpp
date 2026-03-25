#include "SoloReadFeature.h"
#include "ErrorWarning.h"
#include "SoloBinarySpool.h"
#include "streamFuns.h"
#include "SoloFeatureTypes.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>

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

void SoloReadFeature::mergeInlineHash(SoloReadFeature &other)
{
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
} // namespace

void SoloReadFeature::mergePendingAmbiguous(const SoloReadFeature &other)
{
    Parameters &P = other.P;
    const char qsBase = P.pSolo.QSbase;
    const uint32_t qsMax = P.pSolo.QSmax;

    for (const auto &kv : other.pendingAmbiguous_) {
        ReadAlign::AmbigKey key = kv.first;
        const ExtendedAmbiguousEntry &otherEntry = kv.second;

        auto &entry = pendingAmbiguous_[key];
        if (entry.candidateIdx.empty()) {
            entry.candidateIdx = otherEntry.candidateIdx;
            entry.cbSeq = otherEntry.cbSeq;
            entry.cbQual = otherEntry.cbQual;
            entry.umiCounts = otherEntry.umiCounts;
            entry.observations = otherEntry.observations;
            entry.bridgeAmbigUmiGene_ = otherEntry.bridgeAmbigUmiGene_;
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
        } else {
            for (const auto &umiCount : otherEntry.umiCounts) {
                entry.umiCounts[umiCount.first] += umiCount.second;
            }
            entry.observations.insert(entry.observations.end(),
                                      otherEntry.observations.begin(),
                                      otherEntry.observations.end());
            for (const auto &ug : otherEntry.bridgeAmbigUmiGene_) {
                entry.bridgeAmbigUmiGene_[ug.first] += ug.second;
            }
            entry.bridgeAmbigGeneFeatU_ += otherEntry.bridgeAmbigGeneFeatU_;
            entry.bridgeAmbigGeneFeatM_ += otherEntry.bridgeAmbigGeneFeatM_;
            if (!entry.bridgeAmbigGeneHaveSampleU_ && otherEntry.bridgeAmbigGeneHaveSampleU_) {
                entry.bridgeAmbigGeneSampleFlagU_ = otherEntry.bridgeAmbigGeneSampleFlagU_;
                entry.bridgeAmbigGeneHaveSampleU_ = true;
            }
            if (!entry.bridgeAmbigGeneHaveSampleM_ && otherEntry.bridgeAmbigGeneHaveSampleM_) {
                entry.bridgeAmbigGeneSampleFlagM_ = otherEntry.bridgeAmbigGeneSampleFlagM_;
                entry.bridgeAmbigGeneHaveSampleM_ = true;
            }
            entry.bridgeAmbigReadInfoN_ += otherEntry.bridgeAmbigReadInfoN_;
            if (!entry.bridgeAmbigReadInfoHaveSample_ && otherEntry.bridgeAmbigReadInfoHaveSample_) {
                entry.bridgeAmbigReadInfoSampleFlag_ = otherEntry.bridgeAmbigReadInfoSampleFlag_;
                entry.bridgeAmbigReadInfoHaveSample_ = true;
            }

            const int64_t oldScore = bridgeCbQualTotalScore(entry.cbQual, qsBase, qsMax);
            const int64_t newScore = bridgeCbQualTotalScore(otherEntry.cbQual, qsBase, qsMax);
            if (newScore > oldScore
                || (newScore == oldScore && otherEntry.cbQual < entry.cbQual)) {
                entry.cbQual = otherEntry.cbQual;
                entry.cbSeq = otherEntry.cbSeq;
                entry.bridgeAmbigPinCandQuals_ = otherEntry.bridgeAmbigPinCandQuals_;
            } else if (entry.bridgeAmbigPinCandQuals_.empty()
                       && !otherEntry.bridgeAmbigPinCandQuals_.empty()) {
                entry.bridgeAmbigPinCandQuals_ = otherEntry.bridgeAmbigPinCandQuals_;
            }
        }
    }

    for (const auto &okv : other.bridgeAmbigReadInfoOrphan_) {
        ReadAlign::AmbigKey okey = okv.first;
        const SoloReadFeature::BridgeAmbigReadInfoOrphanEntry &oe = okv.second;

        auto pit = pendingAmbiguous_.find(okey);
        if (pit != pendingAmbiguous_.end() && !pit->second.candidateIdx.empty()) {
            auto &pent = pit->second;
            pent.bridgeAmbigReadInfoN_ += oe.readInfoN_;
            if (!pent.bridgeAmbigReadInfoHaveSample_ && oe.haveSample_) {
                pent.bridgeAmbigReadInfoSampleFlag_ = oe.sampleFlag_;
                pent.bridgeAmbigReadInfoHaveSample_ = true;
            }
            const int64_t oldScore = bridgeCbQualTotalScore(pent.cbQual, qsBase, qsMax);
            const int64_t newScore = bridgeCbQualTotalScore(oe.cbQual, qsBase, qsMax);
            if (newScore > oldScore || (newScore == oldScore && oe.cbQual < pent.cbQual)) {
                pent.cbQual = oe.cbQual;
                pent.cbSeq = oe.cbSeq;
                pent.bridgeAmbigPinCandQuals_ = oe.pinCandQuals_;
            }
            continue;
        }

        auto &dest = bridgeAmbigReadInfoOrphan_[okey];
        if (dest.candidateIdx.empty()) {
            dest = oe;
        } else {
            dest.readInfoN_ += oe.readInfoN_;
            if (!dest.haveSample_ && oe.haveSample_) {
                dest.sampleFlag_ = oe.sampleFlag_;
                dest.haveSample_ = true;
            }
            const int64_t oldScore = bridgeCbQualTotalScore(dest.cbQual, qsBase, qsMax);
            const int64_t newScore = bridgeCbQualTotalScore(oe.cbQual, qsBase, qsMax);
            if (newScore > oldScore || (newScore == oldScore && oe.cbQual < dest.cbQual)) {
                dest.cbQual = oe.cbQual;
                dest.cbSeq = oe.cbSeq;
                dest.candidateIdx = oe.candidateIdx;
                dest.pinCandQuals_ = oe.pinCandQuals_;
            }
        }
    }
}

namespace {
void bridgeEvalAmbigKeyPin(const SoloReadFeature::ExtendedAmbiguousEntry &entry,
                           const std::vector<uint32_t> &cbReadCount,
                           const Parameters &P,
                           uint32_t *outPinCb,
                           bool *outPinOk)
{
    const uint32_t kCand = static_cast<uint32_t>(entry.candidateIdx.size());
    const char qsBase = P.pSolo.QSbase;
    const uint32_t qsMax = P.pSolo.QSmax;
    uint32_t cb = 0;
#ifdef MATCH_CellRanger
    double ptot = 0.0, pmax = 0.0, pin;
#else
    float ptot = 0.0, pmax = 0.0, pin;
#endif
    for (uint32_t ii = 0; ii < kCand; ++ii) {
        const uint32_t cbinOneBased = entry.candidateIdx[ii];
        if (cbinOneBased == 0) {
            continue;
        }
        const uint32_t cbin = cbinOneBased - 1u;
        const uint8_t qraw = ii < entry.bridgeAmbigPinCandQuals_.size()
            ? entry.bridgeAmbigPinCandQuals_[ii]
            : static_cast<uint8_t>(qsBase);
        char qin = static_cast<char>(qraw);
        if (cbin < cbReadCount.size() && cbReadCount[cbin] > 0) {
            qin -= qsBase;
            qin = qin < static_cast<char>(qsMax) ? qin : static_cast<char>(qsMax);
#ifdef MATCH_CellRanger
            pin = static_cast<double>(cbReadCount[cbin])
                * std::pow(10.0, -static_cast<double>(qin) / 10.0);
#else
            pin = static_cast<float>(cbReadCount[cbin])
                * std::pow(10.0f, -static_cast<float>(qin) / 10.0f);
#endif
            ptot += pin;
            if (pin > pmax) {
                cb = cbin;
                pmax = pin;
            }
        }
    }
    const bool pinOk = (ptot > 0.0 && pmax >= P.pSolo.cbMinP * ptot);
    *outPinCb = cb;
    *outPinOk = pinOk;
}

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
    uint32_t pinCb0 = 0;
    bool pinOk = false;
    bridgeEvalAmbigKeyPin(entry, cbReadCount, P, &pinCb0, &pinOk);

    const size_t wlN = P.pSolo.cbWLsize;
    if (bridgePinNreadUnique_.size() < wlN) {
        bridgePinNreadUnique_.assign(wlN, 0);
        bridgePinNreadMulti_.assign(wlN, 0);
    }

    const bool readStatsEnabled = P.pSolo.readStatsYes[featureType];
    const uint32_t geneFeatU = entry.bridgeAmbigGeneFeatU_;
    const uint32_t geneFeatM = entry.bridgeAmbigGeneFeatM_;
    const uint32_t readInfoN = entry.bridgeAmbigReadInfoN_;

    const uint32_t accountCb = (bayesResolved && pinOk) ? resolvedCbIdx0 : pinCb0;

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
