#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "Stats.h"
#include "GlobalVariables.h"
#include "SoloReadFeature.h"
#include "ErrorWarning.h"
#include "hash_shims_cpp_compat.h"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {
struct BridgeSourceDigest {
    uint64_t records = 0;
    uint64_t totalCount = 0;
    uint64_t hashSum = 0;
    uint64_t hashXor = 0;
    uint64_t hashSum2 = 0;

    void add(uint64_t h, uint32_t count)
    {
        ++records;
        totalCount += count;
        hashSum += h;
        hashXor ^= h;
        hashSum2 += h * h + 0x9e3779b97f4a7c15ull;
    }

    void merge(const BridgeSourceDigest &other)
    {
        records += other.records;
        totalCount += other.totalCount;
        hashSum += other.hashSum;
        hashXor ^= other.hashXor;
        hashSum2 += other.hashSum2;
    }
};

uint64_t bridgeSourceSplitMix64(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31);
}

uint64_t bridgeSourceTupleHash(uint64_t stageTag,
                               uint32_t a,
                               uint32_t b,
                               uint32_t c,
                               uint32_t d,
                               uint32_t count)
{
    uint64_t h = bridgeSourceSplitMix64(stageTag);
    h ^= bridgeSourceSplitMix64(static_cast<uint64_t>(a) | (static_cast<uint64_t>(b) << 32));
    h ^= bridgeSourceSplitMix64(static_cast<uint64_t>(c) | (static_cast<uint64_t>(d) << 32));
    h ^= bridgeSourceSplitMix64(static_cast<uint64_t>(count) | 0xd1b54a32d192ed03ull);
    return bridgeSourceSplitMix64(h);
}

uint64_t bridgeSourceStringHash(uint64_t stageTag, const std::string &value)
{
    uint64_t h = bridgeSourceSplitMix64(stageTag ^ static_cast<uint64_t>(value.size()));
    for (unsigned char c : value) {
        h = bridgeSourceSplitMix64(h ^ static_cast<uint64_t>(c));
    }
    return h;
}

uint64_t bridgeSourceU8VectorHash(uint64_t stageTag, const std::vector<uint8_t> &value)
{
    uint64_t h = bridgeSourceSplitMix64(stageTag ^ static_cast<uint64_t>(value.size()));
    for (uint8_t v : value) {
        h = bridgeSourceSplitMix64(h ^ static_cast<uint64_t>(v));
    }
    return h;
}

uint64_t bridgeSourceU32VectorHash(uint64_t stageTag, const std::vector<uint32_t> &value)
{
    uint64_t h = bridgeSourceSplitMix64(stageTag ^ static_cast<uint64_t>(value.size()));
    for (uint32_t v : value) {
        h = bridgeSourceSplitMix64(h ^ static_cast<uint64_t>(v));
    }
    return h;
}

uint64_t bridgeSourceDoubleVectorHash(uint64_t stageTag, const std::vector<double> &value)
{
    uint64_t h = bridgeSourceSplitMix64(stageTag ^ static_cast<uint64_t>(value.size()));
    for (double v : value) {
        uint64_t bits = 0;
        std::memcpy(&bits, &v, sizeof(bits));
        h = bridgeSourceSplitMix64(h ^ bits);
    }
    return h;
}

void addBridgeHashDigest(khash_t(cg_agg) *hash, BridgeSourceDigest &digest)
{
    if (hash == nullptr) {
        return;
    }
    for (khiter_t iter = kh_begin(hash); iter != kh_end(hash); ++iter) {
        if (!kh_exist(hash, iter)) {
            continue;
        }
        uint32_t wlCb = 0, umi24 = 0;
        uint16_t gene16 = 0;
        unpackBridgeWlUmiGeneKey(kh_key(hash, iter), &wlCb, &umi24, &gene16);
        const uint32_t count = kh_val(hash, iter);
        digest.add(bridgeSourceTupleHash(0x5352435f48415348ull, wlCb, static_cast<uint32_t>(gene16),
                                         umi24, 0, count),
                   count);
    }
}

void mergeBridgeHashInto(khash_t(cg_agg) *src, khash_t(cg_agg) *dst)
{
    if (src == nullptr || dst == nullptr) {
        return;
    }
    for (khiter_t iter = kh_begin(src); iter != kh_end(src); ++iter) {
        if (!kh_exist(src, iter)) {
            continue;
        }
        int absent = 0;
        khiter_t dest = kh_put(cg_agg, dst, kh_key(src, iter), &absent);
        if (absent) {
            kh_val(dst, dest) = kh_val(src, iter);
        } else {
            kh_val(dst, dest) += kh_val(src, iter);
        }
    }
}

void addBridgePendingDigest(const SoloReadFeature &rf, BridgeSourceDigest &umiGene, BridgeSourceDigest &candidate)
{
    for (const auto &kv : rf.pendingAmbiguous_) {
        const uint64_t ambigKey = kv.first;
        const auto &entry = kv.second;
        for (uint32_t cand : entry.candidateIdx) {
            candidate.add(bridgeSourceTupleHash(0x5352435f43414e44ull,
                                                static_cast<uint32_t>(ambigKey & 0xFFFFFFFFu),
                                                static_cast<uint32_t>(ambigKey >> 32),
                                                cand,
                                                static_cast<uint32_t>(entry.candidateIdx.size()),
                                                1),
                          1);
        }
        for (const auto &ug : entry.bridgeAmbigUmiGene_) {
            const uint64_t ugKey = ug.first;
            const uint32_t umi24 = static_cast<uint32_t>((ugKey >> 16) & 0xFFFFFFu);
            const uint32_t gene16 = static_cast<uint32_t>(ugKey & 0xFFFFu);
            const uint32_t count = ug.second;
            umiGene.add(bridgeSourceTupleHash(0x5352435f50454e44ull,
                                              static_cast<uint32_t>(ambigKey & 0xFFFFFFFFu),
                                              static_cast<uint32_t>(ambigKey >> 32),
                                              umi24,
                                              gene16,
                                              count),
                        count);
        }
    }
}

	void addBridgePendingContextDigests(const SoloReadFeature &rf,
	                                    BridgeSourceDigest &representative,
	                                    BridgeSourceDigest &pinQual,
	                                    BridgeSourceDigest &evidence,
	                                    BridgeSourceDigest &accountingCounts,
	                                    BridgeSourceDigest &accountingFlags)
	{
	    for (const auto &kv : rf.pendingAmbiguous_) {
	        const uint64_t ambigKey = kv.first;
	        const auto &entry = kv.second;
	        const uint64_t keyHash = bridgeSourceSplitMix64(static_cast<uint32_t>(ambigKey & 0xFFFFFFFFu))
	            ^ bridgeSourceSplitMix64(static_cast<uint32_t>(ambigKey >> 32));

	        uint64_t hRep = bridgeSourceSplitMix64(0x5352435f43545250ull) ^ keyHash;
	        hRep ^= bridgeSourceU32VectorHash(0x5352435f4354414eull, entry.candidateIdx);
	        hRep ^= bridgeSourceStringHash(0x5352435f43545345ull, entry.cbSeq);
	        hRep ^= bridgeSourceStringHash(0x5352435f43545155ull, entry.cbQual);
	        representative.add(bridgeSourceSplitMix64(hRep), 1);

	        uint64_t hPin = bridgeSourceSplitMix64(0x5352435f43545051ull) ^ keyHash;
	        hPin ^= bridgeSourceU32VectorHash(0x5352435f4354414eull, entry.candidateIdx);
	        hPin ^= bridgeSourceU8VectorHash(0x5352435f43545051ull, entry.bridgeAmbigPinCandQuals_);
	        pinQual.add(bridgeSourceSplitMix64(hPin), 1);

	        uint64_t hEvidence = bridgeSourceSplitMix64(0x5352435f43544556ull) ^ keyHash;
	        hEvidence ^= bridgeSourceDoubleVectorHash(0x5352435f43544d41ull, entry.cbLogLikMatch);
	        hEvidence ^= bridgeSourceDoubleVectorHash(0x5352435f43544d49ull, entry.cbLogLikMismatch);
	        hEvidence ^= bridgeSourceSplitMix64(static_cast<uint64_t>(entry.cbEvidenceReads)
	                                            | (static_cast<uint64_t>(entry.cbLogLikMatch.size()) << 32));
	        evidence.add(bridgeSourceSplitMix64(hEvidence), 1);

	        uint64_t hAccountingCounts = bridgeSourceSplitMix64(0x5352435f43544143ull) ^ keyHash;
	        hAccountingCounts ^= bridgeSourceSplitMix64(static_cast<uint64_t>(entry.bridgeAmbigGeneFeatU_)
	                                                   | (static_cast<uint64_t>(entry.bridgeAmbigGeneFeatM_) << 32));
	        hAccountingCounts ^= bridgeSourceSplitMix64(static_cast<uint64_t>(entry.bridgeAmbigReadInfoN_));
	        accountingCounts.add(bridgeSourceSplitMix64(hAccountingCounts), 1);

	        uint64_t hAccountingFlags = bridgeSourceSplitMix64(0x5352435f43544146ull) ^ keyHash;
	        hAccountingFlags ^= bridgeSourceSplitMix64(static_cast<uint64_t>(entry.bridgeAmbigGeneHaveSampleU_)
	                                                  | (static_cast<uint64_t>(entry.bridgeAmbigGeneHaveSampleM_) << 1)
	                                                  | (static_cast<uint64_t>(entry.bridgeAmbigReadInfoHaveSample_) << 2));
	        hAccountingFlags ^= bridgeSourceSplitMix64(static_cast<uint64_t>(entry.bridgeAmbigGeneSampleFlagU_)
	                                                  ^ (static_cast<uint64_t>(entry.bridgeAmbigGeneSampleFlagM_) << 21)
	                                                  ^ (static_cast<uint64_t>(entry.bridgeAmbigReadInfoSampleFlag_) << 42));
	        accountingFlags.add(bridgeSourceSplitMix64(hAccountingFlags), 1);
	    }
	}

	void addBridgePinReadDigest(const SoloReadFeature &rf, BridgeSourceDigest &digest)
	{
	    const size_t n = std::max(rf.bridgePinNreadUnique_.size(), rf.bridgePinNreadMulti_.size());
	    for (size_t cb = 0; cb < n; ++cb) {
	        const uint32_t uniqueN = cb < rf.bridgePinNreadUnique_.size() ? rf.bridgePinNreadUnique_[cb] : 0;
	        const uint32_t multiN = cb < rf.bridgePinNreadMulti_.size() ? rf.bridgePinNreadMulti_[cb] : 0;
	        if (uniqueN == 0 && multiN == 0) {
	            continue;
	        }
	        digest.add(bridgeSourceTupleHash(0x5352435f50494e52ull,
	                                         static_cast<uint32_t>(cb),
	                                         uniqueN,
	                                         multiN,
	                                         0,
	                                         uniqueN + multiN),
	                   uniqueN + multiN);
	    }
	}

	void addBridgeFeatureStatsDigest(const SoloReadFeatureStats &stats, BridgeSourceDigest &digest)
	{
	    for (uint32_t ii = 0; ii < stats.nStats; ++ii) {
	        const uint64_t value = stats.V[ii];
	        digest.add(bridgeSourceTupleHash(0x5352435f53544154ull,
	                                         ii,
	                                         static_cast<uint32_t>(value & 0xFFFFFFFFu),
	                                         static_cast<uint32_t>(value >> 32),
	                                         0,
	                                         static_cast<uint32_t>(value & 0xFFFFFFFFu)),
	                   static_cast<uint32_t>(value & 0xFFFFFFFFu));
	    }
	}

void addBridgeCbReadCountDigest(const SoloReadFeature &rf, BridgeSourceDigest &digest)
{
    for (uint32_t cb = 0; cb < rf.cbReadCount.size(); ++cb) {
        const uint32_t count = rf.cbReadCount[cb];
        if (count == 0) {
            continue;
        }
        digest.add(bridgeSourceTupleHash(0x5352435f43425243ull, cb, 0, 0, 0, count), count);
    }
}

void writeBridgeSourceDigest(const std::string &path,
                             const BridgeSourceDigest &threadHash,
	                             const BridgeSourceDigest &threadHashLogical,
	                             const BridgeSourceDigest &pendingUmiGene,
	                             const BridgeSourceDigest &pendingCandidate,
	                             const BridgeSourceDigest &pendingRepresentative,
	                             const BridgeSourceDigest &pendingPinQual,
	                             const BridgeSourceDigest &pendingEvidence,
	                             const BridgeSourceDigest &pendingAccountingCounts,
	                             const BridgeSourceDigest &pendingAccountingFlags,
	                             const BridgeSourceDigest &cbReadCount,
	                             const BridgeSourceDigest &resolvedAmbigHash,
	                             const BridgeSourceDigest &postAccountingPinReads,
	                             const BridgeSourceDigest &postAccountingStats,
	                             const std::vector<BridgeSourceDigest> &threadHashByThread,
                             uint64_t pendingAmbiguousN,
                             uint32_t nCB,
                             Parameters &P)
{
    std::ofstream out(path.c_str());
    if (!out.good()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal OUTPUT FILE error: could not open STAR_SOLO_BRIDGE_SOURCE_DETERMINISM_OUT="
               << path << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    out << "# pending_ambiguous=" << pendingAmbiguousN << "\n";
    out << "# nCB_from_cbReadCount=" << nCB << "\n";
    out << "stage\trecords\ttotal_count\thash_sum\thash_xor\thash_sum2\n";
    auto emit = [&](const char *stage, const BridgeSourceDigest &d) {
        out << stage << '\t'
            << d.records << '\t'
            << d.totalCount << '\t'
            << d.hashSum << '\t'
            << d.hashXor << '\t'
            << d.hashSum2 << '\n';
    };
    emit("thread_hash_pre_resolve", threadHash);
	    emit("thread_hash_logical_pre_resolve", threadHashLogical);
	    emit("pending_ambig_umi_gene_pre_resolve", pendingUmiGene);
	    emit("pending_ambig_candidates_pre_resolve", pendingCandidate);
	    emit("pending_ambig_context_representative_pre_resolve", pendingRepresentative);
	    emit("pending_ambig_context_pin_quals_pre_resolve", pendingPinQual);
	    emit("pending_ambig_context_evidence_pre_resolve", pendingEvidence);
	    emit("pending_ambig_context_accounting_counts_pre_resolve", pendingAccountingCounts);
	    emit("pending_ambig_context_accounting_flags_pre_resolve", pendingAccountingFlags);
	    emit("cb_read_count_pre_resolve", cbReadCount);
	    emit("resolved_ambig_hash_post_resolve", resolvedAmbigHash);
	    emit("post_accounting_pin_reads", postAccountingPinReads);
	    emit("post_accounting_feature_stats", postAccountingStats);
    out.close();

    const std::string threadPath = path + ".threads.tsv";
    std::ofstream tout(threadPath.c_str());
    if (!tout.good()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal OUTPUT FILE error: could not open bridge source per-thread digest file="
               << threadPath << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    tout << "thread\trecords\ttotal_count\thash_sum\thash_xor\thash_sum2\n";
    for (size_t ii = 0; ii < threadHashByThread.size(); ++ii) {
        const BridgeSourceDigest &d = threadHashByThread[ii];
        tout << ii << '\t'
             << d.records << '\t'
             << d.totalCount << '\t'
             << d.hashSum << '\t'
             << d.hashXor << '\t'
             << d.hashSum2 << '\n';
    }
}
} // namespace

void SoloFeature::sumThreads()
{   
    //stats
    nReadsInput=g_statsAll.readN+1; //reserve 1 extra

    auto releaseMergedThreadState = [](SoloReadFeature *rf, bool keepInlineHashAndBridgeMaps) {
        if (rf == nullptr) {
            return;
        }
        if (!keepInlineHashAndBridgeMaps && rf->inlineHash_) {
            kh_destroy(cg_agg, rf->inlineHash_);
            rf->inlineHash_ = nullptr;
        }
        if (rf->readIdTracker_) {
            kh_destroy(readid_cbumi, rf->readIdTracker_);
            rf->readIdTracker_ = nullptr;
        }
        decltype(rf->pendingAmbiguous_)().swap(rf->pendingAmbiguous_);
        decltype(rf->bridgeAmbigReadInfoOrphan_)().swap(rf->bridgeAmbigReadInfoOrphan_);
        decltype(rf->bridgeImmediateReadCounts_)().swap(rf->bridgeImmediateReadCounts_);
        rf->bridgePinNreadUnique_.clear();
        rf->bridgePinNreadMulti_.clear();
        decltype(rf->readFlag.flagCounts)().swap(rf->readFlag.flagCounts);
        rf->readFlag.flagCountsNoCB = {};
        std::vector<uint32_t>().swap(rf->cbReadCount);
        decltype(rf->cbReadCountMap)().swap(rf->cbReadCountMap);
    };

    const bool nonFlexDirectBridge = pSolo.inlineHashMode
        && !pSolo.flexMode
        && std::getenv("STAR_SOLO_NONFLEX_HASH_BRIDGE") != nullptr;
    
    ///////////////////////////// collect per-thread SoloReadFeature objects
    for (int ii=0; ii<P.runThreadN; ii++) {//point to
        if (RAchunk != nullptr) {
            if (RAchunk[ii] == nullptr || RAchunk[ii]->RA == nullptr ||
                RAchunk[ii]->RA->soloRead == nullptr ||
                RAchunk[ii]->RA->soloRead->readFeat[pSolo.featureInd[featureType]] == nullptr) {
                ostringstream errOut;
                errOut << "EXITING because of fatal ERROR: missing per-thread Solo feature for "
                       << SoloFeatureTypes::Names[featureType]
                       << " thread " << ii << "\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            }
            readFeatAll[ii]= RAchunk[ii]->RA->soloRead->readFeat[pSolo.featureInd[featureType]];
        } else if (readFeatAll[ii] == nullptr) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: no-genome Solo aggregation missing thread feature for "
                   << SoloFeatureTypes::Names[featureType]
                   << " thread " << ii << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        }
        readFeatAll[ii]->setOwner(this);
        if (readFeatAll[ii]->streamReads) {
            readFeatAll[ii]->streamReads->flush();
        }
        
        // Merge inline hash if enabled
        if (pSolo.inlineHashMode) {
            if (nonFlexDirectBridge) {
                // Keep the bulky non-ambiguous thread-local hashes for direct draining later;
                // merge only the smaller ambiguous/deferred sidecars into readFeatSum.
                readFeatSum->mergePendingAmbiguous(*readFeatAll[ii]);
            } else {
                readFeatSum->mergeInlineHash(*readFeatAll[ii]);
            }
            readFeatSum->addStats(*readFeatAll[ii]);
        }
        
        readFeatSum->addCounts(*readFeatAll[ii]);

        if (pSolo.inlineHashMode) {
            releaseMergedThreadState(readFeatAll[ii], nonFlexDirectBridge);
        }
    };       
    
    // if WL was not defined
    if (!pSolo.cbWLyes) {//now we can define WL and counts ??? we do not need to do it for every feature???
        pSolo.cbWLsize=readFeatSum->cbReadCountMap.size();
        pSolo.cbWL.resize(pSolo.cbWLsize);
        pSolo.cbWLstr.resize(pSolo.cbWLsize);
        readFeatSum->cbReadCount.resize(pSolo.cbWLsize);
        readBarSum->cbReadCountExact.resize(pSolo.cbWLsize);

        if (pSolo.CBtype.type==1) {//sequence cb
            uint64 icb=0;
            for (auto &cb : readFeatSum->cbReadCountMap) {
                pSolo.cbWL[icb] = cb.first;
                pSolo.cbWLstr[icb] = convertNuclInt64toString(pSolo.cbWL[icb],pSolo.cbL);
                readFeatSum->cbReadCount[icb]=cb.second;
                readBarSum->cbReadCountExact[icb]=cb.second;
                icb++;
            };
        } else if (pSolo.CBtype.type==2) {//string cb
            vector< std::unordered_map<string,uint32>::iterator > cbiter(pSolo.CBtype.strMap.size());
            for (auto cbi=pSolo.CBtype.strMap.begin(); cbi!=pSolo.CBtype.strMap.end(); cbi++)
                cbiter[cbi->second] = cbi;

            uint64 icb=0;
            for (auto &cb : readFeatSum->cbReadCountMap) {
                pSolo.cbWL[icb] = cb.first;
                pSolo.cbWLstr[icb] = cbiter[cb.first]->first;
                readFeatSum->cbReadCount[icb]=cb.second;
                readBarSum->cbReadCountExact[icb]=cb.second;
                icb++;
            };
        };
        pSolo.cbWLstrOut = pSolo.cbWLstr;

        //pseudocounts
        if (pSolo.CBmatchWL.mm1_multi_pc) {
            for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
                readBarSum->cbReadCountExact[ii]++;//add one to exact counts
            };
        };
    };

    // if restarting from _STARtmp/solo* file
    if (P.runRestart.type==1) {//this could happen if the run is restarted. Would be better to save/load cbReadCount, or recalculate it from
        for (int ii=0; ii<P.runThreadN; ii++) {
            if (readFeatAll[ii] != nullptr && readFeatAll[ii]->binarySpool) {
                ostringstream errOut;
                errOut << "EXITING because STAR_SOLO_BINARY_SPOOL does not support restart parsing from existing solo temp files.\n"
                       << "SOLUTION: rerun without restart files or disable STAR_SOLO_BINARY_SPOOL.\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
            }
        }
        if (pSolo.soloFlexMinimalMemory) {
            // Warn and skip restart even if inlineHashMode is false (flag will be ignored but restart still skipped)
            P.inOut->logMain << "WARNING: --soloFlexMinimalMemory is enabled; skipping restart logic that depends on stream files/packed read info" << endl;
        } else {
            // Existing restart logic (with null check for safety)
            for (int ii=0; ii<P.runThreadN; ii++) {
                if (readFeatAll[ii]->streamReads) {
                    readFeatAll[ii]->streamReads->clear(); //just in case EOF was reached in previous reading
                    readFeatAll[ii]->streamReads->seekg(0,ios::beg);
                    string line1;
                    while (std::getline(*readFeatAll[ii]->streamReads, line1)) {
                        istringstream line1stream(line1);
                        uint64 cb1;            
                        line1stream >> cb1 >> cb1 >> cb1;
                        if (featureType==SoloFeatureTypes::SJ)
                            line1stream >> cb1;
                        line1stream >> cb1;
                        //if (cb1>readFeatSum->cbReadCount.size())
                        //    continue;//this should not happen!
                        readFeatSum->cbReadCount[cb1]++;
                    };
                }
            };
        }
    };    
    
    //detected CBs
    nCB=0;nReadsMapped=0;
    for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
        if (readFeatSum->cbReadCount[ii]>0) {
            nCB++;
            nReadsMapped += readFeatSum->cbReadCount[ii];
        };
    };
    
    indCBwl.resize(pSolo.cbWLsize, (uint32) -1);
    indCB.resize(nCB);
    nCB=0;//will count it again below
    for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
        if (readFeatSum->cbReadCount[ii]>0) {
            indCB[nCB]=ii;
            indCBwl[ii]=nCB;
            ++nCB;
        };
    };

    // Non-Flex direct hash bridge: global CB support is complete after addCounts; resolve merged
    // ambiguous CB keys into readFeatSum's small hash before countCBgeneUMI
    // (thread-local bulk hashes stay untouched until collapse drain).
    if (nonFlexDirectBridge && readFeatSum) {
        const char *sourceDetPathEnv = std::getenv("STAR_SOLO_BRIDGE_SOURCE_DETERMINISM_OUT");
        const bool sourceDeterminismTrace = sourceDetPathEnv != nullptr && sourceDetPathEnv[0] != '\0';
        const std::string sourceDeterminismPath = sourceDeterminismTrace ? sourceDetPathEnv : "";
        BridgeSourceDigest sourceThreadHash;
        BridgeSourceDigest sourceThreadHashLogical;
        BridgeSourceDigest sourcePendingUmiGene;
        BridgeSourceDigest sourcePendingCandidate;
        BridgeSourceDigest sourcePendingRepresentative;
        BridgeSourceDigest sourcePendingPinQual;
        BridgeSourceDigest sourcePendingEvidence;
        BridgeSourceDigest sourcePendingAccountingCounts;
        BridgeSourceDigest sourcePendingAccountingFlags;
        BridgeSourceDigest sourceCbReadCount;
        std::vector<BridgeSourceDigest> sourceThreadHashByThread;
        uint64_t sourcePendingAmbiguousN = 0;
        if (sourceDeterminismTrace) {
            sourcePendingAmbiguousN = readFeatSum->pendingAmbiguous_.size();
            sourceThreadHashByThread.resize(static_cast<size_t>(P.runThreadN));
            size_t sourceThreadHashEntryEstimate = 0;
            for (int ii = 0; ii < P.runThreadN; ++ii) {
                if (readFeatAll[ii] != nullptr && readFeatAll[ii]->inlineHash_ != nullptr) {
                    sourceThreadHashEntryEstimate += kh_size(readFeatAll[ii]->inlineHash_);
                }
            }
            khash_t(cg_agg) *sourceMergedDirectHash = kh_init(cg_agg);
            const size_t khMax = static_cast<size_t>(std::numeric_limits<khint_t>::max());
            const size_t resizeTarget = sourceThreadHashEntryEstimate > (khMax - 1) / 2
                ? khMax
                : sourceThreadHashEntryEstimate * 2 + 1;
            if (resizeTarget > 0) {
                kh_resize(cg_agg, sourceMergedDirectHash, static_cast<khint_t>(resizeTarget));
            }
            for (int ii = 0; ii < P.runThreadN; ++ii) {
                if (readFeatAll[ii] == nullptr) {
                    continue;
                }
                addBridgeHashDigest(readFeatAll[ii]->inlineHash_, sourceThreadHashByThread[static_cast<size_t>(ii)]);
                sourceThreadHash.merge(sourceThreadHashByThread[static_cast<size_t>(ii)]);
                mergeBridgeHashInto(readFeatAll[ii]->inlineHash_, sourceMergedDirectHash);
            }
            addBridgeHashDigest(sourceMergedDirectHash, sourceThreadHashLogical);
            kh_destroy(cg_agg, sourceMergedDirectHash);
            addBridgePendingDigest(*readFeatSum, sourcePendingUmiGene, sourcePendingCandidate);
            addBridgePendingContextDigests(*readFeatSum,
                                           sourcePendingRepresentative,
                                           sourcePendingPinQual,
                                           sourcePendingEvidence,
                                           sourcePendingAccountingCounts,
                                           sourcePendingAccountingFlags);
            addBridgeCbReadCountDigest(*readFeatSum, sourceCbReadCount);
        }

        time_t rawTime;
        time(&rawTime);
        P.inOut->logMain << timeMonthDayTime(rawTime)
                         << " ... Non-Flex hash bridge: post-merge ambiguous resolve"
                         << endl;
        resolvePendingAmbiguousToHash(true);

        if (sourceDeterminismTrace) {
            BridgeSourceDigest sourceResolvedAmbigHash;
            BridgeSourceDigest sourcePostAccountingPinReads;
            BridgeSourceDigest sourcePostAccountingStats;
            addBridgeHashDigest(readFeatSum->inlineHash_, sourceResolvedAmbigHash);
            addBridgePinReadDigest(*readFeatSum, sourcePostAccountingPinReads);
            addBridgeFeatureStatsDigest(readFeatSum->stats, sourcePostAccountingStats);
            writeBridgeSourceDigest(sourceDeterminismPath,
                                    sourceThreadHash,
                                    sourceThreadHashLogical,
                                    sourcePendingUmiGene,
                                    sourcePendingCandidate,
                                    sourcePendingRepresentative,
                                    sourcePendingPinQual,
                                    sourcePendingEvidence,
                                    sourcePendingAccountingCounts,
                                    sourcePendingAccountingFlags,
                                    sourceCbReadCount,
                                    sourceResolvedAmbigHash,
                                    sourcePostAccountingPinReads,
                                    sourcePostAccountingStats,
                                    sourceThreadHashByThread,
                                    sourcePendingAmbiguousN,
                                    nCB,
                                    P);
            P.inOut->logMain << "Wrote bridge source determinism digest: "
                             << sourceDeterminismPath << endl;
        }
    }
    
};
    
