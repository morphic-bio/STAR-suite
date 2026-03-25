#ifndef H_SoloReadFeature
#define H_SoloReadFeature
#include <set>
#include <map>
#include <unordered_map>
#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "SoloReadBarcode.h"
#include "SoloCommon.h"
#include "SoloReadFeatureStats.h"
#include "ReadAnnotations.h"
#include "ReadAlign.h"
#include "hash_shims_cpp_compat.h"
#include "SoloBinarySpool.h"
#include <functional>

class SoloFeature;
class ProbeListIndex;
class ReadSoloFeatures;

class SoloReadFeature {
public:
    // Owning SoloFeature set post-construction (not owned)
    void setOwner(SoloFeature *ownerIn) { owner = ownerIn; }
    SoloFeature* getOwner() const { return owner; }

    uint32 homoPolymer[4];//homopolymer constants

    vector<uint32> cbReadCount;
    map <uintCB,uint32> cbReadCountMap;
    
    vector<uint32> transcriptDistCount;
    
    bool readInfoYes ,readIndexYes;
    bool binarySpool;
    bool binarySpoolInMemory;
    uint64_t binarySpoolMemoryLimitBytes;
    string binarySpoolFileName;

    fstream *streamReads;
    SoloBinarySpool::MemoryBuffer binarySpoolBuffer;
    size_t binarySpoolReadPos;

    // Inline hash mode: per-thread hash table (replaces temp stream files)
    khash_t(cg_agg) *inlineHash_; // nullptr if not using inline hash mode
    
    // Parallel readId tracker for sorted BAM CB/UB tag injection (Option C)
    // Maps readId -> packed(cbIdx, umi24, status) for populating packedReadInfo after collapse
    // Only allocated when pSolo.trackReadIdsForTags is true
    khash_t(readid_cbumi) *readIdTracker_; // nullptr if not tracking readIds
    
    // Extended ambiguous entry to store gene/tag info for hash re-keying after resolution
    struct ExtendedAmbiguousEntry : public ReadAlign::AmbiguousEntry {
        // Store (gene, tag, umi24) combinations for each ambiguous CB observation
        // After resolution, we'll create hash entries with resolved CB + stored gene/tag/umi
        struct AmbiguousObservation {
            uint32_t geneIdx;
            uint8_t tagIdx;
            uint32_t umi24;
            uint32_t count;
        };
        std::vector<AmbiguousObservation> observations; // All (gene, tag, umi) observations for this ambiguous CB
    };
    std::unordered_map<ReadAlign::AmbigKey, ExtendedAmbiguousEntry> pendingAmbiguous_; // Ambiguous CB accumulation with gene/tag info

    struct BridgeDeferredReadAccounting {
        SoloReadFlagClass::typeFlag readFlag = 0;
        uint32_t candidateOffset = 0;
        uint16_t candidateCount = 0;
        uint8_t featGood = 0;
        uint8_t multiFeature = 0;
    };
    std::unordered_map<uint32_t, uint64_t> bridgeImmediateReadCounts_; // key: wlCb, value: low32=unique, high32=multi
    std::vector<BridgeDeferredReadAccounting> bridgeDeferredAccounting_;
    std::vector<uint32_t> bridgeDeferredCandidates_; // packed candidate [cb24|qual8]
    std::unordered_map<uint32_t, uint32_t> bridgeCbCompactByWl_;
    std::vector<uint32_t> bridgeCbWlByCompact_;
    std::unordered_map<uint32_t, uint16_t> bridgeGeneCompactByFull_;
    std::vector<uint32_t> bridgeGeneFullByCompact_;

    // STAR_SOLO_NONFLEX_HASH_BRIDGE: kh_val(inlineHash_) = slot id; packed payload below.
    std::vector<uint64_t> bridgePackedSlots_;
    uint64_t bridgeSlotOverflowEvents_ = 0;

    string cbSeq, umiSeq, cbQual, umiQual;

    SoloReadFlagClass readFlag;
    
    SoloReadFeatureStats stats;

    SoloReadFeature (int32 feTy, Parameters &Pin, int iChunk);
    ~SoloReadFeature();
    void record(SoloReadBarcode &soloBar, uint nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot);
    void addCounts(const SoloReadFeature &soloCBin);
    void addStats(const SoloReadFeature &soloCBin);
    void statsOut(ofstream &streamOut);
    void mergeInlineHash(SoloReadFeature &other); // Merge inlineHash_ and pendingAmbiguous_ from other
    void mergePendingAmbiguous(const SoloReadFeature &other);
    void mergeDeferredBridgeAccounting(const SoloReadFeature &other);
    uint32_t getOrCreateBridgeCompactCb(uint32_t wlIdx);
    uint32_t bridgeCompactToWl(uint32_t compactIdx) const;
    uint16_t getOrCreateBridgeCompactGene(uint32_t geneIdx);
    uint32_t bridgeCompactToGene(uint16_t compactIdx) const;
    /** Direct-bridge path: upsert (tupleKey -> slot); slot stores geneFull in 19-bit field. */
    void bridgeDirectTupleAdd(uint64_t tupleKey, uint32_t geneFull, uint32_t umi24, uint32_t delta);
    void maybeSpillBinarySpool(size_t extraBytes);
    // Legacy overload removed
    // Overload that emits read info via a sink (avoids requiring legacy vector storage)
    void inputRecords(uint32 **cbP, uint32 cbPstride, vector<uint32> &cbReadCountTotal, SoloReadFlagClass &readFlagCounts,
                      vector<uint32> &nReadPerCBunique1, vector<uint32> &nReadPerCBmulti1,
                      const std::function<void(uint64, uint32, uint32, uint8)> &recordSink);

private:
    friend class SoloReadInfoLoader;
    friend void record_base(SoloReadFeature *soloReadFeat, SoloReadBarcode &soloBar, uint nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot);
    friend void record_flex(SoloReadFeature *soloReadFeat, SoloReadBarcode &soloBar, uint nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot);
    friend uint32 outputReadCB_base(fstream *streamOut, const uint64 iRead, const int32 featureType, SoloReadBarcode &soloBar,
                                    const ReadSoloFeatures &reFe, const ReadAnnotations &readAnnot, const SoloReadFlagClass &readFlag,
                                    SoloReadFeature *soloReadFeat);
    friend uint32 outputReadCB_flex(fstream *streamOut, const uint64 iRead, const int32 featureType, SoloReadBarcode &soloBar,
                                    const ReadSoloFeatures &reFe, const ReadAnnotations &readAnnot, const SoloReadFlagClass &readFlag,
                                    SoloReadFeature *soloReadFeat);
    friend bool record_flex_hash_screen_keep(SoloReadFeature *soloReadFeat, SoloReadBarcode &soloBar, uint64 iRead, uint16_t geneIdx15, uint8_t cacheClass);
    friend void record_flex_hash_screen_deny(SoloReadFeature *soloReadFeat, SoloReadBarcode &soloBar, uint64 iRead, const char *reason);
    friend const ProbeListIndex* getGlobalProbeIndex(const SoloReadFeature* rf);
    const int32 featureType;

    Parameters &P;
    ParametersSolo &pSolo;
    SoloFeature *owner = nullptr;
};

#endif
