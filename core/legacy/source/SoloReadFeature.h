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
struct FlexGeneInlineResolveResult;

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
        std::vector<AmbiguousObservation> observations; // Flex: per-read observations (unchanged)
        // Non-Flex direct bridge: aggregate (umi24,gene16) -> read counts (no per-read observation vector)
        std::unordered_map<uint64_t, uint32_t> bridgeAmbigUmiGene_;
        std::vector<double> cbLogLikMatch;
        std::vector<double> cbLogLikMismatch;
        uint32_t cbEvidenceReads = 0;
        // Per-key ambiguous read accounting (non-Flex bridge): no per-read replay vectors
        uint32_t bridgeAmbigGeneFeatU_ = 0;
        uint32_t bridgeAmbigGeneFeatM_ = 0;
        bool bridgeAmbigGeneHaveSampleU_ = false;
        bool bridgeAmbigGeneHaveSampleM_ = false;
        SoloReadFlagClass::typeFlag bridgeAmbigGeneSampleFlagU_ = 0;
        SoloReadFlagClass::typeFlag bridgeAmbigGeneSampleFlagM_ = 0;
        uint32_t bridgeAmbigReadInfoN_ = 0;
        bool bridgeAmbigReadInfoHaveSample_ = false;
        SoloReadFlagClass::typeFlag bridgeAmbigReadInfoSampleFlag_ = 0;
        // Per-candidate qualities for one key-level pin evaluation (same source as representative cbQual)
        std::vector<uint8_t> bridgeAmbigPinCandQuals_;
    };
    std::unordered_map<ReadAlign::AmbigKey, ExtendedAmbiguousEntry> pendingAmbiguous_; // Ambiguous CB accumulation with gene/tag info

    // Aggregated readInfo-only ambiguous CB state (not yet in pendingAmbiguous_; merged on first gene hit)
    struct BridgeAmbigReadInfoOrphanEntry {
        std::vector<uint32_t> candidateIdx;
        std::string cbSeq;
        std::string cbQual;
        std::vector<uint8_t> pinCandQuals_;
        uint32_t readInfoN_ = 0;
        bool haveSample_ = false;
        SoloReadFlagClass::typeFlag sampleFlag_ = 0;
    };
    std::unordered_map<ReadAlign::AmbigKey, BridgeAmbigReadInfoOrphanEntry> bridgeAmbigReadInfoOrphan_;

    std::unordered_map<uint32_t, uint64_t> bridgeImmediateReadCounts_; // key: wlCb, value: low32=unique, high32=multi
    // Filled during resolvePendingAmbiguousToHash (aggregated ambiguous accounting)
    std::vector<uint32_t> bridgePinNreadUnique_;
    std::vector<uint32_t> bridgePinNreadMulti_;

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
    void applyBridgeAmbiguousAggregatedReadAccounting(Parameters &P, int32 featureType,
                                                      const ExtendedAmbiguousEntry &entry,
                                                      bool bayesResolved,
                                                      uint32_t resolvedCbIdx0);
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
                                    SoloReadFeature *soloReadFeat,
                                    const FlexGeneInlineResolveResult *preResolved);
    friend bool record_flex_hash_screen_keep(SoloReadFeature *soloReadFeat, SoloReadBarcode &soloBar, uint64 iRead, uint16_t geneIdx15, uint8_t cacheClass);
    friend void record_flex_hash_screen_deny(SoloReadFeature *soloReadFeat, SoloReadBarcode &soloBar, uint64 iRead, const char *reason);
    friend const ProbeListIndex* getGlobalProbeIndex(const SoloReadFeature* rf);
    const int32 featureType;

    Parameters &P;
    ParametersSolo &pSolo;
    SoloFeature *owner = nullptr;
};

#endif
