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
#include <functional>

class SoloFeature;
class ProbeListIndex;

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

    fstream *streamReads;

    // Inline hash mode: per-thread hash table (replaces temp stream files)
    khash_t(cg_agg) *inlineHash_; // nullptr if not using inline hash mode
    // Extended ambiguous entry to store gene/tag info for hash re-keying after resolution
    struct ExtendedAmbiguousEntry : public ReadAlign::AmbiguousEntry {
        // Store (gene, tag, umi24) combinations for each ambiguous CB observation
        // After resolution, we'll create hash entries with resolved CB + stored gene/tag/umi
        struct AmbiguousObservation {
            uint16_t geneIdx;
            uint8_t tagIdx;
            uint32_t umi24;
            uint32_t count;
        };
        std::vector<AmbiguousObservation> observations; // All (gene, tag, umi) observations for this ambiguous CB
    };
    std::unordered_map<ReadAlign::AmbigKey, ExtendedAmbiguousEntry> pendingAmbiguous_; // Ambiguous CB accumulation with gene/tag info

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
    // Legacy overload removed
    // Overload that emits read info via a sink (avoids requiring legacy vector storage)
    void inputRecords(uint32 **cbP, uint32 cbPstride, vector<uint32> &cbReadCountTotal, SoloReadFlagClass &readFlagCounts,
                      vector<uint32> &nReadPerCBunique1, vector<uint32> &nReadPerCBmulti1,
                      const std::function<void(uint64, uint32, uint32, uint8)> &recordSink);

private:
    friend class SoloReadInfoLoader;
    friend const ProbeListIndex* getGlobalProbeIndex(const SoloReadFeature* rf);
    const int32 featureType;

    Parameters &P;
    ParametersSolo &pSolo;
    SoloFeature *owner = nullptr;
};

#endif
