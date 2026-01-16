#ifndef SLAM_COMPAT_H
#define SLAM_COMPAT_H

#include "IncludeDefine.h"
#include <vector>
#include <set>
#include <cstdint>

class Transcriptome;
class Transcript;

struct SlamCompatConfig {
    bool intronic = false;
    bool lenientOverlap = false;
    bool overlapWeight = false;
    bool ignoreOverlap = false;
    int trim5p = 0;
    int trim3p = 0;
};

class SlamCompat {
public:
    // Production constructor - builds caches from Transcriptome
    explicit SlamCompat(const Transcriptome& tr, const SlamCompatConfig& cfg);
    
    // Test constructor - accepts pre-built caches
    SlamCompat(const SlamCompatConfig& cfg,
               std::vector<std::vector<uint32_t>>&& geneToTranscripts,
               std::vector<std::vector<std::pair<uint64_t, uint64_t>>>&& transcriptIntrons);
    
    const SlamCompatConfig& cfg() const { return cfg_; }
    
    // Update trim values (for mid-run auto-trim)
    void updateTrims(int trim5p, int trim3p) {
        cfg_.trim5p = trim5p;
        cfg_.trim3p = trim3p;
    }
    
    // Overlap-gene weighting helper
    double compatOverlapWeight(double baseWeight, size_t readLevelGeneCount) const;
    
    // Position filtering helper (trim only - overlap handled by caller)
    bool compatShouldCountPos(uint32_t mateLocalPos, uint32_t mateLen) const;
    
    // Intronic classification helper
    bool compatIsIntronic(const Transcript& aln,
                          const std::set<uint32_t>& geneCandidates,
                          std::set<uint32_t>& outIntronicGenes) const;
    
    // Overlap computation helper (for lenient overlap)
    uint32_t computeExonOverlap(const Transcript& aln, uint64_t trStart,
                                uint16_t trExN, const uint32_t* trExSE) const;

private:
    SlamCompatConfig cfg_;
    
    // Cached transcript data
    std::vector<std::vector<uint32_t>> geneToTranscripts_;  // geneId -> transcript IDs
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> transcriptIntrons_;  // transcriptId -> intron intervals
    
    // Helper: count transcripts of gene where [alnStart, alnEnd] intersects any intron
    size_t countTranscriptIntronHits(uint32_t geneId, uint64_t alnStart, uint64_t alnEnd) const;
    
    // Build caches from Transcriptome
    void buildCaches(const Transcriptome& tr);
};

#endif
