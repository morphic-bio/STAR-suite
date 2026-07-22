#ifndef SLAM_COMPAT_H
#define SLAM_COMPAT_H

#include "IncludeDefine.h"
#include <array>
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
    std::array<int, 2> trim5p = {{0, 0}};
    std::array<int, 2> trim3p = {{0, 0}};
};

class SlamCompat {
public:
    explicit SlamCompat(const Transcriptome& tr, const SlamCompatConfig& cfg);

    SlamCompat(const SlamCompatConfig& cfg,
               std::vector<std::vector<uint32_t>>&& geneToTranscripts,
               std::vector<std::vector<std::pair<uint64_t, uint64_t>>>&& transcriptIntrons);

    const SlamCompatConfig& cfg() const { return cfg_; }

    void updateTrims(int trim5p_m1, int trim3p_m1, int trim5p_m2, int trim3p_m2) {
        cfg_.trim5p[0] = trim5p_m1;
        cfg_.trim3p[0] = trim3p_m1;
        cfg_.trim5p[1] = trim5p_m2;
        cfg_.trim3p[1] = trim3p_m2;
    }

    void updateTrims(int trim5p, int trim3p) {
        updateTrims(trim5p, trim3p, trim5p, trim3p);
    }

    double compatOverlapWeight(double baseWeight, size_t readLevelGeneCount) const;

    bool compatShouldCountPos(uint32_t mateLocalPos, uint32_t mateLen,
                              uint32_t mateIndex = 0) const;

    bool compatIsIntronic(const Transcript& aln,
                          const std::set<uint32_t>& geneCandidates,
                          std::set<uint32_t>& outIntronicGenes) const;

    uint32_t computeExonOverlap(const Transcript& aln, uint64_t trStart,
                                uint16_t trExN, const uint32_t* trExSE) const;

private:
    SlamCompatConfig cfg_;

    std::vector<std::vector<uint32_t>> geneToTranscripts_;
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> transcriptIntrons_;

    size_t countTranscriptIntronHits(uint32_t geneId, uint64_t alnStart, uint64_t alnEnd) const;

    void buildCaches(const Transcriptome& tr);
};

#endif
