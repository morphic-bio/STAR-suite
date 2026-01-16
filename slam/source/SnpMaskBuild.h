#ifndef SNP_MASK_BUILD_H
#define SNP_MASK_BUILD_H

#include "Parameters.h"
#include "Genome.h"
#include "SlamQuant.h"
#include "libem/slam_snp_em.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>

// Forward declarations
class Transcriptome;
class ReadAlignChunk;

struct SnpMaskBuildStats {
    uint64_t totalSites = 0;
    uint64_t candidateSites = 0;
    uint64_t maskedSites = 0;
    uint64_t coverageOverflowCount = 0;
    double globalBaselineBefore = 0.0;
    double globalBaselineAfter = 0.0;
};

class SnpMaskBuild {
public:
    explicit SnpMaskBuild(Parameters& P, const Genome& genome);
    
    // Parse FOFN file and return list of FASTQ pairs/singles
    bool parseFofn(const std::string& fofnPath, 
                   std::vector<std::pair<std::string, std::string>>& fastqPairs,
                   std::string* err);
    
    // Run mask build: align reads, collect stats, fit EM, write outputs
    bool buildMask(const std::vector<std::pair<std::string, std::string>>& fastqPairs,
                   Transcriptome* transcriptome,
                   SnpMaskBuildStats* stats,
                   std::string* err);
    
    // Get temporary SlamQuant instance for use during alignment
    SlamQuant* getTempSlamQuant() const { return tempSlamQuant_.get(); }
    
    // Extract data from SlamQuant after alignment
    void extractFromSlamQuant(SlamQuant* slamQuant);
    
    // Write BED output (sorted, bgzip-compressed, tabix-indexed)
    bool writeBed(const std::string& bedPath, const Genome& genome, std::string* err);
    
    // Write summary TSV
    bool writeSummary(const std::string& summaryPath, const SnpMaskBuildStats& stats, std::string* err);
    
    // Get EM result (after buildMask completes)
    const SlamSnpEMResult& getEMResult() const { return emResult_; }
    
private:
    Parameters& P_;
    const Genome& genome_;
    
    // Per-position counts: genomic position -> (coverage, mismatches)
    // Packed as uint32_t: (coverage << 16) | mismatches
    // Cap coverage/mismatches at 65535
    // This will be populated from SlamQuant's snpMask_ after alignment
    std::unordered_map<uint64_t, uint32_t> positionCounts_;
    
    // Track which positions have been touched (for candidate filtering)
    std::unordered_set<uint64_t> touchedPositions_;
    
    // Temporary SlamQuant instance for collecting observations during alignment
    std::unique_ptr<SlamQuant> tempSlamQuant_;
    
    // EM model and results
    std::unique_ptr<SlamSnpEM> emModel_;
    SlamSnpEMResult emResult_;
    
    // Masked positions (after thresholding)
    std::unordered_set<uint64_t> maskedPositions_;
    
    // Cache for p-values (binomial model) - avoids recomputation in writeBed
    std::unordered_map<uint64_t, double> positionPvalues_;
    
    // Lookup table for binomial model: min_k required at each coverage level
    // min_k_for_snp_[n] = minimum k such that P[X >= k | n, p_err] < pval_threshold
    // Precomputed up to MAX_LOOKUP_COV, use direct computation above that
    static constexpr uint32_t MAX_LOOKUP_COV = 1024;
    std::vector<uint32_t> min_k_for_snp_;
    double lookup_p_err_ = 0.0;  // p_err used for lookup table
    
    // Helper: build lookup table for given p_err
    void buildBinomialLookup(double p_err);
    
    // Helper: record observation at a position
    void recordObservation(uint64_t pos, bool isMismatch);
    
    // Helper: build histogram from position counts
    SnpHistogram buildHistogram() const;
    
    // Helper: filter candidates and compute posteriors (EM model)
    void filterAndComputePosteriors();
    
    // Helper: filter candidates using binomial p-value (GEDI-style)
    void filterAndComputeBinomial(double p_err);
    
    // Helper: get reference base at position (for BED ref/alt columns)
    char getRefBase(uint64_t pos) const;
    
    // Helper: get alt base (most common mismatch) at position
    char getAltBase(uint64_t pos) const;
    
    // Helper: get component assignment (ERR/HET/HOM) for a site
    std::string getComponent(uint32_t n, uint32_t k) const;
};

#endif // SNP_MASK_BUILD_H
