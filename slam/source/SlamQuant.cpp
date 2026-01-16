#include "SlamQuant.h"
#include "SlamCompat.h"

#include "Genome.h"
#include "Transcriptome.h"
#include "htslib/htslib/bgzf.h"
#include "htslib/htslib/kstring.h"
#include "htslib/htslib/tbx.h"
#include "htslib/htslib/vcf.h"

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <cstdlib>

namespace {
constexpr uint32_t kSnpMinCoverage = 10;
constexpr double kSnpFallbackMismatchFrac = 0.22;  // Fallback threshold (GEDI parity)
constexpr double kSnpMinClamp = 0.10;              // Minimum allowed threshold
constexpr double kSnpMaxClamp = 0.60;              // Maximum allowed threshold
constexpr uint64_t kSnpMinSitesForAuto = 1000;     // Minimum eligible sites for auto-estimation
constexpr double kSnpKneeEpsilon = 0.02;           // Minimum knee strength
constexpr size_t kSnpHistogramBins = 100;          // Number of bins for mismatch fraction histogram
const char* debugDropReasonName(SlamDebugDropReason reason) {
    switch (reason) {
        case SlamDebugDropReason::None:
            return "PASS";
        case SlamDebugDropReason::NoGenes:
            return "DROP_NO_GENES";
        case SlamDebugDropReason::SnpMask:
            return "DROP_SNP_MASK";
        case SlamDebugDropReason::Strandness:
            return "DROP_STRANDNESS";
        case SlamDebugDropReason::ZeroWeight:
            return "DROP_ZERO_WEIGHT";
        default:
            return "UNKNOWN";
    }
}
}

SlamQuant::SlamQuant(uint32_t nGenes, bool snpDetect, double snpMismatchFrac, bool snpObsAnyMismatch)
    : geneStats_(nGenes),
      snpDetectEnabled_(snpDetect),
      snpObsAnyMismatch_(snpObsAnyMismatch),
      snpMismatchFrac_(snpMismatchFrac) {}

SlamQuant::SlamQuant(uint32_t nGenes, std::vector<uint8_t> allowedGenes, bool snpDetect, double snpMismatchFrac, bool snpObsAnyMismatch)
    : geneStats_(nGenes),
      snpDetectEnabled_(snpDetect),
      snpObsAnyMismatch_(snpObsAnyMismatch),
      snpMismatchFrac_(snpMismatchFrac),
      allowedGenes_(std::move(allowedGenes)) {}

void SlamQuant::enableVarianceAnalysis(uint32_t maxReads, uint32_t minReads,
                                       uint32_t smoothWindow, uint32_t minSegLen, uint32_t maxTrim) {
    varianceMaxReads_ = maxReads;
    varianceMinReads_ = minReads;
    varianceAnalyzer_.reset(new SlamVarianceAnalyzer(maxReads, minReads, smoothWindow, minSegLen, maxTrim));
}

uint32_t SlamQuant::getVarianceMaxReads() const {
    return varianceMaxReads_;
}

uint32_t SlamQuant::getVarianceMinReads() const {
    return varianceMinReads_;
}

void SlamQuant::resetVarianceAnalysis() {
    if (varianceAnalyzer_) {
        varianceAnalyzer_->reset();
    }
}

bool SlamQuant::recordVarianceRead() {
    if (varianceAnalyzer_) {
        return varianceAnalyzer_->recordRead();
    }
    return false;
}

void SlamQuant::recordVariancePosition(uint32_t readPos, uint8_t qual, bool isT, bool isTc) {
    if (varianceAnalyzer_) {
        varianceAnalyzer_->recordPosition(readPos, qual, isT, isTc);
    }
}

SlamVarianceTrimResult SlamQuant::computeVarianceTrim(uint32_t readLength) {
    if (!varianceAnalyzer_) {
        SlamVarianceTrimResult result;
        result.success = false;
        result.mode = "disabled";
        return result;
    }
    return varianceAnalyzer_->computeTrim(readLength);
}

void SlamQuant::initDebug(const Transcriptome& tr,
                          const std::unordered_set<std::string>& debugGenes,
                          const std::unordered_set<std::string>& debugReads,
                          size_t maxReads,
                          const std::string& outPrefix) {
    // Enable debug mode if any debug category is requested. SNP-site debug is configured separately
    // but shares the debug output prefix.
    debugEnabled_ = !debugGenes.empty() || !debugReads.empty() || !outPrefix.empty();
    debugMaxReads_ = maxReads;
    debugOutPrefix_ = outPrefix;
    debugReadRecords_.clear();
    debugReadSet_.clear();
    debugGeneMask_.clear();
    debugGeneStats_.clear();
    if (!debugEnabled_) {
        return;
    }
    debugReadSet_ = debugReads;
    if (!debugGenes.empty()) {
        debugGeneMask_.assign(tr.nGe, 0);
        debugGeneStats_.assign(tr.nGe, SlamDebugGeneStats());
        std::unordered_set<uint32_t> numericGenes;
        for (const auto& name : debugGenes) {
            if (name.empty()) {
                continue;
            }
            bool numeric = true;
            for (char c : name) {
                if (c < '0' || c > '9') {
                    numeric = false;
                    break;
                }
            }
            if (numeric) {
                uint32_t gid = static_cast<uint32_t>(std::stoul(name));
                if (gid < tr.nGe) {
                    numericGenes.insert(gid);
                }
            }
        }
        for (uint32_t i = 0; i < tr.nGe; ++i) {
            if (numericGenes.count(i) > 0) {
                debugGeneMask_[i] = 1;
                continue;
            }
            const std::string& gid = tr.geID[i];
            const std::string& gname = tr.geName[i];
            const std::string& gcanon = (i < tr.geIDCanonical.size()) ? tr.geIDCanonical[i] : std::string();
            if ((!gid.empty() && debugGenes.count(gid) > 0) ||
                (!gname.empty() && debugGenes.count(gname) > 0) ||
                (!gcanon.empty() && debugGenes.count(gcanon) > 0)) {
                debugGeneMask_[i] = 1;
            }
        }
    }
}

bool SlamQuant::debugGeneEnabled(uint32_t geneId) const {
    return debugEnabled_ && geneId < debugGeneMask_.size() && debugGeneMask_[geneId];
}

bool SlamQuant::debugReadMatch(const char* readName) const {
    if (!debugEnabled_ || debugReadSet_.empty() || readName == nullptr) {
        return false;
    }
    std::string name(readName);
    size_t end = name.find_first_of(" \t");
    if (end != std::string::npos) {
        name = name.substr(0, end);
    }
    if (!name.empty() && name[0] == '@') {
        name.erase(0, 1);
    }
    return debugReadSet_.count(name) > 0;
}

void SlamQuant::debugCountDrop(uint32_t geneId, SlamDebugDropReason reason) {
    if (!debugGeneEnabled(geneId)) {
        return;
    }
    SlamDebugGeneStats& stats = debugGeneStats_[geneId];
    if (reason == SlamDebugDropReason::SnpMask) {
        stats.dropsSnpMask++;
    } else if (reason == SlamDebugDropReason::Strandness) {
        stats.dropsStrandness++;
    }
}

void SlamQuant::debugAddAssignment(uint32_t geneId, double weight, bool intronic,
                                   bool oppositeStrand, uint16_t nT, uint8_t k) {
    if (!debugGeneEnabled(geneId)) {
        return;
    }
    SlamDebugGeneStats& stats = debugGeneStats_[geneId];
    stats.readsAssigned++;
    if (intronic) {
        stats.readsAssignedIntronic++;
        stats.intronicWeight += weight;
    } else {
        stats.exonicWeight += weight;
        stats.coverage += static_cast<double>(nT) * weight;
        stats.conversions += static_cast<double>(k) * weight;
    }
    if (oppositeStrand) {
        stats.antisenseWeight += weight;
    } else {
        stats.senseWeight += weight;
    }
}

void SlamQuant::debugLogRead(const SlamDebugReadRecord& record) {
    if (!debugEnabled_ || debugMaxReads_ == 0) {
        return;
    }
    if (debugReadRecords_.size() >= debugMaxReads_) {
        return;
    }
    debugReadRecords_.push_back(record);
}

void SlamQuant::enableSnpSiteDebug(uint64_t absPos, int window, const std::string& locString) {
    if (!debugEnabled_) {
        // Debug output prefix is required to write files; caller should have ensured initDebug().
        debugEnabled_ = true;
    }
    debugSnpEnabled_ = true;
    debugSnpAbsPos_ = absPos;
    debugSnpWindow_ = (window < 0) ? 0 : window;
    debugSnpLoc_ = locString;
    size_t n = static_cast<size_t>(debugSnpWindow_ * 2 + 1);
    auto resetVec = [&](std::vector<uint64_t>& v) { v.assign(n, 0); };
    resetVec(debugSnpCov_);
    resetVec(debugSnpAnyMis_);
    resetVec(debugSnpConvMis_);
    resetVec(debugSnpCovPrimary_);
    resetVec(debugSnpAnyMisPrimary_);
    resetVec(debugSnpConvMisPrimary_);
    resetVec(debugSnpCovNh1_);
    resetVec(debugSnpAnyMisNh1_);
    resetVec(debugSnpConvMisNh1_);
    resetVec(debugSnpCovNhGt1_);
    resetVec(debugSnpAnyMisNhGt1_);
    resetVec(debugSnpConvMisNhGt1_);
    resetVec(debugSnpCovMapq20_);
    resetVec(debugSnpAnyMisMapq20_);
    resetVec(debugSnpConvMisMapq20_);
}

void SlamQuant::debugSnpSiteObserve(uint64_t absPos, bool anyMismatch, bool convMismatch,
                                   double weight, bool primaryFlag, int mapq) {
    if (!debugSnpEnabled_) {
        return;
    }
    int64_t delta = static_cast<int64_t>(absPos) - static_cast<int64_t>(debugSnpAbsPos_);
    if (delta < -debugSnpWindow_ || delta > debugSnpWindow_) {
        return;
    }
    size_t idx = static_cast<size_t>(delta + debugSnpWindow_);

    debugSnpCov_[idx] += 1;
    if (anyMismatch) debugSnpAnyMis_[idx] += 1;
    if (convMismatch) debugSnpConvMis_[idx] += 1;

    if (primaryFlag) {
        debugSnpCovPrimary_[idx] += 1;
        if (anyMismatch) debugSnpAnyMisPrimary_[idx] += 1;
        if (convMismatch) debugSnpConvMisPrimary_[idx] += 1;
    }

    // Infer NH from weight (Alignments mode => weight~1/NH; Uniform => weight==1)
    uint32_t nh = 1;
    if (weight > 0.0) {
        double inv = 1.0 / weight;
        if (inv < 1.0) inv = 1.0;
        if (inv > 1e6) inv = 1e6;
        nh = static_cast<uint32_t>(std::llround(inv));
        if (nh == 0) nh = 1;
    }
    if (nh <= 1) {
        debugSnpCovNh1_[idx] += 1;
        if (anyMismatch) debugSnpAnyMisNh1_[idx] += 1;
        if (convMismatch) debugSnpConvMisNh1_[idx] += 1;
    } else {
        debugSnpCovNhGt1_[idx] += 1;
        if (anyMismatch) debugSnpAnyMisNhGt1_[idx] += 1;
        if (convMismatch) debugSnpConvMisNhGt1_[idx] += 1;
    }

    if (mapq >= 20) {
        debugSnpCovMapq20_[idx] += 1;
        if (anyMismatch) debugSnpAnyMisMapq20_[idx] += 1;
        if (convMismatch) debugSnpConvMisMapq20_[idx] += 1;
    }
}

const char* slamMismatchCategoryName(SlamMismatchCategory cat) {
    switch (cat) {
        case SlamMismatchCategory::Exonic:
            return "Exonic";
        case SlamMismatchCategory::ExonicSense:
            return "ExonicSense";
        case SlamMismatchCategory::Intronic:
            return "Intronic";
        case SlamMismatchCategory::IntronicSense:
            return "IntronicSense";
        default:
            return "Unknown";
    }
}

void SlamQuant::addRead(uint32_t geneId, uint16_t nT, uint8_t k, double weight) {
    if (geneId >= geneStats_.size() || weight <= 0.0) {
        return;
    }
    if (!allowedGenes_.empty() && geneId < allowedGenes_.size() && allowedGenes_[geneId] == 0) {
        return;
    }
    if (nT > 511) nT = 511;
    if (k > 255) k = 255;
    uint16_t key = static_cast<uint16_t>((nT << 8) | k);
    SlamGeneStats& stats = geneStats_[geneId];
    stats.histogram[key] += weight;
    stats.readCount += weight;
    stats.conversions += weight * static_cast<double>(k);
    stats.coverage += weight * static_cast<double>(nT);
}

void SlamQuant::recordSnpObservation(uint64_t pos, bool anyMismatch, bool convMismatch) {
    if (!snpDetectEnabled_) {
        return;
    }
    const bool isMismatch = snpObsAnyMismatch_ ? anyMismatch : convMismatch;
    uint32_t &entry = snpMask_[pos];
    uint32_t cov = entry >> 16;
    uint32_t mis = entry & 0xFFFF;
    if (cov < 0xFFFF) {
        ++cov;
    }
    if (isMismatch && mis < 0xFFFF) {
        ++mis;
    }
    entry = (cov << 16) | mis;
}

std::unordered_map<uint64_t, uint32_t> SlamQuant::getSnpMaskData() const {
    return snpMask_;  // Return copy
}

void SlamQuant::bufferSnpRead(uint32_t geneId, uint16_t nT,
                              const std::vector<uint32_t>& mismatchPositions, double weight) {
    if (!snpDetectEnabled_ || mismatchPositions.empty() || weight <= 0.0) {
        return;
    }
    if (!allowedGenes_.empty() && geneId < allowedGenes_.size() && allowedGenes_[geneId] == 0) {
        return;
    }
    uint16_t k = static_cast<uint16_t>(mismatchPositions.size() > 0xFFFF
                                           ? 0xFFFF
                                           : mismatchPositions.size());
    snpReadBuffer_.push_back(geneId);
    snpReadBuffer_.push_back((static_cast<uint32_t>(k) << 16) | nT);
    for (uint16_t i = 0; i < k; ++i) {
        snpReadBuffer_.push_back(mismatchPositions[i]);
    }
    snpReadWeights_.push_back(weight);
}

// Estimate SNP mismatch threshold using knee/elbow detection on cumulative distribution
// Returns: pair<threshold, kneeBin>; threshold=0 indicates fallback
static std::pair<double, uint32_t> estimateSnpMismatchFrac(
    const std::unordered_map<uint64_t, uint32_t>& snpMask, uint64_t* eligibleSites) {
    
    // Build histogram of mismatch fractions for eligible sites
    std::vector<uint64_t> histogram(kSnpHistogramBins, 0);
    uint64_t totalEligible = 0;
    
    for (const auto& kv : snpMask) {
        uint32_t cov = kv.second >> 16;
        uint32_t mis = kv.second & 0xFFFF;
        if (cov >= kSnpMinCoverage && cov > 0) {
            double frac = static_cast<double>(mis) / static_cast<double>(cov);
            // Clamp to [0, 1] and bin
            if (frac < 0.0) frac = 0.0;
            if (frac > 1.0) frac = 1.0;
            size_t bin = static_cast<size_t>(frac * (kSnpHistogramBins - 1));
            if (bin >= kSnpHistogramBins) bin = kSnpHistogramBins - 1;
            histogram[bin]++;
            totalEligible++;
        }
    }
    
    if (eligibleSites) *eligibleSites = totalEligible;
    
    // Check minimum sites requirement
    if (totalEligible < kSnpMinSitesForAuto) {
        return {0.0, 0};  // Fallback
    }
    
    // Compute cumulative counts N(>= f) from high to low bins
    std::vector<double> cumulative(kSnpHistogramBins, 0.0);
    uint64_t runningSum = 0;
    for (int i = static_cast<int>(kSnpHistogramBins) - 1; i >= 0; --i) {
        runningSum += histogram[i];
        cumulative[i] = static_cast<double>(runningSum);
    }
    
    // Apply log1p transformation for stability
    std::vector<double> logCumulative(kSnpHistogramBins);
    for (size_t i = 0; i < kSnpHistogramBins; ++i) {
        logCumulative[i] = std::log1p(cumulative[i]);
    }
    
    // Normalize x (bin centers) and y (log counts) to [0, 1]
    double xMin = 0.0;
    double xMax = 1.0;
    double yMin = logCumulative[kSnpHistogramBins - 1];  // Smallest (high frac end)
    double yMax = logCumulative[0];                       // Largest (low frac end)
    
    if (yMax <= yMin) {
        return {0.0, 0};  // Flat distribution, fallback
    }
    
    // Kneedle algorithm: find max distance from diagonal y = 1 - x
    // In normalized space, diagonal goes from (0,1) to (1,0)
    // Distance formula: maximize (y_norm + x_norm - 1)
    double maxDist = -1.0;
    size_t kneeBin = 0;
    
    for (size_t i = 0; i < kSnpHistogramBins; ++i) {
        double xNorm = (static_cast<double>(i) / (kSnpHistogramBins - 1) - xMin) / (xMax - xMin);
        double yNorm = (logCumulative[i] - yMin) / (yMax - yMin);
        double dist = yNorm + xNorm - 1.0;
        if (dist > maxDist) {
            maxDist = dist;
            kneeBin = i;
        }
    }
    
    // Check knee strength
    if (maxDist < kSnpKneeEpsilon) {
        return {0.0, 0};  // Weak knee, fallback
    }
    
    // Convert bin to mismatch fraction threshold
    double threshold = static_cast<double>(kneeBin) / (kSnpHistogramBins - 1);
    
    // Clamp to valid range
    if (threshold < kSnpMinClamp) threshold = kSnpMinClamp;
    if (threshold > kSnpMaxClamp) threshold = kSnpMaxClamp;
    
    return {threshold, static_cast<uint32_t>(kneeBin)};
}

void SlamQuant::finalizeSnpMask(SlamSnpBufferStats* outStats) {
    if (outStats) {
        *outStats = SlamSnpBufferStats();
        outStats->bufferedReads = snpReadWeights_.size();
        outStats->bufferEntries = snpReadBuffer_.size();
        outStats->maskEntries = snpMask_.size();
        outStats->bufferBytes =
            static_cast<uint64_t>(snpReadBuffer_.size()) * sizeof(uint32_t) +
            static_cast<uint64_t>(snpReadWeights_.size()) * sizeof(double);
    }
    if (!snpDetectEnabled_ || snpFinalized_) {
        return;
    }
    snpFinalized_ = true;
    if (snpReadBuffer_.empty()) {
        return;
    }
    
    // Determine mismatch fraction threshold
    double mismatchFracUsed = 0.0;
    double mismatchFracAuto = 0.0;
    uint32_t kneeBin = 0;
    uint64_t eligibleSites = 0;
    std::string mode;
    
    if (snpMismatchFrac_ > 0.0) {
        // Explicit threshold provided
        mismatchFracUsed = snpMismatchFrac_;
        mode = "explicit";
    } else {
        // Auto-estimate threshold
        auto [autoThreshold, autoBin] = estimateSnpMismatchFrac(snpMask_, &eligibleSites);
        mismatchFracAuto = autoThreshold;
        kneeBin = autoBin;
        
        if (autoThreshold > 0.0) {
            mismatchFracUsed = autoThreshold;
            mode = "auto";
        } else {
            // Fallback to default
            mismatchFracUsed = kSnpFallbackMismatchFrac;
            mode = "auto_fallback";
        }
    }
    
    // Build blacklist using determined threshold
    std::unordered_set<uint64_t> blacklist;
    blacklist.reserve(snpMask_.size());
    for (const auto &kv : snpMask_) {
        uint32_t cov = kv.second >> 16;
        uint32_t mis = kv.second & 0xFFFF;
        if (cov >= kSnpMinCoverage && cov > 0) {
            double rate = static_cast<double>(mis) / static_cast<double>(cov);
            if (rate > mismatchFracUsed) {
                blacklist.insert(kv.first);
            }
        }
    }
    
    if (outStats) {
        outStats->blacklistEntries = blacklist.size();
        outStats->mismatchFracUsed = mismatchFracUsed;
        outStats->mismatchFracAuto = mismatchFracAuto;
        outStats->mismatchFracMode = mode;
        outStats->kneeBin = kneeBin;
        outStats->eligibleSites = eligibleSites;
    }
    
    uint64_t totalMismatches = 0;
    uint64_t totalMismatchesKept = 0;
    size_t idx = 0;
    size_t readIdx = 0;
    while (idx + 1 < snpReadBuffer_.size() && readIdx < snpReadWeights_.size()) {
        uint32_t geneId = snpReadBuffer_[idx++];
        uint32_t header = snpReadBuffer_[idx++];
        uint16_t nT = static_cast<uint16_t>(header & 0xFFFF);
        uint16_t k = static_cast<uint16_t>(header >> 16);
        uint16_t validK = 0;
        totalMismatches += k;
        for (uint16_t i = 0; i < k && idx < snpReadBuffer_.size(); ++i, ++idx) {
            uint32_t pos = snpReadBuffer_[idx];
            if (blacklist.find(pos) == blacklist.end()) {
                ++validK;
            }
        }
        totalMismatchesKept += validK;
        uint8_t k8 = static_cast<uint8_t>(validK > 255 ? 255 : validK);
        addRead(geneId, nT, k8, snpReadWeights_[readIdx]);
        ++readIdx;
    }
    if (outStats) {
        outStats->bufferedMismatches = totalMismatches;
        outStats->bufferedMismatchesKept = totalMismatchesKept;
        uint64_t reads = snpReadWeights_.size();
        if (reads > 0) {
            outStats->avgMismatches = static_cast<double>(totalMismatches) /
                static_cast<double>(reads);
            outStats->avgMismatchesKept = static_cast<double>(totalMismatchesKept) /
                static_cast<double>(reads);
        }
    }
    snpReadBuffer_.clear();
    snpReadWeights_.clear();
    snpMask_.clear();
}

void SlamQuant::addTransitionBase(SlamMismatchCategory category, uint32_t readPos, bool secondMate,
                                  bool overlap, bool opposite, int genomicBase, int readBase, double weight) {
    if (weight <= 0.0 || genomicBase < 0 || genomicBase > 3 || readBase < 0 || readBase > 3) {
        return;
    }
    const size_t cat = static_cast<size_t>(category);
    transitions_[cat].coverage[genomicBase] += weight;
    if (genomicBase != readBase) {
        transitions_[cat].mismatches[genomicBase][readBase] += weight;
    }
    if (!overlap) {
        SlamTransitionStats& orient = secondMate ? transitionsSecond_[cat] : transitionsFirst_[cat];
        orient.coverage[genomicBase] += weight;
        if (genomicBase != readBase) {
            orient.mismatches[genomicBase][readBase] += weight;
        }
    }
    if (readPos >= positionTransitions_[cat].size()) {
        positionTransitions_[cat].resize(readPos + 1);
    }
    SlamPositionStats& pos = positionTransitions_[cat][readPos];
    const size_t ov = overlap ? 1 : 0;
    const size_t opp = opposite ? 1 : 0;
    pos.coverage[ov][opp][genomicBase] += weight;
    if (genomicBase != readBase) {
        pos.mismatches[ov][opp][genomicBase][readBase] += weight;
    }
}

void SlamQuant::merge(const SlamQuant& other) {
    if (other.geneStats_.size() != geneStats_.size()) {
        return;
    }
    if (other.snpDetectEnabled_) {
        snpDetectEnabled_ = true;
    }
    for (size_t i = 0; i < geneStats_.size(); ++i) {
        SlamGeneStats& dst = geneStats_[i];
        const SlamGeneStats& src = other.geneStats_[i];
        dst.readCount += src.readCount;
        dst.conversions += src.conversions;
        dst.coverage += src.coverage;
        for (const auto& kv : src.histogram) {
            dst.histogram[kv.first] += kv.second;
        }
    }
    // Merge variance analyzer stats
    if (other.varianceAnalyzer_ && varianceAnalyzer_) {
        varianceAnalyzer_->merge(*other.varianceAnalyzer_);
    } else if (other.varianceAnalyzer_ && !varianceAnalyzer_) {
        // Copy variance analyzer if we don't have one (use defaults for segmented regression params)
        uint32_t maxReads = other.varianceAnalyzer_->readsAnalyzed() > 0 ? 
            static_cast<uint32_t>(other.varianceAnalyzer_->readsAnalyzed()) : 100000;
        varianceAnalyzer_.reset(new SlamVarianceAnalyzer(maxReads, 1000, 5, 3, 15));
        varianceAnalyzer_->merge(*other.varianceAnalyzer_);
    }
    
    // Merge diagnostics
    diag_.readsDroppedSnpMask += other.diag_.readsDroppedSnpMask;
    diag_.readsDroppedStrandness += other.diag_.readsDroppedStrandness;
    diag_.readsZeroGenes += other.diag_.readsZeroGenes;
    diag_.readsProcessed += other.diag_.readsProcessed;
    diag_.readsNAlignWithGeneZero += other.diag_.readsNAlignWithGeneZero;
    diag_.compatAlignsReclassifiedIntronic += other.diag_.compatAlignsReclassifiedIntronic;
    diag_.compatAlignsLenientAccepted += other.diag_.compatAlignsLenientAccepted;
    diag_.compatAlignsOverlapWeightApplied += other.diag_.compatAlignsOverlapWeightApplied;
    diag_.compatPositionsSkippedOverlap += other.diag_.compatPositionsSkippedOverlap;
    diag_.compatPositionsSkippedTrim += other.diag_.compatPositionsSkippedTrim;
    diag_.readsSumWeightLessThanOne += other.diag_.readsSumWeightLessThanOne;
    for (const auto& kv : other.diag_.nTrDistribution) {
        diag_.nTrDistribution[kv.first] += kv.second;
    }
    for (const auto& kv : other.diag_.geneSetSizeDistribution) {
        diag_.geneSetSizeDistribution[kv.first] += kv.second;
    }
    for (const auto& kv : other.diag_.nAlignWithGeneDistribution) {
        diag_.nAlignWithGeneDistribution[kv.first] += kv.second;
    }
    for (const auto& kv : other.diag_.sumWeightDistribution) {
        diag_.sumWeightDistribution[kv.first] += kv.second;
    }
    for (size_t c = 0; c < kSlamMismatchCategoryCount; ++c) {
        for (int g = 0; g < 4; ++g) {
            transitions_[c].coverage[g] += other.transitions_[c].coverage[g];
            transitionsFirst_[c].coverage[g] += other.transitionsFirst_[c].coverage[g];
            transitionsSecond_[c].coverage[g] += other.transitionsSecond_[c].coverage[g];
            for (int r = 0; r < 4; ++r) {
                transitions_[c].mismatches[g][r] += other.transitions_[c].mismatches[g][r];
                transitionsFirst_[c].mismatches[g][r] += other.transitionsFirst_[c].mismatches[g][r];
                transitionsSecond_[c].mismatches[g][r] += other.transitionsSecond_[c].mismatches[g][r];
            }
        }
        if (other.positionTransitions_[c].size() > positionTransitions_[c].size()) {
            positionTransitions_[c].resize(other.positionTransitions_[c].size());
        }
        for (size_t i = 0; i < other.positionTransitions_[c].size(); ++i) {
            for (int ov = 0; ov < 2; ++ov) {
                for (int opp = 0; opp < 2; ++opp) {
                    for (int g = 0; g < 4; ++g) {
                        positionTransitions_[c][i].coverage[ov][opp][g] +=
                            other.positionTransitions_[c][i].coverage[ov][opp][g];
                        for (int r = 0; r < 4; ++r) {
                            positionTransitions_[c][i].mismatches[ov][opp][g][r] +=
                                other.positionTransitions_[c][i].mismatches[ov][opp][g][r];
                        }
                    }
                }
            }
        }
    }
    if (!other.snpMask_.empty()) {
        for (const auto &kv : other.snpMask_) {
            uint32_t &entry = snpMask_[kv.first];
            uint32_t cov = entry >> 16;
            uint32_t mis = entry & 0xFFFF;
            uint32_t addCov = kv.second >> 16;
            uint32_t addMis = kv.second & 0xFFFF;
            cov = std::min<uint32_t>(0xFFFF, cov + addCov);
            mis = std::min<uint32_t>(0xFFFF, mis + addMis);
            entry = (cov << 16) | mis;
        }
    }
    if (!other.snpReadBuffer_.empty()) {
        snpReadBuffer_.insert(snpReadBuffer_.end(), other.snpReadBuffer_.begin(), other.snpReadBuffer_.end());
        snpReadWeights_.insert(snpReadWeights_.end(), other.snpReadWeights_.begin(), other.snpReadWeights_.end());
    }
    if (debugEnabled_ && other.debugEnabled_) {
        if (debugGeneStats_.size() == other.debugGeneStats_.size()) {
            for (size_t i = 0; i < debugGeneStats_.size(); ++i) {
                SlamDebugGeneStats& dst = debugGeneStats_[i];
                const SlamDebugGeneStats& src = other.debugGeneStats_[i];
                dst.exonicWeight += src.exonicWeight;
                dst.intronicWeight += src.intronicWeight;
                dst.senseWeight += src.senseWeight;
                dst.antisenseWeight += src.antisenseWeight;
                dst.conversions += src.conversions;
                dst.coverage += src.coverage;
                dst.readsAssigned += src.readsAssigned;
                dst.readsAssignedIntronic += src.readsAssignedIntronic;
                dst.dropsSnpMask += src.dropsSnpMask;
                dst.dropsStrandness += src.dropsStrandness;
            }
        }
        if (debugMaxReads_ != 0 && !other.debugReadRecords_.empty()) {
            for (const auto& rec : other.debugReadRecords_) {
                if (debugReadRecords_.size() >= debugMaxReads_) {
                    break;
                }
                debugReadRecords_.push_back(rec);
            }
        }
    }

    if (debugSnpEnabled_ && other.debugSnpEnabled_ &&
        debugSnpAbsPos_ == other.debugSnpAbsPos_ &&
        debugSnpWindow_ == other.debugSnpWindow_) {
        auto addVec = [&](std::vector<uint64_t>& dst, const std::vector<uint64_t>& src) {
            if (dst.size() != src.size()) return;
            for (size_t i = 0; i < dst.size(); ++i) {
                dst[i] += src[i];
            }
        };
        addVec(debugSnpCov_, other.debugSnpCov_);
        addVec(debugSnpAnyMis_, other.debugSnpAnyMis_);
        addVec(debugSnpConvMis_, other.debugSnpConvMis_);
        addVec(debugSnpCovPrimary_, other.debugSnpCovPrimary_);
        addVec(debugSnpAnyMisPrimary_, other.debugSnpAnyMisPrimary_);
        addVec(debugSnpConvMisPrimary_, other.debugSnpConvMisPrimary_);
        addVec(debugSnpCovNh1_, other.debugSnpCovNh1_);
        addVec(debugSnpAnyMisNh1_, other.debugSnpAnyMisNh1_);
        addVec(debugSnpConvMisNh1_, other.debugSnpConvMisNh1_);
        addVec(debugSnpCovNhGt1_, other.debugSnpCovNhGt1_);
        addVec(debugSnpAnyMisNhGt1_, other.debugSnpAnyMisNhGt1_);
        addVec(debugSnpConvMisNhGt1_, other.debugSnpConvMisNhGt1_);
        addVec(debugSnpCovMapq20_, other.debugSnpCovMapq20_);
        addVec(debugSnpAnyMisMapq20_, other.debugSnpAnyMisMapq20_);
        addVec(debugSnpConvMisMapq20_, other.debugSnpConvMisMapq20_);
    }
}

void SlamQuant::write(const Transcriptome& tr, const std::string& outFile,
                      double errorRate, double convRate) const {
    std::ofstream out(outFile.c_str());
    if (!out.good()) {
        return;
    }
    out << "Gene\tSymbol\tReadCount\tConversions\tCoverage\tNTR\tMAP\tSigma\tLogLikelihood\n";

    SlamSolver solver(errorRate, convRate);
    for (size_t i = 0; i < geneStats_.size(); ++i) {
        const SlamGeneStats& stats = geneStats_[i];
        if (stats.readCount <= 0.0) {
            continue;
        }
        SlamResult res = solver.solve(stats.histogram);
        const std::string& geneId = tr.geID[i];
        const std::string& geneName = tr.geName[i].empty() ? tr.geID[i] : tr.geName[i];
        out << geneId << "\t"
            << geneName << "\t"
            << stats.readCount << "\t"
            << stats.conversions << "\t"
            << stats.coverage << "\t"
            << res.ntr << "\t"
            << res.ntr << "\t"
            << res.sigma << "\t"
            << res.log_likelihood << "\n";
    }
}

static std::string cleanGrandSlamPrefix(const std::string& outFileNamePrefix) {
    std::string base = outFileNamePrefix;
    size_t slash = base.find_last_of('/');
    if (slash != std::string::npos) {
        base = base.substr(slash + 1);
    }
    while (!base.empty() && (base.back() == '_' || base.back() == '.')) {
        base.pop_back();
    }
    if (base.empty()) {
        base = "STAR";
    }
    return base;
}

void SlamQuant::writeGrandSlam(const Transcriptome& tr, const std::string& outFile,
                               const std::string& outFileNamePrefix,
                               double errorRate, double convRate) const {
    std::ofstream out(outFile.c_str());
    if (!out.good()) {
        return;
    }
    const std::string prefix = cleanGrandSlamPrefix(outFileNamePrefix);
    out << "Gene\tSymbol\t"
        << prefix << " Readcount\t"
        << prefix << " 0.05 quantile\t"
        << prefix << " Mean\t"
        << prefix << " MAP\t"
        << prefix << " 0.95 quantile\t"
        << prefix << " alpha\t"
        << prefix << " beta\t"
        << prefix << " Conversions\t"
        << prefix << " Coverage\t"
        << prefix << " Double-Hits\t"
        << prefix << " Double-Hit Coverage\t"
        << prefix << " min2\t"
        << "Length\n";

    SlamSolver solver(errorRate, convRate);
    for (size_t i = 0; i < geneStats_.size(); ++i) {
        const SlamGeneStats& stats = geneStats_[i];
        if (stats.readCount <= 0.0) {
            continue;
        }
        SlamResult res = solver.solve(stats.histogram);
        const std::string& geneId = tr.geID[i];
        const std::string& geneName = tr.geName[i].empty() ? tr.geID[i] : tr.geName[i];
        double ntr = res.ntr;
        // Placeholder: GEDI posterior not implemented; use NTR for MAP/mean/quantiles.
        // alpha/beta and double-hit metrics are emitted as 0 until tracked.
        out << geneId << "\t"
            << geneName << "\t"
            << stats.readCount << "\t"
            << ntr << "\t"
            << ntr << "\t"
            << ntr << "\t"
            << ntr << "\t"
            << 0 << "\t"
            << 0 << "\t"
            << stats.conversions << "\t"
            << stats.coverage << "\t"
            << 0 << "\t"
            << 0 << "\t"
            << 0 << "\t"
            << 0 << "\n";
    }
}

void SlamQuant::writeDiagnostics(const std::string& diagFile) const {
    std::ofstream out(diagFile.c_str());
    if (!out.good()) {
        return;
    }
    out << "Metric\tValue\n";
    out << "readsProcessed\t" << diag_.readsProcessed << "\n";
    out << "readsDroppedSnpMask\t" << diag_.readsDroppedSnpMask << "\n";
    out << "readsDroppedStrandness\t" << diag_.readsDroppedStrandness << "\n";
    out << "readsZeroGenes\t" << diag_.readsZeroGenes << "\n";
    out << "readsNAlignWithGeneZero\t" << diag_.readsNAlignWithGeneZero << "\n";
    out << "readsSumWeightLessThanOne\t" << diag_.readsSumWeightLessThanOne << "\n";
    out << "\nnTrDistribution:\n";
    for (const auto& kv : diag_.nTrDistribution) {
        out << "nTr_" << kv.first << "\t" << kv.second << "\n";
    }
    out << "\ngeneSetSizeDistribution:\n";
    for (const auto& kv : diag_.geneSetSizeDistribution) {
        out << "geneSetSize_" << kv.first << "\t" << kv.second << "\n";
    }
    out << "\nnAlignWithGeneDistribution:\n";
    for (const auto& kv : diag_.nAlignWithGeneDistribution) {
        out << "nAlignWithGene_" << kv.first << "\t" << kv.second << "\n";
    }
    out << "\nsumWeightDistribution (bucketed):\n";
    for (const auto& kv : diag_.sumWeightDistribution) {
        double bucketStart = kv.first * 0.1;
        double bucketEnd = (kv.first == 10) ? 999.0 : (kv.first + 1) * 0.1;
        out << "sumWeight_" << bucketStart << "_to_" << bucketEnd << "\t" << kv.second << "\n";
    }
    
    // Compat mode counters
    if (diag_.compatAlignsReclassifiedIntronic > 0 || diag_.compatAlignsLenientAccepted > 0 ||
        diag_.compatAlignsOverlapWeightApplied > 0 || diag_.compatPositionsSkippedOverlap > 0 ||
        diag_.compatPositionsSkippedTrim > 0) {
        out << "\nCompatModeCounters:\n";
        out << "compatAlignsReclassifiedIntronic\t" << diag_.compatAlignsReclassifiedIntronic << "\n";
        out << "compatAlignsLenientAccepted\t" << diag_.compatAlignsLenientAccepted << "\n";
        out << "compatAlignsOverlapWeightApplied\t" << diag_.compatAlignsOverlapWeightApplied << "\n";
        out << "compatPositionsSkippedOverlap\t" << diag_.compatPositionsSkippedOverlap << "\n";
        out << "compatPositionsSkippedTrim\t" << diag_.compatPositionsSkippedTrim << "\n";
    }
}

void SlamQuant::writeDebug(const Transcriptome& tr, double errorRate, double convRate) const {
    if (!debugEnabled_) {
        return;
    }
    std::string base = debugOutPrefix_.empty() ? "SlamQuant.debug" : debugOutPrefix_;
    if (!debugGeneStats_.empty()) {
        std::ofstream out((base + ".gene.tsv").c_str());
        if (out.good()) {
            out << "Gene\tSymbol\tReadCount\tConversions\tCoverage\tNTR\tk_over_nT\t"
                << "DebugExonicWeight\tDebugIntronicWeight\tDebugSenseWeight\tDebugAntisenseWeight\t"
                << "DebugReadsAssigned\tDebugReadsIntronic\tDropsSnpMask\tDropsStrandness\n";
            SlamSolver solver(errorRate, convRate);
            for (size_t i = 0; i < debugGeneStats_.size(); ++i) {
                if (debugGeneMask_.empty() || debugGeneMask_[i] == 0) {
                    continue;
                }
                const SlamGeneStats& stats = geneStats_[i];
                SlamResult res = solver.solve(stats.histogram);
                const std::string& geneId = tr.geID[i];
                const std::string& geneName = tr.geName[i].empty() ? tr.geID[i] : tr.geName[i];
                const SlamDebugGeneStats& dbg = debugGeneStats_[i];
                double kOverNT = (stats.coverage > 0.0) ? (stats.conversions / stats.coverage) : 0.0;
                out << geneId << "\t"
                    << geneName << "\t"
                    << stats.readCount << "\t"
                    << stats.conversions << "\t"
                    << stats.coverage << "\t"
                    << res.ntr << "\t"
                    << kOverNT << "\t"
                    << dbg.exonicWeight << "\t"
                    << dbg.intronicWeight << "\t"
                    << dbg.senseWeight << "\t"
                    << dbg.antisenseWeight << "\t"
                    << dbg.readsAssigned << "\t"
                    << dbg.readsAssignedIntronic << "\t"
                    << dbg.dropsSnpMask << "\t"
                    << dbg.dropsStrandness << "\n";
            }
        }
    }
    if (!debugReadRecords_.empty() && debugMaxReads_ != 0) {
        std::ofstream out((base + ".reads.tsv").c_str());
        if (out.good()) {
            out << "ReadName\tReadLoc\tGene\tSymbol\tStatus\tIntronic\tOppositeStrand\tWeight\t"
                << "nT\tk\tReadLength\tSnpBuffered\tConvReadPos\tConvGenPos\n";
            for (const auto& rec : debugReadRecords_) {
                std::string geneId = "NA";
                std::string geneName = "NA";
                if (rec.geneId < tr.geID.size()) {
                    geneId = tr.geID[rec.geneId];
                    const std::string& gname = tr.geName[rec.geneId];
                    geneName = gname.empty() ? geneId : gname;
                }
                out << rec.readName << "\t"
                    << rec.readLoc << "\t"
                    << geneId << "\t"
                    << geneName << "\t"
                    << debugDropReasonName(rec.status) << "\t"
                    << (rec.intronic ? 1 : 0) << "\t"
                    << (rec.oppositeStrand ? 1 : 0) << "\t"
                    << rec.weight << "\t"
                    << rec.nT << "\t"
                    << rec.k << "\t"
                    << rec.readLength << "\t"
                    << (rec.snpBuffered ? 1 : 0) << "\t"
                    << rec.convReadPos << "\t"
                    << rec.convGenPos << "\n";
            }
        }
    }

    if (debugSnpEnabled_ && !debugSnpCov_.empty()) {
        std::ofstream out((base + ".snp_sites.tsv").c_str());
        if (out.good()) {
            out << "loc\tabsPos0\twindow\toffset\tabsPos0_at_offset\t"
                << "cov\tany_mis\tconv_mis\t"
                << "cov_primary\tany_mis_primary\tconv_mis_primary\t"
                << "cov_nh1\tany_mis_nh1\tconv_mis_nh1\t"
                << "cov_nhgt1\tany_mis_nhgt1\tconv_mis_nhgt1\t"
                << "cov_mapq20\tany_mis_mapq20\tconv_mis_mapq20\n";
            for (int off = -debugSnpWindow_; off <= debugSnpWindow_; ++off) {
                size_t idx = static_cast<size_t>(off + debugSnpWindow_);
                uint64_t pos = static_cast<uint64_t>(static_cast<int64_t>(debugSnpAbsPos_) + off);
                out << debugSnpLoc_ << "\t"
                    << debugSnpAbsPos_ << "\t"
                    << debugSnpWindow_ << "\t"
                    << off << "\t"
                    << pos << "\t"
                    << debugSnpCov_[idx] << "\t"
                    << debugSnpAnyMis_[idx] << "\t"
                    << debugSnpConvMis_[idx] << "\t"
                    << debugSnpCovPrimary_[idx] << "\t"
                    << debugSnpAnyMisPrimary_[idx] << "\t"
                    << debugSnpConvMisPrimary_[idx] << "\t"
                    << debugSnpCovNh1_[idx] << "\t"
                    << debugSnpAnyMisNh1_[idx] << "\t"
                    << debugSnpConvMisNh1_[idx] << "\t"
                    << debugSnpCovNhGt1_[idx] << "\t"
                    << debugSnpAnyMisNhGt1_[idx] << "\t"
                    << debugSnpConvMisNhGt1_[idx] << "\t"
                    << debugSnpCovMapq20_[idx] << "\t"
                    << debugSnpAnyMisMapq20_[idx] << "\t"
                    << debugSnpConvMisMapq20_[idx] << "\n";
            }
        }
    }
}

void SlamQuant::writeTransitions(const std::string& outFile) const {
    std::ofstream out(outFile.c_str());
    if (!out.good()) {
        return;
    }
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    out << "Genomic\tRead\tCoverage\tMismatches\n";
    const SlamTransitionStats& stats = transitions_[static_cast<size_t>(SlamMismatchCategory::Exonic)];
    for (int g = 0; g < 4; ++g) {
        for (int r = 0; r < 4; ++r) {
            if (g == r) {
                continue;
            }
            out << bases[g] << "\t" << bases[r] << "\t"
                << stats.coverage[g] << "\t"
                << stats.mismatches[g][r] << "\n";
        }
    }
}

void SlamQuant::writeMismatches(const std::string& outFile, const std::string& condition) const {
    std::ofstream out(outFile.c_str());
    if (!out.good()) {
        // Try to create parent directory if it doesn't exist
        size_t lastSlash = outFile.find_last_of("/\\");
        if (lastSlash != std::string::npos) {
            std::string dir = outFile.substr(0, lastSlash);
            system(("mkdir -p " + dir).c_str());
            out.open(outFile.c_str());
        }
        if (!out.good()) {
            return;
        }
    }
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    out << "Category\tCondition\tOrientation\tGenomic\tRead\tCoverage\tMismatches\n";
    const char* labels[2] = {"First", "Second"};
    for (size_t c = 0; c < kSlamMismatchCategoryCount; ++c) {
        const char* category = slamMismatchCategoryName(static_cast<SlamMismatchCategory>(c));
        const SlamTransitionStats* orientations[2] = {&transitionsFirst_[c], &transitionsSecond_[c]};
        for (int o = 0; o < 2; ++o) {
            const SlamTransitionStats& stats = *orientations[o];
            for (int g = 0; g < 4; ++g) {
                for (int r = 0; r < 4; ++r) {
                    if (g == r) {
                        continue;
                    }
                    out << category << "\t" << condition << "\t" << labels[o] << "\t"
                        << bases[g] << "\t" << bases[r] << "\t"
                        << stats.coverage[g] << "\t" << stats.mismatches[g][r] << "\n";
                }
            }
        }
    }
}

void SlamQuant::writeMismatchDetails(const std::string& outFile) const {
    std::ofstream out(outFile.c_str(), std::ios::out | std::ios::trunc);
    if (!out.is_open()) {
        return;
    }
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    out << "Category\tGenomic\tRead\tPosition\tOverlap\tOpposite\tCoverage\tMismatches\n";
    for (size_t c = 0; c < kSlamMismatchCategoryCount; ++c) {
        const char* category = slamMismatchCategoryName(static_cast<SlamMismatchCategory>(c));
        const auto& positions = positionTransitions_[c];
        for (size_t pos = 0; pos < positions.size(); ++pos) {
            const SlamPositionStats& stats = positions[pos];
            for (int ov = 0; ov < 2; ++ov) {
                for (int opp = 0; opp < 2; ++opp) {
                    for (int g = 0; g < 4; ++g) {
                        double cov = stats.coverage[ov][opp][g];
                        if (cov <= 0.0) {
                            continue;
                        }
                        for (int r = 0; r < 4; ++r) {
                            if (g == r) {
                                continue;
                            }
                            out << category << "\t" << bases[g] << "\t" << bases[r] << "\t"
                                << pos << "\t" << ov << "\t" << opp << "\t" << cov << "\t"
                                << stats.mismatches[ov][opp][g][r] << "\n";
                        }
                    }
                }
            }
        }
    }
    out.close();
}

std::unordered_map<uint32_t, std::tuple<double, double, double, double>> SlamQuant::getPositionTransitionData() const {
    std::unordered_map<uint32_t, std::tuple<double, double, double, double>> result;
    
    // Extract ExonicSense category (index 1)
    size_t catIdx = static_cast<size_t>(SlamMismatchCategory::ExonicSense);
    const auto& positions = positionTransitions_[catIdx];
    
    // Genomic bases: A=0, C=1, G=2, T=3
    // Transcript T->C corresponds to:
    //   plus-strand: T->C
    //   minus-strand: A->G (complement of T->C)
    const int A_BASE = 0;
    const int C_BASE = 1;
    const int G_BASE = 2;
    const int T_BASE = 3;
    
    for (size_t pos = 0; pos < positions.size(); ++pos) {
        const SlamPositionStats& stats = positions[pos];
        
        double t_cov = 0.0;
        double a_cov = 0.0;
        double tc_mm = 0.0;
        double ta_mm = 0.0;
        
        // Sum across overlap and opposite dimensions
        for (int ov = 0; ov < 2; ++ov) {
            for (int opp = 0; opp < 2; ++opp) {
                t_cov += stats.coverage[ov][opp][T_BASE];
                a_cov += stats.coverage[ov][opp][A_BASE];
                
                // Transcript T→C: genomic T->C plus genomic A->G
                tc_mm += stats.mismatches[ov][opp][T_BASE][C_BASE];
                tc_mm += stats.mismatches[ov][opp][A_BASE][G_BASE];
                
                // Transcript T→A control: genomic T->A plus genomic A->T
                ta_mm += stats.mismatches[ov][opp][T_BASE][A_BASE];
                ta_mm += stats.mismatches[ov][opp][A_BASE][T_BASE];
            }
        }
        
        double total_cov = t_cov + a_cov;
        if (total_cov > 0.0) {
            result[static_cast<uint32_t>(pos)] = std::make_tuple(total_cov, tc_mm, total_cov, ta_mm);
        }
    }
    
    return result;
}

void SlamQuant::writeTopMismatches(const Transcriptome& tr, const std::string& refFile,
                                   const std::string& mismatchFile, size_t topN) const {
    // Parse reference file (handle gzipped files)
    std::unordered_map<std::string, std::map<std::string, double>> refData;
    std::ifstream ref(refFile.c_str(), std::ios::binary);
    if (!ref.good()) {
        return;  // Skip if reference file not found
    }
    
    // Check if gzipped by reading magic bytes
    unsigned char magic[2];
    ref.read(reinterpret_cast<char*>(magic), 2);
    ref.seekg(0);
    bool isGzip = (magic[0] == 0x1f && magic[1] == 0x8b);
    
    if (isGzip) {
        // For gzipped files, use zcat via popen
        std::string cmd = "zcat " + refFile;
        FILE* pipe = popen(cmd.c_str(), "r");
        if (!pipe) {
            return;
        }
        char buffer[4096];
        std::string line;
        bool headerSkipped = false;
        while (fgets(buffer, sizeof(buffer), pipe)) {
            line = buffer;
            if (line.back() == '\n') line.pop_back();
            if (!headerSkipped) {
                headerSkipped = true;
                continue;
            }
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string geneId, symbol, readCount, conversions, coverage, ntr;
            if (iss >> geneId >> symbol >> readCount >> conversions >> coverage >> ntr) {
                try {
                    refData[geneId]["ReadCount"] = std::stod(readCount);
                    refData[geneId]["Conversions"] = std::stod(conversions);
                    refData[geneId]["Coverage"] = std::stod(coverage);
                    refData[geneId]["NTR"] = std::stod(ntr);
                } catch (...) {
                    continue;  // Skip malformed lines
                }
            }
        }
        pclose(pipe);
    } else {
        // Plain text file
        std::string line;
        std::getline(ref, line);  // Skip header
        while (std::getline(ref, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string geneId, symbol, readCount, conversions, coverage, ntr;
            if (iss >> geneId >> symbol >> readCount >> conversions >> coverage >> ntr) {
                try {
                    refData[geneId]["ReadCount"] = std::stod(readCount);
                    refData[geneId]["Conversions"] = std::stod(conversions);
                    refData[geneId]["Coverage"] = std::stod(coverage);
                    refData[geneId]["NTR"] = std::stod(ntr);
                } catch (...) {
                    continue;  // Skip malformed lines
                }
            }
        }
    }

    // Collect mismatches
    std::vector<std::pair<double, size_t>> mismatches;
    for (size_t i = 0; i < geneStats_.size(); ++i) {
        const SlamGeneStats& stats = geneStats_[i];
        const std::string& geneId = tr.geID[i];
        auto it = refData.find(geneId);
        if (it != refData.end()) {
            double delta = std::abs(stats.readCount - it->second.at("ReadCount"));
            if (delta > 0.1) {
                mismatches.push_back({delta, i});
            }
        }
    }
    std::sort(mismatches.rbegin(), mismatches.rend());
    if (mismatches.size() > topN) {
        mismatches.resize(topN);
    }

    // Write mismatch file
    std::ofstream out(mismatchFile.c_str());
    if (!out.good()) {
        return;
    }
    out << "Gene\tSymbol\tReadCount_ref\tReadCount_test\tReadCount_delta\t"
        << "Conversions_ref\tConversions_test\tConversions_delta\t"
        << "Coverage_ref\tCoverage_test\tCoverage_delta\t"
        << "NTR_ref\tNTR_test\tNTR_delta\n";
    
    SlamSolver solver(0.001, 0.05);  // Default rates
    for (const auto& mismatch : mismatches) {
        size_t i = mismatch.second;
        const SlamGeneStats& stats = geneStats_[i];
        const std::string& geneId = tr.geID[i];
        const std::string& geneName = tr.geName[i].empty() ? tr.geID[i] : tr.geName[i];
        auto it = refData.find(geneId);
        if (it == refData.end()) continue;
        
        SlamResult res = solver.solve(stats.histogram);
        const auto& ref = it->second;
        
        out << geneId << "\t" << geneName << "\t"
            << ref.at("ReadCount") << "\t" << stats.readCount << "\t"
            << (stats.readCount - ref.at("ReadCount")) << "\t"
            << ref.at("Conversions") << "\t" << stats.conversions << "\t"
            << (stats.conversions - ref.at("Conversions")) << "\t"
            << ref.at("Coverage") << "\t" << stats.coverage << "\t"
            << (stats.coverage - ref.at("Coverage")) << "\t"
            << ref.at("NTR") << "\t" << res.ntr << "\t"
            << (res.ntr - ref.at("NTR")) << "\n";
    }
}

namespace {
bool isAcmgBase(char base) {
    switch (base) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            return true;
        default:
            return false;
    }
}

bool isSnpAllele(const char* allele) {
    if (allele == nullptr) {
        return false;
    }
    return allele[0] != '\0' && allele[1] == '\0' && isAcmgBase(allele[0]);
}

std::string normalizeContigName(const std::string& contig,
                                const std::unordered_map<std::string, uint32_t>& chrIndex) {
    if (chrIndex.find(contig) != chrIndex.end()) {
        return contig;
    }
    if (contig == "M" || contig == "MT") {
        if (chrIndex.find("chrM") != chrIndex.end()) {
            return "chrM";
        }
        if (chrIndex.find("MT") != chrIndex.end()) {
            return "MT";
        }
    }
    if (contig.size() > 3 && contig.rfind("chr", 0) == 0) {
        std::string stripped = contig.substr(3);
        if (chrIndex.find(stripped) != chrIndex.end()) {
            return stripped;
        }
    } else {
        std::string prefixed = "chr" + contig;
        if (chrIndex.find(prefixed) != chrIndex.end()) {
            return prefixed;
        }
    }
    return "";
}

bool writeSimpleBed(const std::string& bedPath,
                    const std::vector<std::tuple<uint32_t, uint64_t, char, std::string>>& entries,
                    const std::vector<std::string>& chrNames,
                    std::string* err) {
    BGZF* bgzf = bgzf_open(bedPath.c_str(), "w");
    if (!bgzf) {
        if (err) {
            *err = "Cannot open BED output file: " + bedPath;
        }
        return false;
    }
    std::string header = "#chrom\tstart\tend\tref\talt\n";
    bgzf_write(bgzf, header.c_str(), header.size());
    for (const auto& entry : entries) {
        uint32_t chrIdx = std::get<0>(entry);
        uint64_t pos0 = std::get<1>(entry);
        char ref = std::get<2>(entry);
        const std::string& alt = std::get<3>(entry);
        if (chrIdx >= chrNames.size()) {
            continue;
        }
        std::ostringstream oss;
        oss << chrNames[chrIdx] << "\t"
            << pos0 << "\t"
            << (pos0 + 1) << "\t"
            << ref << "\t"
            << (alt.empty() ? "." : alt) << "\n";
        std::string line = oss.str();
        bgzf_write(bgzf, line.c_str(), line.size());
    }
    bgzf_close(bgzf);
    int ret = tbx_index_build(bedPath.c_str(), 0, &tbx_conf_bed);
    if (ret != 0) {
        if (err) {
            *err = "Failed to build tabix index for: " + bedPath;
        }
        return false;
    }
    return true;
}

bool writeVcfSummary(const std::string& summaryPath,
                     const SlamSnpMaskVcfStats& stats,
                     std::string* err) {
    std::ofstream out(summaryPath.c_str());
    if (!out.good()) {
        if (err) {
            *err = "Cannot open summary output file: " + summaryPath;
        }
        return false;
    }
    out << "metric\tvalue\n";
    out << "records_total\t" << stats.recordsTotal << "\n";
    out << "records_filtered\t" << stats.recordsFiltered << "\n";
    out << "records_no_snp_alt\t" << stats.recordsNoSnpAlt << "\n";
    out << "records_non_snp_alt\t" << stats.recordsNonSnpAlt << "\n";
    out << "records_missing_gt\t" << stats.recordsMissingGt << "\n";
    out << "records_no_alt_gt\t" << stats.recordsNoAltGt << "\n";
    out << "records_unknown_contig\t" << stats.recordsUnknownContig << "\n";
    out << "sites_added\t" << stats.sitesAdded << "\n";
    out << "sites_duplicate\t" << stats.sitesDuplicate << "\n";
    return true;
}
template <typename T>
bool loadBedInternal(const std::string& path,
                     const std::vector<std::string>& chrNames,
                     const std::vector<T>& chrStart,
                     std::unordered_set<uint64_t>& positions,
                     std::string* err) {
    positions.clear();
    std::unordered_map<std::string, uint32_t> chrMap;
    chrMap.reserve(chrNames.size());
    for (uint32_t i = 0; i < chrNames.size(); ++i) {
        chrMap[chrNames[i]] = i;
    }

    auto parseLine = [&](const std::string& line) {
        if (line.empty() || line[0] == '#') return;
        std::istringstream iss(line);
        std::string chr;
        uint64_t start = 0;
        uint64_t end = 0;
        if (!(iss >> chr >> start >> end)) {
            return;
        }
        auto it = chrMap.find(chr);
        if (it == chrMap.end()) {
            return;
        }
        if (it->second >= chrStart.size()) {
            return;
        }
        uint64_t cstart = static_cast<uint64_t>(chrStart[it->second]);
        if (end <= start) {
            return;
        }
        uint64_t maxLen = end - start;
        if (maxLen > 1000) {
            return;
        }
        for (uint64_t pos = start; pos < end; ++pos) {
            positions.insert(cstart + pos);
        }
    };

    const bool isGz = (path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0);
    if (isGz) {
        BGZF* fp = bgzf_open(path.c_str(), "r");
        if (!fp) {
            if (err) *err = "Failed to open SNP BED (bgzip/gzip): " + path;
            return false;
        }
        kstring_t ks = {0, 0, nullptr};
        while (bgzf_getline(fp, '\n', &ks) >= 0) {
            parseLine(std::string(ks.s, ks.l));
        }
        if (ks.s) {
            free(ks.s);
        }
        bgzf_close(fp);
    } else {
        std::ifstream in(path.c_str());
        if (!in.good()) {
            if (err) *err = "Failed to open SNP BED: " + path;
            return false;
        }
        std::string line;
        while (std::getline(in, line)) {
            parseLine(line);
        }
    }
    return true;
}
} // namespace

bool SlamSnpMask::loadBed(const std::string& path, const Genome& genome, std::string* err) {
    return loadBedInternal(path, genome.chrName, genome.chrStart, positions_, err);
}

bool SlamSnpMask::loadBedWithChrMap(const std::string& path,
                                    const std::vector<std::string>& chrNames,
                                    const std::vector<uint64_t>& chrStart,
                                    std::string* err) {
    return loadBedInternal(path, chrNames, chrStart, positions_, err);
}

bool SlamSnpMask::loadVcf(const std::string& path,
                          const Genome& genome,
                          const SlamSnpMaskVcfOptions& opts,
                          SlamSnpMaskVcfStats* stats,
                          std::string* err) {
    std::vector<uint64_t> chrStart64(genome.chrStart.begin(), genome.chrStart.end());
    return loadVcfWithChrMap(path, genome.chrName, chrStart64, opts, stats, err);
}

bool SlamSnpMask::loadVcfWithChrMap(const std::string& path,
                                    const std::vector<std::string>& chrNames,
                                    const std::vector<uint64_t>& chrStart,
                                    const SlamSnpMaskVcfOptions& opts,
                                    SlamSnpMaskVcfStats* stats,
                                    std::string* err) {
    if (stats) {
        *stats = SlamSnpMaskVcfStats();
    }
    std::unordered_map<std::string, uint32_t> chrIndex;
    chrIndex.reserve(chrNames.size());
    for (uint32_t i = 0; i < chrNames.size(); ++i) {
        chrIndex[chrNames[i]] = i;
    }

    htsFile* fp = hts_open(path.c_str(), "r");
    if (!fp) {
        if (err) {
            *err = "Failed to open VCF: " + path;
        }
        return false;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        if (err) {
            *err = "Failed to read VCF header: " + path;
        }
        hts_close(fp);
        return false;
    }

    int nsamples = bcf_hdr_nsamples(hdr);
    int sampleIndex = 0;
    if (opts.mode == "gt") {
        if (!opts.sample.empty()) {
            sampleIndex = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, opts.sample.c_str());
            if (sampleIndex < 0) {
                if (err) {
                    *err = "Sample not found in VCF header: " + opts.sample;
                }
                bcf_hdr_destroy(hdr);
                hts_close(fp);
                return false;
            }
        } else {
            if (nsamples == 1) {
                sampleIndex = 0;
            } else {
                if (err) {
                    *err = "VCF has multiple samples; specify --slamSnpMaskVcfSample";
                }
                bcf_hdr_destroy(hdr);
                hts_close(fp);
                return false;
            }
        }
    }

    bcf1_t* rec = bcf_init();
    int* gt = nullptr;
    int ngt = 0;
    std::vector<std::tuple<uint32_t, uint64_t, char, std::string>> bedEntries;

    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_STR | BCF_UN_FLT);
        if (stats) {
            stats->recordsTotal++;
        }
        if (opts.filter == "pass") {
            bool pass = (rec->d.n_flt == 0) || (bcf_has_filter(hdr, rec, const_cast<char*>("PASS")) > 0);
            if (!pass) {
                if (stats) {
                    stats->recordsFiltered++;
                }
                continue;
            }
        }

        if (rec->n_allele <= 1) {
            if (stats) {
                stats->recordsNoSnpAlt++;
            }
            continue;
        }

        const char* ref = rec->d.allele[0];
        if (!isSnpAllele(ref)) {
            if (stats) {
                stats->recordsNonSnpAlt++;
            }
            continue;
        }

        std::vector<int> snpAltIdx;
        std::vector<std::string> snpAltList;
        for (int a = 1; a < rec->n_allele; ++a) {
            const char* alt = rec->d.allele[a];
            if (isSnpAllele(alt) && alt[0] != ref[0]) {
                snpAltIdx.push_back(a);
                snpAltList.emplace_back(alt);
            }
        }
        if (snpAltIdx.empty()) {
            if (stats) {
                stats->recordsNonSnpAlt++;
            }
            continue;
        }

        if (opts.mode == "gt") {
            if (nsamples <= 0) {
                if (stats) {
                    stats->recordsMissingGt++;
                }
                continue;
            }
            int ngt_ret = bcf_get_genotypes(hdr, rec, &gt, &ngt);
            if (ngt_ret <= 0) {
                if (stats) {
                    stats->recordsMissingGt++;
                }
                continue;
            }
            int ploidy = ngt_ret / nsamples;
            const int* sampleGt = gt + sampleIndex * ploidy;
            bool missing = false;
            bool hasAlt = false;
            bool hasSnpAlt = false;
            for (int i = 0; i < ploidy; ++i) {
                if (sampleGt[i] == bcf_gt_missing) {
                    missing = true;
                    continue;
                }
                int allele = bcf_gt_allele(sampleGt[i]);
                if (allele > 0) {
                    hasAlt = true;
                    if (std::find(snpAltIdx.begin(), snpAltIdx.end(), allele) != snpAltIdx.end()) {
                        hasSnpAlt = true;
                    }
                }
            }
            if (missing) {
                if (stats) {
                    stats->recordsMissingGt++;
                }
                continue;
            }
            if (!hasAlt || !hasSnpAlt) {
                if (stats) {
                    stats->recordsNoAltGt++;
                }
                continue;
            }
        }

        const char* contig = bcf_hdr_id2name(hdr, rec->rid);
        std::string normContig = normalizeContigName(contig ? contig : "", chrIndex);
        if (normContig.empty()) {
            if (stats) {
                stats->recordsUnknownContig++;
            }
            continue;
        }
        uint32_t chrIdx = chrIndex[normContig];
        if (chrIdx >= chrStart.size()) {
            if (stats) {
                stats->recordsUnknownContig++;
            }
            continue;
        }
        if (rec->pos < 0) {
            continue;
        }
        uint64_t pos0 = static_cast<uint64_t>(rec->pos);
        uint64_t gpos = chrStart[chrIdx] + pos0;
        auto inserted = positions_.insert(gpos);
        if (!inserted.second) {
            if (stats) {
                stats->sitesDuplicate++;
            }
            continue;
        }
        if (stats) {
            stats->sitesAdded++;
        }
        std::string altJoined;
        for (size_t i = 0; i < snpAltList.size(); ++i) {
            if (i > 0) {
                altJoined += ",";
            }
            altJoined += snpAltList[i];
        }
        bedEntries.emplace_back(chrIdx, pos0, ref[0], altJoined);
    }

    if (gt) {
        free(gt);
    }
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);

    if (!opts.bedOut.empty()) {
        std::sort(bedEntries.begin(), bedEntries.end(),
                  [](const std::tuple<uint32_t, uint64_t, char, std::string>& a,
                     const std::tuple<uint32_t, uint64_t, char, std::string>& b) {
                      if (std::get<0>(a) != std::get<0>(b)) {
                          return std::get<0>(a) < std::get<0>(b);
                      }
                      return std::get<1>(a) < std::get<1>(b);
                  });
        if (!writeSimpleBed(opts.bedOut, bedEntries, chrNames, err)) {
            return false;
        }
    }
    if (!opts.summaryOut.empty() && stats) {
        if (!writeVcfSummary(opts.summaryOut, *stats, err)) {
            return false;
        }
    }
    return true;
}

// Read buffer methods for auto-trim replay
void SlamQuant::enableReadBuffer(uint64_t maxReads) {
    readBuffer_.reset(new SlamReadBuffer(maxReads));
}

bool SlamQuant::bufferRead(SlamBufferedRead&& read) {
    if (!readBuffer_) {
        return false;
    }
    return readBuffer_->addRead(std::move(read));
}

void SlamQuant::clearReadBuffer() {
    if (readBuffer_) {
        readBuffer_->clear();
    }
}

void SlamQuant::enableDumpBuffer(uint64_t maxReads) {
    dumpBuffer_.reset(new SlamReadBuffer(maxReads));
}

bool SlamQuant::bufferDumpRead(SlamBufferedRead&& read) {
    if (!dumpBuffer_) {
        return false;
    }
    return dumpBuffer_->addRead(std::move(read));
}

// Replay buffered reads with trim applied
// This is the core replay logic - processes buffered reads through the same
// counting paths but with trim filtering applied
uint64_t SlamQuant::replayBufferedReads(SlamCompat* compat, const SlamSnpMask* snpMask, int strandness) {
    if (!readBuffer_ || readBuffer_->size() == 0) {
        return 0;
    }
    
    uint64_t replayed = 0;
    
    for (const SlamBufferedRead& read : readBuffer_->reads()) {
        // Skip if no genes assigned
        if (read.geneIds.empty()) {
            diag_.readsZeroGenes++;
            continue;
        }
        
        // Apply strandness filter (same logic as ReadAlign_slamQuant)
        if (strandness == 1 && read.oppositeStrand) {
            diag_.readsDroppedStrandness++;
            continue;
        }
        if (strandness == 2 && !read.oppositeStrand) {
            diag_.readsDroppedStrandness++;
            continue;
        }
        
        // Apply SNP mask filter (check if any buffered genomic positions overlap mask)
        if (snpMask != nullptr) {
            bool hasSnpOverlap = false;
            for (const SlamBufferedPosition& pos : read.positions) {
                if (snpMask->contains(pos.genomicPos)) {
                    hasSnpOverlap = true;
                    break;
                }
            }
            if (hasSnpOverlap) {
                diag_.readsDroppedSnpMask++;
                continue;
            }
        }
        
        // Determine mismatch category
        SlamMismatchCategory category = read.isIntronic ? SlamMismatchCategory::Intronic : SlamMismatchCategory::Exonic;
        SlamMismatchCategory senseCategory = read.isIntronic ? SlamMismatchCategory::IntronicSense : SlamMismatchCategory::ExonicSense;
        
        // Count T bases and conversions with trim filtering
        uint16_t nT = 0;
        uint16_t k = 0;
        std::vector<uint32_t> mismatchPositions;
        bool snpDetect = snpDetectEnabled_ && !read.isIntronic;
        
        for (const SlamBufferedPosition& pos : read.positions) {
            // Apply trim filtering via SlamCompat
            if (compat) {
                // Convert readPos to mate-local coordinates
                uint32_t mateLocalPos;
                uint32_t mateLen;
                if (pos.secondMate) {
                    mateLocalPos = pos.readPos - read.readLength0;
                    mateLen = read.readLength1;
                } else {
                    mateLocalPos = pos.readPos;
                    mateLen = read.readLength0;
                }
                
                // Check overlap skip (if enabled)
                if (compat->cfg().ignoreOverlap && pos.overlap) {
                    diag_.compatPositionsSkippedOverlap++;
                    continue;
                }
                
                // Check trim guards
                if (!compat->compatShouldCountPos(mateLocalPos, mateLen)) {
                    diag_.compatPositionsSkippedTrim++;
                    continue;
                }
            }
            
            // Record transition base (same logic as ReadAlign_slamQuant)
            bool skipMismatch = pos.overlap && pos.secondMate;
            if (!skipMismatch) {
                addTransitionBase(category, pos.readPos, pos.secondMate, pos.overlap, 
                                  read.oppositeStrand, pos.refBase, pos.readBase, read.weight);
                if (!read.oppositeStrand) {
                    addTransitionBase(senseCategory, pos.readPos, pos.secondMate, pos.overlap,
                                      false, pos.refBase, pos.readBase, read.weight);
                }
            }
            
            // Count T bases and conversions (for EM histogram)
            if (!read.isIntronic) {
                bool isT = false;
                bool isTc = false;
                if (!read.isMinus) {
                    isT = (pos.refBase == 3); // T
                    isTc = (pos.refBase == 3 && pos.readBase == 1); // T→C
                } else {
                    isT = (pos.refBase == 0); // A (complement of T)
                    isTc = (pos.refBase == 0 && pos.readBase == 2); // A→G (complement of T→C)
                }
                
                if (isT) {
                    ++nT;
                    if (isTc) {
                        ++k;
                        if (snpDetect) {
                            mismatchPositions.push_back(static_cast<uint32_t>(pos.genomicPos));
                        }
                    }
                }
            }
        }
        
        // Add to gene counts (same logic as ReadAlign_slamQuant)
        if (!read.isIntronic) {
            uint8_t k8 = static_cast<uint8_t>(k > 255 ? 255 : k);
            
            if (snpDetect && !mismatchPositions.empty()) {
                for (uint32_t geneId : read.geneIds) {
                    bufferSnpRead(geneId, nT, mismatchPositions, read.weight);
                }
            } else {
                for (uint32_t geneId : read.geneIds) {
                    addRead(geneId, nT, k8, read.weight);
                }
            }
        }
        
        diag_.readsProcessed++;
        ++replayed;
    }
    
    return replayed;
}
