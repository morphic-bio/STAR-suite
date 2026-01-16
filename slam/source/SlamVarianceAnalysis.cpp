#include "SlamVarianceAnalysis.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

SlamVarianceAnalyzer::SlamVarianceAnalyzer(uint32_t maxReads, uint32_t minReads,
                                           uint32_t smoothWindow, uint32_t minSegLen, uint32_t maxTrim)
    : readsAnalyzed_(0), maxReads_(maxReads), minReads_(minReads),
      smoothWindow_(smoothWindow), minSegLen_(minSegLen), maxTrim_(maxTrim) {
}

bool SlamVarianceAnalyzer::recordRead() {
    // Check if we've reached max reads
    if (maxReads_ > 0 && readsAnalyzed_ >= maxReads_) {
        return false; // Stop collecting
    }
    readsAnalyzed_++;
    return true; // Continue collecting
}

void SlamVarianceAnalyzer::recordPosition(uint32_t readPos, uint8_t qual, bool isT, bool isTc) {
    // Gate on readsAnalyzed_ <= maxReads_ (allow Nth read's positions)
    if (maxReads_ > 0 && readsAnalyzed_ > maxReads_) {
        return; // Stop collecting positions
    }
    
    auto& stats = positionStats_[readPos];
    stats.readCount++;
    
    // Record quality statistics
    double q = static_cast<double>(qual);
    stats.qualSum += q;
    stats.qualCount++;
    stats.qualSumSq += q * q;
    
    // Record T→C statistics
    if (isT) {
        stats.tCount++;
        if (isTc) {
            stats.tcCount++;
        }
        // Track per-T-base T→C indicator (0 or 1) for variance calculation
        double tcRate = isTc ? 1.0 : 0.0;
        stats.tcRateSum += tcRate;
        stats.tcRateSumSq += tcRate * tcRate;  // For variance: sum of squares
        stats.tcRateCount++;
    }
}

void SlamVarianceAnalyzer::reset() {
    positionStats_.clear();
    readsAnalyzed_ = 0;
}

void SlamVarianceAnalyzer::merge(const SlamVarianceAnalyzer& other) {
    // Merge read count (sum across threads, but cap at maxReads_)
    readsAnalyzed_ = std::min(readsAnalyzed_ + other.readsAnalyzed_, 
                              static_cast<uint64_t>(maxReads_ > 0 ? maxReads_ : UINT64_MAX));
    
    // Merge position stats
    for (const auto& kv : other.positionStats_) {
        auto& stats = positionStats_[kv.first];
        const auto& otherStats = kv.second;
        
        stats.readCount += otherStats.readCount;
        stats.tCount += otherStats.tCount;
        stats.tcCount += otherStats.tcCount;
        stats.qualSum += otherStats.qualSum;
        stats.qualCount += otherStats.qualCount;
        stats.qualSumSq += otherStats.qualSumSq;
        stats.tcRateSum += otherStats.tcRateSum;
        stats.tcRateSumSq += otherStats.tcRateSumSq;
        stats.tcRateCount += otherStats.tcRateCount;
    }
}

std::tuple<uint64_t, uint64_t, double> SlamVarianceAnalyzer::computeGlobalTcErrorRate(int trim5p, int trim3p, uint32_t readLength) const {
    uint64_t t_total = 0;
    uint64_t tc_total = 0;
    
    for (const auto& kv : positionStats_) {
        uint32_t pos = kv.first;
        const auto& stats = kv.second;
        
        // Apply trim window if provided
        if (trim5p > 0 && static_cast<int>(pos) < trim5p) {
            continue;  // Before trim5p window
        }
        if (trim3p > 0 && readLength > 0 && static_cast<int>(pos) >= static_cast<int>(readLength) - trim3p) {
            continue;  // After trim3p window
        }
        
        t_total += stats.tCount;
        tc_total += stats.tcCount;
    }
    
    double p_est = 0.0;
    if (t_total > 0) {
        p_est = static_cast<double>(tc_total) / static_cast<double>(t_total);
    }
    
    return std::make_tuple(t_total, tc_total, p_est);
}

std::vector<double> SlamVarianceAnalyzer::smoothMedian(const std::vector<double>& values, uint32_t window) {
    if (window <= 1 || values.empty()) {
        return values;
    }
    
    std::vector<double> result(values.size());
    int half = static_cast<int>(window / 2);
    
    for (size_t i = 0; i < values.size(); ++i) {
        int start = std::max(0, static_cast<int>(i) - half);
        int end = std::min(static_cast<int>(values.size()), static_cast<int>(i) + half + 1);
        
        // Collect non-NaN values in window
        std::vector<double> chunk;
        for (int j = start; j < end; ++j) {
            if (!std::isnan(values[j])) {
                chunk.push_back(values[j]);
            }
        }
        
        if (chunk.empty()) {
            result[i] = std::nan("");
        } else {
            // Compute median
            std::sort(chunk.begin(), chunk.end());
            size_t mid = chunk.size() / 2;
            if (chunk.size() % 2 == 0) {
                result[i] = (chunk[mid - 1] + chunk[mid]) / 2.0;
            } else {
                result[i] = chunk[mid];
            }
        }
    }
    
    return result;
}

void SlamVarianceAnalyzer::interpolateMissing(std::vector<double>& values) {
    if (values.empty()) return;
    
    // Find valid indices
    std::vector<size_t> validIdx;
    for (size_t i = 0; i < values.size(); ++i) {
        if (!std::isnan(values[i])) {
            validIdx.push_back(i);
        }
    }
    
    if (validIdx.size() < 2) {
        // Not enough valid values for interpolation, fill with 0
        for (auto& v : values) {
            if (std::isnan(v)) v = 0.0;
        }
        return;
    }
    
    // Linear interpolation for NaN values
    for (size_t i = 0; i < values.size(); ++i) {
        if (std::isnan(values[i])) {
            // Find surrounding valid indices
            size_t left = 0, right = validIdx.size() - 1;
            for (size_t j = 0; j < validIdx.size(); ++j) {
                if (validIdx[j] < i) left = j;
                if (validIdx[j] > i) { right = j; break; }
            }
            
            size_t li = validIdx[left];
            size_t ri = validIdx[right];
            
            if (li == ri) {
                values[i] = values[li];
            } else {
                // Linear interpolation
                double t = static_cast<double>(i - li) / static_cast<double>(ri - li);
                values[i] = values[li] + t * (values[ri] - values[li]);
            }
        }
    }
}

SegmentFit SlamVarianceAnalyzer::fitSegment(
    const std::vector<double>& prefixN,
    const std::vector<double>& prefixX,
    const std::vector<double>& prefixXX,
    const std::vector<double>& prefixY,
    const std::vector<double>& prefixXY,
    const std::vector<double>& prefixYY,
    const std::vector<double>& y,
    uint32_t start, uint32_t end) {
    
    SegmentFit fit;
    
    if (end < start) {
        return fit;
    }
    
    double n = prefixN[end + 1] - prefixN[start];
    if (n < 2.0) {
        fit.intercept = 0.0;
        fit.slope = 0.0;
        fit.sse = std::numeric_limits<double>::max();
        return fit;
    }
    
    // Use prefix sums: sum from [start, end] = prefix[end+1] - prefix[start]
    double Sx = prefixX[end + 1] - prefixX[start];
    double Sxx = prefixXX[end + 1] - prefixXX[start];
    double Sy = prefixY[end + 1] - prefixY[start];
    double Sxy = prefixXY[end + 1] - prefixXY[start];
    double Syy = prefixYY[end + 1] - prefixYY[start];
    
    double den = n * Sxx - Sx * Sx;
    if (std::abs(den) < 1e-9) {
        fit.slope = 0.0;
    } else {
        fit.slope = (n * Sxy - Sx * Sy) / den;
    }
    fit.intercept = (Sy - fit.slope * Sx) / n;
    
    // SSE = sum(y^2) + m^2*sum(x^2) + n*b^2 + 2*m*b*sum(x) - 2*m*sum(xy) - 2*b*sum(y)
    fit.sse = Syy + fit.slope * fit.slope * Sxx + n * fit.intercept * fit.intercept
              + 2.0 * fit.slope * fit.intercept * Sx
              - 2.0 * fit.slope * Sxy
              - 2.0 * fit.intercept * Sy;
    
    // Clamp negative SSE (numerical precision)
    if (fit.sse < 0.0) fit.sse = 0.0;
    
    return fit;
}

std::tuple<uint32_t, uint32_t, uint32_t, double, SegmentFit, SegmentFit, SegmentFit, SegmentFit, uint32_t>
SlamVarianceAnalyzer::segmentedRegression(const std::vector<double>& y, uint32_t minSegLen, uint32_t maxTrim) {
    uint32_t n = static_cast<uint32_t>(y.size());
    
    // Default return: no breakpoints found
    if (n < minSegLen * 3 + 1) {
        return {0, n - 1, n - 1, std::numeric_limits<double>::max(), {}, {}, {}, {}, 0};
    }
    
    // Build prefix sums for efficient segment fitting
    std::vector<double> prefixN(n + 1, 0.0);
    std::vector<double> prefixX(n + 1, 0.0);
    std::vector<double> prefixXX(n + 1, 0.0);
    std::vector<double> prefixY(n + 1, 0.0);
    std::vector<double> prefixXY(n + 1, 0.0);
    std::vector<double> prefixYY(n + 1, 0.0);
    
    for (uint32_t i = 0; i < n; ++i) {
        double x = static_cast<double>(i);
        double yi = y[i];
        double valid = std::isnan(yi) ? 0.0 : 1.0;
        double yv = std::isnan(yi) ? 0.0 : yi;
        prefixN[i + 1] = prefixN[i] + valid;
        prefixX[i + 1] = prefixX[i] + (valid * x);
        prefixXX[i + 1] = prefixXX[i] + (valid * x * x);
        prefixY[i + 1] = prefixY[i] + yv;
        prefixXY[i + 1] = prefixXY[i] + (valid * x * yv);
        prefixYY[i + 1] = prefixYY[i] + (valid * yv * yv);
    }
    
    // Candidate breakpoint ranges
    // b1: first breakpoint (end of segment 1, segment 2 starts at b1)
    // b2: second breakpoint (end of segment 2, segment 3 starts at b2+1)
    // Segments: [0, b1-1], [b1, b2], [b2+1, n-1]
    
    uint32_t b1Min = minSegLen;
    uint32_t b1Max = std::min(maxTrim, n - 2 * minSegLen - 1);
    uint32_t b2MinFloor = std::max(minSegLen - 1, n - 1 - maxTrim);
    uint32_t b2Max = n - minSegLen - 1;
    
    if (b1Max < b1Min || b2MinFloor > b2Max) {
        return {0, n - 1, n - 1, std::numeric_limits<double>::max(), {}, {}, {}, {}, 0};
    }
    
    double bestSSE = std::numeric_limits<double>::max();
    uint32_t bestB1 = 0, bestB2 = n - 1;
    SegmentFit bestSeg1, bestSeg2, bestSeg3;
    
    for (uint32_t b1 = b1Min; b1 <= b1Max; ++b1) {
        uint32_t b2Min = std::max(b1 + minSegLen - 1, b2MinFloor);
        
        for (uint32_t b2 = b2Min; b2 <= b2Max; ++b2) {
            // Segment 1: [0, b1-1]
            SegmentFit seg1 = fitSegment(prefixN, prefixX, prefixXX, prefixY, prefixXY, prefixYY, y, 0, b1 - 1);
            // Segment 2: [b1, b2]
            SegmentFit seg2 = fitSegment(prefixN, prefixX, prefixXX, prefixY, prefixXY, prefixYY, y, b1, b2);
            // Segment 3: [b2+1, n-1]
            SegmentFit seg3 = fitSegment(prefixN, prefixX, prefixXX, prefixY, prefixXY, prefixYY, y, b2 + 1, n - 1);
            
            double totalSSE = seg1.sse + seg2.sse + seg3.sse;
            
            if (totalSSE < bestSSE) {
                bestSSE = totalSSE;
                bestB1 = b1;
                bestB2 = b2;
                bestSeg1 = seg1;
                bestSeg2 = seg2;
                bestSeg3 = seg3;
            }
        }
    }
    
    // Best 2-segment fit (single breakpoint).
    double bestSSE2 = std::numeric_limits<double>::max();
    uint32_t bestB1_2 = 0;
    SegmentFit bestSeg1_2, bestSeg2_2;
    uint32_t b1Min2 = minSegLen;
    uint32_t b1Max2 = n - minSegLen - 1;
    for (uint32_t b1 = b1Min2; b1 <= b1Max2; ++b1) {
        SegmentFit seg1 = fitSegment(prefixN, prefixX, prefixXX, prefixY, prefixXY, prefixYY, y, 0, b1 - 1);
        SegmentFit seg2 = fitSegment(prefixN, prefixX, prefixXX, prefixY, prefixXY, prefixYY, y, b1, n - 1);
        double totalSSE = seg1.sse + seg2.sse;
        if (totalSSE < bestSSE2) {
            bestSSE2 = totalSSE;
            bestB1_2 = b1;
            bestSeg1_2 = seg1;
            bestSeg2_2 = seg2;
        }
    }

    auto bic = [](double rss, uint32_t nPoints, uint32_t kParams) {
        const double eps = 1e-12;
        double rssAdj = (rss <= 0.0) ? eps : rss;
        return static_cast<double>(nPoints) * std::log(rssAdj / static_cast<double>(nPoints)) +
               static_cast<double>(kParams) * std::log(static_cast<double>(nPoints));
    };

    // Best 4-segment fit (three breakpoints).
    double bestSSE4 = std::numeric_limits<double>::max();
    uint32_t bestB1_4 = 0, bestB2_4 = 0, bestB3_4 = n - 1;
    SegmentFit bestSeg1_4, bestSeg2_4, bestSeg3_4, bestSeg4_4;
    uint32_t b1Max4 = std::min(maxTrim, n - 3 * minSegLen - 1);
    if (b1Max4 >= b1Min) {
        for (uint32_t b1 = b1Min; b1 <= b1Max4; ++b1) {
            uint32_t b2Min = b1 + minSegLen - 1;
            uint32_t b2Max = n - 2 * minSegLen - 1;
            for (uint32_t b2 = b2Min; b2 <= b2Max; ++b2) {
                uint32_t b3Min = std::max(b2 + minSegLen - 1, n - 1 - maxTrim);
                uint32_t b3Max = n - minSegLen - 1;
                for (uint32_t b3 = b3Min; b3 <= b3Max; ++b3) {
                    SegmentFit seg1 = fitSegment(prefixN, prefixX, prefixXX, prefixY, prefixXY, prefixYY, y, 0, b1 - 1);
                    SegmentFit seg2 = fitSegment(prefixN, prefixX, prefixXX, prefixY, prefixXY, prefixYY, y, b1, b2);
                    SegmentFit seg3 = fitSegment(prefixN, prefixX, prefixXX, prefixY, prefixXY, prefixYY, y, b2 + 1, b3);
                    SegmentFit seg4 = fitSegment(prefixN, prefixX, prefixXX, prefixY, prefixXY, prefixYY, y, b3 + 1, n - 1);
                    double totalSSE = seg1.sse + seg2.sse + seg3.sse + seg4.sse;
                    if (totalSSE < bestSSE4) {
                        bestSSE4 = totalSSE;
                        bestB1_4 = b1;
                        bestB2_4 = b2;
                        bestB3_4 = b3;
                        bestSeg1_4 = seg1;
                        bestSeg2_4 = seg2;
                        bestSeg3_4 = seg3;
                        bestSeg4_4 = seg4;
                    }
                }
            }
        }
    }

    // Compare 2-, 3-, 4-segment models with BIC (penalize extra parameters).
    // 2 segments: 2 lines (slope+intercept each) => k=4
    // 3 segments: 3 lines (slope+intercept each) => k=6
    // 4 segments: 4 lines (slope+intercept each) => k=8
    double nValid = prefixN[n];
    if (nValid < 2.0) {
        return {0, n - 1, n - 1, std::numeric_limits<double>::max(), {}, {}, {}, {}, 0};
    }
    double bic2 = bic(bestSSE2, static_cast<uint32_t>(nValid), 4);
    double bic3 = bic(bestSSE, static_cast<uint32_t>(nValid), 6);
    double bic4 = bic(bestSSE4, static_cast<uint32_t>(nValid), 8);
    if (bic2 <= bic3 && bic2 <= bic4) {
        // Use 2 segments; set b2/b3 to end so trim3p=0.
        return {bestB1_2, n - 1, n - 1, bestSSE2, bestSeg1_2, bestSeg2_2, {}, {}, 2};
    }
    if (bic4 <= bic3) {
        return {bestB1_4, bestB2_4, bestB3_4, bestSSE4, bestSeg1_4, bestSeg2_4, bestSeg3_4, bestSeg4_4, 4};
    }

    return {bestB1, bestB2, bestB2, bestSSE, bestSeg1, bestSeg2, bestSeg3, {}, 3};
}

SlamVarianceTrimResult SlamVarianceAnalyzer::computeTrim(uint32_t readLength) {
    SlamVarianceTrimResult result;
    result.readsAnalyzed = readsAnalyzed_;
    
    if (readsAnalyzed_ < minReads_) {
        result.mode = "auto_fallback";
        result.success = false;
        return result;
    }
    
    if (positionStats_.empty()) {
        result.mode = "auto_fallback";
        result.success = false;
        return result;
    }
    
    // Build stdev curve from position statistics
    std::vector<double> stdevCurve(readLength, std::nan(""));
    uint32_t maxPos = 0;
    
    for (const auto& kv : positionStats_) {
        uint32_t pos = kv.first;
        if (pos >= readLength) {
            continue;
        }
        const auto& stats = kv.second;
        
        // Use T→C rate stdev as the signal
        double tcStddev = stats.stddevTcRate();
        stdevCurve[pos] = tcStddev;
        
        if (pos > maxPos) {
            maxPos = pos;
        }
    }
    
    if (maxPos == 0) {
        result.mode = "auto_fallback";
        result.success = false;
        return result;
    }
    
    // Truncate to actual read length observed
    stdevCurve.resize(maxPos + 1);
    
    // Identify first/last non-zero on the UNSMOOTHED curve (ignore soft-trim zeros).
    auto isZeroLike = [](double v) {
        return std::isnan(v) || v == 0.0;
    };
    size_t firstNonZero = 0;
    while (firstNonZero < stdevCurve.size() && isZeroLike(stdevCurve[firstNonZero])) {
        ++firstNonZero;
    }
    size_t lastNonZero = stdevCurve.size();
    while (lastNonZero > firstNonZero && isZeroLike(stdevCurve[lastNonZero - 1])) {
        --lastNonZero;
    }
    if (firstNonZero >= lastNonZero) {
        result.mode = "auto_fallback";
        result.success = false;
        return result;
    }

    // Build a reduced vector for fitting: only points within [firstNonZero, lastNonZero).
    std::vector<double> stdevTrimmed(stdevCurve.begin() + firstNonZero, stdevCurve.begin() + lastNonZero);

    // Do not smooth or count zero points: set zeros to NaN before smoothing.
    for (double &v : stdevTrimmed) {
        if (isZeroLike(v)) {
            v = std::nan("");
        }
    }

    // Smooth with median window (NaNs ignored), within the trimmed window only.
    std::vector<double> smoothedRaw = smoothMedian(stdevTrimmed, smoothWindow_);

    // For QC output, embed a filled smoothed curve back into full-length vector.
    std::vector<double> smoothedOutFull(stdevCurve.size(), 0.0);
    std::vector<double> smoothedFilled = smoothedRaw;
    interpolateMissing(smoothedFilled);
    for (size_t i = 0; i < smoothedFilled.size(); ++i) {
        smoothedOutFull[firstNonZero + i] = smoothedFilled[i];
    }
    result.smoothedCurve = smoothedOutFull;

    // For fitting, use the trimmed smoothed curve with NaNs intact.
    std::vector<double> fitCurve = smoothedRaw;

    // Split interval in two halves and fit each half with 1- vs 2-segment model (BIC).
    auto fitHalf = [&](size_t start, size_t end) {
        struct HalfFit {
            uint32_t breakpoint = 0;
            SegmentFit seg1;
            SegmentFit seg2;
            uint32_t segments = 1;
            double sse = 0.0;
            double nValid = 0.0;
        };
        HalfFit out;
        if (end < start) {
            out.sse = std::numeric_limits<double>::max();
            return out;
        }

        uint32_t n = static_cast<uint32_t>(end - start + 1);
        std::vector<double> prefixN(n + 1, 0.0);
        std::vector<double> prefixX(n + 1, 0.0);
        std::vector<double> prefixXX(n + 1, 0.0);
        std::vector<double> prefixY(n + 1, 0.0);
        std::vector<double> prefixXY(n + 1, 0.0);
        std::vector<double> prefixYY(n + 1, 0.0);

        for (uint32_t i = 0; i < n; ++i) {
            double x = static_cast<double>(i);
            double yi = fitCurve[start + i];
            double valid = std::isnan(yi) ? 0.0 : 1.0;
            double yv = std::isnan(yi) ? 0.0 : yi;
            prefixN[i + 1] = prefixN[i] + valid;
            prefixX[i + 1] = prefixX[i] + (valid * x);
            prefixXX[i + 1] = prefixXX[i] + (valid * x * x);
            prefixY[i + 1] = prefixY[i] + yv;
            prefixXY[i + 1] = prefixXY[i] + (valid * x * yv);
            prefixYY[i + 1] = prefixYY[i] + (valid * yv * yv);
        }

        out.nValid = prefixN[n];
        if (out.nValid < 2.0) {
            out.sse = std::numeric_limits<double>::max();
            return out;
        }

        auto bic = [](double rss, uint32_t nPoints, uint32_t kParams) {
            const double eps = 1e-12;
            double rssAdj = (rss <= 0.0) ? eps : rss;
            return static_cast<double>(nPoints) * std::log(rssAdj / static_cast<double>(nPoints)) +
                   static_cast<double>(kParams) * std::log(static_cast<double>(nPoints));
        };

        // 1-segment fit
        SegmentFit seg1 = fitSegment(prefixN, prefixX, prefixXX, prefixY, prefixXY, prefixYY, fitCurve, 0, n - 1);
        double sse1 = seg1.sse;

        // 2-segment fit
        double bestSSE2 = std::numeric_limits<double>::max();
        uint32_t bestB = 0;
        SegmentFit bestSeg1, bestSeg2;
        if (n >= static_cast<uint32_t>(minSegLen_ * 2 + 1)) {
            uint32_t bMin = minSegLen_;
            uint32_t bMax = n - minSegLen_ - 1;
            for (uint32_t b = bMin; b <= bMax; ++b) {
                SegmentFit s1 = fitSegment(prefixN, prefixX, prefixXX, prefixY, prefixXY, prefixYY, fitCurve, 0, b - 1);
                SegmentFit s2 = fitSegment(prefixN, prefixX, prefixXX, prefixY, prefixXY, prefixYY, fitCurve, b, n - 1);
                double total = s1.sse + s2.sse;
                if (total < bestSSE2) {
                    bestSSE2 = total;
                    bestB = b;
                    bestSeg1 = s1;
                    bestSeg2 = s2;
                }
            }
        }

        double bic1 = bic(sse1, static_cast<uint32_t>(out.nValid), 2);
        double bic2 = bic(bestSSE2, static_cast<uint32_t>(out.nValid), 4);
        if (bestSSE2 < std::numeric_limits<double>::max() && bic2 < bic1) {
            out.segments = 2;
            out.breakpoint = static_cast<uint32_t>(start + bestB);
            out.seg1 = bestSeg1;
            out.seg2 = bestSeg2;
            out.sse = bestSSE2;
        } else {
            out.segments = 1;
            out.breakpoint = static_cast<uint32_t>(end);
            out.seg1 = seg1;
            out.seg2 = {};
            out.sse = sse1;
        }

        return out;
    };

    size_t fitN = fitCurve.size();
    if (fitN < static_cast<size_t>(minSegLen_ * 2 + 1)) {
        result.mode = "auto_fallback";
        result.success = false;
        return result;
    }
    size_t mid = fitN / 2;
    if (mid == 0 || mid >= fitN - 1) {
        result.mode = "auto_fallback";
        result.success = false;
        return result;
    }

    auto leftFit = fitHalf(0, mid - 1);
    auto rightFit = fitHalf(mid, fitN - 1);

    // Map breakpoints back to full read coordinates.
    uint32_t b1 = leftFit.breakpoint;
    uint32_t b2 = static_cast<uint32_t>(mid);
    uint32_t b3 = rightFit.breakpoint;
    double sse = leftFit.sse + rightFit.sse;
    SegmentFit seg1 = leftFit.seg1;
    SegmentFit seg2 = leftFit.seg2;
    SegmentFit seg3 = rightFit.seg1;
    SegmentFit seg4 = rightFit.seg2;
    uint32_t nSegments = (leftFit.segments == 2 ? 2u : 1u) + (rightFit.segments == 2 ? 2u : 1u);
    if (leftFit.nValid < 2.0 || rightFit.nValid < 2.0) {
        result.mode = "auto_fallback";
        result.success = false;
        return result;
    }
    
    // Set trim values
    // b1 = first position of middle segment (trim5p = b1)
    // b2 = last position of middle segment (trim3p = readLen - 1 - b2)
    uint32_t n = static_cast<uint32_t>(stdevCurve.size());
    uint32_t b1Full = static_cast<uint32_t>(firstNonZero + b1);
    uint32_t b2Full = static_cast<uint32_t>(firstNonZero + b2);
    uint32_t b3Full = static_cast<uint32_t>(firstNonZero + b3);
    result.trim5p = std::min(b1Full, maxTrim_);
    result.trim3p = std::min(n - 1 - b3Full, maxTrim_);
    
    // Store breakpoints for backward compatibility and QC
    result.kneeBin5p = b1Full;
    result.kneeBinMid = b2Full;
    result.kneeBin3p = b3Full;
    result.totalSSE = sse;
    result.seg1 = seg1;
    result.seg2 = seg2;
    result.seg3 = seg3;
    result.seg4 = seg4;
    
    result.mode = "auto_segmented_halves_bic2";
    result.success = true;
    
    return result;
}
