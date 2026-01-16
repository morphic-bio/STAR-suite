#ifndef SLAM_QC_OUTPUT_H
#define SLAM_QC_OUTPUT_H

#include "SlamVarianceAnalysis.h"
#include "SlamQuant.h"
#include <string>
#include <cstdint>
#include <vector>

// Write QC JSON output with per-position variance stats and segmented regression info
bool writeSlamQcJson(const SlamVarianceAnalyzer& analyzer, 
                     const std::string& outputPath,
                     uint32_t fileIndex,
                     const std::string& trimScope,
                     int trim5p,
                     int trim3p,
                     uint64_t readsAnalyzed,
                     const SlamVarianceTrimResult* trimResult = nullptr,
                     const std::string& trimSource = "",
                     double snpErrEst = 0.0,
                     double snpErrUsed = 0.0,
                     const std::string& snpErrFallbackReason = "");

// Write comprehensive QC JSON with T→C/T→A rates, PHRED, trim, segmented fits
// Extracts data from SlamQuant positionTransitions_ and varianceAnalyzer
bool writeSlamQcComprehensiveJson(const SlamQuant& slamQuant,
                                 const std::string& outputPath,
                                 int trim5p,
                                 int trim3p,
                                 const SlamVarianceTrimResult* trimResult = nullptr,
                                 const std::vector<double>* varianceStddevTcRate = nullptr);

// Write QC HTML report using Plotly CDN (4-panel comprehensive report)
bool writeSlamQcComprehensiveHtml(const std::string& jsonPath,
                                   const std::string& htmlPath);

// Write QC HTML report using Plotly CDN (legacy, variance-focused)
bool writeSlamQcHtml(const std::string& jsonPath,
                     const std::string& htmlPath,
                     uint32_t fileIndex);

#endif // SLAM_QC_OUTPUT_H
