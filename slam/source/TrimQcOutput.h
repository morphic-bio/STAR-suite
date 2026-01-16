#ifndef TRIM_QC_OUTPUT_H
#define TRIM_QC_OUTPUT_H

#include <string>
#include "TrimQc.h"
#include "Stats.h"

bool writeTrimQcJson(const TrimQcCollector& qc,
                     const Stats& stats,
                     const std::string& outputPath,
                     const std::string& stageLabel);

bool writeTrimQcHtml(const std::string& jsonPath,
                     const std::string& htmlPath);

#endif
