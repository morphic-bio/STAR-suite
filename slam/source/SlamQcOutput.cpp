#include "SlamQcOutput.h"
#include "SlamQuant.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <limits>

static uint32_t maxPosInMap(const std::unordered_map<uint32_t, SlamPositionVarianceStats>& stats) {
    uint32_t m = 0;
    for (const auto& kv : stats) {
        if (kv.first > m) {
            m = kv.first;
        }
    }
    return m;
}

static void writeSmoothedStdevCurve(std::ostream& out, const std::vector<double>& curve) {
    out << "[";
    for (size_t i = 0; i < curve.size(); ++i) {
        if (i > 0) {
            out << ", ";
        }
        if (std::isnan(curve[i])) {
            out << "null";
        } else {
            out << std::fixed << std::setprecision(6) << curve[i];
        }
    }
    out << "]";
}

static void writeSegmentedRegressionLegacy(std::ostream& out, const SlamVarianceTrimResult& tr) {
    out << "  \"segmented_regression\": {\n";
    out << "    \"breakpoint_b1\": " << tr.kneeBin5p << ",\n";
    out << "    \"breakpoint_b2\": " << tr.kneeBinMid << ",\n";
    if (tr.mode == "auto_segmented_bic4" || tr.mode == "auto_segmented_halves_bic2") {
        out << "    \"breakpoint_b3\": " << tr.kneeBin3p << ",\n";
    }
    out << "    \"total_sse\": " << std::fixed << std::setprecision(6) << tr.totalSSE << ",\n";
    out << "    \"mode\": \"" << tr.mode << "\",\n";
    out << "    \"segment1\": {\n";
    out << "      \"slope\": " << std::fixed << std::setprecision(6) << tr.seg1.slope << ",\n";
    out << "      \"intercept\": " << std::fixed << std::setprecision(6) << tr.seg1.intercept << ",\n";
    out << "      \"sse\": " << std::fixed << std::setprecision(6) << tr.seg1.sse << "\n";
    out << "    },\n";
    out << "    \"segment2\": {\n";
    out << "      \"slope\": " << std::fixed << std::setprecision(6) << tr.seg2.slope << ",\n";
    out << "      \"intercept\": " << std::fixed << std::setprecision(6) << tr.seg2.intercept << ",\n";
    out << "      \"sse\": " << std::fixed << std::setprecision(6) << tr.seg2.sse << "\n";
    out << "    },\n";
    out << "    \"segment3\": {\n";
    out << "      \"slope\": " << std::fixed << std::setprecision(6) << tr.seg3.slope << ",\n";
    out << "      \"intercept\": " << std::fixed << std::setprecision(6) << tr.seg3.intercept << ",\n";
    out << "      \"sse\": " << std::fixed << std::setprecision(6) << tr.seg3.sse << "\n";
    out << "    },\n";
    if (tr.mode == "auto_segmented_bic4" || tr.mode == "auto_segmented_halves_bic2") {
        out << "    \"segment4\": {\n";
        out << "      \"slope\": " << std::fixed << std::setprecision(6) << tr.seg4.slope << ",\n";
        out << "      \"intercept\": " << std::fixed << std::setprecision(6) << tr.seg4.intercept << ",\n";
        out << "      \"sse\": " << std::fixed << std::setprecision(6) << tr.seg4.sse << "\n";
        out << "    },\n";
    }
    out << "    \"smoothed_stdev_curve\": ";
    if (!tr.smoothedCurve.empty()) {
        writeSmoothedStdevCurve(out, tr.smoothedCurve);
        out << "\n";
    } else {
        out << "[]\n";
    }
    out << "  }";
}

static void writeMateSliceRegression(std::ostream& out, const SlamVarianceMateSlice& m, const std::string& modeStr) {
    out << "    {\n";
    out << "      \"success\": " << (m.success ? "true" : "false") << ",\n";
    out << "      \"mode\": \"" << modeStr << "\",\n";
    out << "      \"trim5p\": " << m.trim5p << ",\n";
    out << "      \"trim3p\": " << m.trim3p << ",\n";
    out << "      \"breakpoint_b1\": " << m.kneeBin5p << ",\n";
    out << "      \"breakpoint_b2\": " << m.kneeBinMid << ",\n";
    out << "      \"breakpoint_b3\": " << m.kneeBin3p << ",\n";
    out << "      \"total_sse\": " << std::fixed << std::setprecision(6) << m.totalSSE << ",\n";
    out << "      \"segment1\": {\n";
    out << "        \"slope\": " << std::fixed << std::setprecision(6) << m.seg1.slope << ",\n";
    out << "        \"intercept\": " << std::fixed << std::setprecision(6) << m.seg1.intercept << ",\n";
    out << "        \"sse\": " << std::fixed << std::setprecision(6) << m.seg1.sse << "\n";
    out << "      },\n";
    out << "      \"segment2\": {\n";
    out << "        \"slope\": " << std::fixed << std::setprecision(6) << m.seg2.slope << ",\n";
    out << "        \"intercept\": " << std::fixed << std::setprecision(6) << m.seg2.intercept << ",\n";
    out << "        \"sse\": " << std::fixed << std::setprecision(6) << m.seg2.sse << "\n";
    out << "      },\n";
    out << "      \"segment3\": {\n";
    out << "        \"slope\": " << std::fixed << std::setprecision(6) << m.seg3.slope << ",\n";
    out << "        \"intercept\": " << std::fixed << std::setprecision(6) << m.seg3.intercept << ",\n";
    out << "        \"sse\": " << std::fixed << std::setprecision(6) << m.seg3.sse << "\n";
    out << "      },\n";
    out << "      \"segment4\": {\n";
    out << "        \"slope\": " << std::fixed << std::setprecision(6) << m.seg4.slope << ",\n";
    out << "        \"intercept\": " << std::fixed << std::setprecision(6) << m.seg4.intercept << ",\n";
    out << "        \"sse\": " << std::fixed << std::setprecision(6) << m.seg4.sse << "\n";
    out << "      },\n";
    out << "      \"smoothed_stdev_curve\": ";
    if (!m.smoothedCurve.empty()) {
        writeSmoothedStdevCurve(out, m.smoothedCurve);
    } else {
        out << "[]";
    }
    out << "\n    }";
}

static void writePositionsFromMap(std::ostream& out, const std::unordered_map<uint32_t, SlamPositionVarianceStats>& stats,
                                  const char* indentTuple) {
    bool first = true;
    uint32_t maxPos = maxPosInMap(stats);
    for (uint32_t pos = 0; pos <= maxPos; ++pos) {
        auto it = stats.find(pos);
        if (it == stats.end()) {
            continue;
        }
        if (!first) {
            out << ",\n";
        }
        first = false;
        const auto& s = it->second;
        double meanQual = s.meanQual();
        double stddevQual = s.stddevQual();
        double meanTcRate = s.meanTcRate();
        double stddevTcRate = s.stddevTcRate();
        out << indentTuple << "{\n";
        out << indentTuple << "  \"position\": " << pos << ",\n";
        out << indentTuple << "  \"read_count\": " << s.readCount << ",\n";
        out << indentTuple << "  \"t_count\": " << s.tCount << ",\n";
        out << indentTuple << "  \"tc_count\": " << s.tcCount << ",\n";
        out << indentTuple << "  \"mean_qual\": " << std::fixed << std::setprecision(2) << meanQual << ",\n";
        out << indentTuple << "  \"stddev_qual\": " << std::fixed << std::setprecision(2) << stddevQual << ",\n";
        out << indentTuple << "  \"variance_qual\": " << std::fixed << std::setprecision(4) << s.varianceQual()
            << ",\n";
        out << indentTuple << "  \"mean_tc_rate\": " << std::fixed << std::setprecision(6) << meanTcRate << ",\n";
        out << indentTuple << "  \"stddev_tc_rate\": " << std::fixed << std::setprecision(6) << stddevTcRate << ",\n";
        out << indentTuple << "  \"variance_tc_rate\": " << std::fixed << std::setprecision(6) << s.varianceTcRate()
            << "\n";
        out << indentTuple << "}";
    }
}

// Mate-local positions synthesized from cached per-index stddev(T→C rate) curves (detection pass),
// used when merged SlamQuant has no per-mate histogram maps but STAR cached curves for --slamQcReport.
static void writePositionsFromStddevCurve(std::ostream& out, const std::vector<double>* curve,
                                          const char* indentTuple) {
    if (curve == nullptr || curve->empty()) {
        return;
    }
    bool first = true;
    for (size_t i = 0; i < curve->size(); ++i) {
        double v = (*curve)[i];
        if (!first) {
            out << ",\n";
        }
        first = false;
        out << indentTuple << "{\n";
        out << indentTuple << "  \"position\": " << i << ",\n";
        out << indentTuple << "  \"read_count\": 0,\n";
        out << indentTuple << "  \"t_count\": 0,\n";
        out << indentTuple << "  \"tc_count\": 0,\n";
        out << indentTuple << "  \"mean_qual\": null,\n";
        out << indentTuple << "  \"stddev_qual\": null,\n";
        out << indentTuple << "  \"variance_qual\": null,\n";
        out << indentTuple << "  \"mean_tc_rate\": null,\n";
        if (std::isnan(v)) {
            out << indentTuple << "  \"stddev_tc_rate\": null,\n";
        } else {
            out << indentTuple << "  \"stddev_tc_rate\": " << std::fixed << std::setprecision(6) << v << ",\n";
        }
        out << indentTuple << "  \"variance_tc_rate\": null\n";
        out << indentTuple << "}";
    }
}

static bool readTextFile(const std::string& path, std::string* out, std::string* err) {
    if (out == nullptr) {
        if (err) *err = "Null output buffer";
        return false;
    }
    std::ifstream in(path.c_str(), std::ios::in | std::ios::binary);
    if (!in.good()) {
        if (err) *err = "Could not open file: " + path;
        return false;
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    *out = ss.str();
    return true;
}

bool writeSlamQcJson(const SlamVarianceAnalyzer& analyzer,
                     const std::string& outputPath,
                     uint32_t fileIndex,
                     const std::string& trimScope,
                     int trim5p_mate1,
                     int trim3p_mate1,
                     int trim5p_mate2,
                     int trim3p_mate2,
                     uint64_t readsAnalyzed,
                     const SlamVarianceTrimResult* trimResult,
                     const std::string& trimSource,
                     double snpErrEst,
                     double snpErrUsed,
                     const std::string& snpErrFallbackReason) {
    std::ofstream out(outputPath.c_str());
    if (!out.good()) {
        return false;
    }

    const bool perMate = analyzer.separateMateHistograms();

    out << "{\n";
    out << "  \"version\": \"1.3\",\n";
    out << "  \"algorithm\": \"segmented_regression\",\n";
    out << "  \"file_index\": " << fileIndex << ",\n";
    out << "  \"variance_histogram_mode\": \"" << (perMate ? "per_mate_separate" : "combined") << "\",\n";
    out << "  \"trim_scope\": \"" << trimScope << "\",\n";
    if (!trimSource.empty()) {
        out << "  \"trim_source\": \"" << trimSource << "\",\n";
    }
    out << "  \"trim5p\": " << trim5p_mate1 << ",\n";
    out << "  \"trim3p\": " << trim3p_mate1 << ",\n";
    out << "  \"trim5p_mate1\": " << trim5p_mate1 << ",\n";
    out << "  \"trim3p_mate1\": " << trim3p_mate1 << ",\n";
    out << "  \"trim5p_mate2\": " << trim5p_mate2 << ",\n";
    out << "  \"trim3p_mate2\": " << trim3p_mate2 << ",\n";
    out << "  \"reads_analyzed\": " << readsAnalyzed << ",\n";

    if (snpErrEst > 0.0 || snpErrUsed > 0.0) {
        out << "  \"snp_err_est\": " << std::fixed << std::setprecision(6) << snpErrEst << ",\n";
        out << "  \"snp_err_used\": " << std::fixed << std::setprecision(6) << snpErrUsed << ",\n";
        if (!snpErrFallbackReason.empty()) {
            out << "  \"snp_err_fallback_reason\": \"" << snpErrFallbackReason << "\",\n";
        }
    }

    if (trimResult != nullptr) {
        writeSegmentedRegressionLegacy(out, *trimResult);
        out << ",\n";
    }

    if (perMate) {
        out << "  \"mates\": [\n";
        for (size_t mi = 0; mi < 2; ++mi) {
            const auto& st = analyzer.getStats(mi);
            int t5 = (mi == 0) ? trim5p_mate1 : trim5p_mate2;
            int t3 = (mi == 0) ? trim3p_mate1 : trim3p_mate2;
            out << "    {\n";
            out << "      \"mate_index\": " << mi << ",\n";
            out << "      \"trim5p\": " << t5 << ",\n";
            out << "      \"trim3p\": " << t3 << ",\n";
            out << "      \"positions\": [\n";
            writePositionsFromMap(out, st, "        ");
            out << "\n      ]";
            if (trimResult != nullptr) {
                out << ",\n";
                out << "      \"segmented_regression\": ";
                writeMateSliceRegression(out, trimResult->mates[mi], trimResult->mode);
            }
            out << "\n    }";
            if (mi < 1) {
                out << ",";
            }
            out << "\n";
        }
        out << "  ],\n";

        out << "  \"positions\": []\n";
    } else {
        const auto& stats = analyzer.getStats();
        uint32_t maxPos = maxPosInMap(stats);

        out << "  \"positions\": [\n";

        bool first = true;
        for (uint32_t pos = 0; pos <= maxPos; ++pos) {
            auto it = stats.find(pos);
            if (it == stats.end()) {
                continue;
            }

            if (!first) {
                out << ",\n";
            }
            first = false;

            const auto& s = it->second;
            double meanQual = s.meanQual();
            double stddevQual = s.stddevQual();
            double meanTcRate = s.meanTcRate();
            double stddevTcRate = s.stddevTcRate();

            out << "    {\n";
            out << "      \"position\": " << pos << ",\n";
            out << "      \"read_count\": " << s.readCount << ",\n";
            out << "      \"t_count\": " << s.tCount << ",\n";
            out << "      \"tc_count\": " << s.tcCount << ",\n";
            out << "      \"mean_qual\": " << std::fixed << std::setprecision(2) << meanQual << ",\n";
            out << "      \"stddev_qual\": " << std::fixed << std::setprecision(2) << stddevQual << ",\n";
            out << "      \"variance_qual\": " << std::fixed << std::setprecision(4) << s.varianceQual() << ",\n";
            out << "      \"mean_tc_rate\": " << std::fixed << std::setprecision(6) << meanTcRate << ",\n";
            out << "      \"stddev_tc_rate\": " << std::fixed << std::setprecision(6) << stddevTcRate << ",\n";
            out << "      \"variance_tc_rate\": " << std::fixed << std::setprecision(6) << s.varianceTcRate() << "\n";
            out << "    }";
        }

        out << "\n  ]\n";
    }

    out << "}\n";

    return out.good();
}

bool writeSlamQcHtml(const std::string& jsonPath,
                     const std::string& htmlPath,
                     uint32_t fileIndex) {
    std::ofstream out(htmlPath.c_str());
    if (!out.good()) {
        return false;
    }

    std::string jsonText;
    std::string readErr;
    if (!readTextFile(jsonPath, &jsonText, &readErr)) {
        return false;
    }
    
    out << "<!DOCTYPE html>\n";
    out << "<html>\n";
    out << "<head>\n";
    out << "  <title>SLAM QC Report - File " << fileIndex << "</title>\n";
    out << "  <script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>\n";
    out << "  <style>\n";
    out << "    body { font-family: 'Segoe UI', Arial, sans-serif; margin: 20px; background: #f5f5f5; }\n";
    out << "    .container { max-width: 1200px; margin: 0 auto; }\n";
    out << "    h1 { color: #333; }\n";
    out << "    .summary { background: white; padding: 15px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }\n";
    out << "    .summary-item { display: inline-block; margin-right: 30px; }\n";
    out << "    .summary-label { color: #666; font-size: 12px; }\n";
    out << "    .summary-value { font-size: 24px; font-weight: bold; color: #2196F3; }\n";
    out << "    .plot-container { background: white; padding: 15px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }\n";
    out << "  </style>\n";
    out << "</head>\n";
    out << "<body>\n";
    out << "  <div class=\"container\">\n";
    out << "    <h1>SLAM Auto-Trim QC Report</h1>\n";
    out << "    <div class=\"summary\" id=\"summary\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"tcStddevPlot\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"qualityPlot\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"tcStddevPlot_m0\" style=\"display:none\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"qualityPlot_m0\" style=\"display:none\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"tcStddevPlot_m1\" style=\"display:none\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"qualityPlot_m1\" style=\"display:none\"></div>\n";
    out << "  </div>\n";
    // Embed JSON so the HTML works when opened directly via file:// (no fetch/CORS issues).
    out << "  <script type=\"application/json\" id=\"qc_json\">\n";
    out << jsonText;
    if (jsonText.empty() || jsonText.back() != '\n') {
        out << "\n";
    }
    out << "  </script>\n";
    out << "  <script>\n";
    out << "    function qcFail(msg) {\n";
    out << "      const el = document.getElementById('summary');\n";
    out << "      if (el) el.innerHTML = '<div style=\"color:#b00020;font-weight:bold\">' + msg + '</div>';\n";
    out << "    }\n";
    out << "    try {\n";
    out << "      if (typeof Plotly === 'undefined') { throw new Error('Plotly failed to load. If you are offline or the CDN is blocked, the plots cannot render.'); }\n";
    out << "      const data = JSON.parse(document.getElementById('qc_json').textContent);\n";
    out << "      const mode = data.variance_histogram_mode || 'combined';\n";
    out << "\n";
    out << "      const mateSummaryVals = (d) => {\n";
    out << "        const t51 = (d.trim5p_mate1 !== undefined && d.trim5p_mate1 !== null) ? d.trim5p_mate1 : d.trim5p;\n";
    out << "        const t31 = (d.trim3p_mate1 !== undefined && d.trim3p_mate1 !== null) ? d.trim3p_mate1 : d.trim3p;\n";
    out << "        const t52 = d.trim5p_mate2;\n";
    out << "        const t32 = d.trim3p_mate2;\n";
    out << "        return { t51, t31, t52, t32 };\n";
    out << "      };\n";
    out << "\n";
    out << "      const plotOneMate = (mate, tcDiv, qualDiv, titlePrefix) => {\n";
    out << "        const seg = (mate && mate.segmented_regression) ? mate.segmented_regression : {};\n";
    out << "        const rows = mate && mate.positions ? mate.positions : [];\n";
    out << "        const positions = rows.map(p => p.position);\n";
    out << "        const stddevTcRate = rows.map(p => p.stddev_tc_rate);\n";
    out << "        const meanQual = rows.map(p => p.mean_qual);\n";
    out << "        const stddevQual = rows.map(p => p.stddev_qual);\n";
    out << "\n";
    out << "        const traces1 = [{\n";
    out << "          x: positions,\n";
    out << "          y: stddevTcRate,\n";
    out << "          type: 'scatter',\n";
    out << "          mode: 'lines+markers',\n";
    out << "          name: 'T→C Stdev (raw)',\n";
    out << "          marker: { size: 4, color: '#2196F3' },\n";
    out << "          line: { width: 1 }\n";
    out << "        }];\n";
    out << "\n";
    out << "        if (seg && seg.smoothed_stdev_curve) {\n";
    out << "          const smoothed = seg.smoothed_stdev_curve;\n";
    out << "          traces1.push({\n";
    out << "            x: positions.slice(0, smoothed.length),\n";
    out << "            y: smoothed,\n";
    out << "            type: 'scatter',\n";
    out << "            mode: 'lines',\n";
    out << "            name: 'Smoothed',\n";
    out << "            line: { width: 2, color: '#333', dash: 'dash' }\n";
    out << "          });\n";
    out << "\n";
    out << "          const b1 = seg.breakpoint_b1;\n";
    out << "          const b2 = seg.breakpoint_b2;\n";
    out << "          const b3 = seg.breakpoint_b3;\n";
    out << "          const seg1 = seg.segment1;\n";
    out << "          const seg2 = seg.segment2;\n";
    out << "          const seg3 = seg.segment3;\n";
    out << "          const seg4 = seg.segment4;\n";
    out << "\n";
    out << "          if (b1 > 0) {\n";
    out << "            const x1 = Array.from({length: b1}, (_, i) => i);\n";
    out << "            const y1 = x1.map(x => seg1.slope * x + seg1.intercept);\n";
    out << "            traces1.push({ x: x1, y: y1, type: 'scatter', mode: 'lines', name: 'Seg1 (5\\' artifact)', line: { width: 2, color: '#f44336' } });\n";
    out << "          }\n";
    out << "\n";
    out << "          const x2 = Array.from({length: b2 - b1 + 1}, (_, i) => b1 + i);\n";
    out << "          const y2 = x2.map(x => seg2.slope * x + seg2.intercept);\n";
    out << "          traces1.push({ x: x2, y: y2, type: 'scatter', mode: 'lines', name: 'Seg2 (signal)', line: { width: 2, color: '#4CAF50' } });\n";
    out << "\n";
    out << "          if (b3 !== undefined && b3 !== null && b3 > 0) {\n";
    out << "            if (b3 > b2 + 1) {\n";
    out << "              const x3 = Array.from({length: b3 - b2 - 1}, (_, i) => b2 + 1 + i);\n";
    out << "              const y3 = x3.map(x => seg3.slope * (x - b2) + seg3.intercept);\n";
    out << "              traces1.push({ x: x3, y: y3, type: 'scatter', mode: 'lines', name: 'Seg3 (mid)', line: { width: 2, color: '#ff9800' } });\n";
    out << "            }\n";
    out << "            if (seg4 && b3 < smoothed.length - 1) {\n";
    out << "              const x4 = Array.from({length: smoothed.length - b3 - 1}, (_, i) => b3 + 1 + i);\n";
    out << "              const y4 = x4.map(x => seg4.slope * (x - b2) + seg4.intercept);\n";
    out << "              traces1.push({ x: x4, y: y4, type: 'scatter', mode: 'lines', name: 'Seg4 (3\\' artifact)', line: { width: 2, color: '#9c27b0' } });\n";
    out << "            }\n";
    out << "          } else if (b2 < smoothed.length - 1) {\n";
    out << "            const x3 = Array.from({length: smoothed.length - b2 - 1}, (_, i) => b2 + 1 + i);\n";
    out << "            const y3 = x3.map(x => seg3.slope * (x - b2) + seg3.intercept);\n";
    out << "            traces1.push({ x: x3, y: y3, type: 'scatter', mode: 'lines', name: 'Seg3 (3\\' artifact)', line: { width: 2, color: '#ff9800' } });\n";
    out << "          }\n";
    out << "        }\n";
    out << "\n";
    out << "        const t5 = mate.trim5p;\n";
    out << "        const t3 = mate.trim3p;\n";
    out << "        const shapes = [];\n";
    out << "        if (t5 > 0) {\n";
    out << "          shapes.push({ type: 'line', x0: t5, x1: t5, y0: 0, y1: 1, yref: 'paper', line: { color: 'purple', width: 2, dash: 'dot' } });\n";
    out << "        }\n";
    out << "        if (t3 > 0 && positions.length > 0) {\n";
    out << "          const trim3pPos = positions[positions.length - 1] - t3;\n";
    out << "          shapes.push({ type: 'line', x0: trim3pPos, x1: trim3pPos, y0: 0, y1: 1, yref: 'paper', line: { color: 'purple', width: 2, dash: 'dot' } });\n";
    out << "        }\n";
    out << "\n";
    out << "        const ann = [];\n";
    out << "        if (positions.length > 0) {\n";
    out << "          ann.push({ x: t5, y: 1, yref: 'paper', text: 'trim5p=' + t5, showarrow: false, yanchor: 'bottom' });\n";
    out << "          ann.push({ x: positions[positions.length - 1] - t3, y: 1, yref: 'paper', text: 'trim3p=' + t3, showarrow: false, yanchor: 'bottom' });\n";
    out << "        }\n";
    out << "\n";
    out << "        Plotly.newPlot(tcDiv, traces1, {\n";
    out << "          title: titlePrefix + ' — T→C stdev (segmented regression)',\n";
    out << "          xaxis: { title: 'Read position (mate-local)' },\n";
    out << "          yaxis: { title: 'T→C Stdev' },\n";
    out << "          shapes: shapes,\n";
    out << "          hovermode: 'closest',\n";
    out << "          annotations: ann\n";
    out << "        });\n";
    out << "\n";
    out << "        Plotly.newPlot(qualDiv, [\n";
    out << "          { x: positions, y: meanQual, type: 'scatter', mode: 'lines+markers', name: 'Mean Quality', marker: { size: 4 } },\n";
    out << "          { x: positions, y: stddevQual, type: 'scatter', mode: 'lines+markers', name: 'Quality Stdev', yaxis: 'y2', marker: { size: 4 } }\n";
    out << "        ], {\n";
    out << "          title: titlePrefix + ' — Quality',\n";
    out << "          xaxis: { title: 'Read position (mate-local)' },\n";
    out << "          yaxis: { title: 'Mean Quality', side: 'left' },\n";
    out << "          yaxis2: { title: 'Quality Stdev', side: 'right', overlaying: 'y' },\n";
    out << "          shapes: shapes,\n";
    out << "          hovermode: 'closest'\n";
    out << "        });\n";
    out << "      };\n";
    out << "\n";
    out << "      if (mode === 'per_mate_separate' && data.mates && data.mates.length >= 2) {\n";
    out << "        const mv = mateSummaryVals(data);\n";
    out << "        const summaryHtml = `\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Histogram mode</div><div class=\"summary-value\">${data.variance_histogram_mode}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Algorithm</div><div class=\"summary-value\">${data.algorithm || 'segmented_regression'}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Mate1 trim 5'/3'</div><div class=\"summary-value\">${mv.t51} / ${mv.t31}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Mate2 trim 5'/3'</div><div class=\"summary-value\">${mv.t52} / ${mv.t32}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Reads analyzed</div><div class=\"summary-value\">${data.reads_analyzed.toLocaleString()}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">File index</div><div class=\"summary-value\">${data.file_index}</div></div>\n";
    out << "        `;\n";
    out << "        document.getElementById('summary').innerHTML = summaryHtml;\n";
    out << "        document.getElementById('tcStddevPlot').style.display = 'none';\n";
    out << "        document.getElementById('qualityPlot').style.display = 'none';\n";
    out << "        document.getElementById('tcStddevPlot_m0').style.display = 'block';\n";
    out << "        document.getElementById('qualityPlot_m0').style.display = 'block';\n";
    out << "        document.getElementById('tcStddevPlot_m1').style.display = 'block';\n";
    out << "        document.getElementById('qualityPlot_m1').style.display = 'block';\n";
    out << "        plotOneMate(data.mates[0], 'tcStddevPlot_m0', 'qualityPlot_m0', 'Mate 1');\n";
    out << "        plotOneMate(data.mates[1], 'tcStddevPlot_m1', 'qualityPlot_m1', 'Mate 2');\n";
    out << "      } else {\n";
    out << "        const mv = mateSummaryVals(data);\n";
    out << "        const summaryHtml = `\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Algorithm</div><div class=\"summary-value\">${data.algorithm || 'segmented_regression'}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Trim 5' (mate1 legacy)</div><div class=\"summary-value\">${mv.t51}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Trim 3' (mate1 legacy)</div><div class=\"summary-value\">${mv.t31}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Mate2 trim 5'/3'</div><div class=\"summary-value\">${mv.t52} / ${mv.t32}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Reads analyzed</div><div class=\"summary-value\">${data.reads_analyzed.toLocaleString()}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">File index</div><div class=\"summary-value\">${data.file_index}</div></div>\n";
    out << "        `;\n";
    out << "        document.getElementById('summary').innerHTML = summaryHtml;\n";
    out << "\n";
    out << "        const positions = (data.positions || []).map(p => p.position);\n";
    out << "        const stddevTcRate = (data.positions || []).map(p => p.stddev_tc_rate);\n";
    out << "        const meanQual = (data.positions || []).map(p => p.mean_qual);\n";
    out << "        const stddevQual = (data.positions || []).map(p => p.stddev_qual);\n";
    out << "\n";
    out << "        const traces1 = [{\n";
    out << "          x: positions,\n";
    out << "          y: stddevTcRate,\n";
    out << "          type: 'scatter',\n";
    out << "          mode: 'lines+markers',\n";
    out << "          name: 'T→C Stdev (raw)',\n";
    out << "          marker: { size: 4, color: '#2196F3' },\n";
    out << "          line: { width: 1 }\n";
    out << "        }];\n";
    out << "\n";
    out << "        if (data.segmented_regression && data.segmented_regression.smoothed_stdev_curve) {\n";
    out << "          const smoothed = data.segmented_regression.smoothed_stdev_curve;\n";
    out << "          traces1.push({\n";
    out << "            x: positions.slice(0, smoothed.length),\n";
    out << "            y: smoothed,\n";
    out << "            type: 'scatter',\n";
    out << "            mode: 'lines',\n";
    out << "            name: 'Smoothed',\n";
    out << "            line: { width: 2, color: '#333', dash: 'dash' }\n";
    out << "          });\n";
    out << "\n";
    out << "          const b1 = data.segmented_regression.breakpoint_b1;\n";
    out << "          const b2 = data.segmented_regression.breakpoint_b2;\n";
    out << "          const b3 = data.segmented_regression.breakpoint_b3;\n";
    out << "          const seg1 = data.segmented_regression.segment1;\n";
    out << "          const seg2 = data.segmented_regression.segment2;\n";
    out << "          const seg3 = data.segmented_regression.segment3;\n";
    out << "          const seg4 = data.segmented_regression.segment4;\n";
    out << "\n";
    out << "          if (b1 > 0) {\n";
    out << "            const x1 = Array.from({length: b1}, (_, i) => i);\n";
    out << "            const y1 = x1.map(x => seg1.slope * x + seg1.intercept);\n";
    out << "            traces1.push({ x: x1, y: y1, type: 'scatter', mode: 'lines', name: 'Seg1 (5\\' artifact)', line: { width: 2, color: '#f44336' } });\n";
    out << "          }\n";
    out << "\n";
    out << "          const x2 = Array.from({length: b2 - b1 + 1}, (_, i) => b1 + i);\n";
    out << "          const y2 = x2.map(x => seg2.slope * x + seg2.intercept);\n";
    out << "          traces1.push({ x: x2, y: y2, type: 'scatter', mode: 'lines', name: 'Seg2 (signal)', line: { width: 2, color: '#4CAF50' } });\n";
    out << "\n";
    out << "          if (b3 !== undefined && b3 !== null && b3 > 0) {\n";
    out << "            if (b3 > b2 + 1) {\n";
    out << "              const x3 = Array.from({length: b3 - b2 - 1}, (_, i) => b2 + 1 + i);\n";
    out << "              const y3 = x3.map(x => seg3.slope * (x - b2) + seg3.intercept);\n";
    out << "              traces1.push({ x: x3, y: y3, type: 'scatter', mode: 'lines', name: 'Seg3 (mid)', line: { width: 2, color: '#ff9800' } });\n";
    out << "            }\n";
    out << "            if (seg4 && b3 < smoothed.length - 1) {\n";
    out << "              const x4 = Array.from({length: smoothed.length - b3 - 1}, (_, i) => b3 + 1 + i);\n";
    out << "              const y4 = x4.map(x => seg4.slope * (x - b2) + seg4.intercept);\n";
    out << "              traces1.push({ x: x4, y: y4, type: 'scatter', mode: 'lines', name: 'Seg4 (3\\' artifact)', line: { width: 2, color: '#9c27b0' } });\n";
    out << "            }\n";
    out << "          } else if (b2 < smoothed.length - 1) {\n";
    out << "            const x3 = Array.from({length: smoothed.length - b2 - 1}, (_, i) => b2 + 1 + i);\n";
    out << "            const y3 = x3.map(x => seg3.slope * (x - b2) + seg3.intercept);\n";
    out << "            traces1.push({ x: x3, y: y3, type: 'scatter', mode: 'lines', name: 'Seg3 (3\\' artifact)', line: { width: 2, color: '#ff9800' } });\n";
    out << "          }\n";
    out << "        }\n";
    out << "\n";
    out << "        const shapes = [];\n";
    out << "        if (data.trim5p > 0) {\n";
    out << "          shapes.push({ type: 'line', x0: data.trim5p, x1: data.trim5p, y0: 0, y1: 1, yref: 'paper', line: { color: 'purple', width: 2, dash: 'dot' } });\n";
    out << "        }\n";
    out << "        if (data.trim3p > 0 && positions.length > 0) {\n";
    out << "          const trim3pPos = positions[positions.length - 1] - data.trim3p;\n";
    out << "          shapes.push({ type: 'line', x0: trim3pPos, x1: trim3pPos, y0: 0, y1: 1, yref: 'paper', line: { color: 'purple', width: 2, dash: 'dot' } });\n";
    out << "        }\n";
    out << "\n";
    out << "        Plotly.newPlot('tcStddevPlot', traces1, {\n";
    out << "          title: 'T→C rate standard deviation by position (segmented regression)',\n";
    out << "          xaxis: { title: 'Read position' },\n";
    out << "          yaxis: { title: 'T→C Stdev' },\n";
    out << "          shapes: shapes,\n";
    out << "          hovermode: 'closest',\n";
    out << "          annotations: [\n";
    out << "            { x: data.trim5p, y: 1, yref: 'paper', text: 'trim5p=' + data.trim5p, showarrow: false, yanchor: 'bottom' },\n";
    out << "            { x: positions.length ? positions[positions.length - 1] - data.trim3p : data.trim3p, y: 1, yref: 'paper', text: 'trim3p=' + data.trim3p, showarrow: false, yanchor: 'bottom' }\n";
    out << "          ]\n";
    out << "        });\n";
    out << "\n";
    out << "        Plotly.newPlot('qualityPlot', [\n";
    out << "          { x: positions, y: meanQual, type: 'scatter', mode: 'lines+markers', name: 'Mean Quality', marker: { size: 4 } },\n";
    out << "          { x: positions, y: stddevQual, type: 'scatter', mode: 'lines+markers', name: 'Quality Stdev', yaxis: 'y2', marker: { size: 4 } }\n";
    out << "        ], {\n";
    out << "          title: 'Quality statistics by position',\n";
    out << "          xaxis: { title: 'Read position' },\n";
    out << "          yaxis: { title: 'Mean Quality', side: 'left' },\n";
    out << "          yaxis2: { title: 'Quality Stdev', side: 'right', overlaying: 'y' },\n";
    out << "          shapes: shapes,\n";
    out << "          hovermode: 'closest'\n";
    out << "        });\n";
    out << "      }\n";
    out << "\n";
    out << "    } catch (e) {\n";
    out << "      qcFail('Failed to render plots: ' + (e && e.message ? e.message : String(e)));\n";
    out << "    }\n";
    out << "  </script>\n";
    out << "</body>\n";
    out << "</html>\n";
    
    return out.good();
}

bool writeSlamQcComprehensiveJson(const SlamQuant& slamQuant,
                                   const std::string& outputPath,
                                   int trim5p_mate1,
                                   int trim3p_mate1,
                                   int trim5p_mate2,
                                   int trim3p_mate2,
                                   const SlamVarianceTrimResult* trimResult,
                                   const std::vector<double>* varianceStddevTcRate,
                                   const std::vector<double>* varianceStddevTcRateMate2,
                                   uint32_t mateLen0,
                                   uint32_t mateLen1) {
    std::ofstream out(outputPath.c_str());
    if (!out.good()) {
        return false;
    }

    auto transitionData = slamQuant.getPositionTransitionData();

    const SlamVarianceAnalyzer* analyzer = slamQuant.varianceAnalyzer();
    const bool perMateFromAnalyzer = analyzer != nullptr && analyzer->separateMateHistograms();
    const bool perMateFromCurves = mateLen1 > 0 && varianceStddevTcRateMate2 != nullptr &&
                                   !varianceStddevTcRateMate2->empty();
    const bool perMate = perMateFromAnalyzer || perMateFromCurves;
    std::unordered_map<uint32_t, SlamPositionVarianceStats> var0;
    std::unordered_map<uint32_t, SlamPositionVarianceStats> var1;
    if (analyzer != nullptr) {
        if (perMate) {
            var0 = analyzer->getStats(0);
            var1 = analyzer->getStats(1);
        } else {
            var0 = analyzer->getStats();
        }
    }

    uint32_t maxPos = 0;
    for (const auto& kv : transitionData) {
        if (kv.first > maxPos) {
            maxPos = kv.first;
        }
    }
    for (const auto& kv : var0) {
        if (kv.first > maxPos) {
            maxPos = kv.first;
        }
    }
    if (perMate && mateLen1 > 0) {
        for (const auto& kv : var1) {
            const uint32_t concatPos = mateLen0 + kv.first;
            if (concatPos > maxPos) {
                maxPos = concatPos;
            }
        }
    }
    if (perMateFromCurves) {
        if (varianceStddevTcRate != nullptr && !varianceStddevTcRate->empty()) {
            const uint32_t last0 = static_cast<uint32_t>(varianceStddevTcRate->size() - 1);
            if (last0 > maxPos) {
                maxPos = last0;
            }
        }
        if (varianceStddevTcRateMate2 != nullptr && !varianceStddevTcRateMate2->empty() && mateLen0 > 0 &&
            mateLen1 > 0) {
            const uint32_t lastConcat =
                mateLen0 + static_cast<uint32_t>(varianceStddevTcRateMate2->size() - 1);
            if (lastConcat > maxPos) {
                maxPos = lastConcat;
            }
        }
    }

    out << "{\n";
    out << "  \"version\": \"2.1\",\n";
    out << "  \"type\": \"comprehensive_qc\",\n";
    out << "  \"variance_histogram_mode\": \"" << (perMate ? "per_mate_separate" : "combined") << "\",\n";
    out << "  \"mate_len_mate1\": " << mateLen0 << ",\n";
    out << "  \"mate_len_mate2\": " << mateLen1 << ",\n";
    out << "  \"trim5p\": " << trim5p_mate1 << ",\n";
    out << "  \"trim3p\": " << trim3p_mate1 << ",\n";
    out << "  \"trim5p_mate1\": " << trim5p_mate1 << ",\n";
    out << "  \"trim3p_mate1\": " << trim3p_mate1 << ",\n";
    out << "  \"trim5p_mate2\": " << trim5p_mate2 << ",\n";
    out << "  \"trim3p_mate2\": " << trim3p_mate2 << ",\n";

    if (trimResult != nullptr) {
        writeSegmentedRegressionLegacy(out, *trimResult);
        out << ",\n";
    } else {
        out << "  \"segmented_regression\": null,\n";
    }

    out << "  \"positions\": [\n";

    bool first = true;
    for (uint32_t pos = 0; pos <= maxPos; ++pos) {
        const bool hasTransition = transitionData.count(pos) > 0;

        bool hasVariance = false;
        const SlamPositionVarianceStats* vp = nullptr;
        if (analyzer != nullptr) {
            if (!perMate) {
                if (var0.count(pos) > 0) {
                    hasVariance = true;
                    vp = &var0.at(pos);
                }
            } else {
                if (pos < mateLen0 && var0.count(pos) > 0) {
                    hasVariance = true;
                    vp = &var0.at(pos);
                } else if (mateLen1 > 0 && pos >= mateLen0) {
                    const uint32_t lp = pos - mateLen0;
                    if (var1.count(lp) > 0) {
                        hasVariance = true;
                        vp = &var1.at(lp);
                    }
                }
            }
        }

        double curveV = std::numeric_limits<double>::quiet_NaN();
        bool haveCurveSample = false;
        if (perMate && varianceStddevTcRate != nullptr && varianceStddevTcRateMate2 != nullptr) {
            if (mateLen0 > 0 && pos < mateLen0 && pos < varianceStddevTcRate->size()) {
                curveV = (*varianceStddevTcRate)[pos];
                haveCurveSample = true;
            } else if (mateLen1 > 0 && pos >= mateLen0) {
                const size_t li = static_cast<size_t>(pos - mateLen0);
                if (li < varianceStddevTcRateMate2->size()) {
                    curveV = (*varianceStddevTcRateMate2)[li];
                    haveCurveSample = true;
                }
            }
        } else if (!perMate && varianceStddevTcRate != nullptr && pos < varianceStddevTcRate->size()) {
            curveV = (*varianceStddevTcRate)[pos];
            haveCurveSample = true;
        }
        const bool useCurveNumeric = haveCurveSample && !std::isnan(curveV);

        if (!hasTransition && !hasVariance && !useCurveNumeric) {
            continue;
        }

        if (!first) {
            out << ",\n";
        }
        first = false;

        out << "    {\n";
        out << "      \"position\": " << (pos + 1) << ",\n";

        if (hasTransition) {
            const auto& trans = transitionData.at(pos);
            double tc_cov = std::get<0>(trans);
            double tc_mm = std::get<1>(trans);
            double ta_cov = std::get<2>(trans);
            double ta_mm = std::get<3>(trans);

            double tc_rate = (tc_cov > 0.0) ? (tc_mm / tc_cov) * 100.0 : 0.0;
            double ta_rate = (ta_cov > 0.0) ? (ta_mm / ta_cov) * 100.0 : 0.0;
            double tc_stdev = (tc_cov > 0.0) ? std::sqrt((tc_mm / tc_cov) * (1.0 - tc_mm / tc_cov)) * 100.0 : 0.0;

            out << "      \"star_tc_cov\": " << std::fixed << std::setprecision(2) << tc_cov << ",\n";
            out << "      \"star_tc_mm\": " << std::fixed << std::setprecision(2) << tc_mm << ",\n";
            out << "      \"star_tc_rate\": " << std::fixed << std::setprecision(6) << tc_rate << ",\n";
            out << "      \"star_tc_stdev\": " << std::fixed << std::setprecision(6) << tc_stdev << ",\n";
            out << "      \"star_ta_cov\": " << std::fixed << std::setprecision(2) << ta_cov << ",\n";
            out << "      \"star_ta_mm\": " << std::fixed << std::setprecision(2) << ta_mm << ",\n";
            out << "      \"star_ta_rate\": " << std::fixed << std::setprecision(6) << ta_rate << ",\n";
        } else {
            out << "      \"star_tc_cov\": 0,\n";
            out << "      \"star_tc_mm\": 0,\n";
            out << "      \"star_tc_rate\": 0,\n";
            out << "      \"star_tc_stdev\": 0,\n";
            out << "      \"star_ta_cov\": 0,\n";
            out << "      \"star_ta_mm\": 0,\n";
            out << "      \"star_ta_rate\": 0,\n";
        }

        if (hasVariance && vp != nullptr) {
            double meanQual = vp->meanQual();
            double stddevQual = vp->stddevQual();
            double meanTcRate = vp->meanTcRate();
            double stddevTcRate = vp->stddevTcRate();

            out << "      \"mean_qual\": " << std::fixed << std::setprecision(2) << meanQual << ",\n";
            out << "      \"stddev_qual\": " << std::fixed << std::setprecision(2) << stddevQual << ",\n";
            out << "      \"t_count\": " << vp->tCount << ",\n";
            out << "      \"mean_tc_rate\": " << std::fixed << std::setprecision(6) << meanTcRate << ",\n";
            out << "      \"stddev_tc_rate\": " << std::fixed << std::setprecision(6) << stddevTcRate << "\n";
        } else {
            out << "      \"mean_qual\": null,\n";
            out << "      \"stddev_qual\": null,\n";
            out << "      \"t_count\": 0,\n";
            out << "      \"mean_tc_rate\": null,\n";

            if (useCurveNumeric) {
                out << "      \"stddev_tc_rate\": " << std::fixed << std::setprecision(6) << curveV << "\n";
            } else {
                out << "      \"stddev_tc_rate\": null\n";
            }
        }

        out << "    }";
    }

    out << "\n  ]";

    if (perMate) {
        out << ",\n";
        out << "  \"mates\": [\n";
        for (size_t mi = 0; mi < 2; ++mi) {
            const auto& stMap = (mi == 0) ? var0 : var1;
            int t5 = (mi == 0) ? trim5p_mate1 : trim5p_mate2;
            int t3 = (mi == 0) ? trim3p_mate1 : trim3p_mate2;
            out << "    {\n";
            out << "      \"mate_index\": " << mi << ",\n";
            out << "      \"trim5p\": " << t5 << ",\n";
            out << "      \"trim3p\": " << t3 << ",\n";
            out << "      \"positions\": [\n";
            if (!stMap.empty()) {
                writePositionsFromMap(out, stMap, "        ");
            } else if (perMateFromCurves) {
                const std::vector<double>* crv =
                    (mi == 0) ? varianceStddevTcRate : varianceStddevTcRateMate2;
                writePositionsFromStddevCurve(out, crv, "        ");
            }
            out << "\n      ]";
            if (trimResult != nullptr) {
                out << ",\n";
                out << "      \"segmented_regression\": ";
                writeMateSliceRegression(out, trimResult->mates[mi], trimResult->mode);
            }
            out << "\n    }";
            if (mi < 1) {
                out << ",";
            }
            out << "\n";
        }
        out << "  ]\n";
    } else {
        out << "\n";
    }

    out << "}\n";

    return out.good();
}


bool writeSlamQcComprehensiveHtml(const std::string& jsonPath,
                                   const std::string& htmlPath) {
    std::ofstream out(htmlPath.c_str());
    if (!out.good()) {
        return false;
    }

    std::string jsonText;
    std::string readErr;
    if (!readTextFile(jsonPath, &jsonText, &readErr)) {
        return false;
    }
    
    // Write comprehensive HTML with 4 plots matching analyze_phred_by_position.py
    out << "<!DOCTYPE html>\n";
    out << "<html>\n";
    out << "<head>\n";
    out << "  <title>SLAM QC Report</title>\n";
    out << "  <script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>\n";
    out << "  <style>\n";
    out << "    body { font-family: 'Segoe UI', Arial, sans-serif; margin: 20px; background: #f5f5f5; }\n";
    out << "    .container { max-width: 1400px; margin: 0 auto; }\n";
    out << "    h1 { color: #333; }\n";
    out << "    .summary { background: white; padding: 15px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }\n";
    out << "    .summary-item { display: inline-block; margin-right: 30px; }\n";
    out << "    .summary-label { color: #666; font-size: 12px; }\n";
    out << "    .summary-value { font-size: 24px; font-weight: bold; color: #2196F3; }\n";
    out << "    .plot-container { background: white; padding: 15px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }\n";
    out << "  </style>\n";
    out << "</head>\n";
    out << "<body>\n";
    out << "  <div class=\"container\">\n";
    out << "    <h1>SLAM QC Report</h1>\n";
    out << "    <div class=\"summary\" id=\"summary\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"plot1\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"plot2\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"plot3\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"plot4\"></div>\n";
    out << "  </div>\n";
    // Embed JSON so the HTML works when opened directly via file:// (no fetch/CORS issues).
    out << "  <script type=\"application/json\" id=\"qc_json\">\n";
    out << jsonText;
    if (jsonText.empty() || jsonText.back() != '\n') {
        out << "\n";
    }
    out << "  </script>\n";
    out << "  <script>\n";
    out << "    function qcFail(msg) {\n";
    out << "      const el = document.getElementById('summary');\n";
    out << "      if (el) el.innerHTML = '<div style=\"color:#b00020;font-weight:bold\">' + msg + '</div>';\n";
    out << "    }\n";
    out << "    try {\n";
    out << "      if (typeof Plotly === 'undefined') { throw new Error('Plotly failed to load. If you are offline or the CDN is blocked, the plots cannot render.'); }\n";
    out << "      const data = JSON.parse(document.getElementById('qc_json').textContent);\n";
    out << "        const positions = data.positions.map(p => p.position);\n";
    out << "        const starTcRate = data.positions.map(p => p.star_tc_rate || 0);\n";
    out << "        const starTaRate = data.positions.map(p => p.star_ta_rate || 0);\n";
    out << "        const tcStdev = data.positions.map(p => {\n";
    out << "          if (p.stddev_tc_rate !== null && p.stddev_tc_rate !== undefined) return p.stddev_tc_rate;\n";
    out << "          return null;\n";
    out << "        });\n";
    out << "        const meanQual = data.positions.map(p => p.mean_qual);\n";
    out << "        const tCount = data.positions.map(p => p.t_count || 0);\n";
    out << "        const seg = data.segmented_regression || null;\n";
    out << "        const hasVarianceRaw = tcStdev.some(v => v !== null && !isNaN(v));\n";
    out << "        const hasVarianceSmooth = !!(seg && seg.smoothed_stdev_curve && seg.smoothed_stdev_curve.some(v => v !== null && !isNaN(v)));\n";
    out << "        const hasVariance = hasVarianceRaw || hasVarianceSmooth;\n";
    out << "        \n";
    out << "        // Trim bars\n";
    out << "        const shapes = [];\n";
    out << "        if (data.trim5p > 0) {\n";
    out << "          shapes.push({ type: 'line', x0: data.trim5p + 1, x1: data.trim5p + 1, y0: 0, y1: 1, yref: 'paper', line: { color: 'purple', width: 2, dash: 'dash' } });\n";
    out << "        }\n";
    out << "        if (data.trim3p > 0 && positions.length > 0) {\n";
    out << "          const trim3pPos = positions[positions.length - 1] - data.trim3p;\n";
    out << "          shapes.push({ type: 'line', x0: trim3pPos, x1: trim3pPos, y0: 0, y1: 1, yref: 'paper', line: { color: 'purple', width: 2, dash: 'dot' } });\n";
    out << "        }\n";
    out << "        \n";
    out << "        // Summary\n";
    out << "        const mv = (d) => ({\n";
    out << "          t51: (d.trim5p_mate1 !== undefined && d.trim5p_mate1 !== null) ? d.trim5p_mate1 : d.trim5p,\n";
    out << "          t31: (d.trim3p_mate1 !== undefined && d.trim3p_mate1 !== null) ? d.trim3p_mate1 : d.trim3p,\n";
    out << "          t52: d.trim5p_mate2,\n";
    out << "          t32: d.trim3p_mate2\n";
    out << "        });\n";
    out << "        const mv0 = mv(data);\n";
    out << "        document.getElementById('summary').innerHTML = `\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Mode</div><div class=\"summary-value\">${data.variance_histogram_mode || 'combined'}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Trim M1 5'/3'</div><div class=\"summary-value\">${mv0.t51} / ${mv0.t31}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Trim M2 5'/3'</div><div class=\"summary-value\">${mv0.t52} / ${mv0.t32}</div></div>\n";
    out << "        `;\n";
    out << "        \n";
    out << "        // Plot 1: T transition rates (T→C and T→A)\n";
    out << "        Plotly.newPlot('plot1', [\n";
    out << "          { x: positions, y: starTcRate, type: 'scatter', mode: 'lines', name: 'STAR T→C', line: { width: 2, color: '#2196F3' } },\n";
    out << "          { x: positions, y: starTaRate, type: 'scatter', mode: 'lines', name: 'STAR T→A', line: { width: 1, color: '#2196F3', dash: 'dash' } }\n";
    out << "        ], {\n";
    out << "          title: 'Raw Mismatch Rates by Position',\n";
    out << "          xaxis: { title: 'Read Position (1-based)' },\n";
    out << "          yaxis: { title: 'Mismatch Rate (%)', type: 'log' },\n";
    out << "          shapes: shapes,\n";
    out << "          hovermode: 'closest',\n";
    out << "          legend: { x: 0.7, y: 0.9 }\n";
    out << "        });\n";
    out << "        \n";
    out << "        // Plot 2: T→C stdev + smoothed + segmented fit (only if variance data exists)\n";
    out << "        if (hasVariance) {\n";
    out << "          const traces2 = [\n";
    out << "            { x: positions, y: tcStdev, type: 'scatter', mode: 'lines+markers', name: 'T→C Stdev (variance)', marker: { size: 4, color: '#2196F3' }, line: { width: 1 } }\n";
    out << "          ];\n";
    out << "          \n";
    out << "          if (seg) {\n";
    out << "            if (seg.smoothed_stdev_curve && seg.smoothed_stdev_curve.length > 0) {\n";
    out << "              const smoothed = seg.smoothed_stdev_curve;\n";
    out << "              const smoothedPos = Array.from({length: smoothed.length}, (_, i) => i + 1);\n";
    out << "              traces2.push({ x: smoothedPos, y: smoothed, type: 'scatter', mode: 'lines', name: 'Smoothed', line: { width: 2, color: '#333', dash: 'dash' } });\n";
    out << "              \n";
    out << "              // Segmented regression lines\n";
    out << "              const b1 = seg.breakpoint_b1;\n";
    out << "              const b2 = seg.breakpoint_b2;\n";
    out << "              const b3 = seg.breakpoint_b3;\n";
    out << "              const mode = seg.mode || '';\n";
    out << "              if (b1 > 0) {\n";
    out << "                const x1 = Array.from({length: b1}, (_, i) => i + 1);\n";
    out << "                const y1 = x1.map(x => seg.segment1.slope * (x - 1) + seg.segment1.intercept);\n";
    out << "                traces2.push({ x: x1, y: y1, type: 'scatter', mode: 'lines', name: 'Seg1', line: { width: 2, color: '#f44336' } });\n";
    out << "              }\n";
    out << "              const x2 = Array.from({length: b2 - b1 + 1}, (_, i) => b1 + i + 1);\n";
    out << "              const y2 = x2.map(x => seg.segment2.slope * (x - 1) + seg.segment2.intercept);\n";
    out << "              traces2.push({ x: x2, y: y2, type: 'scatter', mode: 'lines', name: 'Seg2', line: { width: 2, color: '#4CAF50' } });\n";
    out << "              if (b3 !== undefined && b3 !== null && b3 > 0) {\n";
    out << "                if (b3 > b2 + 1) {\n";
    out << "                  const x3 = Array.from({length: b3 - b2 - 1}, (_, i) => b2 + i + 2);\n";
    out << "                  const y3 = x3.map(x => {\n";
    out << "                    const xv = (mode === 'auto_segmented_halves_bic2') ? (x - b2) : (x - 1);\n";
    out << "                    return seg.segment3.slope * xv + seg.segment3.intercept;\n";
    out << "                  });\n";
    out << "                  traces2.push({ x: x3, y: y3, type: 'scatter', mode: 'lines', name: 'Seg3', line: { width: 2, color: '#ff9800' } });\n";
    out << "                }\n";
    out << "                if (seg.segment4 && b3 < smoothed.length - 1) {\n";
    out << "                  const x4 = Array.from({length: smoothed.length - b3 - 1}, (_, i) => b3 + i + 2);\n";
    out << "                  const y4 = x4.map(x => {\n";
    out << "                    const xv = (mode === 'auto_segmented_halves_bic2') ? (x - b2) : (x - 1);\n";
    out << "                    return seg.segment4.slope * xv + seg.segment4.intercept;\n";
    out << "                  });\n";
    out << "                  traces2.push({ x: x4, y: y4, type: 'scatter', mode: 'lines', name: 'Seg4', line: { width: 2, color: '#9c27b0' } });\n";
    out << "                }\n";
    out << "              } else if (b2 < smoothed.length - 1) {\n";
    out << "                const x3 = Array.from({length: smoothed.length - b2 - 1}, (_, i) => b2 + i + 2);\n";
    out << "                const y3 = x3.map(x => {\n";
    out << "                  const xv = (mode === 'auto_segmented_halves_bic2') ? (x - b2) : (x - 1);\n";
    out << "                  return seg.segment3.slope * xv + seg.segment3.intercept;\n";
    out << "                });\n";
    out << "                traces2.push({ x: x3, y: y3, type: 'scatter', mode: 'lines', name: 'Seg3', line: { width: 2, color: '#ff9800' } });\n";
    out << "              }\n";
    out << "            }\n";
    out << "          }\n";
    out << "          \n";
    out << "          Plotly.newPlot('plot2', traces2, {\n";
    out << "            title: 'T→C Standard Deviation by Position (segmented fit)',\n";
    out << "            xaxis: { title: 'Read Position (1-based)' },\n";
    out << "            yaxis: { title: 'T→C Stdev (%)' },\n";
    out << "            shapes: shapes,\n";
    out << "            hovermode: 'closest',\n";
    out << "            legend: { x: 0.7, y: 0.9 }\n";
    out << "          });\n";
    out << "        } else {\n";
    out << "          const plot2 = document.getElementById('plot2');\n";
    out << "          if (plot2) plot2.style.display = 'none';\n";
    out << "        }\n";
    out << "        \n";
    out << "        // Plot 3: PHRED + T assignment rate\n";
    out << "        Plotly.newPlot('plot3', [\n";
    out << "          { x: positions, y: meanQual, type: 'scatter', mode: 'lines+markers', name: 'Avg PHRED (T positions)', yaxis: 'y', marker: { size: 4, color: '#2196F3' }, line: { width: 2 } },\n";
    out << "          { x: positions, y: tCount, type: 'scatter', mode: 'lines+markers', name: 'T Count', yaxis: 'y2', marker: { size: 4, color: '#4CAF50' }, line: { width: 2 } }\n";
    out << "        ], {\n";
    out << "          title: 'PHRED Quality and T Assignment Rate by Position',\n";
    out << "          xaxis: { title: 'Read Position (1-based)' },\n";
    out << "          yaxis: { title: 'Average PHRED Quality', side: 'left', range: [0, 45] },\n";
    out << "          yaxis2: { title: 'T Count', side: 'right', overlaying: 'y' },\n";
    out << "          shapes: shapes,\n";
    out << "          hovermode: 'closest',\n";
    out << "          legend: { x: 0.7, y: 0.9 }\n";
    out << "        });\n";
    out << "        \n";
    out << "        // Plot 4: Trim overlay (duplicate of plot 1 with trim emphasis)\n";
    out << "        Plotly.newPlot('plot4', [\n";
    out << "          { x: positions, y: starTcRate, type: 'scatter', mode: 'lines', name: 'STAR T→C', line: { width: 2, color: '#2196F3' } },\n";
    out << "          { x: positions, y: starTaRate, type: 'scatter', mode: 'lines', name: 'STAR T→A', line: { width: 1, color: '#2196F3', dash: 'dash' } }\n";
    out << "        ], {\n";
    out << "          title: 'Trim Overlay on T Transition Rates',\n";
    out << "          xaxis: { title: 'Read Position (1-based)' },\n";
    out << "          yaxis: { title: 'Mismatch Rate (%)', type: 'log' },\n";
    out << "          shapes: shapes,\n";
    out << "          hovermode: 'closest',\n";
    out << "          annotations: [\n";
    out << "            { x: data.trim5p + 1, y: 1, yref: 'paper', text: 'trim5p=' + data.trim5p, showarrow: true, arrowhead: 2, ax: 0, ay: -20 },\n";
    out << "            { x: positions[positions.length - 1] - data.trim3p, y: 1, yref: 'paper', text: 'trim3p=' + data.trim3p, showarrow: true, arrowhead: 2, ax: 0, ay: -20 }\n";
    out << "          ],\n";
    out << "          legend: { x: 0.7, y: 0.9 }\n";
    out << "        });\n";
    out << "    } catch (e) {\n";
    out << "      qcFail('Failed to render plots: ' + (e && e.message ? e.message : String(e)));\n";
    out << "    }\n";
    out << "  </script>\n";
    out << "</body>\n";
    out << "</html>\n";
    
    return out.good();
}
