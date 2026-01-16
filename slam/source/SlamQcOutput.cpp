#include "SlamQcOutput.h"
#include "SlamQuant.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <tuple>
#include <algorithm>

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
                     int trim5p,
                     int trim3p,
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
    
    const auto& stats = analyzer.getStats();
    
    // Find max position
    uint32_t maxPos = 0;
    for (const auto& kv : stats) {
        if (kv.first > maxPos) {
            maxPos = kv.first;
        }
    }
    
    out << "{\n";
    out << "  \"version\": \"1.2\",\n";
    out << "  \"algorithm\": \"segmented_regression\",\n";
    out << "  \"file_index\": " << fileIndex << ",\n";
    out << "  \"trim_scope\": \"" << trimScope << "\",\n";
    if (!trimSource.empty()) {
        out << "  \"trim_source\": \"" << trimSource << "\",\n";
    }
    out << "  \"trim5p\": " << trim5p << ",\n";
    out << "  \"trim3p\": " << trim3p << ",\n";
    out << "  \"reads_analyzed\": " << readsAnalyzed << ",\n";
    
    // Include SNP error rate estimation
    if (snpErrEst > 0.0 || snpErrUsed > 0.0) {
        out << "  \"snp_err_est\": " << std::fixed << std::setprecision(6) << snpErrEst << ",\n";
        out << "  \"snp_err_used\": " << std::fixed << std::setprecision(6) << snpErrUsed << ",\n";
        if (!snpErrFallbackReason.empty()) {
            out << "  \"snp_err_fallback_reason\": \"" << snpErrFallbackReason << "\",\n";
        }
    }
    
    // Include segmented regression info if available
    if (trimResult != nullptr) {
        out << "  \"segmented_regression\": {\n";
        out << "    \"breakpoint_b1\": " << trimResult->kneeBin5p << ",\n";
        out << "    \"breakpoint_b2\": " << trimResult->kneeBinMid << ",\n";
        if (trimResult->mode == "auto_segmented_bic4" || trimResult->mode == "auto_segmented_halves_bic2") {
            out << "    \"breakpoint_b3\": " << trimResult->kneeBin3p << ",\n";
        }
        out << "    \"total_sse\": " << std::fixed << std::setprecision(6) << trimResult->totalSSE << ",\n";
        out << "    \"mode\": \"" << trimResult->mode << "\",\n";
        out << "    \"segment1\": {\n";
        out << "      \"slope\": " << std::fixed << std::setprecision(6) << trimResult->seg1.slope << ",\n";
        out << "      \"intercept\": " << std::fixed << std::setprecision(6) << trimResult->seg1.intercept << ",\n";
        out << "      \"sse\": " << std::fixed << std::setprecision(6) << trimResult->seg1.sse << "\n";
        out << "    },\n";
        out << "    \"segment2\": {\n";
        out << "      \"slope\": " << std::fixed << std::setprecision(6) << trimResult->seg2.slope << ",\n";
        out << "      \"intercept\": " << std::fixed << std::setprecision(6) << trimResult->seg2.intercept << ",\n";
        out << "      \"sse\": " << std::fixed << std::setprecision(6) << trimResult->seg2.sse << "\n";
        out << "    },\n";
        out << "    \"segment3\": {\n";
        out << "      \"slope\": " << std::fixed << std::setprecision(6) << trimResult->seg3.slope << ",\n";
        out << "      \"intercept\": " << std::fixed << std::setprecision(6) << trimResult->seg3.intercept << ",\n";
        out << "      \"sse\": " << std::fixed << std::setprecision(6) << trimResult->seg3.sse << "\n";
        out << "    },\n";
        if (trimResult->mode == "auto_segmented_bic4" || trimResult->mode == "auto_segmented_halves_bic2") {
            out << "    \"segment4\": {\n";
            out << "      \"slope\": " << std::fixed << std::setprecision(6) << trimResult->seg4.slope << ",\n";
            out << "      \"intercept\": " << std::fixed << std::setprecision(6) << trimResult->seg4.intercept << ",\n";
            out << "      \"sse\": " << std::fixed << std::setprecision(6) << trimResult->seg4.sse << "\n";
            out << "    },\n";
        }
        
        // Include smoothed curve if available
        if (!trimResult->smoothedCurve.empty()) {
            out << "    \"smoothed_stdev_curve\": [";
            for (size_t i = 0; i < trimResult->smoothedCurve.size(); ++i) {
                if (i > 0) out << ", ";
                if (std::isnan(trimResult->smoothedCurve[i])) {
                    out << "null";
                } else {
                    out << std::fixed << std::setprecision(6) << trimResult->smoothedCurve[i];
                }
            }
            out << "]\n";
        } else {
            out << "    \"smoothed_stdev_curve\": []\n";
        }
        out << "  },\n";
    }
    
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
    out << "        // Summary\n";
    out << "        const summaryHtml = `\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Algorithm</div><div class=\"summary-value\">${data.algorithm || 'segmented_regression'}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Trim 5'</div><div class=\"summary-value\">${data.trim5p}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Trim 3'</div><div class=\"summary-value\">${data.trim3p}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Reads Analyzed</div><div class=\"summary-value\">${data.reads_analyzed.toLocaleString()}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">File Index</div><div class=\"summary-value\">${data.file_index}</div></div>\n";
    out << "        `;\n";
    out << "        document.getElementById('summary').innerHTML = summaryHtml;\n";
    out << "        \n";
    out << "        const positions = data.positions.map(p => p.position);\n";
    out << "        const stddevTcRate = data.positions.map(p => p.stddev_tc_rate);\n";
    out << "        const meanQual = data.positions.map(p => p.mean_qual);\n";
    out << "        const stddevQual = data.positions.map(p => p.stddev_qual);\n";
    out << "        \n";
    out << "        // T→C Stdev plot with segmented regression\n";
    out << "        const traces1 = [{\n";
    out << "          x: positions,\n";
    out << "          y: stddevTcRate,\n";
    out << "          type: 'scatter',\n";
    out << "          mode: 'lines+markers',\n";
    out << "          name: 'T→C Stdev (raw)',\n";
    out << "          marker: { size: 4, color: '#2196F3' },\n";
    out << "          line: { width: 1 }\n";
    out << "        }];\n";
    out << "        \n";
    out << "        // Add smoothed curve and segment fits if available\n";
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
    out << "          \n";
    out << "          // Segment fits\n";
    out << "          const b1 = data.segmented_regression.breakpoint_b1;\n";
    out << "          const b2 = data.segmented_regression.breakpoint_b2;\n";
    out << "          const b3 = data.segmented_regression.breakpoint_b3;\n";
    out << "          const seg1 = data.segmented_regression.segment1;\n";
    out << "          const seg2 = data.segmented_regression.segment2;\n";
    out << "          const seg3 = data.segmented_regression.segment3;\n";
    out << "          const seg4 = data.segmented_regression.segment4;\n";
    out << "          \n";
    out << "          if (b1 > 0) {\n";
    out << "            const x1 = Array.from({length: b1}, (_, i) => i);\n";
    out << "            const y1 = x1.map(x => seg1.slope * x + seg1.intercept);\n";
    out << "            traces1.push({ x: x1, y: y1, type: 'scatter', mode: 'lines', name: 'Seg1 (5\\' artifact)', line: { width: 2, color: '#f44336' } });\n";
    out << "          }\n";
    out << "          \n";
    out << "          const x2 = Array.from({length: b2 - b1 + 1}, (_, i) => b1 + i);\n";
    out << "          const y2 = x2.map(x => seg2.slope * x + seg2.intercept);\n";
    out << "          traces1.push({ x: x2, y: y2, type: 'scatter', mode: 'lines', name: 'Seg2 (signal)', line: { width: 2, color: '#4CAF50' } });\n";
    out << "          \n";
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
    out << "        \n";
    out << "        const shapes = [];\n";
    out << "        if (data.trim5p > 0) {\n";
    out << "          shapes.push({ type: 'line', x0: data.trim5p, x1: data.trim5p, y0: 0, y1: 1, yref: 'paper', line: { color: 'purple', width: 2, dash: 'dot' } });\n";
    out << "        }\n";
    out << "        if (data.trim3p > 0) {\n";
    out << "          const trim3pPos = positions[positions.length - 1] - data.trim3p;\n";
    out << "          shapes.push({ type: 'line', x0: trim3pPos, x1: trim3pPos, y0: 0, y1: 1, yref: 'paper', line: { color: 'purple', width: 2, dash: 'dot' } });\n";
    out << "        }\n";
    out << "        \n";
    out << "        Plotly.newPlot('tcStddevPlot', traces1, {\n";
    out << "          title: 'T→C Rate Standard Deviation by Position (Segmented Regression)',\n";
    out << "          xaxis: { title: 'Read Position' },\n";
    out << "          yaxis: { title: 'T→C Stdev' },\n";
    out << "          shapes: shapes,\n";
    out << "          hovermode: 'closest',\n";
    out << "          annotations: [\n";
    out << "            { x: data.trim5p, y: 1, yref: 'paper', text: 'trim5p=' + data.trim5p, showarrow: false, yanchor: 'bottom' },\n";
    out << "            { x: positions[positions.length - 1] - data.trim3p, y: 1, yref: 'paper', text: 'trim3p=' + data.trim3p, showarrow: false, yanchor: 'bottom' }\n";
    out << "          ]\n";
    out << "        });\n";
    out << "        \n";
    out << "        // Quality plot\n";
    out << "        Plotly.newPlot('qualityPlot', [\n";
    out << "          { x: positions, y: meanQual, type: 'scatter', mode: 'lines+markers', name: 'Mean Quality', marker: { size: 4 } },\n";
    out << "          { x: positions, y: stddevQual, type: 'scatter', mode: 'lines+markers', name: 'Quality Stdev', yaxis: 'y2', marker: { size: 4 } }\n";
    out << "        ], {\n";
    out << "          title: 'Quality Statistics by Position',\n";
    out << "          xaxis: { title: 'Read Position' },\n";
    out << "          yaxis: { title: 'Mean Quality', side: 'left' },\n";
    out << "          yaxis2: { title: 'Quality Stdev', side: 'right', overlaying: 'y' },\n";
    out << "          shapes: shapes,\n";
    out << "          hovermode: 'closest'\n";
    out << "        });\n";
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
                                   int trim5p,
                                   int trim3p,
                                   const SlamVarianceTrimResult* trimResult,
                                   const std::vector<double>* varianceStddevTcRate) {
    std::ofstream out(outputPath.c_str());
    if (!out.good()) {
        return false;
    }
    
    // Get position transition data (T→C and T→A)
    auto transitionData = slamQuant.getPositionTransitionData();
    
    // Get variance analyzer stats (PHRED quality, T counts)
    const SlamVarianceAnalyzer* analyzer = slamQuant.varianceAnalyzer();
    std::unordered_map<uint32_t, SlamPositionVarianceStats> varianceStats;
    if (analyzer != nullptr) {
        varianceStats = analyzer->getStats();
    }
    
    // Find max position
    uint32_t maxPos = 0;
    for (const auto& kv : transitionData) {
        if (kv.first > maxPos) maxPos = kv.first;
    }
    for (const auto& kv : varianceStats) {
        if (kv.first > maxPos) maxPos = kv.first;
    }
    
    out << "{\n";
    out << "  \"version\": \"2.0\",\n";
    out << "  \"type\": \"comprehensive_qc\",\n";
    out << "  \"trim5p\": " << trim5p << ",\n";
    out << "  \"trim3p\": " << trim3p << ",\n";
    
    // Segmented regression info
    if (trimResult != nullptr) {
        out << "  \"segmented_regression\": {\n";
        out << "    \"breakpoint_b1\": " << trimResult->kneeBin5p << ",\n";
        out << "    \"breakpoint_b2\": " << trimResult->kneeBinMid << ",\n";
        if (trimResult->mode == "auto_segmented_bic4" || trimResult->mode == "auto_segmented_halves_bic2") {
            out << "    \"breakpoint_b3\": " << trimResult->kneeBin3p << ",\n";
        }
        out << "    \"total_sse\": " << std::fixed << std::setprecision(6) << trimResult->totalSSE << ",\n";
        out << "    \"mode\": \"" << trimResult->mode << "\",\n";
        out << "    \"segment1\": {\n";
        out << "      \"slope\": " << std::fixed << std::setprecision(6) << trimResult->seg1.slope << ",\n";
        out << "      \"intercept\": " << std::fixed << std::setprecision(6) << trimResult->seg1.intercept << "\n";
        out << "    },\n";
        out << "    \"segment2\": {\n";
        out << "      \"slope\": " << std::fixed << std::setprecision(6) << trimResult->seg2.slope << ",\n";
        out << "      \"intercept\": " << std::fixed << std::setprecision(6) << trimResult->seg2.intercept << "\n";
        out << "    },\n";
        out << "    \"segment3\": {\n";
        out << "      \"slope\": " << std::fixed << std::setprecision(6) << trimResult->seg3.slope << ",\n";
        out << "      \"intercept\": " << std::fixed << std::setprecision(6) << trimResult->seg3.intercept << "\n";
        out << "    },\n";
        if (trimResult->mode == "auto_segmented_bic4" || trimResult->mode == "auto_segmented_halves_bic2") {
            out << "    \"segment4\": {\n";
            out << "      \"slope\": " << std::fixed << std::setprecision(6) << trimResult->seg4.slope << ",\n";
            out << "      \"intercept\": " << std::fixed << std::setprecision(6) << trimResult->seg4.intercept << "\n";
            out << "    },\n";
        }
        if (!trimResult->smoothedCurve.empty()) {
            out << "    \"smoothed_stdev_curve\": [";
            for (size_t i = 0; i < trimResult->smoothedCurve.size(); ++i) {
                if (i > 0) out << ", ";
                if (std::isnan(trimResult->smoothedCurve[i])) {
                    out << "null";
                } else {
                    out << std::fixed << std::setprecision(6) << trimResult->smoothedCurve[i];
                }
            }
            out << "]\n";
        } else {
            out << "    \"smoothed_stdev_curve\": []\n";
        }
        out << "  },\n";
    } else {
        out << "  \"segmented_regression\": null,\n";
    }
    
    out << "  \"positions\": [\n";
    
    bool first = true;
    for (uint32_t pos = 0; pos <= maxPos; ++pos) {
        // Skip positions with no data
        bool hasTransition = transitionData.count(pos) > 0;
        bool hasVariance = varianceStats.count(pos) > 0;
        if (!hasTransition && !hasVariance) {
            continue;
        }
        
        if (!first) {
            out << ",\n";
        }
        first = false;
        
        out << "    {\n";
        out << "      \"position\": " << (pos + 1) << ",\n";  // 1-based for plotting
        
        // T→C and T→A rates
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
        
        // PHRED quality and T counts from variance analyzer
        if (hasVariance) {
            const auto& vstats = varianceStats.at(pos);
            double meanQual = vstats.meanQual();
            double stddevQual = vstats.stddevQual();
            double meanTcRate = vstats.meanTcRate();
            double stddevTcRate = vstats.stddevTcRate();
            
            out << "      \"mean_qual\": " << std::fixed << std::setprecision(2) << meanQual << ",\n";
            out << "      \"stddev_qual\": " << std::fixed << std::setprecision(2) << stddevQual << ",\n";
            out << "      \"t_count\": " << vstats.tCount << ",\n";
            out << "      \"mean_tc_rate\": " << std::fixed << std::setprecision(6) << meanTcRate << ",\n";
            out << "      \"stddev_tc_rate\": " << std::fixed << std::setprecision(6) << stddevTcRate << "\n";
        } else {
            out << "      \"mean_qual\": null,\n";
            out << "      \"stddev_qual\": null,\n";
            out << "      \"t_count\": 0,\n";
            out << "      \"mean_tc_rate\": null,\n";
            if (varianceStddevTcRate && pos < varianceStddevTcRate->size()) {
                double v = (*varianceStddevTcRate)[pos];
                if (std::isnan(v)) {
                    out << "      \"stddev_tc_rate\": null\n";
                } else {
                    out << "      \"stddev_tc_rate\": " << std::fixed << std::setprecision(6) << v << "\n";
                }
            } else {
                out << "      \"stddev_tc_rate\": null\n";
            }
        }
        
        out << "    }";
    }
    
    out << "\n  ]\n";
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
    out << "        document.getElementById('summary').innerHTML = `\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Trim 5'</div><div class=\"summary-value\">${data.trim5p}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Trim 3'</div><div class=\"summary-value\">${data.trim3p}</div></div>\n";
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
