#include "TrimQcOutput.h"
#include <fstream>
#include <iomanip>
#include <sstream>

static std::string basenameFromPath(const std::string& path) {
    size_t lastSlash = path.find_last_of("/\\");
    if (lastSlash == std::string::npos) {
        return path;
    }
    return path.substr(lastSlash + 1);
}

bool writeTrimQcJson(const TrimQcCollector& qc,
                     const Stats& stats,
                     const std::string& outputPath,
                     const std::string& stageLabel) {
    std::ofstream out(outputPath.c_str());
    if (!out.good()) {
        return false;
    }

    out << "{\n";
    out << "  \"version\": \"1.0\",\n";
    out << "  \"stage\": \"" << stageLabel << "\",\n";
    out << "  \"read_mates\": " << qc.mateCount() << ",\n";
    out << "  \"total_reads\": " << qc.totalReads() << ",\n";
    out << "  \"quality_base\": " << qc.qualityBase() << ",\n";
    out << "  \"trim_stats\": {\n";
    out << "    \"reads_processed\": " << stats.trimReadsProcessed << ",\n";
    out << "    \"reads_trimmed\": " << stats.trimReadsTrimmed << ",\n";
    out << "    \"reads_too_short\": " << stats.trimReadsTooShort << ",\n";
    out << "    \"bases_quality_trimmed\": " << stats.trimBasesQualityTrimmed << ",\n";
    out << "    \"bases_adapter_trimmed\": " << stats.trimBasesAdapterTrimmed << ",\n";
    out << "    \"pairs_processed\": " << stats.trimPairsProcessed << ",\n";
    out << "    \"pairs_dropped\": " << stats.trimPairsDropped << ",\n";
    out << "    \"pairs_kept\": " << stats.trimPairsKept << "\n";
    out << "  },\n";
    out << "  \"mates\": [\n";

    for (size_t m = 0; m < qc.mateCount(); ++m) {
        const auto& mate = qc.mate(m);
        if (m > 0) {
            out << ",\n";
        }
        out << "    {\n";
        out << "      \"mate\": " << (m + 1) << ",\n";
        out << "      \"reads\": " << mate.reads << ",\n";
        out << "      \"length_hist\": [";
        for (size_t i = 0; i < mate.lengthHist.size(); ++i) {
            if (i > 0) out << ", ";
            out << mate.lengthHist[i];
        }
        out << "],\n";
        out << "      \"gc_hist\": [";
        for (size_t i = 0; i < mate.gcHist.size(); ++i) {
            if (i > 0) out << ", ";
            out << mate.gcHist[i];
        }
        out << "],\n";
        out << "      \"positions\": [\n";
        bool firstPos = true;
        for (size_t i = 0; i < mate.posCount.size(); ++i) {
            if (mate.posCount[i] == 0) {
                continue;
            }
            if (!firstPos) {
                out << ",\n";
            }
            firstPos = false;
            double meanQual = mate.posCount[i] > 0 ?
                static_cast<double>(mate.posQualSum[i]) / static_cast<double>(mate.posCount[i]) : 0.0;
            out << "        {\n";
            out << "          \"position\": " << (i + 1) << ",\n";
            out << "          \"count\": " << mate.posCount[i] << ",\n";
            out << "          \"mean_qual\": " << std::fixed << std::setprecision(2) << meanQual << ",\n";
            out << "          \"base_counts\": {";
            out << "\"A\": " << mate.posBaseCounts[i][0] << ", ";
            out << "\"C\": " << mate.posBaseCounts[i][1] << ", ";
            out << "\"G\": " << mate.posBaseCounts[i][2] << ", ";
            out << "\"T\": " << mate.posBaseCounts[i][3] << ", ";
            out << "\"N\": " << mate.posBaseCounts[i][4];
            out << "}\n";
            out << "        }";
        }
        out << "\n      ]\n";
        out << "    }";
    }

    out << "\n  ]\n";
    out << "}\n";

    return out.good();
}

bool writeTrimQcHtml(const std::string& jsonPath,
                     const std::string& htmlPath) {
    std::ofstream out(htmlPath.c_str());
    if (!out.good()) {
        return false;
    }

    std::string jsonText;
    {
        std::ifstream in(jsonPath.c_str(), std::ios::in | std::ios::binary);
        if (!in.good()) {
            return false;
        }
        std::ostringstream ss;
        ss << in.rdbuf();
        jsonText = ss.str();
    }

    out << "<!DOCTYPE html>\n";
    out << "<html>\n";
    out << "<head>\n";
    out << "  <title>Trim QC Report</title>\n";
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
    out << "    <h1>Trim QC Report</h1>\n";
    out << "    <div class=\"summary\" id=\"summary\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"meanQualPlot\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"baseCompPlot\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"lengthPlot\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"gcPlot\"></div>\n";
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
    out << "        const summaryHtml = `\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Stage</div><div class=\"summary-value\">${data.stage}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Total Reads</div><div class=\"summary-value\">${data.total_reads.toLocaleString()}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Mates</div><div class=\"summary-value\">${data.read_mates}</div></div>\n";
    out << "        `;\n";
    out << "        document.getElementById('summary').innerHTML = summaryHtml;\n";
    out << "\n";
    out << "        const meanQualTraces = [];\n";
    out << "        const baseCompTraces = [];\n";
    out << "        const lengthTraces = [];\n";
    out << "        const gcTraces = [];\n";
    out << "        const baseKeys = ['A','C','G','T','N'];\n";
    out << "        const baseColors = {A:'#4caf50', C:'#2196f3', G:'#ff9800', T:'#f44336', N:'#9e9e9e'};\n";
    out << "\n";
    out << "        data.mates.forEach((mate, idx) => {\n";
    out << "          const positions = mate.positions.map(p => p.position);\n";
    out << "          const meanQual = mate.positions.map(p => p.mean_qual);\n";
    out << "          meanQualTraces.push({\n";
    out << "            x: positions,\n";
    out << "            y: meanQual,\n";
    out << "            type: 'scatter',\n";
    out << "            mode: 'lines',\n";
    out << "            name: `Mate ${mate.mate}`\n";
    out << "          });\n";
    out << "\n";
    out << "          baseKeys.forEach(base => {\n";
    out << "            const basePct = mate.positions.map(p => {\n";
    out << "              const total = p.count || 1;\n";
    out << "              return 100.0 * (p.base_counts[base] || 0) / total;\n";
    out << "            });\n";
    out << "            baseCompTraces.push({\n";
    out << "              x: positions,\n";
    out << "              y: basePct,\n";
    out << "              type: 'scatter',\n";
    out << "              mode: 'lines',\n";
    out << "              name: `Mate ${mate.mate} ${base}`,\n";
    out << "              line: {color: baseColors[base]}\n";
    out << "            });\n";
    out << "          });\n";
    out << "\n";
    out << "          const lengths = mate.length_hist.map((_, i) => i);\n";
    out << "          lengthTraces.push({\n";
    out << "            x: lengths,\n";
    out << "            y: mate.length_hist,\n";
    out << "            type: 'bar',\n";
    out << "            name: `Mate ${mate.mate}`\n";
    out << "          });\n";
    out << "\n";
    out << "          const gcBins = mate.gc_hist.map((_, i) => i);\n";
    out << "          gcTraces.push({\n";
    out << "            x: gcBins,\n";
    out << "            y: mate.gc_hist,\n";
    out << "            type: 'bar',\n";
    out << "            name: `Mate ${mate.mate}`\n";
    out << "          });\n";
    out << "        });\n";
    out << "\n";
    out << "        Plotly.newPlot('meanQualPlot', meanQualTraces, {\n";
    out << "          title: 'Per-base mean quality',\n";
    out << "          xaxis: {title: 'Position'},\n";
    out << "          yaxis: {title: 'Mean Phred'}\n";
    out << "        });\n";
    out << "\n";
    out << "        Plotly.newPlot('baseCompPlot', baseCompTraces, {\n";
    out << "          title: 'Per-base base composition',\n";
    out << "          xaxis: {title: 'Position'},\n";
    out << "          yaxis: {title: 'Percent'}\n";
    out << "        });\n";
    out << "\n";
    out << "        Plotly.newPlot('lengthPlot', lengthTraces, {\n";
    out << "          title: 'Read length distribution',\n";
    out << "          xaxis: {title: 'Length'},\n";
    out << "          yaxis: {title: 'Reads'}\n";
    out << "        });\n";
    out << "\n";
    out << "        Plotly.newPlot('gcPlot', gcTraces, {\n";
    out << "          title: 'GC content distribution',\n";
    out << "          xaxis: {title: 'GC%'} ,\n";
    out << "          yaxis: {title: 'Reads'}\n";
    out << "        });\n";
    out << "    } catch (e) {\n";
    out << "      qcFail('Failed to render plots: ' + (e && e.message ? e.message : String(e)));\n";
    out << "    }\n";
    out << "  </script>\n";
    out << "</body>\n";
    out << "</html>\n";

    return out.good();
}
