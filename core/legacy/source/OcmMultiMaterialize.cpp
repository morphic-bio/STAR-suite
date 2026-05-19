#include "OcmMultiMaterialize.h"
#include "OcmMultiConfig.h"
#include "PfMultiMerge.h"
#include "VelocytoMexWriter.h"
#include "streamFuns.h"
#include "VERSION"
#include "SoloFeature.h"
#include "TimeFunctions.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unordered_map>
#include <zlib.h>

namespace {

static const int kGzLevel = 3;

static bool isUnsetTokenLocal(const string& input) {
    return input.empty() || input == "-";
}

static string normalizeEnableMode(string mode) {
    std::transform(mode.begin(), mode.end(), mode.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return mode;
}

static bool hasMexFiles(const string& dirPath) {
    struct stat st;
    string features = dirPath + "/features.tsv";
    string featuresGz = dirPath + "/features.tsv.gz";
    return (stat(features.c_str(), &st) == 0) || (stat(featuresGz.c_str(), &st) == 0);
}

static void resolveGexMexDirs(const Parameters& P,
                              const string& soloOut,
                              string& rawOut,
                              string& filteredOut,
                              string& gexFeatureLabel) {
    const string geneOut = soloOut + "/Gene";
    const string geneFullOut = soloOut + "/GeneFull";
    rawOut = geneOut + "/raw";
    filteredOut = geneOut + "/filtered";
    gexFeatureLabel = "Gene";

    if (P.pSolo.crGexFeature == ParametersSolo::CrGexGeneFull) {
        rawOut = geneFullOut + "/raw";
        filteredOut = geneFullOut + "/filtered";
        gexFeatureLabel = "GeneFull";
    } else if (P.pSolo.crGexFeature == ParametersSolo::CrGexGene) {
        return;
    } else if (hasMexFiles(geneFullOut + "/raw") || hasMexFiles(geneFullOut + "/filtered")) {
        rawOut = geneFullOut + "/raw";
        filteredOut = geneFullOut + "/filtered";
        gexFeatureLabel = "GeneFull";
    }
}

static PfMultiMerge::MexData subsetMexColumns(const PfMultiMerge::MexData& data,
                                              const vector<uint32_t>& colIndices) {
    PfMultiMerge::MexData out;
    out.features = data.features;
    out.featureNames = data.featureNames;
    out.featureTypes = data.featureTypes;
    vector<uint32_t> oldToNew(data.barcodes.size(), UINT32_MAX);
    out.barcodes.reserve(colIndices.size());
    for (uint32_t oldIdx : colIndices) {
        if (oldIdx >= data.barcodes.size()) {
            continue;
        }
        oldToNew[oldIdx] = static_cast<uint32_t>(out.barcodes.size());
        out.barcodes.push_back(data.barcodes[oldIdx]);
    }
    out.triplets.reserve(data.triplets.size());
    for (const auto& t : data.triplets) {
        if (t.cell_idx < oldToNew.size() && oldToNew[t.cell_idx] != UINT32_MAX) {
            MexWriter::Triplet nt;
            nt.cell_idx = oldToNew[t.cell_idx];
            nt.gene_idx = t.gene_idx;
            nt.count = t.count;
            out.triplets.push_back(nt);
        }
    }
    return out;
}

static void validateFeatureAxesMatch(const PfMultiMerge::MexData& raw,
                                     const PfMultiMerge::MexData& filtered) {
    if (raw.features != filtered.features) {
        throw std::runtime_error("OCM materializer: raw/filtered feature IDs do not match");
    }
    if (raw.featureNames != filtered.featureNames) {
        throw std::runtime_error("OCM materializer: raw/filtered feature names do not match");
    }
    if (raw.featureTypes != filtered.featureTypes) {
        throw std::runtime_error("OCM materializer: raw/filtered feature types do not match");
    }
}

static map<string, vector<uint32_t>> buildTagColumnIndices(const vector<string>& barcodes,
                                                           uint64_t& unknownCount) {
    map<string, vector<uint32_t>> indicesPerTag;
    unknownCount = 0;
    for (size_t i = 0; i < barcodes.size(); ++i) {
        const string tag = OcmMultiMaterialize::classifyBarcodeTag(barcodes[i]);
        if (tag.empty()) {
            unknownCount++;
            continue;
        }
        indicesPerTag[tag].push_back(static_cast<uint32_t>(i));
    }
    return indicesPerTag;
}

static map<string, vector<string>> buildCellsPerTag(const vector<string>& barcodes) {
    map<string, vector<string>> cellsPerTag;
    for (const auto& barcode : barcodes) {
        const string tag = OcmMultiMaterialize::classifyBarcodeTag(barcode);
        if (!tag.empty()) {
            cellsPerTag[tag].push_back(barcode);
        }
    }
    for (auto& kv : cellsPerTag) {
        std::sort(kv.second.begin(), kv.second.end());
    }
    return cellsPerTag;
}

static vector<uint32_t> unionColumnIndices(const map<string, vector<uint32_t>>& indicesPerTag,
                                           const vector<string>& ocmIds) {
    vector<uint32_t> merged;
    std::set<uint32_t> seen;
    for (const auto& ocmId : ocmIds) {
        auto it = indicesPerTag.find(ocmId);
        if (it == indicesPerTag.end()) {
            continue;
        }
        for (uint32_t idx : it->second) {
            if (seen.insert(idx).second) {
                merged.push_back(idx);
            }
        }
    }
    std::sort(merged.begin(), merged.end());
    return merged;
}

static string jsonEscape(const string& value) {
    string out;
    out.reserve(value.size() + 8);
    for (char c : value) {
        if (c == '\\' || c == '"') {
            out.push_back('\\');
        }
        out.push_back(c);
    }
    return out;
}

static void writeCellsPerTagJson(const string& path, const map<string, vector<string>>& cellsPerTag) {
    ofstream out(path.c_str());
    if (!out.is_open()) {
        throw std::runtime_error("Failed to write cells_per_tag.json: " + path);
    }
    out << "{\n";
    bool firstTag = true;
    for (const char* tag : OcmMultiConfig::kValidOcmIds) {
        if (!firstTag) {
            out << ",\n";
        }
        firstTag = false;
        out << "  \"" << tag << "\": [";
        auto it = cellsPerTag.find(tag);
        if (it != cellsPerTag.end()) {
            for (size_t i = 0; i < it->second.size(); ++i) {
                if (i > 0) {
                    out << ", ";
                }
                out << "\"" << jsonEscape(it->second[i]) << "\"";
            }
        }
        out << "]";
    }
    out << "\n}\n";
}

static string resolveFilteredBarcodeGenomeLabel(const PfMultiConfig::Config& config) {
    // Cell Ranger sample_filtered_barcodes.csv uses a genome prefix. Derive from config when
    // possible; JAX human OCM runs use GRCh38.
    if (config.referencePath.find("GRCh38") != string::npos) {
        return "GRCh38";
    }
    return "GRCh38";
}

static void writeSampleFilteredBarcodesCsv(const string& path,
                                           const vector<string>& barcodes,
                                           const string& genomeLabel) {
    ofstream out(path.c_str());
    if (!out.is_open()) {
        throw std::runtime_error("Failed to write sample_filtered_barcodes.csv: " + path);
    }
    for (const auto& barcode : barcodes) {
        out << genomeLabel << "," << barcode << "\n";
    }
}

static string resolveProjectRoot(const string& outPrefix) {
    string prefix = outPrefix;
    if (!prefix.empty() && prefix.back() != '/') {
        prefix.push_back('/');
    }
    const string runSuffix = "/run/";
    size_t runPos = prefix.rfind(runSuffix);
    if (runPos != string::npos) {
        const string samplesToken = "/samples/";
        size_t samplesPos = prefix.rfind(samplesToken, runPos);
        if (samplesPos != string::npos) {
            return prefix.substr(0, samplesPos + 1);
        }
        if (runPos == 0) {
            return "./";
        }
        return prefix.substr(0, runPos + 1);
    }
    size_t lastSlash = prefix.find_last_of('/');
    if (lastSlash == string::npos || lastSlash == 0) {
        return "./";
    }
    return prefix.substr(0, lastSlash + 1);
}

static bool writeGzLines(const string& path, const vector<string>& lines) {
    gzFile gz = gzopen(path.c_str(), "wb");
    if (gz == nullptr) {
        return false;
    }
    gzbuffer(gz, 1 << 20);
    gzsetparams(gz, kGzLevel, Z_DEFAULT_STRATEGY);
    for (const auto& line : lines) {
        if (gzwrite(gz, line.data(), line.size()) <= 0) {
            gzclose(gz);
            return false;
        }
        if (gzwrite(gz, "\n", 1) <= 0) {
            gzclose(gz);
            return false;
        }
    }
    return gzclose(gz) == Z_OK;
}

static int writeEmptyGexMexGz(const string& outputDir,
                            const PfMultiMerge::MexData& featureAxis,
                            Parameters& P) {
    createDirectory(outputDir + "/", P.runDirPerm, "OCM MEX output", P);
    vector<string> featureLines;
    featureLines.reserve(featureAxis.features.size());
    for (size_t i = 0; i < featureAxis.features.size(); ++i) {
        const string name =
            (i < featureAxis.featureNames.size()) ? featureAxis.featureNames[i] : featureAxis.features[i];
        const string type =
            (i < featureAxis.featureTypes.size()) ? featureAxis.featureTypes[i] : "Gene Expression";
        featureLines.push_back(featureAxis.features[i] + "\t" + name + "\t" + type);
    }
    {
        gzFile gz = gzopen((outputDir + "/features.tsv.gz").c_str(), "wb");
        if (gz == nullptr) {
            return -1;
        }
        gzbuffer(gz, 1 << 20);
        gzsetparams(gz, kGzLevel, Z_DEFAULT_STRATEGY);
        for (const auto& line : featureLines) {
            if (gzwrite(gz, line.data(), line.size()) <= 0 || gzwrite(gz, "\n", 1) <= 0) {
                gzclose(gz);
                return -1;
            }
        }
        if (gzclose(gz) != Z_OK) {
            return -1;
        }
    }
    if (!writeGzLines(outputDir + "/barcodes.tsv.gz", vector<string>())) {
        return -1;
    }
    {
        gzFile gz = gzopen((outputDir + "/matrix.mtx.gz").c_str(), "wb");
        if (gz == nullptr) {
            return -1;
        }
        gzbuffer(gz, 1 << 20);
        gzsetparams(gz, kGzLevel, Z_DEFAULT_STRATEGY);
        ostringstream hdr;
        hdr << "%%MatrixMarket matrix coordinate integer general\n%\n"
            << featureAxis.features.size() << " 0 0\n";
        const string hdrStr = hdr.str();
        if (gzwrite(gz, hdrStr.data(), hdrStr.size()) <= 0 || gzclose(gz) != Z_OK) {
            return -1;
        }
    }
    return 0;
}

static int writeMexGzDir(const string& outputDir,
                        const PfMultiMerge::MexData& data,
                        ofstream& logStream,
                        Parameters& P) {
    createDirectory(outputDir + "/", P.runDirPerm, "OCM MEX output", P);
    if (data.barcodes.empty()) {
        return writeEmptyGexMexGz(outputDir, data, P);
    }
    return PfMultiMerge::writeCombinedMex(outputDir, data, "1", logStream, data.barcodes);
}

static string barcodeJoinKey(const string& barcode) {
    const size_t dashPos = barcode.find_last_of('-');
    if (dashPos == string::npos || dashPos + 1 >= barcode.size()) {
        return barcode;
    }
    for (size_t i = dashPos + 1; i < barcode.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(barcode[i]))) {
            return barcode;
        }
    }
    return barcode.substr(0, dashPos);
}

static vector<uint32_t> mapGexColumnIndicesToVelocyto(const vector<uint32_t>& gexColIndices,
                                                      const vector<string>& gexBarcodes,
                                                      const vector<string>& velocytoBarcodes) {
    unordered_map<string, uint32_t> veloIndex;
    veloIndex.reserve(velocytoBarcodes.size() * 2);
    for (size_t i = 0; i < velocytoBarcodes.size(); ++i) {
        veloIndex.emplace(barcodeJoinKey(velocytoBarcodes[i]), static_cast<uint32_t>(i));
        veloIndex.emplace(velocytoBarcodes[i], static_cast<uint32_t>(i));
    }

    vector<uint32_t> veloCols;
    veloCols.reserve(gexColIndices.size());
    for (uint32_t gexIdx : gexColIndices) {
        if (gexIdx >= gexBarcodes.size()) {
            ostringstream err;
            err << "OCM materializer: GEX column index " << gexIdx << " out of range";
            throw std::runtime_error(err.str());
        }
        const string& gexBarcode = gexBarcodes[gexIdx];
        auto it = veloIndex.find(barcodeJoinKey(gexBarcode));
        if (it == veloIndex.end()) {
            it = veloIndex.find(gexBarcode);
        }
        if (it == veloIndex.end()) {
            ostringstream err;
            err << "OCM materializer: GEX barcode not present in Velocyto matrix: " << gexBarcode;
            throw std::runtime_error(err.str());
        }
        veloCols.push_back(it->second);
    }
    return veloCols;
}

static void validateVelocytoBarcodes(const vector<string>& velocytoBarcodes,
                                     const vector<string>& gexBarcodes) {
    unordered_map<string, uint32_t> gexIndex;
    gexIndex.reserve(gexBarcodes.size() * 2);
    for (size_t i = 0; i < gexBarcodes.size(); ++i) {
        gexIndex.emplace(barcodeJoinKey(gexBarcodes[i]), static_cast<uint32_t>(i));
        gexIndex.emplace(gexBarcodes[i], static_cast<uint32_t>(i));
    }
    size_t missing = 0;
    for (const auto& bc : velocytoBarcodes) {
        if (gexIndex.find(barcodeJoinKey(bc)) == gexIndex.end() && gexIndex.find(bc) == gexIndex.end()) {
            missing++;
        }
    }
    if (missing > 0) {
        ostringstream err;
        err << "OCM materializer: " << missing << " of " << velocytoBarcodes.size()
            << " Velocyto barcodes could not be joined to the GeneFull barcode namespace "
            << "(gex_barcodes=" << gexBarcodes.size() << ")";
        throw std::runtime_error(err.str());
    }
}

static void writeMaterializationSummary(const string& path,
                                        const string& configPath,
                                        const string& gexFeatureLabel,
                                        const map<string, vector<uint32_t>>& rawTagIndices,
                                        const map<string, vector<uint32_t>>& filteredTagIndices,
                                        uint64_t rawUnknown,
                                        uint64_t filteredUnknown,
                                        const PfMultiConfig::Config& config,
                                        const map<string, size_t>& perSampleRawCounts,
                                        const map<string, size_t>& perSampleFilteredCounts) {
    ofstream out(path.c_str());
    if (!out.is_open()) {
        throw std::runtime_error("Failed to write ocm_materialization_summary.json: " + path);
    }
    out << "{\n";
    out << "  \"config_path\": \"" << jsonEscape(configPath) << "\",\n";
    out << "  \"gex_feature_surface\": \"" << jsonEscape(gexFeatureLabel) << "\",\n";
    out << "  \"star_version\": \"" << jsonEscape(STAR_VERSION) << "\",\n";
    out << "  \"raw_unknown_overhangs\": " << rawUnknown << ",\n";
    out << "  \"filtered_unknown_overhangs\": " << filteredUnknown << ",\n";
    out << "  \"per_tag\": {\n";
    bool firstTag = true;
    for (const char* tag : OcmMultiConfig::kValidOcmIds) {
        if (!firstTag) {
            out << ",\n";
        }
        firstTag = false;
        size_t rawCount = 0;
        size_t filteredCount = 0;
        auto rawIt = rawTagIndices.find(tag);
        auto filtIt = filteredTagIndices.find(tag);
        if (rawIt != rawTagIndices.end()) {
            rawCount = rawIt->second.size();
        }
        if (filtIt != filteredTagIndices.end()) {
            filteredCount = filtIt->second.size();
        }
        out << "    \"" << tag << "\": {\"raw\": " << rawCount << ", \"filtered\": " << filteredCount << "}";
    }
    out << "\n  },\n  \"per_sample\": {\n";
    bool firstSample = true;
    for (const auto& sample : config.samples) {
        if (!firstSample) {
            out << ",\n";
        }
        firstSample = false;
        size_t rawCount = 0;
        size_t filteredCount = 0;
        auto rawIt = perSampleRawCounts.find(sample.sample_id);
        auto filtIt = perSampleFilteredCounts.find(sample.sample_id);
        if (rawIt != perSampleRawCounts.end()) {
            rawCount = rawIt->second;
        }
        if (filtIt != perSampleFilteredCounts.end()) {
            filteredCount = filtIt->second;
        }
        out << "    \"" << jsonEscape(sample.sample_id) << "\": {\"ocm_ids\": \""
            << jsonEscape(sample.ocm_barcode_ids) << "\", \"raw\": " << rawCount
            << ", \"filtered\": " << filteredCount << "}";
    }
    out << "\n  }\n}\n";
}

} // namespace

namespace OcmMultiMaterialize {

static string stripBarcodeSuffixForClassification(const string& barcode) {
    return barcodeJoinKey(barcode);
}

string classifyBarcodeTag(const string& barcode) {
    const string core = stripBarcodeSuffixForClassification(barcode);
    if (core.size() < 9) {
        return string();
    }
    const string overhang = core.substr(7, 2);
    if (overhang == "GT") {
        return "OB1";
    }
    if (overhang == "CA") {
        return "OB2";
    }
    if (overhang == "TC") {
        return "OB3";
    }
    if (overhang == "AG") {
        return "OB4";
    }
    return string();
}

} // namespace OcmMultiMaterialize

static bool isCellRangerOutputCompat(const string& mode) {
    string normalized = mode;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return normalized.empty() || normalized == "cellranger";
}

int runOcmMultiMaterialize(Parameters& P) {
    const string enableMode = normalizeEnableMode(P.pfMulti.ocmMultiEnable);
    if (enableMode.empty() || enableMode == "no") {
        return 0;
    }
    if (!isCellRangerOutputCompat(P.pfMulti.ocmMultiOutputCompat)) {
        P.inOut->logMain << "EXITING because of fatal PARAMETERS error: unsupported --ocmMultiOutputCompat="
                         << P.pfMulti.ocmMultiOutputCompat << " (only cellranger is implemented)\n";
        return 1;
    }
    if (enableMode != "yes" && enableMode != "auto") {
        P.inOut->logMain << "EXITING because of fatal PARAMETERS error: unrecognized --ocmMultiEnable="
                         << P.pfMulti.ocmMultiEnable << " (expected no|yes|auto)\n";
        return 1;
    }

    string configPath;
    if (!isUnsetTokenLocal(P.pfMulti.ocmMultiConfig)) {
        configPath = P.pfMulti.ocmMultiConfig;
    } else if (!isUnsetTokenLocal(P.pfMulti.pfMultiConfig)) {
        configPath = P.pfMulti.pfMultiConfig;
        P.inOut->logMain << "NOTICE: --ocmMultiConfig unset; reusing --pfMultiConfig path: "
                         << configPath << "\n";
    } else if (enableMode == "yes") {
        P.inOut->logMain << "EXITING because of fatal PARAMETERS error: --ocmMultiEnable yes requires "
                         << "--ocmMultiConfig (or --pfMultiConfig fallback)\n";
        return 1;
    } else {
        P.inOut->logMain << "WARNING: --ocmMultiEnable auto with no config path; skipping OCM materialization\n";
        return 0;
    }

    PfMultiConfig::Config config;
    try {
        config = OcmMultiConfig::parseAndValidate(configPath, P.inOut->logMain);
    } catch (const std::exception& ex) {
        if (enableMode == "auto") {
            P.inOut->logMain << "WARNING: OCM materialization skipped in auto mode: " << ex.what() << "\n";
            return 0;
        }
        P.inOut->logMain << "EXITING because of fatal OCM CONFIG error: " << ex.what() << "\n";
        return 1;
    }

    const string outPrefix = P.outFileNamePrefix;
    const string soloOut = outPrefix + "Solo.out";
    string rawOut;
    string filteredOut;
    string gexFeatureLabel;
    resolveGexMexDirs(P, soloOut, rawOut, filteredOut, gexFeatureLabel);

    if (!hasMexFiles(rawOut) || !hasMexFiles(filteredOut)) {
        const string err = "OCM materializer: missing GeneFull/Gene MEX at raw=" + rawOut
                           + " filtered=" + filteredOut;
        if (enableMode == "auto") {
            P.inOut->logMain << "WARNING: " << err << "\n";
            return 0;
        }
        P.inOut->logMain << "EXITING because of fatal OCM MATERIALIZATION error: " << err << "\n";
        return 1;
    }

    P.inOut->logMain << timeMonthDayTime() << " ..... started OCM multi MEX materialization\n";
    P.inOut->logMain << "  config=" << configPath << " gex_surface=" << gexFeatureLabel << "\n";

    try {
        PfMultiMerge::MexData rawData = PfMultiMerge::readMex(rawOut);
        PfMultiMerge::MexData filteredData = PfMultiMerge::readMex(filteredOut);
        validateFeatureAxesMatch(rawData, filteredData);

        uint64_t rawUnknown = 0;
        uint64_t filteredUnknown = 0;
        const map<string, vector<uint32_t>> rawTagIndices =
            buildTagColumnIndices(rawData.barcodes, rawUnknown);
        const map<string, vector<uint32_t>> filteredTagIndices =
            buildTagColumnIndices(filteredData.barcodes, filteredUnknown);
        const map<string, vector<string>> cellsPerTag = buildCellsPerTag(filteredData.barcodes);

        const string outsDir = outPrefix + "outs";
        createDirectory(outsDir + "/", P.runDirPerm, "OCM outs", P);
        const string multiRawDir = outsDir + "/multi/count/raw_feature_bc_matrix";
        const string multiMuxDir = outsDir + "/multi/multiplexing_analysis";
        createDirectory(multiMuxDir + "/", P.runDirPerm, "OCM multi multiplexing_analysis", P);

        const string filteredBarcodeGenomeLabel = resolveFilteredBarcodeGenomeLabel(config);

        if (writeMexGzDir(multiRawDir, rawData, P.inOut->logMain, P) != 0) {
            throw std::runtime_error("Failed to write outs/multi/count/raw_feature_bc_matrix");
        }
        writeCellsPerTagJson(multiMuxDir + "/cells_per_tag.json", cellsPerTag);

        const string projectRoot = resolveProjectRoot(outPrefix);
        const bool velocytoRequested = P.pSolo.featureYes[SoloFeatureTypes::Velocyto];
        VelocytoMexWriter::VelocytoMexData velocytoRaw;
        VelocytoMexWriter::VelocytoMexData velocytoFiltered;
        bool haveVelocyto = false;
        if (velocytoRequested) {
            if (VelocytoMexWriter::soloVelocytoRawReady(soloOut)) {
                const VelocytoMexWriter::VelocytoRunMex velocytoRun =
                    VelocytoMexWriter::loadSoloVelocytoMex(soloOut);
                velocytoRaw = velocytoRun.raw;
                velocytoFiltered = velocytoRun.filtered;
                validateVelocytoBarcodes(velocytoRaw.barcodes, rawData.barcodes);
                validateVelocytoBarcodes(velocytoFiltered.barcodes, filteredData.barcodes);
                haveVelocyto = true;
            } else {
                P.inOut->logMain << "WARNING: Velocyto requested but Solo.out/Velocyto/raw MEX not found; "
                                 << "skipping per-sample Velocyto mirrors\n";
            }
        }

        map<string, size_t> perSampleRawCounts;
        map<string, size_t> perSampleFilteredCounts;
        for (const auto& sample : config.samples) {
            const vector<string> ocmIds = sample.resolvedOcmIds();
            const vector<uint32_t> rawCols = unionColumnIndices(rawTagIndices, ocmIds);
            const vector<uint32_t> filteredCols = unionColumnIndices(filteredTagIndices, ocmIds);

            PfMultiMerge::MexData sampleRaw = subsetMexColumns(rawData, rawCols);
            PfMultiMerge::MexData sampleFiltered = subsetMexColumns(filteredData, filteredCols);
            perSampleRawCounts[sample.sample_id] = sampleRaw.barcodes.size();
            perSampleFilteredCounts[sample.sample_id] = sampleFiltered.barcodes.size();

            const string sampleCountDir =
                outsDir + "/per_sample_outs/" + sample.sample_id + "/count";
            const string sampleRawDir = sampleCountDir + "/sample_raw_feature_bc_matrix";
            const string sampleFilteredDir = sampleCountDir + "/sample_filtered_feature_bc_matrix";
            const string sampleFilteredCsv = sampleCountDir + "/sample_filtered_barcodes.csv";

            if (writeMexGzDir(sampleRawDir, sampleRaw, P.inOut->logMain, P) != 0) {
                throw std::runtime_error("Failed to write per-sample raw MEX for " + sample.sample_id);
            }
            if (writeMexGzDir(sampleFilteredDir, sampleFiltered, P.inOut->logMain, P) != 0) {
                throw std::runtime_error("Failed to write per-sample filtered MEX for " + sample.sample_id);
            }
            if (!sampleFiltered.barcodes.empty()) {
                writeSampleFilteredBarcodesCsv(sampleFilteredCsv, sampleFiltered.barcodes,
                                             filteredBarcodeGenomeLabel);
            } else {
                ofstream emptyCsv(sampleFilteredCsv.c_str());
                if (!emptyCsv.is_open()) {
                    throw std::runtime_error("Failed to write empty sample_filtered_barcodes.csv for "
                                             + sample.sample_id);
                }
            }

            const string downstreamOuts =
                projectRoot + "samples/" + sample.sample_id + "/run/outs";
            createDirectory(downstreamOuts + "/", P.runDirPerm, "OCM downstream outs", P);
            if (writeMexGzDir(downstreamOuts + "/raw_feature_bc_matrix", sampleRaw, P.inOut->logMain, P) != 0) {
                throw std::runtime_error("Failed downstream raw_feature_bc_matrix for " + sample.sample_id);
            }
            if (writeMexGzDir(downstreamOuts + "/filtered_feature_bc_matrix", sampleFiltered,
                              P.inOut->logMain, P) != 0) {
                throw std::runtime_error("Failed downstream filtered_feature_bc_matrix for "
                                         + sample.sample_id);
            }
            createDirectory(downstreamOuts + "/multiplexing_analysis/", P.runDirPerm,
                            "OCM downstream multiplexing_analysis", P);
            writeCellsPerTagJson(downstreamOuts + "/multiplexing_analysis/cells_per_tag.json",
                                 buildCellsPerTag(sampleFiltered.barcodes));

            if (haveVelocyto) {
                const vector<uint32_t> rawVeloCols =
                    mapGexColumnIndicesToVelocyto(rawCols, rawData.barcodes, velocytoRaw.barcodes);
                const vector<uint32_t> filteredVeloCols = mapGexColumnIndicesToVelocyto(
                    filteredCols, filteredData.barcodes, velocytoFiltered.barcodes);
                VelocytoMexWriter::VelocytoMexData sampleVeloRaw =
                    VelocytoMexWriter::subsetVelocytoColumns(velocytoRaw, rawVeloCols);
                VelocytoMexWriter::VelocytoMexData sampleVeloFiltered =
                    VelocytoMexWriter::subsetVelocytoColumns(velocytoFiltered, filteredVeloCols);
                if (VelocytoMexWriter::writeVelocytoGzDir(downstreamOuts + "/raw_velocyto_feature_bc_matrix",
                                                          sampleVeloRaw) != 0) {
                    throw std::runtime_error("Failed downstream raw_velocyto for " + sample.sample_id);
                }
                if (VelocytoMexWriter::writeVelocytoGzDir(downstreamOuts + "/filtered_velocyto_feature_bc_matrix",
                                                          sampleVeloFiltered) != 0) {
                    throw std::runtime_error("Failed downstream filtered_velocyto for " + sample.sample_id);
                }
            }

            P.inOut->logMain << "  OCM sample " << sample.sample_id << " tags="
                             << sample.ocm_barcode_ids << " raw_cells=" << sampleRaw.barcodes.size()
                             << " filtered_cells=" << sampleFiltered.barcodes.size() << "\n";
        }

        writeMaterializationSummary(multiMuxDir + "/ocm_materialization_summary.json",
                                  configPath,
                                  gexFeatureLabel,
                                  rawTagIndices,
                                  filteredTagIndices,
                                  rawUnknown,
                                  filteredUnknown,
                                  config,
                                  perSampleRawCounts,
                                  perSampleFilteredCounts);

        P.inOut->logMain << "OCM materialization summary: raw_unknown_overhangs=" << rawUnknown
                         << " filtered_unknown_overhangs=" << filteredUnknown << "\n";
        P.inOut->logMain << timeMonthDayTime() << " ..... finished OCM multi MEX materialization\n";
    } catch (const std::exception& ex) {
        P.inOut->logMain << "EXITING because of fatal OCM MATERIALIZATION error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
