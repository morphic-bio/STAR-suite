#include "OcmMultiMaterialize.h"
#include "OcmMultiConfig.h"
#include "PfMultiMerge.h"
#include "VelocytoMexWriter.h"
#include "SoloMemoryProfile.h"
#include "streamFuns.h"
#include "VERSION"
#include "SoloFeature.h"
#include "TimeFunctions.h"
#include "scrna_api.h"
#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
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

static void validateFeatureAxesMatch(const PfMultiMerge::MexAxes& raw,
                                     const PfMultiMerge::MexAxes& filtered) {
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

static bool writeGzLinesLocal(const string& path, const vector<string>& lines) {
    gzFile gz = gzopen(path.c_str(), "wb");
    if (gz == nullptr) {
        return false;
    }
    gzbuffer(gz, 1 << 20);
    gzsetparams(gz, kGzLevel, Z_DEFAULT_STRATEGY);
    for (const auto& line : lines) {
        if (gzwrite(gz, line.data(), line.size()) <= 0 ||
            gzwrite(gz, "\n", 1) <= 0) {
            gzclose(gz);
            return false;
        }
    }
    return gzclose(gz) == Z_OK;
}

static vector<string> featureLinesFromAxes(const PfMultiMerge::MexAxes& axes) {
    vector<string> lines;
    lines.reserve(axes.features.size());
    for (size_t i = 0; i < axes.features.size(); ++i) {
        const string name = (i < axes.featureNames.size()) ? axes.featureNames[i] : axes.features[i];
        const string type = (i < axes.featureTypes.size()) ? axes.featureTypes[i] : "Gene Expression";
        lines.push_back(axes.features[i] + "\t" + name + "\t" + type);
    }
    return lines;
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

static PfMultiMerge::MexAxes axesWithOcmFlexTagsStripped(const PfMultiMerge::MexAxes& axes) {
    PfMultiMerge::MexAxes out = axes;
    for (string& bc : out.barcodes) {
        bc = OcmMultiMaterialize::stripFlexTag8(bc);
    }
    return out;
}

static bool gzGetLineLocal(gzFile file, string& lineOut) {
    lineOut.clear();
    if (file == nullptr) {
        return false;
    }
    char buffer[8192];
    while (true) {
        char* ret = gzgets(file, buffer, static_cast<int>(sizeof(buffer)));
        if (ret == nullptr) {
            return !lineOut.empty();
        }
        lineOut.append(buffer);
        if (!lineOut.empty() && lineOut.back() == '\n') {
            while (!lineOut.empty() && (lineOut.back() == '\n' || lineOut.back() == '\r')) {
                lineOut.pop_back();
            }
            return true;
        }
    }
}

template <typename ShapeFn, typename EntryFn>
static void streamMtxEntries(const string& matrixPath, ShapeFn onShape, EntryFn onEntry) {
    const bool isGz = matrixPath.size() > 3 && matrixPath.substr(matrixPath.size() - 3) == ".gz";
    bool shapeSeen = false;
    auto processLine = [&](const string& line) {
        if (line.empty() || line[0] == '%') {
            return;
        }
        if (!shapeSeen) {
            uint64_t nRows = 0;
            uint64_t nCols = 0;
            uint64_t nnz = 0;
            istringstream ss(line);
            if (!(ss >> nRows >> nCols >> nnz)) {
                throw std::runtime_error("Invalid MatrixMarket dimensions in " + matrixPath);
            }
            onShape(nRows, nCols, nnz);
            shapeSeen = true;
            return;
        }
        uint32_t row = 0;
        uint32_t col = 0;
        double value = 0.0;
        istringstream ss(line);
        if (ss >> row >> col >> value) {
            onEntry(row, col, value);
        }
    };

    if (isGz) {
        gzFile gz = gzopen(matrixPath.c_str(), "rb");
        if (gz == nullptr) {
            throw std::runtime_error("Failed to open MatrixMarket file: " + matrixPath);
        }
        string line;
        while (gzGetLineLocal(gz, line)) {
            processLine(line);
        }
        const int rc = gzclose(gz);
        if (rc != Z_OK) {
            throw std::runtime_error("Failed while reading MatrixMarket file: " + matrixPath);
        }
    } else {
        ifstream in(matrixPath.c_str());
        if (!in.is_open()) {
            throw std::runtime_error("Failed to open MatrixMarket file: " + matrixPath);
        }
        string line;
        while (getline(in, line)) {
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
                line.pop_back();
            }
            processLine(line);
        }
    }
    if (!shapeSeen) {
        throw std::runtime_error("Missing MatrixMarket dimensions in " + matrixPath);
    }
}

struct RoutedMexGroup {
    string name;
    string outputDir;
    PfMultiMerge::CrBarcodeLayout layout;
};

static void writeMexMetadataForGroups(const vector<RoutedMexGroup>& groups,
                                      const vector<string>& featureLines,
                                      Parameters& P) {
    for (const auto& group : groups) {
        createDirectory(group.outputDir + "/", P.runDirPerm, "OCM routed MEX output", P);
        if (!writeGzLinesLocal(group.outputDir + "/features.tsv.gz", featureLines) ||
            !writeGzLinesLocal(group.outputDir + "/barcodes.tsv.gz", group.layout.sortedBarcodes)) {
            throw std::runtime_error("Failed writing MEX axis files for " + group.outputDir);
        }
    }
}

static void finalizeMatrixBodyToGz(const string& bodyPath,
                                   const string& outputPath,
                                   uint64_t nRows,
                                   uint64_t nCols,
                                   uint64_t nnz) {
    gzFile gz = gzopen(outputPath.c_str(), "wb");
    if (gz == nullptr) {
        throw std::runtime_error("Failed opening output matrix: " + outputPath);
    }
    gzbuffer(gz, 1 << 20);
    gzsetparams(gz, kGzLevel, Z_DEFAULT_STRATEGY);
    ostringstream header;
    header << "%%MatrixMarket matrix coordinate integer general\n"
           << "%\n"
           << nRows << " " << nCols << " " << nnz << "\n";
    const string headerStr = header.str();
    if (gzwrite(gz, headerStr.data(), headerStr.size()) <= 0) {
        gzclose(gz);
        throw std::runtime_error("Failed writing output matrix header: " + outputPath);
    }

    ifstream body(bodyPath.c_str(), ios::binary);
    if (!body.is_open()) {
        gzclose(gz);
        throw std::runtime_error("Failed opening matrix body: " + bodyPath);
    }
    vector<char> buffer(1 << 20);
    while (body.good()) {
        body.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const std::streamsize n = body.gcount();
        if (n > 0 && gzwrite(gz, buffer.data(), static_cast<unsigned int>(n)) <= 0) {
            gzclose(gz);
            throw std::runtime_error("Failed writing output matrix body: " + outputPath);
        }
    }
    if (gzclose(gz) != Z_OK) {
        throw std::runtime_error("Failed closing output matrix: " + outputPath);
    }
    std::remove(bodyPath.c_str());
}

static void streamMexMatrixToGroups(const string& inputMexDir,
                                    const string& matrixName,
                                    const vector<RoutedMexGroup>& groups,
                                    const string& tempDir,
                                    const string& label,
                                    Parameters& P) {
    if (groups.empty()) {
        return;
    }
    createDirectory(tempDir + "/", P.runDirPerm, "OCM materializer temp", P);
    const string matrixPath = PfMultiMerge::resolveMexFile(inputMexDir, matrixName);
    vector<ofstream> bodies(groups.size());
    vector<string> bodyPaths(groups.size());
    vector<uint64_t> nnz(groups.size(), 0);
    for (size_t i = 0; i < groups.size(); ++i) {
        string safeMatrix = matrixName;
        std::replace(safeMatrix.begin(), safeMatrix.end(), '/', '_');
        bodyPaths[i] = tempDir + "/" + label + "." + groups[i].name + "." + safeMatrix + ".body";
        bodies[i].open(bodyPaths[i].c_str(), ios::binary);
        if (!bodies[i].is_open()) {
            throw std::runtime_error("Failed opening temp matrix body: " + bodyPaths[i]);
        }
    }

    uint64_t nRows = 0;
    uint64_t nCols = 0;
    streamMtxEntries(
        matrixPath,
        [&](uint64_t rows, uint64_t cols, uint64_t) {
            nRows = rows;
            nCols = cols;
            for (const auto& group : groups) {
                if (group.layout.sourceColToSorted.size() != nCols) {
                    ostringstream err;
                    err << "OCM routed MEX map size mismatch for " << group.outputDir
                        << ": map=" << group.layout.sourceColToSorted.size()
                        << " matrix_cols=" << nCols;
                    throw std::runtime_error(err.str());
                }
            }
        },
        [&](uint32_t row, uint32_t col, double value) {
            if (col == 0 || col > nCols) {
                return;
            }
            const uint32_t oldCol = col - 1;
            const uint32_t intValue = value < 0.0 ? 0 : static_cast<uint32_t>(value + 0.5);
            if (intValue == 0) {
                return;
            }
            char line[96];
            for (size_t i = 0; i < groups.size(); ++i) {
                const vector<uint32_t>& remap = groups[i].layout.sourceColToSorted;
                const uint32_t newCol = remap[oldCol];
                if (newCol == UINT32_MAX) {
                    continue;
                }
                const int n = snprintf(line, sizeof(line), "%u %u %u\n", row, newCol + 1, intValue);
                if (n <= 0) {
                    throw std::runtime_error("Failed formatting routed matrix line");
                }
                bodies[i].write(line, n);
                nnz[i]++;
            }
        });
    for (auto& body : bodies) {
        body.close();
    }
    for (size_t i = 0; i < groups.size(); ++i) {
        finalizeMatrixBodyToGz(bodyPaths[i],
                               groups[i].outputDir + "/" + matrixName + ".gz",
                               nRows,
                               groups[i].layout.sortedBarcodes.size(),
                               nnz[i]);
    }
}

static void streamMexToGroups(const string& inputMexDir,
                              const PfMultiMerge::MexAxes& axes,
                              const vector<RoutedMexGroup>& groups,
                              const string& tempDir,
                              const string& label,
                              Parameters& P) {
    writeMexMetadataForGroups(groups, featureLinesFromAxes(axes), P);
    streamMexMatrixToGroups(inputMexDir, "matrix.mtx", groups, tempDir, label, P);
}

static PfMultiMerge::CrBarcodeLayout buildLayoutFromDesiredBarcodes(
    const PfMultiMerge::CrBarcodeLayout& rawLayout,
    const vector<uint32_t>& rawCols,
    const vector<string>& desiredBarcodes,
    size_t sourceColumnCount) {
    unordered_map<string, uint32_t> outputBarcodeToSource;
    outputBarcodeToSource.reserve(rawCols.size() * 2 + 1);
    for (uint32_t sourceCol : rawCols) {
        if (sourceCol >= rawLayout.sourceColToSorted.size()) {
            continue;
        }
        const uint32_t outCol = rawLayout.sourceColToSorted[sourceCol];
        if (outCol == UINT32_MAX || outCol >= rawLayout.sortedBarcodes.size()) {
            continue;
        }
        const string& barcode = rawLayout.sortedBarcodes[outCol];
        outputBarcodeToSource.emplace(barcode, sourceCol);
        outputBarcodeToSource.emplace(barcodeJoinKey(barcode), sourceCol);
    }

    PfMultiMerge::CrBarcodeLayout layout;
    layout.sortedBarcodes.reserve(desiredBarcodes.size());
    layout.sourceColToSorted.assign(sourceColumnCount, UINT32_MAX);
    for (const string& barcode : desiredBarcodes) {
        auto it = outputBarcodeToSource.find(barcode);
        if (it == outputBarcodeToSource.end()) {
            it = outputBarcodeToSource.find(barcodeJoinKey(barcode));
        }
        if (it == outputBarcodeToSource.end()) {
            throw std::runtime_error("OCM EmptyDrops barcode not found in raw sample layout: " + barcode);
        }
        const uint32_t sourceCol = it->second;
        if (layout.sourceColToSorted[sourceCol] != UINT32_MAX) {
            continue;
        }
        layout.sourceColToSorted[sourceCol] = static_cast<uint32_t>(layout.sortedBarcodes.size());
        layout.sortedBarcodes.push_back(barcode);
    }
    return layout;
}

static vector<uint32_t> columnsFromLayout(const PfMultiMerge::CrBarcodeLayout& layout) {
    vector<uint32_t> cols;
    for (size_t i = 0; i < layout.sourceColToSorted.size(); ++i) {
        if (layout.sourceColToSorted[i] != UINT32_MAX) {
            cols.push_back(static_cast<uint32_t>(i));
        }
    }
    return cols;
}

static void addLayoutCellsPerTag(const PfMultiMerge::CrBarcodeLayout& layout,
                                 const vector<string>& sourceBarcodes,
                                 map<string, vector<string>>& cellsPerTag) {
    vector<pair<uint32_t, uint32_t>> ordered;
    ordered.reserve(layout.sortedBarcodes.size());
    for (uint32_t sourceCol = 0; sourceCol < layout.sourceColToSorted.size(); ++sourceCol) {
        const uint32_t outCol = layout.sourceColToSorted[sourceCol];
        if (outCol != UINT32_MAX) {
            ordered.push_back({outCol, sourceCol});
        }
    }
    std::sort(ordered.begin(), ordered.end());
    for (const auto& item : ordered) {
        const uint32_t outCol = item.first;
        const uint32_t sourceCol = item.second;
        if (sourceCol >= sourceBarcodes.size() || outCol >= layout.sortedBarcodes.size()) {
            continue;
        }
        const string tag = OcmMultiMaterialize::classifyBarcodeTag(sourceBarcodes[sourceCol]);
        if (!tag.empty()) {
            cellsPerTag[tag].push_back(layout.sortedBarcodes[outCol]);
        }
    }
}

static void readSparseMatrixForEmptyDrops(const string& matrixPath,
                                          vector<uint32_t>& umiCounts,
                                          vector<uint32_t>& sparseGeneIds,
                                          vector<uint32_t>& sparseCounts,
                                          vector<uint32_t>& sparseCellIndex,
                                          vector<uint32_t>& nGenesPerCell,
                                          uint32_t& nRowsOut,
                                          uint32_t& nColsOut,
                                          size_t& nnzOut) {
    uint64_t nRows64 = 0;
    uint64_t nCols64 = 0;
    uint64_t declaredNnz = 0;
    streamMtxEntries(
        matrixPath,
        [&](uint64_t rows, uint64_t cols, uint64_t nnz) {
            nRows64 = rows;
            nCols64 = cols;
            declaredNnz = nnz;
            if (rows > std::numeric_limits<uint32_t>::max() ||
                cols > std::numeric_limits<uint32_t>::max()) {
                throw std::runtime_error("EmptyDrops matrix dimensions exceed uint32 range");
            }
            nGenesPerCell.assign(static_cast<size_t>(cols), 0);
        },
        [&](uint32_t, uint32_t col, double) {
            if (col > 0 && col <= nCols64) {
                nGenesPerCell[col - 1]++;
            }
        });

    nRowsOut = static_cast<uint32_t>(nRows64);
    nColsOut = static_cast<uint32_t>(nCols64);
    nnzOut = static_cast<size_t>(declaredNnz);
    sparseCellIndex.assign(static_cast<size_t>(nColsOut) + 1, 0);
    for (uint32_t i = 0; i < nColsOut; ++i) {
        sparseCellIndex[i + 1] = sparseCellIndex[i] + nGenesPerCell[i];
    }
    sparseGeneIds.assign(nnzOut, 0);
    sparseCounts.assign(nnzOut, 0);
    umiCounts.assign(nColsOut, 0);
    vector<uint32_t> offsets(nColsOut, 0);

    streamMtxEntries(
        matrixPath,
        [&](uint64_t, uint64_t, uint64_t) {},
        [&](uint32_t row, uint32_t col, double value) {
            if (row == 0 || col == 0 || row > nRowsOut || col > nColsOut) {
                return;
            }
            const uint32_t cell = col - 1;
            const uint32_t pos = sparseCellIndex[cell] + offsets[cell];
            if (pos >= sparseGeneIds.size()) {
                throw std::runtime_error("EmptyDrops sparse fill overflow");
            }
            const uint32_t count = value < 0.0 ? 0 : static_cast<uint32_t>(value + 0.5);
            sparseGeneIds[pos] = row - 1;
            sparseCounts[pos] = count;
            umiCounts[cell] += count;
            offsets[cell]++;
        });
}

static void dropZeroUmiCellsForEmptyDrops(vector<string>& barcodes,
                                          vector<uint32_t>& umiCounts,
                                          vector<uint32_t>& sparseGeneIds,
                                          vector<uint32_t>& sparseCounts,
                                          vector<uint32_t>& sparseCellIndex,
                                          vector<uint32_t>& nGenesPerCell,
                                          size_t& nnz) {
    const uint32_t oldCols = static_cast<uint32_t>(barcodes.size());
    vector<uint32_t> oldToNew(oldCols, UINT32_MAX);
    vector<string> barcodesOut;
    vector<uint32_t> umiOut;
    barcodesOut.reserve(barcodes.size());
    umiOut.reserve(umiCounts.size());
    for (uint32_t old = 0; old < oldCols; ++old) {
        if (umiCounts[old] == 0) {
            continue;
        }
        oldToNew[old] = static_cast<uint32_t>(barcodesOut.size());
        barcodesOut.push_back(barcodes[old]);
        umiOut.push_back(umiCounts[old]);
    }
    if (barcodesOut.empty()) {
        throw std::runtime_error("All barcodes have zero UMIs; cannot run OCM EmptyDrops");
    }

    vector<uint32_t> sparseGeneOut;
    vector<uint32_t> sparseCountOut;
    vector<uint32_t> sparseIndexOut(barcodesOut.size() + 1, 0);
    vector<uint32_t> nGenesOut(barcodesOut.size(), 0);
    sparseGeneOut.reserve(sparseGeneIds.size());
    sparseCountOut.reserve(sparseCounts.size());
    size_t outPos = 0;
    for (uint32_t old = 0; old < oldCols; ++old) {
        const uint32_t next = oldToNew[old];
        if (next == UINT32_MAX) {
            continue;
        }
        sparseIndexOut[next] = static_cast<uint32_t>(outPos);
        const uint32_t start = sparseCellIndex[old];
        const uint32_t nGenes = nGenesPerCell[old];
        nGenesOut[next] = nGenes;
        for (uint32_t k = 0; k < nGenes; ++k) {
            const size_t pos = static_cast<size_t>(start + k);
            sparseGeneOut.push_back(sparseGeneIds[pos]);
            sparseCountOut.push_back(sparseCounts[pos]);
            outPos++;
        }
    }
    sparseIndexOut[barcodesOut.size()] = static_cast<uint32_t>(outPos);
    nnz = outPos;
    barcodes.swap(barcodesOut);
    umiCounts.swap(umiOut);
    sparseGeneIds.swap(sparseGeneOut);
    sparseCounts.swap(sparseCountOut);
    sparseCellIndex.swap(sparseIndexOut);
    nGenesPerCell.swap(nGenesOut);
}

static void configureOcmEmptyDrops(const Parameters& P, scrna_ed_config* config) {
    config->n_expected_cells = 0;
    config->max_percentile = 0.99;
    config->max_min_ratio = 10.0;
    config->ind_min = 45000;
    config->ind_max = 90000;
    config->umi_min = 100;
    config->umi_min_frac_median = 0.01;
    config->cand_max_n = 20000;
    config->fdr = 0.01;
    config->sim_n = 100000;
    config->raw_pvalue_threshold = 0.05;
    config->seed = 1;
    config->lower_testing_bound = 500;
    config->ambient_umi_max = 100;
    config->mc_threads = 0;
    config->disable_occupancy_filter = 1;
    config->ed_retain_count = 0;
    config->use_fdr_gate = 1;
    config->apply_bh_correction = 1;
    config->use_bootstrap = P.pSolo.emptyDropsLegacyKnee ? 0 : 1;

    if (P.pSolo.cellFilter.eDcr.indMin > 0) config->ind_min = P.pSolo.cellFilter.eDcr.indMin;
    if (P.pSolo.cellFilter.eDcr.indMax > 0) config->ind_max = P.pSolo.cellFilter.eDcr.indMax;
    if (P.pSolo.cellFilter.eDcr.umiMinFracMedian > 0.0) {
        config->umi_min_frac_median = P.pSolo.cellFilter.eDcr.umiMinFracMedian;
    }
    if (P.pSolo.cellFilter.eDcr.candMaxN > 0) config->cand_max_n = P.pSolo.cellFilter.eDcr.candMaxN;
    if (P.pSolo.cellFilter.eDcr.FDR > 0.0) config->fdr = P.pSolo.cellFilter.eDcr.FDR;
    if (P.pSolo.flexFilterEdNiters > 0) config->sim_n = P.pSolo.flexFilterEdNiters;
    if (P.pSolo.flexFilterEdFdrThreshold > 0.0) config->fdr = P.pSolo.flexFilterEdFdrThreshold;
}

static vector<string> runOcmEmptyDropsOnRawMex(const string& rawMexDir,
                                               const string& edOutDir,
                                               const string& sampleId,
                                               Parameters& P) {
    const string matrixPath = PfMultiMerge::resolveMexFile(rawMexDir, "matrix.mtx");
    vector<string> barcodes = PfMultiMerge::readLines(PfMultiMerge::resolveMexFile(rawMexDir, "barcodes.tsv"));
    vector<uint32_t> umiCounts;
    vector<uint32_t> sparseGeneIds;
    vector<uint32_t> sparseCounts;
    vector<uint32_t> sparseCellIndex;
    vector<uint32_t> nGenesPerCell;
    uint32_t nRows = 0;
    uint32_t nCols = 0;
    size_t nnz = 0;
    readSparseMatrixForEmptyDrops(matrixPath,
                                  umiCounts,
                                  sparseGeneIds,
                                  sparseCounts,
                                  sparseCellIndex,
                                  nGenesPerCell,
                                  nRows,
                                  nCols,
                                  nnz);
    if (nCols != barcodes.size()) {
        throw std::runtime_error("OCM EmptyDrops barcode count mismatch for " + sampleId);
    }
    const size_t beforeCells = barcodes.size();
    dropZeroUmiCellsForEmptyDrops(barcodes,
                                  umiCounts,
                                  sparseGeneIds,
                                  sparseCounts,
                                  sparseCellIndex,
                                  nGenesPerCell,
                                  nnz);

    vector<char*> barcodePtrs;
    barcodePtrs.reserve(barcodes.size());
    for (auto& barcode : barcodes) {
        barcodePtrs.push_back(const_cast<char*>(barcode.c_str()));
    }
    scrna_matrix_input input;
    std::memset(&input, 0, sizeof(input));
    input.umi_counts = umiCounts.data();
    input.barcodes = barcodePtrs.data();
    input.features = nullptr;
    input.n_cells = static_cast<uint32_t>(barcodes.size());
    input.n_features = nRows;
    input.sparse_gene_ids = sparseGeneIds.data();
    input.sparse_counts = sparseCounts.data();
    input.sparse_cell_index = sparseCellIndex.data();
    input.n_genes_per_cell = nGenesPerCell.data();
    input.sparse_nnz = nnz;

    scrna_ed_config* config = scrna_ed_config_create();
    if (config == nullptr) {
        throw std::runtime_error("OCM EmptyDrops failed to allocate config");
    }
    configureOcmEmptyDrops(P, config);
    P.inOut->logMain << "OCM EmptyDrops sample=" << sampleId
                     << " nonzero_barcodes=" << barcodes.size()
                     << " raw_barcodes=" << beforeCells
                     << " nnz=" << nnz
                     << " simN=" << config->sim_n
                     << " bootstrap=" << (config->use_bootstrap ? "yes" : "no")
                     << "\n";

    scrna_ed_result result;
    std::memset(&result, 0, sizeof(result));
    const int rc = scrna_emptydrops_run(&input, config, &result);
    if (rc != 0) {
        string message = result.error_message ? result.error_message : "unknown error";
        scrna_ed_result_free(&result);
        scrna_ed_config_destroy(config);
        throw std::runtime_error("OCM EmptyDrops failed for " + sampleId + ": " + message);
    }

    vector<string> filtered;
    filtered.reserve(result.n_barcodes);
    for (size_t i = 0; i < result.n_barcodes; ++i) {
        if (result.barcodes[i] != nullptr) {
            filtered.push_back(result.barcodes[i]);
        }
    }
    createDirectory(edOutDir + "/", P.runDirPerm, "OCM EmptyDrops output", P);
    if (scrna_emptydrops_write_outputs(&result, edOutDir.c_str()) != 0) {
        P.inOut->logMain << "WARNING: OCM EmptyDrops detailed output failed for " << sampleId
                         << " at " << edOutDir << "\n";
    }
    scrna_ed_result_free(&result);
    scrna_ed_config_destroy(config);
    return filtered;
}

static bool optionalMexFileExists(const string& mexDir, const string& basename) {
    struct stat st;
    const string plain = mexDir + "/" + basename;
    const string gz = plain + ".gz";
    return stat(plain.c_str(), &st) == 0 || stat(gz.c_str(), &st) == 0;
}

static vector<string> availableVelocytoLayers(const string& veloRawDir) {
    vector<string> layers;
    const char* names[] = {"matrix.mtx", "spliced.mtx", "unspliced.mtx", "ambiguous.mtx"};
    for (const char* name : names) {
        if (optionalMexFileExists(veloRawDir, name)) {
            layers.push_back(name);
        }
    }
    return layers;
}

static void streamVelocytoGroups(const string& veloRawDir,
                                 const VelocytoMexWriter::VelocytoAxes& veloAxes,
                                 const vector<RoutedMexGroup>& groups,
                                 const string& tempDir,
                                 Parameters& P) {
    if (groups.empty()) {
        return;
    }
    writeMexMetadataForGroups(groups, [&]() {
        vector<string> lines;
        lines.reserve(veloAxes.features.size());
        for (size_t i = 0; i < veloAxes.features.size(); ++i) {
            const string name = (i < veloAxes.featureNames.size()) ? veloAxes.featureNames[i] : veloAxes.features[i];
            const string type = (i < veloAxes.featureTypes.size()) ? veloAxes.featureTypes[i] : "Gene Expression";
            lines.push_back(veloAxes.features[i] + "\t" + name + "\t" + type);
        }
        return lines;
    }(), P);

    vector<string> layers = availableVelocytoLayers(veloRawDir);
    const bool haveMatrix = std::find(layers.begin(), layers.end(), "matrix.mtx") != layers.end();
    vector<string> countLayers;
    for (const string& layer : layers) {
        if (layer != "matrix.mtx") {
            countLayers.push_back(layer);
        }
    }
    if (countLayers.empty()) {
        throw std::runtime_error("Velocyto raw MEX has no count layers under " + veloRawDir);
    }
    for (const string& layer : layers) {
        streamMexMatrixToGroups(veloRawDir, layer, groups, tempDir, "Velocyto", P);
    }
    if (!haveMatrix) {
        vector<ofstream> bodies(groups.size());
        vector<string> bodyPaths(groups.size());
        vector<uint64_t> nnz(groups.size(), 0);
        createDirectory(tempDir + "/", P.runDirPerm, "OCM Velocyto total temp", P);
        for (size_t i = 0; i < groups.size(); ++i) {
            bodyPaths[i] = tempDir + "/Velocyto." + groups[i].name + ".matrix_total.body";
            bodies[i].open(bodyPaths[i].c_str(), ios::binary);
            if (!bodies[i].is_open()) {
                throw std::runtime_error("Failed opening Velocyto total temp body");
            }
        }
        uint64_t nRows = 0;
        for (const string& layer : countLayers) {
            const string matrixPath = PfMultiMerge::resolveMexFile(veloRawDir, layer);
            streamMtxEntries(
                matrixPath,
                [&](uint64_t rows, uint64_t cols, uint64_t) {
                    if (nRows == 0) {
                        nRows = rows;
                    } else if (nRows != rows) {
                        throw std::runtime_error("Velocyto layer row mismatch while synthesizing matrix.mtx");
                    }
                    for (const auto& group : groups) {
                        if (group.layout.sourceColToSorted.size() != cols) {
                            throw std::runtime_error("Velocyto layout column mismatch while synthesizing matrix.mtx");
                        }
                    }
                },
                [&](uint32_t row, uint32_t col, double value) {
                    if (col == 0) {
                        return;
                    }
                    const uint32_t oldCol = col - 1;
                    const uint32_t intValue = value < 0.0 ? 0 : static_cast<uint32_t>(value + 0.5);
                    if (intValue == 0) {
                        return;
                    }
                    char line[96];
                    for (size_t i = 0; i < groups.size(); ++i) {
                        if (oldCol >= groups[i].layout.sourceColToSorted.size()) {
                            continue;
                        }
                        const uint32_t newCol = groups[i].layout.sourceColToSorted[oldCol];
                        if (newCol == UINT32_MAX) {
                            continue;
                        }
                        const int n = snprintf(line, sizeof(line), "%u %u %u\n", row, newCol + 1, intValue);
                        bodies[i].write(line, n);
                        nnz[i]++;
                    }
                });
        }
        for (auto& body : bodies) {
            body.close();
        }
        for (size_t i = 0; i < groups.size(); ++i) {
            finalizeMatrixBodyToGz(bodyPaths[i],
                                   groups[i].outputDir + "/matrix.mtx.gz",
                                   nRows,
                                   groups[i].layout.sortedBarcodes.size(),
                                   nnz[i]);
        }
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

static PfMultiMerge::CrBarcodeLayout buildVeloLayoutFromGexLayout(
    const vector<uint32_t>& gexCols,
    const vector<uint32_t>& veloCols,
    const PfMultiMerge::CrBarcodeLayout& gexLayout,
    size_t veloColumnCount) {
    PfMultiMerge::CrBarcodeLayout layout;
    layout.sortedBarcodes = gexLayout.sortedBarcodes;
    layout.sourceColToSorted.assign(veloColumnCount, UINT32_MAX);
    for (size_t i = 0; i < gexCols.size() && i < veloCols.size(); ++i) {
        const uint32_t gexCol = gexCols[i];
        const uint32_t veloCol = veloCols[i];
        if (gexCol >= gexLayout.sourceColToSorted.size() || veloCol >= layout.sourceColToSorted.size()) {
            continue;
        }
        layout.sourceColToSorted[veloCol] = gexLayout.sourceColToSorted[gexCol];
    }
    return layout;
}

static void addFilteredTagIndices(const PfMultiMerge::CrBarcodeLayout& layout,
                                  const vector<string>& sourceBarcodes,
                                  map<string, vector<uint32_t>>& filteredTagIndices) {
    for (uint32_t sourceCol = 0; sourceCol < layout.sourceColToSorted.size(); ++sourceCol) {
        if (layout.sourceColToSorted[sourceCol] == UINT32_MAX || sourceCol >= sourceBarcodes.size()) {
            continue;
        }
        const string tag = OcmMultiMaterialize::classifyBarcodeTag(sourceBarcodes[sourceCol]);
        if (!tag.empty()) {
            filteredTagIndices[tag].push_back(sourceCol);
        }
    }
}

static int runOcmMultiMaterializeNativeEmptyDrops(Parameters& P,
                                                  const string& configPath,
                                                  const PfMultiConfig::Config& config,
                                                  const string& rawOut,
                                                  const string& gexFeatureLabel,
                                                  const PfMultiMerge::MexAxes& rawAxes,
                                                  const PfMultiMerge::MexAxes& rawOutputAxes,
                                                  const map<string, vector<uint32_t>>& rawTagIndices,
                                                  uint64_t rawUnknown,
                                                  const string& soloOut) {
    P.inOut->logMain << "OCM materializer: no pool filtered MEX present; "
                     << "running native per-sample EmptyDrops from OCM raw MEX\n";
    soloMemoryProfileCheckpoint(P.inOut->logMain, "ocm_native_ed_begin");

    const string outPrefix = P.outFileNamePrefix;
    const string projectRoot = resolveProjectRoot(outPrefix);
    const string outsDir = projectRoot + "outs";
    const string multiRawDir = outsDir + "/multi/count/raw_feature_bc_matrix";
    const string multiMuxDir = outsDir + "/multi/multiplexing_analysis";
    const string tempDir = outsDir + "/multi/count/.ocm_native_tmp";
    createDirectory(outsDir + "/", P.runDirPerm, "OCM outs", P);
    createDirectory(multiMuxDir + "/", P.runDirPerm, "OCM multi multiplexing_analysis", P);
    createDirectory(tempDir + "/", P.runDirPerm, "OCM native materialization temp", P);

    if (PfMultiMerge::writeStreamedPoolMexGzCrCompat(rawOut,
                                                     multiRawDir,
                                                     rawOutputAxes,
                                                     P,
                                                     P.inOut->logMain) != 0) {
        P.inOut->logMain << "EXITING because of fatal OCM MATERIALIZATION error: "
                         << "failed writing pooled raw_feature_bc_matrix\n";
        return 1;
    }

    vector<RoutedMexGroup> rawGroups;
    vector<vector<uint32_t>> rawColsBySample;
    map<string, size_t> perSampleRawCounts;
    rawGroups.reserve(config.samples.size());
    rawColsBySample.reserve(config.samples.size());
    for (const auto& sample : config.samples) {
        const vector<uint32_t> rawCols = unionColumnIndices(rawTagIndices, sample.resolvedOcmIds());
        perSampleRawCounts[sample.sample_id] = rawCols.size();
        rawColsBySample.push_back(rawCols);

        RoutedMexGroup group;
        group.name = "raw." + sample.sample_id;
        group.outputDir = outsDir + "/per_sample_outs/" + sample.sample_id
                        + "/count/sample_raw_feature_bc_matrix";
        group.layout = PfMultiMerge::buildCrBarcodeLayoutForColumns(rawOutputAxes.barcodes,
                                                                    rawCols,
                                                                    "1",
                                                                    "TRU",
                                                                    "TRU",
                                                                    P.inOut->logMain);
        rawGroups.push_back(group);
    }

    streamMexToGroups(rawOut, rawOutputAxes, rawGroups, tempDir, "GeneFull.raw", P);
    soloMemoryProfileCheckpoint(P.inOut->logMain, "ocm_native_ed_raw_written");

    vector<RoutedMexGroup> filteredGroups;
    vector<vector<uint32_t>> filteredColsBySample;
    map<string, size_t> perSampleFilteredCounts;
    map<string, vector<uint32_t>> filteredTagIndices;
    map<string, vector<string>> globalCellsPerTag;
    const string filteredBarcodeGenomeLabel = resolveFilteredBarcodeGenomeLabel(config);

    for (size_t i = 0; i < config.samples.size(); ++i) {
        const auto& sample = config.samples[i];
        const string sampleCountDir = outsDir + "/per_sample_outs/" + sample.sample_id + "/count";
        const string edOutDir = sampleCountDir + "/emptydrops";
        vector<string> filteredBarcodes = runOcmEmptyDropsOnRawMex(rawGroups[i].outputDir,
                                                                   edOutDir,
                                                                   sample.sample_id,
                                                                   P);
        {
            ofstream filteredOutFile((edOutDir + "/filtered_barcodes.tsv").c_str());
            if (!filteredOutFile.is_open()) {
                throw std::runtime_error("Failed writing OCM filtered_barcodes.tsv for " + sample.sample_id);
            }
            for (const string& bc : filteredBarcodes) {
                filteredOutFile << bc << "\n";
            }
        }

        RoutedMexGroup filteredGroup;
        filteredGroup.name = "filtered." + sample.sample_id;
        filteredGroup.outputDir = sampleCountDir + "/sample_filtered_feature_bc_matrix";
        filteredGroup.layout = buildLayoutFromDesiredBarcodes(rawGroups[i].layout,
                                                              rawColsBySample[i],
                                                              filteredBarcodes,
                                                              rawAxes.barcodes.size());
        perSampleFilteredCounts[sample.sample_id] = filteredGroup.layout.sortedBarcodes.size();
        filteredColsBySample.push_back(columnsFromLayout(filteredGroup.layout));
        addFilteredTagIndices(filteredGroup.layout, rawAxes.barcodes, filteredTagIndices);
        addLayoutCellsPerTag(filteredGroup.layout, rawAxes.barcodes, globalCellsPerTag);

        writeSampleFilteredBarcodesCsv(sampleCountDir + "/sample_filtered_barcodes.csv",
                                       filteredGroup.layout.sortedBarcodes,
                                       filteredBarcodeGenomeLabel);
        filteredGroups.push_back(filteredGroup);
        P.inOut->logMain << "  OCM native ED sample " << sample.sample_id
                         << " raw_cells=" << rawColsBySample[i].size()
                         << " filtered_cells=" << filteredGroup.layout.sortedBarcodes.size()
                         << "\n";
    }
    soloMemoryProfileCheckpoint(P.inOut->logMain, "ocm_native_ed_cell_calls_done");

    streamMexToGroups(rawOut, rawOutputAxes, filteredGroups, tempDir, "GeneFull.filtered", P);
    writeCellsPerTagJson(multiMuxDir + "/cells_per_tag.json", globalCellsPerTag);

    for (size_t i = 0; i < config.samples.size(); ++i) {
        const auto& sample = config.samples[i];
        const string downstreamOuts = projectRoot + "samples/" + sample.sample_id + "/run/outs";
        createDirectory(downstreamOuts + "/", P.runDirPerm, "OCM downstream outs", P);
        if (PfMultiMerge::copyMexGzDir(rawGroups[i].outputDir,
                                       downstreamOuts + "/raw_feature_bc_matrix",
                                       P) != 0) {
            throw std::runtime_error("Failed downstream raw_feature_bc_matrix for " + sample.sample_id);
        }
        if (PfMultiMerge::copyMexGzDir(filteredGroups[i].outputDir,
                                       downstreamOuts + "/filtered_feature_bc_matrix",
                                       P) != 0) {
            throw std::runtime_error("Failed downstream filtered_feature_bc_matrix for " + sample.sample_id);
        }
        map<string, vector<string>> sampleCellsPerTag;
        addLayoutCellsPerTag(filteredGroups[i].layout, rawAxes.barcodes, sampleCellsPerTag);
        createDirectory(downstreamOuts + "/multiplexing_analysis/", P.runDirPerm,
                        "OCM downstream multiplexing_analysis", P);
        writeCellsPerTagJson(downstreamOuts + "/multiplexing_analysis/cells_per_tag.json",
                             sampleCellsPerTag);
    }

    const bool velocytoRequested = P.pSolo.featureYes[SoloFeatureTypes::Velocyto];
    if (velocytoRequested && VelocytoMexWriter::soloVelocytoRawReady(soloOut)) {
        VelocytoMexWriter::VelocytoRunAxes velocytoAxes;
        velocytoAxes.rawDir = soloOut + "/Velocyto/raw";
        velocytoAxes.filteredDir = soloOut + "/Velocyto/filtered";
        velocytoAxes.raw = VelocytoMexWriter::readVelocytoAxes(velocytoAxes.rawDir);
        validateVelocytoBarcodes(velocytoAxes.raw.barcodes, rawAxes.barcodes);

        vector<RoutedMexGroup> veloGroups;
        veloGroups.reserve(config.samples.size() * 2);
        for (size_t i = 0; i < config.samples.size(); ++i) {
            const auto& sample = config.samples[i];
            const string downstreamOuts = projectRoot + "samples/" + sample.sample_id + "/run/outs";
            const vector<uint32_t> rawVeloCols =
                mapGexColumnIndicesToVelocyto(rawColsBySample[i], rawAxes.barcodes, velocytoAxes.raw.barcodes);
            RoutedMexGroup rawVelo;
            rawVelo.name = "raw_velo." + sample.sample_id;
            rawVelo.outputDir = downstreamOuts + "/raw_velocyto_feature_bc_matrix";
            rawVelo.layout = buildVeloLayoutFromGexLayout(rawColsBySample[i],
                                                          rawVeloCols,
                                                          rawGroups[i].layout,
                                                          velocytoAxes.raw.barcodes.size());
            veloGroups.push_back(rawVelo);

            const vector<uint32_t> filteredVeloCols = mapGexColumnIndicesToVelocyto(
                filteredColsBySample[i], rawAxes.barcodes, velocytoAxes.raw.barcodes);
            RoutedMexGroup filteredVelo;
            filteredVelo.name = "filtered_velo." + sample.sample_id;
            filteredVelo.outputDir = downstreamOuts + "/filtered_velocyto_feature_bc_matrix";
            filteredVelo.layout = buildVeloLayoutFromGexLayout(filteredColsBySample[i],
                                                               filteredVeloCols,
                                                               filteredGroups[i].layout,
                                                               velocytoAxes.raw.barcodes.size());
            veloGroups.push_back(filteredVelo);
        }
        P.inOut->logMain << "OCM native materializer: streaming raw+filtered Velocyto for "
                         << config.samples.size() << " samples\n";
        streamVelocytoGroups(velocytoAxes.rawDir, velocytoAxes.raw, veloGroups, tempDir, P);
    } else if (velocytoRequested) {
        P.inOut->logMain << "WARNING: Velocyto requested but Solo.out/Velocyto/raw MEX not found; "
                         << "skipping per-sample Velocyto mirrors\n";
    }

    writeMaterializationSummary(multiMuxDir + "/ocm_materialization_summary.json",
                                configPath,
                                gexFeatureLabel,
                                rawTagIndices,
                                filteredTagIndices,
                                rawUnknown,
                                0,
                                config,
                                perSampleRawCounts,
                                perSampleFilteredCounts);
    soloMemoryProfileCheckpoint(P.inOut->logMain, "ocm_native_ed_done");
    P.inOut->logMain << "OCM native EmptyDrops materialization summary: raw_unknown_overhangs="
                     << rawUnknown << "\n";
    return 0;
}

} // namespace

namespace OcmMultiMaterialize {

static string stripBarcodeSuffixForClassification(const string& barcode) {
    return barcodeJoinKey(barcode);
}

bool isFlexBarcodeMode(const Parameters& P) {
    string mode = P.pfMulti.ocmMultiBarcodeMode;
    std::transform(mode.begin(), mode.end(), mode.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return mode == "flex";
}

string tag8ForOcmId(const string& ocmId) {
    if (ocmId == "OB1") {
        return "GTGTGTGT";
    }
    if (ocmId == "OB2") {
        return "CACACACA";
    }
    if (ocmId == "OB3") {
        return "TCTCTCTC";
    }
    if (ocmId == "OB4") {
        return "AGAGAGAG";
    }
    return string();
}

static string classifyBarcodeTagByOverhang(const string& core) {
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

static string classifyBarcodeTagByFlexTag8(const string& core) {
    if (core.size() < 24) {
        return string();
    }
    const string tag8 = core.substr(core.size() - 8);
    for (const char* ocmId : OcmMultiConfig::kValidOcmIds) {
        if (tag8 == tag8ForOcmId(ocmId)) {
            return ocmId;
        }
    }
    return string();
}

string classifyBarcodeTag(const string& barcode) {
    const string core = stripBarcodeSuffixForClassification(barcode);
    const string tagFromSuffix = classifyBarcodeTagByFlexTag8(core);
    if (!tagFromSuffix.empty()) {
        return tagFromSuffix;
    }
    return classifyBarcodeTagByOverhang(core);
}

string tag8ForBarcode(const string& barcode) {
    return tag8ForOcmId(classifyBarcodeTagByOverhang(stripBarcodeSuffixForClassification(barcode)));
}

string appendFlexTag8(const string& barcode) {
    const string core = stripBarcodeSuffixForClassification(barcode);
    const string tag8 = tag8ForBarcode(core);
    if (tag8.empty()) {
        return string();
    }
    return core + tag8;
}

string stripFlexTag8(const string& barcode) {
    const size_t dashPos = barcode.find_last_of('-');
    string suffix;
    string core = barcode;
    if (dashPos != string::npos && dashPos + 1 < barcode.size()) {
        bool numericSuffix = true;
        for (size_t i = dashPos + 1; i < barcode.size(); ++i) {
            if (!std::isdigit(static_cast<unsigned char>(barcode[i]))) {
                numericSuffix = false;
                break;
            }
        }
        if (numericSuffix) {
            suffix = barcode.substr(dashPos);
            core = barcode.substr(0, dashPos);
        }
    }
    if (!classifyBarcodeTagByFlexTag8(core).empty()) {
        core = core.substr(0, core.size() - 8);
    }
    return core + suffix;
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

    const bool haveFilteredMex = hasMexFiles(filteredOut);
    if (!hasMexFiles(rawOut)) {
        const string err = "OCM materializer: missing GeneFull/Gene raw MEX at raw=" + rawOut;
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
        soloMemoryProfileCheckpoint(P.inOut->logMain, "ocm_materialize_begin");
        const PfMultiMerge::MexAxes rawAxes = PfMultiMerge::readMexAxes(rawOut);
        const bool stripFlexTag8ForCrOutput = OcmMultiMaterialize::isFlexBarcodeMode(P);
        const PfMultiMerge::MexAxes rawOutputAxes =
            stripFlexTag8ForCrOutput ? axesWithOcmFlexTagsStripped(rawAxes) : rawAxes;
        uint64_t rawUnknown = 0;
        const map<string, vector<uint32_t>> rawTagIndices =
            buildTagColumnIndices(rawAxes.barcodes, rawUnknown);

        if (!haveFilteredMex) {
            const int ret = runOcmMultiMaterializeNativeEmptyDrops(P,
                                                                   configPath,
                                                                   config,
                                                                   rawOut,
                                                                   gexFeatureLabel,
                                                                   rawAxes,
                                                                   rawOutputAxes,
                                                                   rawTagIndices,
                                                                   rawUnknown,
                                                                   soloOut);
            if (ret != 0) {
                return ret;
            }
            P.inOut->logMain << timeMonthDayTime() << " ..... finished OCM multi MEX materialization\n";
            return 0;
        }

        const PfMultiMerge::MexAxes filteredAxes = PfMultiMerge::readMexAxes(filteredOut);
        validateFeatureAxesMatch(rawAxes, filteredAxes);
        const PfMultiMerge::MexAxes filteredOutputAxes =
            stripFlexTag8ForCrOutput ? axesWithOcmFlexTagsStripped(filteredAxes) : filteredAxes;
        soloMemoryProfileCheckpoint(P.inOut->logMain, "ocm_materialize_axes_loaded");

        uint64_t filteredUnknown = 0;
        const map<string, vector<uint32_t>> filteredTagIndices =
            buildTagColumnIndices(filteredAxes.barcodes, filteredUnknown);
        const map<string, vector<string>> cellsPerTag = buildCellsPerTag(filteredAxes.barcodes);

        const string projectRoot = resolveProjectRoot(outPrefix);
        const string outsDir = projectRoot + "outs";
        createDirectory(outsDir + "/", P.runDirPerm, "OCM outs", P);
        const string multiRawDir = outsDir + "/multi/count/raw_feature_bc_matrix";
        const string multiMuxDir = outsDir + "/multi/multiplexing_analysis";
        createDirectory(multiMuxDir + "/", P.runDirPerm, "OCM multi multiplexing_analysis", P);

        const string filteredBarcodeGenomeLabel = resolveFilteredBarcodeGenomeLabel(config);

        if (PfMultiMerge::writeStreamedPoolMexGzCrCompat(rawOut,
                                                         multiRawDir,
                                                         rawOutputAxes,
                                                         P,
                                                         P.inOut->logMain) != 0) {
            throw std::runtime_error("Failed to write outs/multi/count/raw_feature_bc_matrix");
        }
        writeCellsPerTagJson(multiMuxDir + "/cells_per_tag.json", cellsPerTag);
        soloMemoryProfileCheckpoint(P.inOut->logMain, "ocm_materialize_multi_raw_copied");

        const bool velocytoRequested = P.pSolo.featureYes[SoloFeatureTypes::Velocyto];
        VelocytoMexWriter::VelocytoRunAxes velocytoAxes;
        bool haveVelocyto = false;
        if (velocytoRequested) {
            if (VelocytoMexWriter::soloVelocytoRawReady(soloOut)) {
                velocytoAxes = VelocytoMexWriter::loadSoloVelocytoAxes(soloOut);
                validateVelocytoBarcodes(velocytoAxes.raw.barcodes, rawAxes.barcodes);
                validateVelocytoBarcodes(velocytoAxes.filtered.barcodes, filteredAxes.barcodes);
                haveVelocyto = true;
                soloMemoryProfileCheckpoint(P.inOut->logMain, "ocm_materialize_velocyto_axes_loaded");
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
            perSampleRawCounts[sample.sample_id] = rawCols.size();
            perSampleFilteredCounts[sample.sample_id] = filteredCols.size();

            const string sampleCountDir =
                outsDir + "/per_sample_outs/" + sample.sample_id + "/count";
            const string sampleRawDir = sampleCountDir + "/sample_raw_feature_bc_matrix";
            const string sampleFilteredDir = sampleCountDir + "/sample_filtered_feature_bc_matrix";
            const string sampleFilteredCsv = sampleCountDir + "/sample_filtered_barcodes.csv";

            soloMemoryProfileCheckpoint(P.inOut->logMain,
                                        "ocm_materialize_sample_begin:" + sample.sample_id);
            const PfMultiMerge::CrBarcodeLayout rawLayout =
                PfMultiMerge::buildCrBarcodeLayoutForColumns(rawOutputAxes.barcodes,
                                                             rawCols,
                                                             "1",
                                                             "TRU",
                                                             "TRU",
                                                             P.inOut->logMain);
            const PfMultiMerge::CrBarcodeLayout filteredLayout =
                PfMultiMerge::buildCrBarcodeLayoutForColumns(filteredOutputAxes.barcodes,
                                                             filteredCols,
                                                             "1",
                                                             "TRU",
                                                             "TRU",
                                                             P.inOut->logMain);
            if (PfMultiMerge::writeColumnSubsetMexGz(rawOut,
                                                   sampleRawDir,
                                                   rawOutputAxes,
                                                   rawLayout,
                                                   P,
                                                   P.inOut->logMain) != 0) {
                throw std::runtime_error("Failed to stream per-sample raw MEX for " + sample.sample_id);
            }
            if (PfMultiMerge::writeColumnSubsetMexGz(filteredOut,
                                                   sampleFilteredDir,
                                                   filteredOutputAxes,
                                                   filteredLayout,
                                                   P,
                                                   P.inOut->logMain) != 0) {
                throw std::runtime_error("Failed to stream per-sample filtered MEX for " + sample.sample_id);
            }

            const vector<string>& filteredBarcodeLines = filteredLayout.sortedBarcodes;
            if (!filteredBarcodeLines.empty()) {
                writeSampleFilteredBarcodesCsv(sampleFilteredCsv, filteredBarcodeLines,
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
            if (PfMultiMerge::copyMexGzDir(sampleRawDir,
                                           downstreamOuts + "/raw_feature_bc_matrix",
                                           P) != 0) {
                throw std::runtime_error("Failed downstream raw_feature_bc_matrix for " + sample.sample_id);
            }
            if (PfMultiMerge::copyMexGzDir(sampleFilteredDir,
                                           downstreamOuts + "/filtered_feature_bc_matrix",
                                           P) != 0) {
                throw std::runtime_error("Failed downstream filtered_feature_bc_matrix for "
                                         + sample.sample_id);
            }
            createDirectory(downstreamOuts + "/multiplexing_analysis/", P.runDirPerm,
                            "OCM downstream multiplexing_analysis", P);
            writeCellsPerTagJson(downstreamOuts + "/multiplexing_analysis/cells_per_tag.json",
                                 buildCellsPerTag(filteredBarcodeLines));

            if (haveVelocyto) {
                const vector<uint32_t> rawVeloCols =
                    mapGexColumnIndicesToVelocyto(rawCols, rawAxes.barcodes, velocytoAxes.raw.barcodes);
                const vector<uint32_t> filteredVeloCols = mapGexColumnIndicesToVelocyto(
                    filteredCols, filteredAxes.barcodes, velocytoAxes.raw.barcodes);
                auto buildVeloColToOut = [&](const vector<uint32_t>& gexCols,
                                             const vector<uint32_t>& veloCols,
                                             const PfMultiMerge::CrBarcodeLayout& gexLayout) {
                    vector<uint32_t> veloToOut(velocytoAxes.raw.barcodes.size(), UINT32_MAX);
                    for (size_t i = 0; i < gexCols.size() && i < veloCols.size(); ++i) {
                        const uint32_t gexIdx = gexCols[i];
                        const uint32_t veloIdx = veloCols[i];
                        if (gexIdx < gexLayout.sourceColToSorted.size() &&
                            veloIdx < veloToOut.size()) {
                            veloToOut[veloIdx] = gexLayout.sourceColToSorted[gexIdx];
                        }
                    }
                    return veloToOut;
                };
                const vector<uint32_t> rawVeloToOut =
                    buildVeloColToOut(rawCols, rawVeloCols, rawLayout);
                const vector<uint32_t> filteredVeloToOut =
                    buildVeloColToOut(filteredCols, filteredVeloCols, filteredLayout);
                const string sampleVeloRawDir = downstreamOuts + "/raw_velocyto_feature_bc_matrix";
                const string sampleVeloFilteredDir =
                    downstreamOuts + "/filtered_velocyto_feature_bc_matrix";
                if (VelocytoMexWriter::streamVelocytoColumnSubsetToDir(
                        velocytoAxes.rawDir,
                        sampleVeloRawDir,
                        rawVeloCols,
                        velocytoAxes.raw,
                        rawVeloToOut,
                        rawLayout.sortedBarcodes,
                        P) != 0) {
                    throw std::runtime_error("Failed downstream raw_velocyto for " + sample.sample_id);
                }
                if (VelocytoMexWriter::streamVelocytoColumnSubsetToDir(
                        velocytoAxes.rawDir,
                        sampleVeloFilteredDir,
                        filteredVeloCols,
                        velocytoAxes.raw,
                        filteredVeloToOut,
                        filteredLayout.sortedBarcodes,
                        P) != 0) {
                    throw std::runtime_error("Failed downstream filtered_velocyto for "
                                             + sample.sample_id);
                }
            }

            soloMemoryProfileCheckpoint(P.inOut->logMain,
                                        "ocm_materialize_sample_done:" + sample.sample_id);
            P.inOut->logMain << "  OCM sample " << sample.sample_id << " tags="
                             << sample.ocm_barcode_ids << " raw_cells=" << rawCols.size()
                             << " filtered_cells=" << filteredCols.size() << "\n";
        }
        soloMemoryProfileCheckpoint(P.inOut->logMain, "ocm_materialize_all_samples_done");

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
