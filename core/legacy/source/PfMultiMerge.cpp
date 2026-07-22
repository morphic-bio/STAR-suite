#include "PfMultiMerge.h"
#include "ErrorWarning.h"
#include "Parameters.h"
#include "streamFuns.h"
#include "serviceFuns.cpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <sys/stat.h>
#include <stdexcept>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <climits>
#include <zlib.h>
using std::cerr;
using std::endl;

namespace PfMultiMerge {

namespace {

static const int kGzLevel = 3;

static string normalizeChemistry(const string& input) {
    string out = input;
    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return out;
}

static char complementBase(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    }
}

static void translateNxtMiddleTwoBasesInplace(string& barcode) {
    if (barcode.empty()) {
        return;
    }
    size_t coreEnd = barcode.size();
    size_t dashPos = barcode.find_last_of('-');
    if (dashPos != string::npos && dashPos < barcode.size() - 1) {
        bool digits = true;
        for (size_t i = dashPos + 1; i < barcode.size(); ++i) {
            if (!std::isdigit(static_cast<unsigned char>(barcode[i]))) {
                digits = false;
                break;
            }
        }
        if (digits) {
            coreEnd = dashPos;
        }
    }

    if (coreEnd >= 9) {
        barcode[7] = complementBase(barcode[7]);
        barcode[8] = complementBase(barcode[8]);
    }
}

static bool needsNamespaceTranslation(const string& inputChemistry, const string& outputChemistry) {
    string inNorm = normalizeChemistry(inputChemistry);
    string outNorm = normalizeChemistry(outputChemistry);
    if ((inNorm != "NXT" && inNorm != "TRU") || (outNorm != "NXT" && outNorm != "TRU")) {
        return false;
    }
    return inNorm != outNorm;
}

static bool gzGetLine(gzFile file, string& lineOut) {
    lineOut.clear();
    if (file == nullptr) {
        return false;
    }

    char buffer[8192];
    while (true) {
        char* ret = gzgets(file, buffer, static_cast<int>(sizeof(buffer)));
        if (ret == nullptr) {
            if (lineOut.empty()) {
                return false;
            }
            return true;
        }

        lineOut.append(buffer);
        size_t n = lineOut.size();
        if (n > 0 && lineOut[n - 1] == '\n') {
            while (!lineOut.empty() &&
                   (lineOut.back() == '\n' || lineOut.back() == '\r')) {
                lineOut.pop_back();
            }
            return true;
        }
    }
}

static bool writeGzLines(const string& path, const vector<string>& lines) {
    gzFile gz = gzopen(path.c_str(), "wb");
    if (gz == nullptr) return false;
    gzbuffer(gz, 1 << 20);
    gzsetparams(gz, kGzLevel, Z_DEFAULT_STRATEGY);
    for (const auto& line : lines) {
        if (gzwrite(gz, line.data(), line.size()) <= 0) { gzclose(gz); return false; }
        if (gzwrite(gz, "\n", 1) <= 0) { gzclose(gz); return false; }
    }
    return gzclose(gz) == Z_OK;
}

static string stripBarcodeSuffix(const string& bc) {
    size_t dashPos = bc.find_last_of('-');
    if (dashPos != string::npos && dashPos < bc.size() - 1) {
        bool allDigits = true;
        for (size_t i = dashPos + 1; i < bc.size(); ++i) {
            if (!std::isdigit(static_cast<unsigned char>(bc[i]))) {
                allDigits = false;
                break;
            }
        }
        if (allDigits) {
            return bc.substr(0, dashPos);
        }
    }
    return bc;
}

static bool barcodeHasGemSuffix(const string& bc) {
    if (bc.size() < 2) {
        return false;
    }
    size_t dashPos = bc.find_last_of('-');
    if (dashPos == string::npos || dashPos == bc.size() - 1) {
        return false;
    }
    for (size_t i = dashPos + 1; i < bc.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(bc[i]))) {
            return false;
        }
    }
    return true;
}

static string applyGemWellSuffix(const string& bc, const string& gemWell, ostream& logStream) {
    if (barcodeHasGemSuffix(bc)) {
        const string existing = bc.substr(bc.find_last_of('-') + 1);
        if (existing != gemWell) {
            logStream << "WARNING: barcode " << bc << " suffix " << existing
                      << " differs from gem_well=" << gemWell << ", keeping existing suffix\n";
        }
        return bc;
    }
    return bc + "-" + gemWell;
}

} // namespace

CrBarcodeLayout buildCrBarcodeLayoutForColumns(const vector<string>& sourceBarcodes,
                                                const vector<uint32_t>& colIndices,
                                                const string& gemWell,
                                                const string& inputChemistry,
                                                const string& outputChemistry,
                                                ostream& logStream) {
    CrBarcodeLayout layout;
    layout.sourceColToSorted.assign(sourceBarcodes.size(), UINT32_MAX);

    vector<pair<string, uint32_t>> labeled;
    labeled.reserve(colIndices.size());
    for (uint32_t idx : colIndices) {
        if (idx >= sourceBarcodes.size()) {
            continue;
        }
        string bc = sourceBarcodes[idx];
        if (needsNamespaceTranslation(inputChemistry, outputChemistry)) {
            translateNxtMiddleTwoBasesInplace(bc);
        }
        bc = applyGemWellSuffix(bc, gemWell, logStream);
        labeled.emplace_back(bc, idx);
    }

    std::sort(labeled.begin(), labeled.end(),
              [](const pair<string, uint32_t>& a, const pair<string, uint32_t>& b) {
                  return a.first < b.first;
              });

    std::unordered_map<string, size_t> dup;
    dup.reserve(labeled.size() * 2);
    for (const auto& entry : labeled) {
        dup[entry.first]++;
    }
    for (const auto& kv : dup) {
        if (kv.second > 1) {
            ostringstream err;
            err << "ERROR: duplicate barcodes after CR-compat formatting: " << kv.first;
            throw runtime_error(err.str());
        }
    }

    layout.sortedBarcodes.reserve(labeled.size());
    for (size_t out = 0; out < labeled.size(); ++out) {
        layout.sourceColToSorted[labeled[out].second] = static_cast<uint32_t>(out);
        layout.sortedBarcodes.push_back(labeled[out].first);
    }
    return layout;
}

string resolveMexFile(const string& mexDir, const string& basename) {
    string plain = mexDir + "/" + basename;
    string gz = plain + ".gz";
    
    struct stat st;
    if (stat(gz.c_str(), &st) == 0) {
        return gz;
    }
    if (stat(plain.c_str(), &st) == 0) {
        return plain;
    }
    
    ostringstream err;
    err << "Missing " << basename << "(.gz) in " << mexDir;
    throw runtime_error(err.str());
}

vector<string> readLines(const string& path) {
    vector<string> result;
    bool isGz = (path.length() > 3 && path.substr(path.length() - 3) == ".gz");
    
    if (isGz) {
        gzFile file = gzopen(path.c_str(), "rb");
        if (file == nullptr) {
            ostringstream err;
            err << "Failed to open gzipped file: " << path;
            throw runtime_error(err.str());
        }

        string line;
        while (gzGetLine(file, line)) {
            if (!line.empty()) {
                result.push_back(line);
            }
        }
        int rc = gzclose(file);
        if (rc != Z_OK) {
            ostringstream err;
            err << "Failed while reading gzipped file: " << path;
            throw runtime_error(err.str());
        }
    } else {
        ifstream file(path);
        if (!file.is_open()) {
            ostringstream err;
            err << "Failed to open file: " << path;
            throw runtime_error(err.str());
        }
        
        string line;
        while (getline(file, line)) {
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
                line.pop_back();
            }
            if (!line.empty()) {
                result.push_back(line);
            }
        }
    }
    
    return result;
}

MexData readMex(const string& mexDir) {
    MexData data;
    
    string featuresPath = resolveMexFile(mexDir, "features.tsv");
    vector<string> featureLines = readLines(featuresPath);
    for (const auto& line : featureLines) {
        istringstream ss(line);
        string id, name, type;
        if (getline(ss, id, '\t')) {
            data.features.push_back(id);
            if (getline(ss, name, '\t')) {
                data.featureNames.push_back(name);
                if (getline(ss, type, '\t')) {
                    data.featureTypes.push_back(type);
                } else {
                    data.featureTypes.push_back("Gene Expression");
                }
            } else {
                data.featureNames.push_back(id);
                data.featureTypes.push_back("Gene Expression");
            }
        }
    }
    
    string barcodesPath = resolveMexFile(mexDir, "barcodes.tsv");
    data.barcodes = readLines(barcodesPath);
    
    string matrixPath = resolveMexFile(mexDir, "matrix.mtx");
    bool isGz = (matrixPath.length() > 3 && matrixPath.substr(matrixPath.length() - 3) == ".gz");
    uint32_t matrixNrows = 0;
    uint32_t matrixNcols = 0;
    
    if (isGz) {
        gzFile file = gzopen(matrixPath.c_str(), "rb");
        if (file == nullptr) {
            ostringstream err;
            err << "Failed to open matrix.mtx.gz: " << matrixPath;
            throw runtime_error(err.str());
        }

        string line;
        bool headerDone = false;
        uint64_t nnz = 0;
        
        while (gzGetLine(file, line)) {
            if (line.empty()) continue;
            
            if (line[0] == '%') {
                continue;
            }
            
            if (!headerDone) {
                istringstream ss(line);
                ss >> matrixNrows >> matrixNcols >> nnz;
                headerDone = true;
                continue;
            }
            
            istringstream ss(line);
            uint32_t row, col;
            double val;
            if (ss >> row >> col >> val) {
                if (row > 0 && col > 0 && row <= data.features.size() && col <= data.barcodes.size()) {
                    MexWriter::Triplet t;
                    t.gene_idx = row - 1;
                    t.cell_idx = col - 1;
                    t.count = static_cast<uint32_t>(val);
                    data.triplets.push_back(t);
                }
            }
        }
        int rc = gzclose(file);
        if (rc != Z_OK) {
            ostringstream err;
            err << "Failed while reading matrix.mtx.gz: " << matrixPath;
            throw runtime_error(err.str());
        }
    } else {
        ifstream file(matrixPath);
        if (!file.is_open()) {
            ostringstream err;
            err << "Failed to open matrix.mtx: " << matrixPath;
            throw runtime_error(err.str());
        }
        
        string line;
        bool headerDone = false;
        uint64_t nnz = 0;
        
        while (getline(file, line)) {
            if (line.empty()) continue;
            if (line[0] == '%') continue;
            
            if (!headerDone) {
                istringstream ss(line);
                ss >> matrixNrows >> matrixNcols >> nnz;
                headerDone = true;
                continue;
            }
            
            istringstream ss(line);
            uint32_t row, col;
            double val;
            if (ss >> row >> col >> val) {
                if (row > 0 && col > 0 && row <= data.features.size() && col <= data.barcodes.size()) {
                    MexWriter::Triplet t;
                    t.gene_idx = row - 1;
                    t.cell_idx = col - 1;
                    t.count = static_cast<uint32_t>(val);
                    data.triplets.push_back(t);
                }
            }
        }
    }
    
    if (data.features.size() != data.featureNames.size() || 
        data.features.size() != data.featureTypes.size()) {
        throw runtime_error("Features array size mismatch");
    }

    if (matrixNrows != data.features.size() || matrixNcols != data.barcodes.size()) {
        ostringstream err;
        err << "Matrix Market dimension mismatch in " << matrixPath
            << ": header rows=" << matrixNrows << " cols=" << matrixNcols
            << " but features.tsv has " << data.features.size()
            << " rows and barcodes.tsv has " << data.barcodes.size() << " columns";
        throw runtime_error(err.str());
    }
    
    return data;
}

MexData filterByFeatureType(const MexData& data, const string& featureType) {
    MexData filtered;
    map<uint32_t, uint32_t> rowMap;
    
    for (size_t i = 0; i < data.features.size(); ++i) {
        if (data.featureTypes[i] == featureType) {
            rowMap[i] = filtered.features.size();
            filtered.features.push_back(data.features[i]);
            filtered.featureNames.push_back(data.featureNames[i]);
            filtered.featureTypes.push_back(data.featureTypes[i]);
        }
    }
    
    filtered.barcodes = data.barcodes;
    
    for (const auto& t : data.triplets) {
        auto it = rowMap.find(t.gene_idx);
        if (it != rowMap.end()) {
            MexWriter::Triplet newT;
            newT.gene_idx = it->second;
            newT.cell_idx = t.cell_idx;
            newT.count = t.count;
            filtered.triplets.push_back(newT);
        }
    }
    
    return filtered;
}

MexData mergeMex(const MexData& gexData, const vector<MexData>& featureDataVec) {
    MexData merged;
    
    merged.features = gexData.features;
    merged.featureNames = gexData.featureNames;
    merged.featureTypes = gexData.featureTypes;
    merged.barcodes = gexData.barcodes;
    merged.triplets = gexData.triplets;
    
    std::unordered_map<string, uint32_t> barcodeMap;
    barcodeMap.reserve(merged.barcodes.size() * 2);
    for (size_t i = 0; i < merged.barcodes.size(); ++i) {
        barcodeMap[merged.barcodes[i]] = i;
    }
    
    uint32_t rowOffset = merged.features.size();
    uint64_t missingCount = 0;
    
    for (const auto& featData : featureDataVec) {
        for (size_t i = 0; i < featData.features.size(); ++i) {
            merged.features.push_back(featData.features[i]);
            merged.featureNames.push_back(featData.featureNames[i]);
            merged.featureTypes.push_back(featData.featureTypes[i]);
        }
        
        for (const auto& t : featData.triplets) {
            if (t.cell_idx >= featData.barcodes.size()) {
                continue;
            }
            const string& bc = featData.barcodes[t.cell_idx];
            auto it = barcodeMap.find(bc);
            if (it != barcodeMap.end()) {
                MexWriter::Triplet newT;
                newT.gene_idx = rowOffset + t.gene_idx;
                newT.cell_idx = it->second;
                newT.count = t.count;
                merged.triplets.push_back(newT);
            } else {
                missingCount++;
            }
        }
        
        rowOffset += featData.features.size();
    }
    
    if (missingCount > 0) {
        cerr << "WARNING: " << missingCount << " feature entries dropped (barcode not in GEX)" << endl;
    }
    
    return merged;
}

vector<string> computeObservedGexBarcodes(const MexData& gexData) {
    vector<uint64_t> barcodeCounts(gexData.barcodes.size(), 0);
    for (const auto& t : gexData.triplets) {
        if (t.cell_idx < barcodeCounts.size()) {
            barcodeCounts[t.cell_idx] += t.count;
        }
    }
    
    vector<string> observedBarcodes;
    for (size_t i = 0; i < barcodeCounts.size(); ++i) {
        if (barcodeCounts[i] > 0) {
            observedBarcodes.push_back(gexData.barcodes[i]);
        }
    }
    
    return observedBarcodes;
}

void pruneZeroCountFeatures(MexData& data, ofstream& logStream) {
    vector<bool> hasCount(data.features.size(), false);
    for (const auto& t : data.triplets) {
        if (t.gene_idx < hasCount.size()) {
            hasCount[t.gene_idx] = true;
        }
    }

    vector<uint32_t> oldToNew(data.features.size(), UINT32_MAX);
    vector<string> keptFeatures, keptNames, keptTypes;
    size_t prunedNonGex = 0;

    for (size_t i = 0; i < data.features.size(); ++i) {
        bool isGex = (data.featureTypes[i] == "Gene Expression");
        if (isGex || hasCount[i]) {
            oldToNew[i] = keptFeatures.size();
            keptFeatures.push_back(data.features[i]);
            keptNames.push_back(data.featureNames[i]);
            keptTypes.push_back(data.featureTypes[i]);
        } else {
            prunedNonGex++;
        }
    }

    if (prunedNonGex == 0) {
        logStream << "pruneZeroCountFeatures: nothing to prune (" << data.features.size() << " features all have counts)\n";
        return;
    }

    vector<MexWriter::Triplet> keptTriplets;
    keptTriplets.reserve(data.triplets.size());
    for (const auto& t : data.triplets) {
        if (t.gene_idx < oldToNew.size() && oldToNew[t.gene_idx] != UINT32_MAX) {
            MexWriter::Triplet newT;
            newT.gene_idx = oldToNew[t.gene_idx];
            newT.cell_idx = t.cell_idx;
            newT.count = t.count;
            keptTriplets.push_back(newT);
        }
    }

    logStream << "pruneZeroCountFeatures: " << data.features.size() << " -> " << keptFeatures.size()
              << " features (pruned " << prunedNonGex << " zero-count non-GEX features)\n";

    data.features = std::move(keptFeatures);
    data.featureNames = std::move(keptNames);
    data.featureTypes = std::move(keptTypes);
    data.triplets = std::move(keptTriplets);
}

size_t pruneZeroCountFeatures(MexData& data) {
    vector<bool> hasCount(data.features.size(), false);
    for (const auto& t : data.triplets) {
        if (t.gene_idx < hasCount.size()) {
            hasCount[t.gene_idx] = true;
        }
    }

    vector<uint32_t> oldToNew(data.features.size(), UINT32_MAX);
    vector<string> keptFeatures, keptNames, keptTypes;
    size_t pruned = 0;

    for (size_t i = 0; i < data.features.size(); ++i) {
        bool isGex = (data.featureTypes[i] == "Gene Expression");
        if (isGex || hasCount[i]) {
            oldToNew[i] = keptFeatures.size();
            keptFeatures.push_back(data.features[i]);
            keptNames.push_back(data.featureNames[i]);
            keptTypes.push_back(data.featureTypes[i]);
        } else {
            pruned++;
        }
    }

    if (pruned == 0) return 0;

    vector<MexWriter::Triplet> keptTriplets;
    keptTriplets.reserve(data.triplets.size());
    for (const auto& t : data.triplets) {
        if (t.gene_idx < oldToNew.size() && oldToNew[t.gene_idx] != UINT32_MAX) {
            MexWriter::Triplet newT;
            newT.gene_idx = oldToNew[t.gene_idx];
            newT.cell_idx = t.cell_idx;
            newT.count = t.count;
            keptTriplets.push_back(newT);
        }
    }

    data.features = std::move(keptFeatures);
    data.featureNames = std::move(keptNames);
    data.featureTypes = std::move(keptTypes);
    data.triplets = std::move(keptTriplets);
    return pruned;
}

int writeCombinedMex(const string& outputDir,
                     const MexData& data,
                     const string& gemWell,
                     ofstream& logStream,
                     const vector<string>& gexBarcodes,
                     const string& inputChemistry,
                     const string& outputChemistry) {
    TimerScope timer("writeCombinedMex", logStream);

    string cmd = "mkdir -p \"" + outputDir + "\"";
    int ret = system(cmd.c_str());
    if (ret != 0) {
        cerr << "ERROR: Failed to create output directory: " << outputDir << endl;
        return -1;
    }
    
    if (data.features.empty()) {
        cerr << "ERROR: No features in MEX data" << endl;
        return -1;
    }
    
    size_t originalBarcodeCount = data.barcodes.size();
    bool useGexFilter = !gexBarcodes.empty();
    
    auto stripSuffix = [](const string& bc) -> string {
        size_t dashPos = bc.find_last_of('-');
        if (dashPos != string::npos && dashPos < bc.size() - 1) {
            bool allDigits = true;
            for (size_t i = dashPos + 1; i < bc.size(); ++i) {
                if (!std::isdigit(static_cast<unsigned char>(bc[i]))) {
                    allDigits = false;
                    break;
                }
            }
            if (allDigits) {
                return bc.substr(0, dashPos);
            }
        }
        return bc;
    };
    
    std::unordered_set<string> gexBarcodeSet;
    if (useGexFilter) {
        gexBarcodeSet.reserve(gexBarcodes.size() * 2);
        for (const auto& bc : gexBarcodes) {
            gexBarcodeSet.insert(stripSuffix(bc));
        }
        logStream << "CR-compat MEX filtering: Using GEX barcodes only (" << gexBarcodeSet.size() << " barcodes)\n";
    } else {
        logStream << "CR-compat MEX filtering: Using observed barcodes (count > 0)\n";
    }
    
    // Barcode count accumulation (only needed for non-GEX-filter mode)
    vector<uint64_t> barcodeCounts;
    if (!useGexFilter) {
        barcodeCounts.resize(data.barcodes.size(), 0);
        for (const auto& t : data.triplets) {
            if (t.cell_idx < barcodeCounts.size()) {
                barcodeCounts[t.cell_idx] += t.count;
            }
        }
    }
    
    // Vector-based O(1) remap: old_idx -> compact_idx
    vector<uint32_t> oldToCompact(data.barcodes.size(), UINT32_MAX);
    vector<string> filteredBarcodes;
    
    if (useGexFilter) {
        filteredBarcodes.reserve(gexBarcodeSet.size());
        for (size_t i = 0; i < data.barcodes.size(); ++i) {
            string baseBc = stripSuffix(data.barcodes[i]);
            if (gexBarcodeSet.count(baseBc)) {
                oldToCompact[i] = filteredBarcodes.size();
                filteredBarcodes.push_back(data.barcodes[i]);
            }
        }
    } else {
        size_t nObserved = 0;
        for (size_t i = 0; i < barcodeCounts.size(); ++i) {
            if (barcodeCounts[i] > 0) nObserved++;
        }
        filteredBarcodes.reserve(nObserved);
        for (size_t i = 0; i < barcodeCounts.size(); ++i) {
            if (barcodeCounts[i] > 0) {
                oldToCompact[i] = filteredBarcodes.size();
                filteredBarcodes.push_back(data.barcodes[i]);
            }
        }
    }
    
    size_t observedBarcodeCount = filteredBarcodes.size();
    
    if (filteredBarcodes.empty()) {
        if (useGexFilter) {
            cerr << "ERROR: No barcodes match GEX barcode list" << endl;
        } else {
            cerr << "ERROR: No observed barcodes (all have zero counts)" << endl;
        }
        return -1;
    }
    
    // Namespace translation
    vector<string> namespaceAdjustedBarcodes = filteredBarcodes;
    uint64_t namespaceTranslatedCount = 0;
    if (needsNamespaceTranslation(inputChemistry, outputChemistry)) {
        for (size_t i = 0; i < namespaceAdjustedBarcodes.size(); ++i) {
            string translated = namespaceAdjustedBarcodes[i];
            translateNxtMiddleTwoBasesInplace(translated);
            if (translated != namespaceAdjustedBarcodes[i]) {
                namespaceTranslatedCount++;
            }
            namespaceAdjustedBarcodes[i] = translated;
        }
    }

    // GEM suffix handling
    auto hasSuffix = [](const string& bc) -> bool {
        if (bc.size() < 2) return false;
        size_t dashPos = bc.find_last_of('-');
        if (dashPos == string::npos || dashPos == bc.size() - 1) return false;
        for (size_t i = dashPos + 1; i < bc.size(); ++i) {
            if (!std::isdigit(static_cast<unsigned char>(bc[i]))) {
                return false;
            }
        }
        return true;
    };
    
    auto extractSuffix = [](const string& bc) -> string {
        size_t dashPos = bc.find_last_of('-');
        if (dashPos != string::npos && dashPos < bc.size() - 1) {
            return bc.substr(dashPos + 1);
        }
        return "";
    };
    
    vector<string> suffixedBarcodes;
    suffixedBarcodes.reserve(namespaceAdjustedBarcodes.size());
    size_t suffixWarnings = 0;
    
    for (const auto& bc : namespaceAdjustedBarcodes) {
        string newBc = bc;
        if (hasSuffix(bc)) {
            string existingSuffix = extractSuffix(bc);
            if (existingSuffix != gemWell) {
                suffixWarnings++;
            }
        } else {
            newBc += "-" + gemWell;
        }
        suffixedBarcodes.push_back(newBc);
    }
    
    if (suffixWarnings > 0) {
        logStream << "WARNING: " << suffixWarnings << " barcodes already have a different suffix than gem_well=" 
                  << gemWell << ", keeping existing suffix\n";
    }
    
    // Duplicate detection with unordered_map
    {
        std::unordered_map<string, size_t> barcodeDupCounts;
        barcodeDupCounts.reserve(suffixedBarcodes.size() * 2);
        for (const auto& bc : suffixedBarcodes) {
            barcodeDupCounts[bc]++;
        }
        vector<string> duplicates;
        for (const auto& pair : barcodeDupCounts) {
            if (pair.second > 1) {
                duplicates.push_back(pair.first);
            }
        }
        if (!duplicates.empty()) {
            ostringstream err;
            err << "ERROR: Duplicate barcodes after suffixing (e.g., mixed suffixed/unsuffixed input):\n";
            for (size_t i = 0; i < duplicates.size() && i < 10; ++i) {
                err << "  " << duplicates[i] << " (appears " << barcodeDupCounts[duplicates[i]] << " times)\n";
            }
            if (duplicates.size() > 10) {
                err << "  ... and " << (duplicates.size() - 10) << " more\n";
            }
            cerr << err.str();
            return -1;
        }
    }
    
    const bool emitNamespaceArtifacts = needsNamespaceTranslation(inputChemistry, outputChemistry);
    vector<string> nativeSuffixedBarcodes = suffixedBarcodes;
    if (emitNamespaceArtifacts) {
        for (auto& bc : nativeSuffixedBarcodes) {
            translateNxtMiddleTwoBasesInplace(bc);
        }
    }

    // Sort barcodes lexicographically
    vector<size_t> sortIndices(suffixedBarcodes.size());
    for (size_t i = 0; i < sortIndices.size(); ++i) {
        sortIndices[i] = i;
    }
    
    std::sort(sortIndices.begin(), sortIndices.end(), 
              [&](size_t a, size_t b) {
                  return suffixedBarcodes[a] < suffixedBarcodes[b];
              });
    
    // Vector-based O(1) compact_idx -> sorted_idx
    vector<uint32_t> compactToSorted(suffixedBarcodes.size(), UINT32_MAX);
    vector<string> sortedBarcodes;
    vector<string> sortedNativeBarcodes;
    sortedBarcodes.reserve(suffixedBarcodes.size());
    sortedNativeBarcodes.reserve(nativeSuffixedBarcodes.size());
    
    for (size_t i = 0; i < sortIndices.size(); ++i) {
        size_t oldCompactIdx = sortIndices[i];
        compactToSorted[oldCompactIdx] = i;
        sortedBarcodes.push_back(suffixedBarcodes[oldCompactIdx]);
        sortedNativeBarcodes.push_back(nativeSuffixedBarcodes[oldCompactIdx]);
    }
    
    // Fused single-pass remap: old_idx -> sorted_idx via two O(1) vector lookups
    vector<MexWriter::Triplet> remappedTriplets;
    remappedTriplets.reserve(data.triplets.size());
    uint64_t tripletsRetained = 0;
    
    for (const auto& t : data.triplets) {
        if (t.cell_idx >= oldToCompact.size()) continue;
        uint32_t compactIdx = oldToCompact[t.cell_idx];
        if (compactIdx == UINT32_MAX) continue;
        uint32_t sortedIdx = compactToSorted[compactIdx];
        if (sortedIdx == UINT32_MAX) continue;
        MexWriter::Triplet newT;
        newT.gene_idx = t.gene_idx;
        newT.cell_idx = sortedIdx;
        newT.count = t.count;
        remappedTriplets.push_back(newT);
        tripletsRetained++;
    }
    
    std::sort(remappedTriplets.begin(), remappedTriplets.end(),
              [](const MexWriter::Triplet& a, const MexWriter::Triplet& b) {
                  if (a.cell_idx != b.cell_idx) return a.cell_idx < b.cell_idx;
                  return a.gene_idx < b.gene_idx;
              });
    
    // Streaming gzip writes — write directly to .gz files

    // features.tsv.gz
    {
        string featPath = outputDir + "/features.tsv.gz";
        gzFile gz = gzopen(featPath.c_str(), "wb");
        if (gz == nullptr) {
            cerr << "ERROR: Failed to open " << featPath << endl;
            return -1;
        }
        gzbuffer(gz, 1 << 20);
        gzsetparams(gz, kGzLevel, Z_DEFAULT_STRATEGY);
        for (size_t i = 0; i < data.features.size(); ++i) {
            string name = (i < data.featureNames.size()) ? data.featureNames[i] : data.features[i];
            string type = (i < data.featureTypes.size()) ? data.featureTypes[i] : "Gene Expression";
            string line = data.features[i] + "\t" + name + "\t" + type + "\n";
            if (gzwrite(gz, line.data(), line.size()) <= 0) {
                gzclose(gz);
                cerr << "ERROR: Failed writing features.tsv.gz" << endl;
                return -1;
            }
        }
        if (gzclose(gz) != Z_OK) {
            cerr << "ERROR: Failed closing features.tsv.gz" << endl;
            return -1;
        }
    }

    // barcodes.tsv.gz
    {
        string bcPath = outputDir + "/barcodes.tsv.gz";
        if (!writeGzLines(bcPath, sortedBarcodes)) {
            cerr << "ERROR: Failed writing barcodes.tsv.gz" << endl;
            return -1;
        }
    }

    // matrix.mtx.gz — batched triplet writes
    {
        string mtxPath = outputDir + "/matrix.mtx.gz";
        gzFile gz = gzopen(mtxPath.c_str(), "wb");
        if (gz == nullptr) {
            cerr << "ERROR: Failed to open " << mtxPath << endl;
            return -1;
        }
        gzbuffer(gz, 1 << 20);
        gzsetparams(gz, kGzLevel, Z_DEFAULT_STRATEGY);

        ostringstream hdr;
        hdr << "%%MatrixMarket matrix coordinate integer general\n"
            << "%\n"
            << data.features.size() << " " << sortedBarcodes.size() << " " << remappedTriplets.size() << "\n";
        string hdrStr = hdr.str();
        if (gzwrite(gz, hdrStr.data(), hdrStr.size()) <= 0) {
            gzclose(gz);
            return -1;
        }

        const size_t kBatchSize = 4096;
        const size_t kBufCapacity = kBatchSize * 40;
        vector<char> buf(kBufCapacity);
        size_t bufUsed = 0;

        for (size_t i = 0; i < remappedTriplets.size(); ++i) {
            const auto& t = remappedTriplets[i];
            int n = snprintf(buf.data() + bufUsed, kBufCapacity - bufUsed,
                             "%u %u %u\n",
                             t.gene_idx + 1, t.cell_idx + 1, t.count);
            if (n < 0 || static_cast<size_t>(n) >= kBufCapacity - bufUsed) {
                if (bufUsed > 0) {
                    if (gzwrite(gz, buf.data(), bufUsed) <= 0) { gzclose(gz); return -1; }
                    bufUsed = 0;
                }
                n = snprintf(buf.data(), kBufCapacity, "%u %u %u\n",
                             t.gene_idx + 1, t.cell_idx + 1, t.count);
            }
            bufUsed += n;

            if (bufUsed >= kBufCapacity - 40 || ((i + 1) % kBatchSize == 0)) {
                if (gzwrite(gz, buf.data(), bufUsed) <= 0) { gzclose(gz); return -1; }
                bufUsed = 0;
            }
        }
        if (bufUsed > 0) {
            if (gzwrite(gz, buf.data(), bufUsed) <= 0) { gzclose(gz); return -1; }
        }

        if (gzclose(gz) != Z_OK) {
            cerr << "ERROR: Failed closing matrix.mtx.gz" << endl;
            return -1;
        }
    }

    // Namespace artifacts (plain-then-compress since these are small)
    if (emitNamespaceArtifacts) {
        const string nativePath = outputDir + "/barcodes.native.tsv.gz";
        if (!writeGzLines(nativePath, sortedNativeBarcodes)) {
            logStream << "WARNING: Failed to write native barcode file: " << nativePath << "\n";
        }

        const string namespaceMapPath = outputDir + "/barcodes.namespace_map.tsv.gz";
        {
            gzFile gz = gzopen(namespaceMapPath.c_str(), "wb");
            if (gz != nullptr) {
                gzbuffer(gz, 1 << 20);
                gzsetparams(gz, kGzLevel, Z_DEFAULT_STRATEGY);
                string hdrLine = "native_barcode\toutput_barcode\n";
                gzwrite(gz, hdrLine.data(), hdrLine.size());
                for (size_t i = 0; i < sortedNativeBarcodes.size(); ++i) {
                    string line = sortedNativeBarcodes[i] + "\t" + sortedBarcodes[i] + "\n";
                    gzwrite(gz, line.data(), line.size());
                }
                gzclose(gz);
            } else {
                logStream << "WARNING: Failed to write barcode namespace map file: "
                          << namespaceMapPath << "\n";
            }
        }
    }

    // Log metrics
    logStream << "CR-compat MEX formatting:\n";
    logStream << "  Filter mode: " << (useGexFilter ? "GEX barcodes only" : "Observed barcodes (count > 0)") << "\n";
    logStream << "  Original barcode count: " << originalBarcodeCount << "\n";
    logStream << "  Filtered barcode count: " << observedBarcodeCount << "\n";
    logStream << "  Barcode namespace input->output: " << normalizeChemistry(inputChemistry)
              << " -> " << normalizeChemistry(outputChemistry) << "\n";
    logStream << "  Barcode namespace translated: " << namespaceTranslatedCount << "\n";
    if (emitNamespaceArtifacts) {
        logStream << "  Added barcode file: barcodes.native.tsv.gz\n";
        logStream << "  Added barcode file: barcodes.namespace_map.tsv.gz\n";
    } else {
        logStream << "  Namespace conversion artifacts: skipped (input/output namespaces identical)\n";
    }
    logStream << "  Triplets retained: " << tripletsRetained << "\n";
    logStream << "  GEM well used: " << gemWell << "\n";
    
    return 0;
}

MexAxes readMexAxes(const string& mexDir) {
    MexAxes axes;
    const string featuresPath = resolveMexFile(mexDir, "features.tsv");
    const vector<string> featureLines = readLines(featuresPath);
    for (const auto& line : featureLines) {
        istringstream ss(line);
        string id, name, type;
        if (getline(ss, id, '\t')) {
            axes.features.push_back(id);
            if (getline(ss, name, '\t')) {
                axes.featureNames.push_back(name);
                if (getline(ss, type, '\t')) {
                    axes.featureTypes.push_back(type);
                } else {
                    axes.featureTypes.push_back("Gene Expression");
                }
            } else {
                axes.featureNames.push_back(id);
                axes.featureTypes.push_back("Gene Expression");
            }
        }
    }
    axes.barcodes = readLines(resolveMexFile(mexDir, "barcodes.tsv"));
    if (axes.features.size() != axes.featureNames.size() ||
        axes.features.size() != axes.featureTypes.size()) {
        throw runtime_error("readMexAxes: feature axis size mismatch in " + mexDir);
    }
    return axes;
}

vector<uint32_t> buildColumnRemap(const vector<uint32_t>& colIndices, size_t sourceColumnCount) {
    vector<uint32_t> oldToNew(sourceColumnCount, UINT32_MAX);
    uint32_t next = 0;
    for (uint32_t oldIdx : colIndices) {
        if (oldIdx < sourceColumnCount && oldToNew[oldIdx] == UINT32_MAX) {
            oldToNew[oldIdx] = next++;
        }
    }
    return oldToNew;
}

static uint64_t countSubsetMatrixNnz(const string& matrixPath, const vector<uint32_t>& oldToNew) {
    const bool isGz =
        (matrixPath.length() > 3 && matrixPath.substr(matrixPath.length() - 3) == ".gz");

    if (isGz) {
        gzFile file = gzopen(matrixPath.c_str(), "rb");
        if (file == nullptr) {
            throw runtime_error("countSubsetMatrixNnz: failed to open " + matrixPath);
        }
        string line;
        uint64_t localNnz = 0;
        bool headerDone = false;
        while (gzGetLine(file, line)) {
            if (line.empty()) {
                continue;
            }
            if (line[0] == '%') {
                continue;
            }
            if (!headerDone) {
                headerDone = true;
                continue;
            }
            istringstream ss(line);
            uint32_t row = 0;
            uint32_t col = 0;
            double val = 0;
            if (ss >> row >> col >> val) {
                const uint32_t oldCol = col - 1;
                if (oldCol < oldToNew.size() && oldToNew[oldCol] != UINT32_MAX) {
                    localNnz++;
                }
            }
        }
        if (gzclose(file) != Z_OK) {
            throw runtime_error("countSubsetMatrixNnz: gzclose failed for " + matrixPath);
        }
        return localNnz;
    }

    ifstream file(matrixPath.c_str());
    if (!file.is_open()) {
        throw runtime_error("countSubsetMatrixNnz: failed to open " + matrixPath);
    }
    string line;
    uint64_t localNnz = 0;
    bool headerDone = false;
    while (getline(file, line)) {
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        if (line[0] == '%') {
            continue;
        }
        if (!headerDone) {
            headerDone = true;
            continue;
        }
        istringstream ss(line);
        uint32_t row = 0;
        uint32_t col = 0;
        double val = 0;
        if (ss >> row >> col >> val) {
            const uint32_t oldCol = col - 1;
            if (oldCol < oldToNew.size() && oldToNew[oldCol] != UINT32_MAX) {
                localNnz++;
            }
        }
    }
    return localNnz;
}

uint64_t streamMatrixColumnSubset(const string& inputMatrixPath,
                                  const string& outputMatrixGzPath,
                                  const vector<uint32_t>& oldToNew,
                                  size_t nRows,
                                  size_t nColsOut) {
    const uint64_t nnz = countSubsetMatrixNnz(inputMatrixPath, oldToNew);
    const bool isGz =
        (inputMatrixPath.length() > 3 && inputMatrixPath.substr(inputMatrixPath.length() - 3) == ".gz");

    gzFile outGz = gzopen(outputMatrixGzPath.c_str(), "wb");
    if (outGz == nullptr) {
        throw runtime_error("streamMatrixColumnSubset: failed to open output " + outputMatrixGzPath);
    }
    gzbuffer(outGz, 1 << 20);
    gzsetparams(outGz, kGzLevel, Z_DEFAULT_STRATEGY);

    ostringstream hdr;
    hdr << "%%MatrixMarket matrix coordinate integer general\n%\n"
        << nRows << " " << nColsOut << " " << nnz << "\n";
    const string hdrStr = hdr.str();
    if (gzwrite(outGz, hdrStr.data(), hdrStr.size()) <= 0) {
        gzclose(outGz);
        throw runtime_error("streamMatrixColumnSubset: failed to write header");
    }

    auto writeDataLine = [&](const string& line, bool& headerDone) {
        if (line.empty() || line[0] == '%') {
            return;
        }
        if (!headerDone) {
            headerDone = true;
            return;
        }
        istringstream ss(line);
        uint32_t row = 0;
        uint32_t col = 0;
        double val = 0;
        if (ss >> row >> col >> val) {
            const uint32_t oldCol = col - 1;
            if (oldCol < oldToNew.size() && oldToNew[oldCol] != UINT32_MAX) {
                char buf[64];
                const int n = snprintf(buf, sizeof(buf), "%u %u %u\n",
                                       row,
                                       oldToNew[oldCol] + 1,
                                       static_cast<uint32_t>(val));
                if (n <= 0 || gzwrite(outGz, buf, n) <= 0) {
                    throw runtime_error("streamMatrixColumnSubset: write failed");
                }
            }
        }
    };

    if (isGz) {
        gzFile inGz = gzopen(inputMatrixPath.c_str(), "rb");
        if (inGz == nullptr) {
            gzclose(outGz);
            throw runtime_error("streamMatrixColumnSubset: failed to open input " + inputMatrixPath);
        }
        string line;
        bool headerDone = false;
        try {
            while (gzGetLine(inGz, line)) {
                writeDataLine(line, headerDone);
            }
        } catch (...) {
            gzclose(inGz);
            gzclose(outGz);
            throw;
        }
        if (gzclose(inGz) != Z_OK) {
            gzclose(outGz);
            throw runtime_error("streamMatrixColumnSubset: input gzclose failed");
        }
    } else {
        ifstream inFile(inputMatrixPath.c_str());
        if (!inFile.is_open()) {
            gzclose(outGz);
            throw runtime_error("streamMatrixColumnSubset: failed to open input " + inputMatrixPath);
        }
        string line;
        bool headerDone = false;
        while (getline(inFile, line)) {
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
                line.pop_back();
            }
            writeDataLine(line, headerDone);
        }
    }

    if (gzclose(outGz) != Z_OK) {
        throw runtime_error("streamMatrixColumnSubset: output gzclose failed");
    }
    return nnz;
}

static vector<string> featureAxisLines(const MexAxes& axes) {
    vector<string> lines;
    lines.reserve(axes.features.size());
    for (size_t i = 0; i < axes.features.size(); ++i) {
        const string name = (i < axes.featureNames.size()) ? axes.featureNames[i] : axes.features[i];
        const string type = (i < axes.featureTypes.size()) ? axes.featureTypes[i] : "Gene Expression";
        lines.push_back(axes.features[i] + "\t" + name + "\t" + type);
    }
    return lines;
}

int writeColumnSubsetMexGz(const string& inputMexDir,
                           const string& outputDir,
                           const MexAxes& sourceAxes,
                           const CrBarcodeLayout& layout,
                           Parameters& P,
                           ostream& logStream) {
    (void)logStream;
    createDirectory(outputDir + "/", P.runDirPerm, "MEX output", P);

    const vector<string>& sortedBarcodes = layout.sortedBarcodes;
    if (!writeGzLines(outputDir + "/features.tsv.gz", featureAxisLines(sourceAxes))) {
        return -1;
    }
    if (!writeGzLines(outputDir + "/barcodes.tsv.gz", sortedBarcodes)) {
        return -1;
    }
    if (sortedBarcodes.empty()) {
        gzFile gz = gzopen((outputDir + "/matrix.mtx.gz").c_str(), "wb");
        if (gz == nullptr) {
            return -1;
        }
        gzbuffer(gz, 1 << 20);
        gzsetparams(gz, kGzLevel, Z_DEFAULT_STRATEGY);
        ostringstream hdr;
        hdr << "%%MatrixMarket matrix coordinate integer general\n%\n"
            << sourceAxes.features.size() << " 0 0\n";
        const string hdrStr = hdr.str();
        if (gzwrite(gz, hdrStr.data(), hdrStr.size()) <= 0 || gzclose(gz) != Z_OK) {
            return -1;
        }
        return 0;
    }

    const string matrixIn = resolveMexFile(inputMexDir, "matrix.mtx");
    try {
        streamMatrixColumnSubset(matrixIn,
                                 outputDir + "/matrix.mtx.gz",
                                 layout.sourceColToSorted,
                                 sourceAxes.features.size(),
                                 sortedBarcodes.size());
    } catch (const exception& ex) {
        cerr << "writeColumnSubsetMexGz: " << ex.what() << endl;
        return -1;
    }
    return 0;
}

int writeStreamedPoolMexGzCrCompat(const string& inputMexDir,
                                   const string& outputDir,
                                   const MexAxes& sourceAxes,
                                   Parameters& P,
                                   ostream& logStream,
                                   const string& gemWell,
                                   const vector<string>& gexBarcodeFilter,
                                   const string& inputChemistry,
                                   const string& outputChemistry) {
    createDirectory(outputDir + "/", P.runDirPerm, "MEX output", P);

    vector<uint32_t> colIndices;

    if (!gexBarcodeFilter.empty()) {
        std::unordered_set<string> gexSet;
        gexSet.reserve(gexBarcodeFilter.size() * 2);
        for (const auto& bc : gexBarcodeFilter) {
            gexSet.insert(stripBarcodeSuffix(bc));
        }
        colIndices.reserve(gexSet.size());
        for (size_t idx = 0; idx < sourceAxes.barcodes.size(); ++idx) {
            if (gexSet.count(stripBarcodeSuffix(sourceAxes.barcodes[idx]))) {
                colIndices.push_back(static_cast<uint32_t>(idx));
            }
        }
        logStream << "CR-compat pool MEX: GEX barcode filter (" << colIndices.size()
                  << " columns)\n";
    } else {
        colIndices.reserve(sourceAxes.barcodes.size());
        for (size_t idx = 0; idx < sourceAxes.barcodes.size(); ++idx) {
            colIndices.push_back(static_cast<uint32_t>(idx));
        }
        logStream << "CR-compat pool MEX: full raw barcode axis (" << colIndices.size()
                  << " columns)\n";
    }

    if (colIndices.empty()) {
        cerr << "ERROR: No columns for pool MEX output" << endl;
        return -1;
    }

    const CrBarcodeLayout layout = buildCrBarcodeLayoutForColumns(sourceAxes.barcodes,
                                                                  colIndices,
                                                                  gemWell,
                                                                  inputChemistry,
                                                                  outputChemistry,
                                                                  logStream);
    return writeColumnSubsetMexGz(inputMexDir, outputDir, sourceAxes, layout, P, logStream);
}

static bool copyBinaryFile(const string& src, const string& dst) {
    ifstream in(src.c_str(), ios::binary);
    ofstream out(dst.c_str(), ios::binary);
    if (!in.is_open() || !out.is_open()) {
        return false;
    }
    out << in.rdbuf();
    return in.good() || in.eof();
}

int copyMexGzDir(const string& inputMexDir, const string& outputDir, Parameters& P) {
    createDirectory(outputDir + "/", P.runDirPerm, "MEX output", P);
    const char* basenames[] = {"features.tsv", "barcodes.tsv", "matrix.mtx"};
    for (const char* basename : basenames) {
        const string src = resolveMexFile(inputMexDir, basename);
        const string dst = outputDir + "/" + string(basename) + ".gz";
        const bool srcGz = (src.length() > 3 && src.substr(src.length() - 3) == ".gz");
        if (srcGz) {
            if (!copyBinaryFile(src, dst)) {
                cerr << "copyMexGzDir: failed to copy " << src << " -> " << dst << endl;
                return -1;
            }
            continue;
        }
        if (string(basename).find(".tsv") != string::npos) {
            if (!writeGzLines(dst, readLines(src))) {
                return -1;
            }
            continue;
        }
        gzFile outGz = gzopen(dst.c_str(), "wb");
        if (outGz == nullptr) {
            return -1;
        }
        gzbuffer(outGz, 1 << 20);
        gzsetparams(outGz, kGzLevel, Z_DEFAULT_STRATEGY);
        ifstream inFile(src.c_str());
        if (!inFile.is_open()) {
            gzclose(outGz);
            return -1;
        }
        string line;
        while (getline(inFile, line)) {
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            if (gzwrite(outGz, line.data(), line.size()) <= 0 || gzwrite(outGz, "\n", 1) <= 0) {
                gzclose(outGz);
                return -1;
            }
        }
        if (gzclose(outGz) != Z_OK) {
            return -1;
        }
    }
    return 0;
}

} // namespace PfMultiMerge
