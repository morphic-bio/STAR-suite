#include "CrMultiMerge.h"
#include "ErrorWarning.h"
#include "serviceFuns.cpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <sys/stat.h>
#include <stdexcept>
#include <cstdlib>
#include <cstdio>
#include <cctype>
using std::cerr;
using std::endl;

namespace CrMultiMerge {

string resolveMexFile(const string& mexDir, const string& basename) {
    string plain = mexDir + "/" + basename;
    string gz = plain + ".gz";
    
    struct stat st;
    if (stat(plain.c_str(), &st) == 0) {
        return plain;
    }
    if (stat(gz.c_str(), &st) == 0) {
        return gz;
    }
    
    ostringstream err;
    err << "Missing " << basename << "(.gz) in " << mexDir;
    throw runtime_error(err.str());
}

vector<string> readLines(const string& path) {
    vector<string> result;
    bool isGz = (path.length() > 3 && path.substr(path.length() - 3) == ".gz");
    
    if (isGz) {
        // Use zcat/gunzip via pipe for simplicity
        string cmd = "zcat \"" + path + "\" 2>/dev/null || gunzip -c \"" + path + "\" 2>/dev/null";
        FILE* pipe = popen(cmd.c_str(), "r");
        if (!pipe) {
            ostringstream err;
            err << "Failed to open gzipped file: " << path;
            throw runtime_error(err.str());
        }
        
        char buffer[8192];
        string line;
        while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            line = buffer;
            // Remove trailing newline
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
                line.pop_back();
            }
            if (!line.empty()) {
                result.push_back(line);
            }
        }
        pclose(pipe);
    } else {
        ifstream file(path);
        if (!file.is_open()) {
            ostringstream err;
            err << "Failed to open file: " << path;
            throw runtime_error(err.str());
        }
        
        string line;
        while (getline(file, line)) {
            // Remove trailing newline
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
    
    // Read features.tsv
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
    
    // Read barcodes.tsv
    string barcodesPath = resolveMexFile(mexDir, "barcodes.tsv");
    data.barcodes = readLines(barcodesPath);
    
    // Read matrix.mtx
    string matrixPath = resolveMexFile(mexDir, "matrix.mtx");
    bool isGz = (matrixPath.length() > 3 && matrixPath.substr(matrixPath.length() - 3) == ".gz");
    
    if (isGz) {
        // Use zcat/gunzip via pipe
        string cmd = "zcat \"" + matrixPath + "\" 2>/dev/null || gunzip -c \"" + matrixPath + "\" 2>/dev/null";
        FILE* pipe = popen(cmd.c_str(), "r");
        if (!pipe) {
            ostringstream err;
            err << "Failed to open matrix.mtx.gz: " << matrixPath;
            throw runtime_error(err.str());
        }
        
        char buffer[8192];
        string line;
        bool headerDone = false;
        uint32_t nrows = 0, ncols = 0;
        uint64_t nnz = 0;
        
        while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            line = buffer;
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
                line.pop_back();
            }
            
            if (line.empty()) continue;
            
            if (line[0] == '%') {
                continue; // Skip comments
            }
            
            if (!headerDone) {
                istringstream ss(line);
                ss >> nrows >> ncols >> nnz;
                headerDone = true;
                continue;
            }
            
            istringstream ss(line);
            uint32_t row, col;
            double val;
            if (ss >> row >> col >> val) {
                // Convert from 1-based to 0-based
                if (row > 0 && col > 0 && row <= data.features.size() && col <= data.barcodes.size()) {
                    MexWriter::Triplet t;
                    t.gene_idx = row - 1;
                    t.cell_idx = col - 1;
                    t.count = static_cast<uint32_t>(val);
                    data.triplets.push_back(t);
                }
            }
        }
        pclose(pipe);
    } else {
        ifstream file(matrixPath);
        if (!file.is_open()) {
            ostringstream err;
            err << "Failed to open matrix.mtx: " << matrixPath;
            throw runtime_error(err.str());
        }
        
        string line;
        bool headerDone = false;
        uint32_t nrows = 0, ncols = 0;
        uint64_t nnz = 0;
        
        while (getline(file, line)) {
            if (line.empty()) continue;
            if (line[0] == '%') continue;
            
            if (!headerDone) {
                istringstream ss(line);
                ss >> nrows >> ncols >> nnz;
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
    
    // Validate dimensions
    if (data.features.size() != data.featureNames.size() || 
        data.features.size() != data.featureTypes.size()) {
        throw runtime_error("Features array size mismatch");
    }
    
    return data;
}

MexData filterByFeatureType(const MexData& data, const string& featureType) {
    MexData filtered;
    map<uint32_t, uint32_t> rowMap; // old row -> new row
    
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
    
    // Start with GEX features and barcodes
    merged.features = gexData.features;
    merged.featureNames = gexData.featureNames;
    merged.featureTypes = gexData.featureTypes;
    merged.barcodes = gexData.barcodes;
    merged.triplets = gexData.triplets;
    
    // Create barcode map (barcode -> column index in GEX)
    map<string, uint32_t> barcodeMap;
    for (size_t i = 0; i < merged.barcodes.size(); ++i) {
        barcodeMap[merged.barcodes[i]] = i;
    }
    
    // Merge each feature MEX
    uint32_t rowOffset = merged.features.size();
    uint32_t missingCount = 0;
    
    for (const auto& featData : featureDataVec) {
        // Add features
        for (size_t i = 0; i < featData.features.size(); ++i) {
            merged.features.push_back(featData.features[i]);
            merged.featureNames.push_back(featData.featureNames[i]);
            merged.featureTypes.push_back(featData.featureTypes[i]);
        }
        
        // Merge triplets
        for (const auto& t : featData.triplets) {
            if (t.cell_idx >= featData.barcodes.size()) {
                continue; // Skip invalid indices
            }
            string bc = featData.barcodes[t.cell_idx];
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
    // Compute per-barcode counts from GEX triplets
    map<uint32_t, uint64_t> barcodeCounts;
    for (const auto& t : gexData.triplets) {
        if (t.cell_idx < gexData.barcodes.size()) {
            barcodeCounts[t.cell_idx] += t.count;
        }
    }
    
    // Return barcodes with counts > 0
    vector<string> observedBarcodes;
    observedBarcodes.reserve(barcodeCounts.size());
    for (const auto& pair : barcodeCounts) {
        if (pair.second > 0) {
            observedBarcodes.push_back(gexData.barcodes[pair.first]);
        }
    }
    
    return observedBarcodes;
}

int writeCombinedMex(const string& outputDir, const MexData& data, const string& gemWell, ofstream& logStream, const vector<string>& gexBarcodes) {
    // Create directory
    string cmd = "mkdir -p \"" + outputDir + "\"";
    int ret = system(cmd.c_str());
    if (ret != 0) {
        cerr << "ERROR: Failed to create output directory: " << outputDir << endl;
        return -1;
    }
    
    // Validate inputs
    if (data.features.empty()) {
        cerr << "ERROR: No features in MEX data" << endl;
        return -1;
    }
    
    size_t originalBarcodeCount = data.barcodes.size();
    bool useGexFilter = !gexBarcodes.empty();
    
    // Helper function to strip potential suffix from barcode for comparison
    auto stripSuffix = [](const string& bc) -> string {
        size_t dashPos = bc.find_last_of('-');
        if (dashPos != string::npos && dashPos < bc.size() - 1) {
            // Check if everything after '-' is digits
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
    
    // Build GEX barcode set (strip suffixes for comparison)
    map<string, bool> gexBarcodeSet;
    if (useGexFilter) {
        for (const auto& bc : gexBarcodes) {
            string baseBc = stripSuffix(bc);
            gexBarcodeSet[baseBc] = true;
        }
        logStream << "CR-compat MEX filtering: Using GEX barcodes only (" << gexBarcodeSet.size() << " barcodes)\n";
    } else {
        logStream << "CR-compat MEX filtering: Using observed barcodes (count > 0)\n";
    }
    
    // Step 1: Filter barcodes (GEX-only or observed only)
    map<uint32_t, uint64_t> barcodeCounts; // barcode index -> total count
    for (const auto& t : data.triplets) {
        if (t.cell_idx < data.barcodes.size()) {
            barcodeCounts[t.cell_idx] += t.count;
        }
    }
    
    // Build mapping: old_idx -> compact_idx
    map<uint32_t, uint32_t> oldToCompact;
    vector<string> filteredBarcodes;
    
    if (useGexFilter) {
        // Filter to GEX barcodes only
        filteredBarcodes.reserve(gexBarcodeSet.size());
        for (size_t i = 0; i < data.barcodes.size(); ++i) {
            string baseBc = stripSuffix(data.barcodes[i]);
            if (gexBarcodeSet.find(baseBc) != gexBarcodeSet.end()) {
                // Barcode is in GEX set, include it (even if count is 0)
                oldToCompact[i] = filteredBarcodes.size();
                filteredBarcodes.push_back(data.barcodes[i]);
            }
        }
    } else {
        // Filter to observed barcodes only (count > 0)
        filteredBarcodes.reserve(barcodeCounts.size());
        for (const auto& pair : barcodeCounts) {
            if (pair.second > 0) {
                oldToCompact[pair.first] = filteredBarcodes.size();
                filteredBarcodes.push_back(data.barcodes[pair.first]);
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
    
    // Step 2: Append GEM suffix with detection of existing -[0-9]+ pattern
    // Helper function to detect and extract existing suffix
    auto hasSuffix = [](const string& bc) -> bool {
        if (bc.size() < 2) return false;
        size_t dashPos = bc.find_last_of('-');
        if (dashPos == string::npos || dashPos == bc.size() - 1) return false;
        // Check if everything after '-' is digits
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
    suffixedBarcodes.reserve(filteredBarcodes.size());
    size_t suffixWarnings = 0;
    
    for (const auto& bc : filteredBarcodes) {
        string newBc = bc;
        if (hasSuffix(bc)) {
            string existingSuffix = extractSuffix(bc);
            if (existingSuffix != gemWell) {
                suffixWarnings++;
                // Keep existing suffix, don't append
            }
            // Already has suffix, use as-is
        } else {
            // Append GEM well suffix
            newBc += "-" + gemWell;
        }
        suffixedBarcodes.push_back(newBc);
    }
    
    if (suffixWarnings > 0) {
        logStream << "WARNING: " << suffixWarnings << " barcodes already have a different suffix than gem_well=" 
                  << gemWell << ", keeping existing suffix\n";
    }
    
    // Detect duplicate barcodes after suffixing
    map<string, size_t> barcodeDupCounts;
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
    
    // Step 3: Sort barcodes lexicographically after suffix
    vector<size_t> sortIndices(suffixedBarcodes.size());
    for (size_t i = 0; i < sortIndices.size(); ++i) {
        sortIndices[i] = i;
    }
    
    std::sort(sortIndices.begin(), sortIndices.end(), 
              [&](size_t a, size_t b) {
                  return suffixedBarcodes[a] < suffixedBarcodes[b];
              });
    
    // Build sorted barcode list and remap: compact_idx -> sorted_idx
    vector<string> sortedBarcodes;
    sortedBarcodes.reserve(suffixedBarcodes.size());
    map<uint32_t, uint32_t> compactToSorted;
    
    for (size_t i = 0; i < sortIndices.size(); ++i) {
        size_t oldCompactIdx = sortIndices[i];
        compactToSorted[oldCompactIdx] = i;
        sortedBarcodes.push_back(suffixedBarcodes[oldCompactIdx]);
    }
    
    // Step 4: Remap triplet cell_idx to sorted indices
    vector<MexWriter::Triplet> remappedTriplets;
    remappedTriplets.reserve(data.triplets.size());
    uint64_t tripletsRetained = 0;
    
    for (const auto& t : data.triplets) {
        // Map: old_idx -> compact_idx -> sorted_idx
        auto it1 = oldToCompact.find(t.cell_idx);
        if (it1 != oldToCompact.end()) {
            auto it2 = compactToSorted.find(it1->second);
            if (it2 != compactToSorted.end()) {
                MexWriter::Triplet newT;
                newT.gene_idx = t.gene_idx;
                newT.cell_idx = it2->second;
                newT.count = t.count;
                remappedTriplets.push_back(newT);
                tripletsRetained++;
            }
        }
    }
    
    // Sort triplets by (cell_idx, gene_idx) for exact parity with Cell Ranger
    std::sort(remappedTriplets.begin(), remappedTriplets.end(),
              [](const MexWriter::Triplet& a, const MexWriter::Triplet& b) {
                  if (a.cell_idx != b.cell_idx) {
                      return a.cell_idx < b.cell_idx;
                  }
                  return a.gene_idx < b.gene_idx;
              });
    
    // Step 5: Write MEX files
    vector<MexWriter::Feature> features;
    for (size_t i = 0; i < data.features.size(); ++i) {
        string name = (i < data.featureNames.size()) ? data.featureNames[i] : data.features[i];
        string type = (i < data.featureTypes.size()) ? data.featureTypes[i] : "Gene Expression";
        features.emplace_back(data.features[i], name, type);
    }
    
    string outputPrefix = outputDir + "/";
    int result = MexWriter::writeMex(outputPrefix, sortedBarcodes, features, remappedTriplets, -1);
    
    if (result != 0) {
        cerr << "ERROR: Failed to write combined MEX" << endl;
        return -1;
    }
    
    // Step 6: Gzip output files
    vector<string> filesToGzip = {"matrix.mtx", "barcodes.tsv", "features.tsv"};
    for (const auto& filename : filesToGzip) {
        string filePath = outputDir + "/" + filename;
        string gzipCmd = "gzip -f \"" + filePath + "\" 2>&1";
        FILE* gzipPipe = popen(gzipCmd.c_str(), "r");
        if (gzipPipe) {
            int gzipRet = pclose(gzipPipe);
            if (gzipRet != 0) {
                logStream << "WARNING: Failed to gzip " << filename 
                          << " (exit code " << gzipRet << "), leaving uncompressed\n";
            }
        } else {
            logStream << "WARNING: Failed to execute gzip for " << filename << ", leaving uncompressed\n";
        }
    }
    
    // Step 7: Log metrics
    logStream << "CR-compat MEX formatting:\n";
    logStream << "  Filter mode: " << (useGexFilter ? "GEX barcodes only" : "Observed barcodes (count > 0)") << "\n";
    logStream << "  Original barcode count: " << originalBarcodeCount << "\n";
    logStream << "  Filtered barcode count: " << observedBarcodeCount << "\n";
    logStream << "  Triplets retained: " << tripletsRetained << "\n";
    logStream << "  GEM well used: " << gemWell << "\n";
    
    return 0;
}

} // namespace CrMultiMerge
