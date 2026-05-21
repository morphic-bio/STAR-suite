#include "PfMultiConfig.h"
#include "ErrorWarning.h"
#include "serviceFuns.cpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <sys/stat.h>
#include <set>
#include <stdexcept>
#include <cstdlib>

namespace PfMultiConfig {

static void trimInPlace(string& s) {
    size_t first = s.find_first_not_of(" \t\r\n");
    if (first == string::npos) {
        s.clear();
        return;
    }
    size_t last = s.find_last_not_of(" \t\r\n");
    s = s.substr(first, last - first + 1);
}

static vector<string> parseCsvFields(const string& row) {
    vector<string> fields;
    string field;
    bool inQuotes = false;
    for (char c : row) {
        if (c == '"') {
            inQuotes = !inQuotes;
        } else if (c == ',' && !inQuotes) {
            fields.push_back(field);
            field.clear();
        } else {
            field += c;
        }
    }
    if (!field.empty() || (!row.empty() && row.back() == ',')) {
        fields.push_back(field);
    }
    for (auto& value : fields) {
        trimInPlace(value);
        if (!value.empty() && value.front() == '"' && value.back() == '"') {
            value = value.substr(1, value.length() - 2);
            trimInPlace(value);
        }
    }
    return fields;
}

static void parseSamplesRows(const vector<string>& samplesHeader,
                             const vector<string>& rows,
                             vector<SampleEntry>& samplesOut) {
    if (samplesHeader.empty()) {
        throw runtime_error("samples section has no header");
    }
    for (const auto& row : rows) {
        if (row.empty()) {
            continue;
        }
        vector<string> fields = parseCsvFields(row);
        while (fields.size() < samplesHeader.size()) {
            fields.push_back(string());
        }
        SampleEntry entry;
        for (size_t i = 0; i < samplesHeader.size() && i < fields.size(); ++i) {
            string header = samplesHeader[i];
            std::transform(header.begin(), header.end(), header.begin(), ::tolower);
            const string& value = fields[i];
            if (header == "sample_id") {
                entry.sample_id = value;
            } else if (header == "ocm_barcode_ids") {
                entry.ocm_barcode_ids = value;
            } else if (header == "description") {
                entry.description = value;
            }
        }
        if (!entry.sample_id.empty() || !entry.ocm_barcode_ids.empty()) {
            samplesOut.push_back(entry);
        }
    }
}

vector<string> SampleEntry::resolvedOcmIds() const {
    vector<string> ids;
    string field = ocm_barcode_ids;
    trimInPlace(field);
    if (field.empty()) {
        return ids;
    }
    string token;
    for (char c : field) {
        if (c == '|') {
            trimInPlace(token);
            if (!token.empty()) {
                ids.push_back(token);
            }
            token.clear();
        } else {
            token.push_back(c);
        }
    }
    trimInPlace(token);
    if (!token.empty()) {
        ids.push_back(token);
    }
    return ids;
}

string LibraryEntry::normalizedFeatureType() const {
    string normalized;
    normalized.reserve(feature_types.size());
    for (unsigned char c : feature_types) {
        if (std::isalnum(c)) {
            normalized.push_back(static_cast<char>(std::tolower(c)));
        }
    }
    return normalized;
}

vector<LibraryEntry> Config::getGexLibraries() const {
    vector<LibraryEntry> result;
    for (const auto& lib : libraries) {
        string norm = lib.normalizedFeatureType();
        if (norm == "geneexpression" || norm == "gex") {
            result.push_back(lib);
        }
    }
    return result;
}

vector<LibraryEntry> Config::getFeatureLibraries(const string& featureType) const {
    vector<LibraryEntry> result;
    // Normalize the target the same way as normalizedFeatureType():
    // strip all non-alphanumeric, lowercase.
    string targetNorm;
    targetNorm.reserve(featureType.size());
    for (unsigned char c : featureType) {
        if (std::isalnum(c)) {
            targetNorm.push_back(static_cast<char>(std::tolower(c)));
        }
    }

    for (const auto& lib : libraries) {
        if (lib.normalizedFeatureType() == targetNorm) {
            result.push_back(lib);
        }
    }
    return result;
}

Config parseConfig(const string& configPath) {
    Config config;
    ifstream file(configPath);
    if (!file.is_open()) {
        ostringstream err;
        err << "Failed to open multi config file: " << configPath;
        throw runtime_error(err.str());
    }
    
    string currentSection;
    vector<string> currentLines;
    bool inLibraries = false;
    bool inSamples = false;
    vector<string> librariesHeader;
    vector<string> samplesHeader;
    
    string line;
    while (getline(file, line)) {
        // Remove comments and whitespace
        size_t commentPos = line.find_first_of("#;");
        if (commentPos != string::npos) {
            line = line.substr(0, commentPos);
        }
        
        // Trim whitespace safely
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == string::npos) {
            continue;
        }
        size_t last = line.find_last_not_of(" \t\r\n");
        line = line.substr(first, last - first + 1);
        
        if (line.empty()) {
            continue;
        }
        
        // Check for section header
        if (line.front() == '[' && line.back() == ']') {
            // Process previous section
            if (currentSection == "libraries" && inLibraries) {
                // Parse CSV
                if (librariesHeader.empty()) {
                    ostringstream err;
                    err << "libraries section has no header";
                    throw runtime_error(err.str());
                }
                
                for (const auto& row : currentLines) {
                    if (row.empty()) continue;
                    istringstream ss(row);
                    vector<string> fields;
                    string field;
                    bool inQuotes = false;
                    
                    // Simple CSV parsing (handles quoted fields)
                    for (char c : row) {
                        if (c == '"') {
                            inQuotes = !inQuotes;
                        } else if (c == ',' && !inQuotes) {
                            fields.push_back(field);
                            field.clear();
                        } else {
                            field += c;
                        }
                    }
                    if (!field.empty() || row.back() == ',') {
                        fields.push_back(field);
                    }
                    
                    // Pad missing trailing fields with empty strings
                    while (fields.size() < librariesHeader.size()) {
                        fields.push_back(string());
                    }
                    
                    LibraryEntry entry;
                    entry.gem_well = "1"; // Default to "1" if not specified
                    for (size_t i = 0; i < librariesHeader.size() && i < fields.size(); ++i) {
                        string header = librariesHeader[i];
                        std::transform(header.begin(), header.end(), header.begin(), ::tolower);
                        string value = fields[i];
                        // Remove quotes if present
                        if (!value.empty() && value.front() == '"' && value.back() == '"') {
                            value = value.substr(1, value.length() - 2);
                        }
                        trimInPlace(value);
                        
                        if (header == "fastqs") {
                            entry.fastqs = value;
                        } else if (header == "feature_types" || header == "feature_type") {
                            entry.feature_types = value;
                        } else if (header == "sample") {
                            entry.sample = value;
                        } else if (header == "library_type" || header == "librarytype") {
                            entry.library_type = value;
                        } else if (header == "gem_well") {
                            entry.gem_well = value.empty() ? "1" : value;
                        } else if (header == "star_chemistry" || header == "starchemistry") {
                            if (!value.empty()) {
                                string lower = value;
                                std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                                if (lower != "tru" && lower != "nxt" && lower != "auto") {
                                    throw runtime_error("Invalid star_chemistry value '" + value
                                        + "'; must be TRU, NXT, auto, or empty");
                                }
                                entry.starChemistry = lower;
                            }
                        } else if (header == "star_whitelist" || header == "starwhitelist") {
                            entry.starWhitelist = value;
                        } else if (header == "star_feature_ref" || header == "starfeatureref") {
                            entry.starFeatureRef = value;
                        } else if (header == "star_library_id" || header == "starlibraryid") {
                            entry.starLibraryId = value;
                        } else if (header == "star_max_hamming" || header == "starmaxhamming") {
                            if (!value.empty()) {
                                entry.starMaxHamming = std::atoi(value.c_str());
                            }
                        }
                    }
                    if (!entry.fastqs.empty()) {
                        config.libraries.push_back(entry);
                    }
                }
            } else if (currentSection == "feature") {
                // Parse key=value pairs
                for (const auto& kvLine : currentLines) {
                    size_t commaPos = kvLine.find(',');
                    if (commaPos == string::npos) continue;
                    string key = kvLine.substr(0, commaPos);
                    string value = kvLine.substr(commaPos + 1);
                    // Trim
                    trimInPlace(key);
                    trimInPlace(value);
                    
                    string keyLower = key;
                    std::transform(keyLower.begin(), keyLower.end(), keyLower.begin(), ::tolower);
                    if (keyLower == "ref" || keyLower == "reference") {
                        config.featureRef = value;
                    }
                }
            } else if (currentSection == "reference") {
                // Parse key=value pairs
                for (const auto& kvLine : currentLines) {
                    size_t commaPos = kvLine.find(',');
                    if (commaPos == string::npos) continue;
                    string key = kvLine.substr(0, commaPos);
                    string value = kvLine.substr(commaPos + 1);
                    // Trim
                    trimInPlace(key);
                    trimInPlace(value);
                    
                    string keyLower = key;
                    std::transform(keyLower.begin(), keyLower.end(), keyLower.begin(), ::tolower);
                    if (keyLower == "path" || keyLower == "genome") {
                        config.referencePath = value;
                    }
                }
            } else if (currentSection == "samples" && inSamples) {
                parseSamplesRows(samplesHeader, currentLines, config.samples);
            }
            
            // Start new section
            currentSection = line.substr(1, line.length() - 2);
            std::transform(currentSection.begin(), currentSection.end(), currentSection.begin(), ::tolower);
            currentLines.clear();
            inLibraries = false;
            inSamples = false;
            samplesHeader.clear();
            
            if (currentSection == "libraries") {
                inLibraries = true;
                // Read header
                if (!getline(file, line)) {
                    break;
                }
                // Remove comments
                size_t commentPos = line.find_first_of("#;");
                if (commentPos != string::npos) {
                    line = line.substr(0, commentPos);
                }
                // Parse header
                istringstream headerStream(line);
                string headerField;
                while (getline(headerStream, headerField, ',')) {
                    headerField.erase(0, headerField.find_first_not_of(" \t"));
                    headerField.erase(headerField.find_last_not_of(" \t") + 1);
                    std::transform(headerField.begin(), headerField.end(), headerField.begin(), ::tolower);
                    librariesHeader.push_back(headerField);
                }
            } else if (currentSection == "samples") {
                inSamples = true;
                if (!getline(file, line)) {
                    break;
                }
                size_t commentPosSamples = line.find_first_of("#;");
                if (commentPosSamples != string::npos) {
                    line = line.substr(0, commentPosSamples);
                }
                trimInPlace(line);
                istringstream headerStream(line);
                string headerField;
                while (getline(headerStream, headerField, ',')) {
                    trimInPlace(headerField);
                    std::transform(headerField.begin(), headerField.end(), headerField.begin(), ::tolower);
                    samplesHeader.push_back(headerField);
                }
            }
        } else {
            currentLines.push_back(line);
        }
    }
    
    // Process final section
    if (currentSection == "libraries" && inLibraries) {
        // Same parsing as above
        if (!librariesHeader.empty()) {
            for (const auto& row : currentLines) {
                if (row.empty()) continue;
                istringstream ss(row);
                vector<string> fields;
                string field;
                bool inQuotes = false;
                
                for (char c : row) {
                    if (c == '"') {
                        inQuotes = !inQuotes;
                    } else if (c == ',' && !inQuotes) {
                        fields.push_back(field);
                        field.clear();
                    } else {
                        field += c;
                    }
                }
                if (!field.empty() || row.back() == ',') {
                    fields.push_back(field);
                }
                
                // Pad missing trailing fields with empty strings
                while (fields.size() < librariesHeader.size()) {
                    fields.push_back(string());
                }
                
                LibraryEntry entry;
                entry.gem_well = "1"; // Default to "1" if not specified
                for (size_t i = 0; i < librariesHeader.size() && i < fields.size(); ++i) {
                    string header = librariesHeader[i];
                    std::transform(header.begin(), header.end(), header.begin(), ::tolower);
                    string value = fields[i];
                    if (!value.empty() && value.front() == '"' && value.back() == '"') {
                        value = value.substr(1, value.length() - 2);
                    }
                    trimInPlace(value);
                    
                    if (header == "fastqs") {
                        entry.fastqs = value;
                    } else if (header == "feature_types" || header == "feature_type") {
                        entry.feature_types = value;
                    } else if (header == "sample") {
                        entry.sample = value;
                    } else if (header == "library_type" || header == "librarytype") {
                        entry.library_type = value;
                    } else if (header == "gem_well") {
                        entry.gem_well = value.empty() ? "1" : value;
                    } else if (header == "star_chemistry" || header == "starchemistry") {
                        if (!value.empty()) {
                            string lower = value;
                            std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                            if (lower != "tru" && lower != "nxt" && lower != "auto") {
                                throw runtime_error("Invalid star_chemistry value '" + value
                                    + "'; must be TRU, NXT, auto, or empty");
                            }
                            entry.starChemistry = lower;
                        }
                    } else if (header == "star_whitelist" || header == "starwhitelist") {
                        entry.starWhitelist = value;
                    } else if (header == "star_feature_ref" || header == "starfeatureref") {
                        entry.starFeatureRef = value;
                    } else if (header == "star_library_id" || header == "starlibraryid") {
                        entry.starLibraryId = value;
                    } else if (header == "star_max_hamming" || header == "starmaxhamming") {
                        if (!value.empty()) {
                            entry.starMaxHamming = std::atoi(value.c_str());
                        }
                    }
                }
                if (!entry.fastqs.empty()) {
                    config.libraries.push_back(entry);
                }
            }
        }
    } else if (currentSection == "feature") {
        for (const auto& kvLine : currentLines) {
            size_t commaPos = kvLine.find(',');
            if (commaPos == string::npos) continue;
            string key = kvLine.substr(0, commaPos);
            string value = kvLine.substr(commaPos + 1);
            trimInPlace(key);
            trimInPlace(value);
            
            string keyLower = key;
            std::transform(keyLower.begin(), keyLower.end(), keyLower.begin(), ::tolower);
            if (keyLower == "ref" || keyLower == "reference") {
                config.featureRef = value;
            }
        }
    } else if (currentSection == "reference") {
        for (const auto& kvLine : currentLines) {
            size_t commaPos = kvLine.find(',');
            if (commaPos == string::npos) continue;
            string key = kvLine.substr(0, commaPos);
            string value = kvLine.substr(commaPos + 1);
            trimInPlace(key);
            trimInPlace(value);
            
            string keyLower = key;
            std::transform(keyLower.begin(), keyLower.end(), keyLower.begin(), ::tolower);
            if (keyLower == "path" || keyLower == "genome") {
                config.referencePath = value;
            }
        }
    } else if (currentSection == "samples" && inSamples) {
        parseSamplesRows(samplesHeader, currentLines, config.samples);
    }
    
    // Auto-generate star_library_id when absent; check for duplicates
    // on both raw IDs and filesystem-sanitized IDs (non-alnum chars → '_').
    {
        set<string> seenIds;
        set<string> seenSanitized;
        for (size_t i = 0; i < config.libraries.size(); ++i) {
            auto& lib = config.libraries[i];
            if (lib.starLibraryId.empty()) {
                string samplePart = lib.sample.empty() ? "lib" : lib.sample;
                string ftNorm = lib.feature_types;
                std::replace(ftNorm.begin(), ftNorm.end(), ' ', '_');
                lib.starLibraryId = samplePart + "_" + ftNorm + "_" + std::to_string(i);
            }
            if (!seenIds.insert(lib.starLibraryId).second) {
                throw runtime_error("Duplicate star_library_id '" + lib.starLibraryId
                    + "'; each library must have a unique identifier");
            }
            // Sanitize the same way as PfMultiProcess::sanitizeDirName
            string sanitized;
            sanitized.reserve(lib.starLibraryId.size());
            for (unsigned char c : lib.starLibraryId) {
                sanitized.push_back(
                    (std::isalnum(c) || c == '-' || c == '_') ? static_cast<char>(c) : '_');
            }
            if (!seenSanitized.insert(sanitized).second) {
                throw runtime_error("star_library_id '" + lib.starLibraryId
                    + "' collides with another ID after path sanitization (both become '"
                    + sanitized + "'); use IDs that differ in alphanumeric/dash/underscore characters");
            }
        }
    }

    // Resolve and validate star_feature_ref paths.
    // Relative paths are resolved against the config file's directory.
    string configDir;
    {
        size_t lastSlash = configPath.find_last_of("/\\");
        configDir = (lastSlash == string::npos) ? "." : configPath.substr(0, lastSlash);
    }
    for (auto& lib : config.libraries) {
        if (!lib.starWhitelist.empty()) {
            if (lib.starWhitelist[0] != '/') {
                lib.starWhitelist = configDir + "/" + lib.starWhitelist;
            }
            struct stat st;
            if (stat(lib.starWhitelist.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) {
                throw runtime_error("star_whitelist path does not exist or is not a file: "
                    + lib.starWhitelist);
            }
        }
        if (!lib.starFeatureRef.empty()) {
            // Resolve relative paths against config directory
            if (lib.starFeatureRef[0] != '/') {
                lib.starFeatureRef = configDir + "/" + lib.starFeatureRef;
            }
            struct stat st;
            if (stat(lib.starFeatureRef.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) {
                throw runtime_error("star_feature_ref path does not exist or is not a file: "
                    + lib.starFeatureRef);
            }
        }
    }

    return config;
}

string resolveFastqDir(const string& configPath, const string& fastqRoot, 
                       const map<string, string>& fastqMap) {
    // Check fastq_map first
    auto it = fastqMap.find(configPath);
    if (it != fastqMap.end()) {
        return it->second;
    }
    
    // Check if path exists as-is
    struct stat st;
    if (stat(configPath.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
        return configPath;
    }
    
    // Try fastq_root + basename
    if (!fastqRoot.empty()) {
        size_t lastSlash = configPath.find_last_of("/\\");
        string basename = (lastSlash == string::npos) ? configPath : configPath.substr(lastSlash + 1);
        // Remove trailing slash if present
        while (!basename.empty() && (basename.back() == '/' || basename.back() == '\\')) {
            basename.pop_back();
        }
        string candidate = fastqRoot;
        if (candidate.back() != '/' && candidate.back() != '\\') {
            candidate += "/";
        }
        candidate += basename;
        if (stat(candidate.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
            return candidate;
        }
    }
    
    // Return original path (will fail later if invalid)
    return configPath;
}

map<string, string> parseFastqMap(const vector<string>& fastqMapVec) {
    map<string, string> result;
    for (const auto& entry : fastqMapVec) {
        size_t eqPos = entry.find('=');
        if (eqPos == string::npos) {
            continue; // Skip invalid entries
        }
        string key = entry.substr(0, eqPos);
        string value = entry.substr(eqPos + 1);
        // Trim whitespace
        trimInPlace(key);
        trimInPlace(value);
        result[key] = value;
    }
    return result;
}

} // namespace PfMultiConfig
