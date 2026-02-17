#include "PfMultiProcess.h"
#include "Solo.h"
#include "PfMultiConfig.h"
#include "PfMultiAssign.h"
#include "PfMultiMexStub.h"
#include "PfMultiMerge.h"
#include "ErrorWarning.h"
#include "serviceFuns.cpp"
#include "TimeFunctions.h"
#include "call_features.h"
#include "MexWriter.h"
#include <sys/stat.h>
#include <dirent.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <cctype>
using std::cerr;
using std::endl;

namespace {

struct FeatureSpec {
    string libraryType;
    string featureRefType;
};

static string normalizeType(const string& input) {
    string out;
    out.reserve(input.size());
    for (unsigned char c : input) {
        if (std::isalnum(c)) {
            out.push_back(static_cast<char>(std::tolower(c)));
        }
    }
    return out;
}

static string trimCopy(const string& input) {
    size_t first = input.find_first_not_of(" \t\r\n");
    if (first == string::npos) {
        return "";
    }
    size_t last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}

static string lowerCopy(string input) {
    std::transform(input.begin(), input.end(), input.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return input;
}

static bool isUnsetToken(const string& input) {
    string token = lowerCopy(trimCopy(input));
    return token.empty() || token == "-" || token == "none";
}

static string basenameOf(const string& path);

static string upperCopy(string input) {
    std::transform(input.begin(), input.end(), input.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return input;
}

static bool isValidBarcodeSeq(const string& seq) {
    if (seq.empty()) {
        return false;
    }
    for (unsigned char c : seq) {
        unsigned char u = static_cast<unsigned char>(std::toupper(c));
        if (!(u == 'A' || u == 'C' || u == 'G' || u == 'T' || u == 'N')) {
            return false;
        }
    }
    return true;
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

static string translateNxtMiddleTwoBases(string barcode) {
    barcode = upperCopy(barcode);
    if (barcode.size() >= 9) {
        barcode[7] = complementBase(barcode[7]);
        barcode[8] = complementBase(barcode[8]);
    }
    return barcode;
}

static vector<string> splitWhitelistColumns(const string& line) {
    vector<string> fields;
    string current;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (c == '\t' || c == ',' || std::isspace(static_cast<unsigned char>(c))) {
            if (!current.empty()) {
                fields.push_back(current);
                current.clear();
                if (fields.size() >= 2) {
                    break;
                }
            }
        } else {
            current.push_back(c);
        }
    }
    if (!current.empty() && fields.size() < 2) {
        fields.push_back(current);
    }
    return fields;
}

static string detectChemistryFromWhitelistPath(const string& whitelistPath, string& reason) {
    if (isUnsetToken(whitelistPath)) {
        reason = "whitelist path is unset";
        return "unknown";
    }

    // Prefer content-based detection first for 2-column NXT/TRU translation lists.
    std::ifstream in(whitelistPath.c_str());
    if (in.is_open()) {
        const size_t kMaxScanRows = 10000;
        size_t scanned = 0;
        size_t twoColumnRows = 0;
        size_t nxtRuleRows = 0;
        string line;

        while (scanned < kMaxScanRows && std::getline(in, line)) {
            string trimmed = trimCopy(line);
            if (trimmed.empty() || trimmed[0] == '#') {
                continue;
            }

            vector<string> fields = splitWhitelistColumns(trimmed);
            if (fields.empty()) {
                continue;
            }
            scanned++;

            if (fields.size() >= 2 && isValidBarcodeSeq(fields[0]) && isValidBarcodeSeq(fields[1])) {
                twoColumnRows++;
                const string col1 = upperCopy(fields[0]);
                const string col2 = upperCopy(fields[1]);
                if (translateNxtMiddleTwoBases(col1) == col2) {
                    nxtRuleRows++;
                }
            }
        }

        if (twoColumnRows > 0) {
            std::ostringstream msg;
            msg << "whitelist has two barcode columns (" << twoColumnRows
                << " sampled rows)";
            if (nxtRuleRows > 0) {
                msg << "; center-2bp complement rule matched " << nxtRuleRows << "/" << twoColumnRows;
            }
            reason = msg.str();
            return "NXT";
        }
    }

    // Fallback to filename hints.
    string base = lowerCopy(basenameOf(whitelistPath));
    if (base.find("nxt") != string::npos) {
        reason = "whitelist filename contains 'nxt'";
        return "NXT";
    }
    if (base.find("tru") != string::npos) {
        reason = "whitelist filename contains 'tru'";
        return "TRU";
    }

    reason = "no chemistry marker found in whitelist content or filename; defaulting to TRU";
    return "TRU";
}

static string parseChemistryToken(const string& rawValue, const string& flagName) {
    string value = lowerCopy(trimCopy(rawValue));
    if (value.empty() || value == "-" || value == "none") {
        value = "auto";
    }
    if (value != "auto" && value != "nxt" && value != "tru") {
        throw runtime_error("Invalid " + flagName + " value '" + rawValue +
                            "' (allowed: auto|NXT|TRU)");
    }
    return value;
}

static string sanitizeDirName(const string& input) {
    string out = input;
    for (char& c : out) {
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '-' || c == '_')) {
            c = '_';
        }
    }
    return out;
}

static string basenameOf(const string& path) {
    size_t pos = path.find_last_of("/\\");
    return (pos == string::npos) ? path : path.substr(pos + 1);
}

static string findAssignOutputDir(const string& baseDir) {
    struct stat st;
    string matrixPath = baseDir + "/matrix.mtx";
    string matrixGz = baseDir + "/matrix.mtx.gz";
    if (stat(matrixPath.c_str(), &st) == 0 || stat(matrixGz.c_str(), &st) == 0) {
        return baseDir;
    }

    DIR* dir = opendir(baseDir.c_str());
    if (!dir) {
        return baseDir;
    }
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        if (entry->d_name[0] == '.') {
            continue;
        }
        string sub = baseDir + "/" + entry->d_name;
        if (stat(sub.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
            string subMatrix = sub + "/matrix.mtx";
            string subMatrixGz = sub + "/matrix.mtx.gz";
            if (stat(subMatrix.c_str(), &st) == 0 || stat(subMatrixGz.c_str(), &st) == 0) {
                closedir(dir);
                return sub;
            }
        }
    }
    closedir(dir);
    return baseDir;
}

static bool hasMexFiles(const string& dirPath) {
    struct stat st;
    string features = dirPath + "/features.tsv";
    string featuresGz = dirPath + "/features.tsv.gz";
    return (stat(features.c_str(), &st) == 0) || (stat(featuresGz.c_str(), &st) == 0);
}

static bool filterFeatureRefCsv(const string& inputPath, const string& featureType,
                                const string& outputPath) {
    ifstream in(inputPath);
    if (!in.is_open()) {
        return false;
    }

    string headerLine;
    if (!getline(in, headerLine)) {
        return false;
    }

    vector<string> headers;
    {
        string field;
        bool inQuotes = false;
        for (char c : headerLine) {
            if (c == '"') {
                inQuotes = !inQuotes;
            } else if (c == ',' && !inQuotes) {
                headers.push_back(field);
                field.clear();
            } else {
                field += c;
            }
        }
        headers.push_back(field);
    }

    int typeIdx = -1;
    for (size_t i = 0; i < headers.size(); ++i) {
        string h = headers[i];
        std::transform(h.begin(), h.end(), h.begin(), ::tolower);
        if (h == "feature_type" || h == "type") {
            typeIdx = static_cast<int>(i);
            break;
        }
    }

    if (typeIdx < 0) {
        return false;
    }

    string normalizedTarget = normalizeType(featureType);
    ofstream out(outputPath);
    if (!out.is_open()) {
        return false;
    }
    out << headerLine << "\n";

    string line;
    bool wroteAny = false;
    while (getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        vector<string> fields;
        string field;
        bool inQuotes = false;
        for (char c : line) {
            if (c == '"') {
                inQuotes = !inQuotes;
            } else if (c == ',' && !inQuotes) {
                fields.push_back(field);
                field.clear();
            } else {
                field += c;
            }
        }
        fields.push_back(field);
        if (static_cast<int>(fields.size()) <= typeIdx) {
            continue;
        }
        string value = fields[typeIdx];
        if (!value.empty() && value.front() == '"' && value.back() == '"') {
            value = value.substr(1, value.size() - 2);
        }
        string normalizedValue = normalizeType(value);
        if (normalizedValue == normalizedTarget) {
            out << line << "\n";
            wroteAny = true;
        }
    }
    return wroteAny;
}

static bool getFilteredBarcodesFromSolo(const Solo* solo, const Parameters& P, vector<string>& out) {
    if (!solo || !solo->soloFeat) {
        return false;
    }
    int32 featIdx = -1;
    if (P.pSolo.crGexFeature == ParametersSolo::CrGexGene) {
        featIdx = P.pSolo.featureInd[SoloFeatureTypes::Gene];
    } else if (P.pSolo.crGexFeature == ParametersSolo::CrGexGeneFull) {
        featIdx = P.pSolo.featureInd[SoloFeatureTypes::GeneFull];
    } else {
        // Auto: prefer Gene, fallback to GeneFull
        featIdx = P.pSolo.featureInd[SoloFeatureTypes::Gene];
        if (featIdx < 0) {
            featIdx = P.pSolo.featureInd[SoloFeatureTypes::GeneFull];
        }
    }
    if (featIdx < 0 || solo->soloFeat[featIdx] == nullptr) {
        return false;
    }

    const SoloFeature* gex = solo->soloFeat[featIdx];
    if (gex->filteredCells.filtVecBool.empty()) {
        return false;
    }

    out.clear();
    out.reserve(gex->filteredCells.nCells);
    for (uint32 icb = 0; icb < gex->nCB; icb++) {
        if (gex->filteredCells.filtVecBool[icb]) {
            uint32 wlIdx = gex->indCB[icb];
            if (wlIdx < gex->pSolo.cbWLstr.size()) {
                out.push_back(gex->pSolo.cbWLstr[wlIdx]);
            }
        }
    }
    return !out.empty();
}

/**
 * Run CRISPR feature calling on filtered MEX with CRISPR Guide Capture features.
 * 
 * @param filteredMexDir Directory containing filtered_feature_bc_matrix
 * @param outputDir Output directory for crispr_analysis/
 * @param minUmi Minimum UMI threshold for GMM calling
 * @param logStream Log output stream
 * @return 0 on success, -1 on failure
 */
static int runCrisprFeatureCalling(const string& filteredMexDir, const string& outputDir,
                                    int minUmi, ostream& logStream) {
    logStream << timeMonthDayTime() << " ..... starting CRISPR feature calling\n";
    
    // Step 1: Read the filtered MEX
    PfMultiMerge::MexData mexData;
    try {
        mexData = PfMultiMerge::readMex(filteredMexDir);
    } catch (const exception& e) {
        logStream << "ERROR: Failed to read filtered MEX for CRISPR calling: " << e.what() << "\n";
        return -1;
    }
    
    // Step 2: Filter to CRISPR Guide Capture features only
    PfMultiMerge::MexData crisprData = PfMultiMerge::filterByFeatureType(mexData, "CRISPR Guide Capture");
    
    if (crisprData.features.empty()) {
        logStream << "NOTICE: No CRISPR Guide Capture features found, skipping feature calling\n";
        return 0;
    }
    
    logStream << "  CRISPR features found: " << crisprData.features.size() << "\n";
    logStream << "  Barcodes: " << crisprData.barcodes.size() << "\n";
    logStream << "  Non-zero entries: " << crisprData.triplets.size() << "\n";
    
    // Step 3: Write CRISPR-only MEX to temporary directory
    string tempMexDir = outputDir + "/.crispr_mex_tmp";
    string mkdirCmd = "mkdir -p \"" + tempMexDir + "\"";
    if (system(mkdirCmd.c_str()) != 0) {
        logStream << "ERROR: Failed to create temp directory: " << tempMexDir << "\n";
        return -1;
    }
    
    // Convert MexData to MexWriter format
    vector<MexWriter::Feature> features;
    for (size_t i = 0; i < crisprData.features.size(); ++i) {
        string name = (i < crisprData.featureNames.size()) ? crisprData.featureNames[i] : crisprData.features[i];
        string type = (i < crisprData.featureTypes.size()) ? crisprData.featureTypes[i] : "CRISPR Guide Capture";
        features.emplace_back(crisprData.features[i], name, type);
    }
    
    // Write MEX (uncompressed for call_features compatibility)
    string mexPrefix = tempMexDir + "/";
    int ret = MexWriter::writeMex(mexPrefix, crisprData.barcodes, features, crisprData.triplets, -1);
    if (ret != 0) {
        logStream << "ERROR: Failed to write temporary CRISPR MEX\n";
        return -1;
    }
    
    // Step 4: Run GMM feature calling with min_umi=10 (CR-compatible default)
    string crisprAnalysisDir = outputDir + "/crispr_analysis";
    mkdirCmd = "mkdir -p \"" + crisprAnalysisDir + "\"";
    if (system(mkdirCmd.c_str()) != 0) {
        logStream << "ERROR: Failed to create crispr_analysis directory\n";
        string rmCmd = "rm -rf \"" + tempMexDir + "\"";
        system(rmCmd.c_str());
        return -1;
    }
    
    cf_gmm_config *gmm_cfg = cf_gmm_config_create();
    if (!gmm_cfg) {
        logStream << "ERROR: Failed to create GMM config\n";
        return -1;
    }
    gmm_cfg->min_umi_threshold = minUmi;
    
    logStream << "  Calling mode: GMM (CR9-compatible)\n";
    logStream << "  min_umi: " << gmm_cfg->min_umi_threshold << "\n";
    logStream << "  Output: " << crisprAnalysisDir << "\n";
    
    ret = cf_process_mex_dir_gmm(tempMexDir.c_str(), crisprAnalysisDir.c_str(), gmm_cfg);
    cf_gmm_config_destroy(gmm_cfg);
    
    if (ret != 0) {
        logStream << "ERROR: CRISPR feature calling failed\n";
        // Cleanup temp dir
        string rmCmd = "rm -rf \"" + tempMexDir + "\"";
        system(rmCmd.c_str());
        return -1;
    }
    
    // Step 5: Cleanup temporary MEX directory
    string rmCmd = "rm -rf \"" + tempMexDir + "\"";
    system(rmCmd.c_str());
    
    logStream << timeMonthDayTime() << " ..... finished CRISPR feature calling\n";
    logStream << "  Output files:\n";
    logStream << "    " << crisprAnalysisDir << "/protospacer_calls_per_cell.csv\n";
    logStream << "    " << crisprAnalysisDir << "/protospacer_calls_summary.csv\n";
    logStream << "    " << crisprAnalysisDir << "/protospacer_umi_thresholds.csv\n";
    logStream << "    " << crisprAnalysisDir << "/protospacer_umi_thresholds.json\n";
    
    return 0;
}

} // namespace

int processPfMultiConfig(Parameters& P, const Solo* solo) {
    if (isUnsetToken(P.pfMulti.pfMultiConfig)) {
        return 0; // Not enabled
    }
    
    P.inOut->logMain << timeMonthDayTime() << " ..... started pf-multi config processing\n";
    
    try {
        // Parse multi config
        PfMultiConfig::Config config = PfMultiConfig::parseConfig(P.pfMulti.pfMultiConfig);
        
        if (config.libraries.empty()) {
            ostringstream err;
            err << "No libraries found in multi config";
            throw runtime_error(err.str());
        }
        
        // Get feature reference
        string featureRef = P.pfMulti.crFeatureRef;
        if (isUnsetToken(featureRef)) {
            featureRef = config.featureRef;
        }
        if (isUnsetToken(featureRef)) {
            throw runtime_error("Feature reference not provided (use --crFeatureRef or set in config)");
        }
        
        // Get whitelist
        string whitelist = P.pfMulti.crWhitelist;
        if (isUnsetToken(whitelist)) {
            // Try to get from solo whitelist
            if (!P.pSolo.soloCBwhitelist.empty() && !isUnsetToken(P.pSolo.soloCBwhitelist[0])) {
                whitelist = P.pSolo.soloCBwhitelist[0];
            }
        }
        if (isUnsetToken(whitelist)) {
            throw runtime_error("Whitelist not provided (use --crWhitelist or --soloCBwhitelist)");
        }

        // Resolve chemistry precedence: default auto-detect, explicit chemistry override.
        const string requestedChem = parseChemistryToken(P.pfMulti.crChemistry, "--crChemistry");
        const string requestedOutputChem = parseChemistryToken(P.pfMulti.crOutputChemistry, "--crOutputChemistry");

        string inferredReason;
        string inferredChem = detectChemistryFromWhitelistPath(whitelist, inferredReason);
        string effectiveChem = inferredChem;
        if (requestedChem == "nxt") {
            effectiveChem = "NXT";
        } else if (requestedChem == "tru") {
            effectiveChem = "TRU";
        }

        string outputChem = effectiveChem;
        if (requestedOutputChem == "nxt") {
            outputChem = "NXT";
        } else if (requestedOutputChem == "tru") {
            outputChem = "TRU";
        }

        P.inOut->logMain << "pf-multi chemistry: requested=" << requestedChem
                         << " inferred=" << inferredChem
                         << " (" << inferredReason << ")"
                         << " effective=" << effectiveChem
                         << " output_requested=" << requestedOutputChem
                         << " output_effective=" << outputChem << "\n";
        
        // Parse FASTQ map
        map<string, string> fastqMap = PfMultiConfig::parseFastqMap(P.pfMulti.crFastqMap);
        
        // Get output prefix
        string outPrefix = P.outFileNamePrefix;
        // Remove trailing slash if present
        while (!outPrefix.empty() && outPrefix.back() == '/') {
            outPrefix.pop_back();
        }

        string crAssignRoot = outPrefix + "/cr_assign";
        string mkCrAssignRootCmd = "mkdir -p \"" + crAssignRoot + "\"";
        if (system(mkCrAssignRootCmd.c_str()) != 0) {
            throw runtime_error("Failed to create CR assign root directory: " + crAssignRoot);
        }

        PfMultiAssign::AssignOptions assignOpts;
        assignOpts.maxHammingDistance = P.pfMulti.crAssignMaxHamming;
        assignOpts.featureConstantOffset = P.pfMulti.crAssignFeatureOffset;
        assignOpts.limitSearch = P.pfMulti.crAssignLimitSearch;
        assignOpts.minCounts = P.pfMulti.crAssignMinCounts;
        assignOpts.maxBarcodeMismatches = P.pfMulti.crAssignMaxBarcodeMismatches;
        assignOpts.featureN = P.pfMulti.crAssignFeatureN;
        assignOpts.barcodeN = P.pfMulti.crAssignBarcodeN;
        assignOpts.consumerThreadsPerSet = P.pfMulti.crAssignConsumerThreads;
        assignOpts.searchThreads = P.pfMulti.crAssignSearchThreads;
        assignOpts.minPosterior = P.pfMulti.crAssignMinPosterior;

        if (!P.pfMulti.crAssignFilteredBarcodes.empty() && P.pfMulti.crAssignFilteredBarcodes != "-") {
            struct stat stFiltered;
            if (stat(P.pfMulti.crAssignFilteredBarcodes.c_str(), &stFiltered) != 0) {
                throw runtime_error("crAssignFilteredBarcodes file not found: " + P.pfMulti.crAssignFilteredBarcodes);
            }
            assignOpts.filteredBarcodesPath = P.pfMulti.crAssignFilteredBarcodes;
            P.inOut->logMain << "NOTICE: Using explicit filtered barcodes for assignBarcodes: "
                             << assignOpts.filteredBarcodesPath << "\n";
        } else {
            vector<string> filteredAssignBarcodes;
            if (getFilteredBarcodesFromSolo(solo, P, filteredAssignBarcodes)) {
                string filteredPath = crAssignRoot + "/filtered_gex_barcodes.txt";
                std::ofstream fbc(filteredPath.c_str());
                if (!fbc.is_open()) {
                    throw runtime_error("Failed to write filtered barcode list for assignBarcodes: " + filteredPath);
                }
                for (const auto& bc : filteredAssignBarcodes) {
                    fbc << bc << "\n";
                }
                fbc.close();
                assignOpts.filteredBarcodesPath = filteredPath;
                P.inOut->logMain << "NOTICE: Using Solo filtered barcodes for assignBarcodes (" 
                                 << filteredAssignBarcodes.size() << " barcodes)\n";
            }
        }
        
        vector<FeatureSpec> featureSpecs = {
            {"CRISPR Guide Capture", "CRISPR Guide Capture"},
            {"Antibody Capture", "Antibody Capture"},
            {"CellPlex (CMO)", "Multiplexing Capture"},
            {"Multiplexing Capture", "Multiplexing Capture"}
        };

        struct FeatureRun {
            string featureType;
            string assignOut;
            string featureRefPath;
        };
        vector<FeatureRun> featureRuns;

        for (const auto& spec : featureSpecs) {
            vector<PfMultiConfig::LibraryEntry> libs = config.getFeatureLibraries(spec.libraryType);
            if (libs.empty()) {
                continue;
            }
            string featureDir = sanitizeDirName(spec.libraryType);
            string assignBase = outPrefix + "/cr_assign/" + featureDir;

            for (const auto& lib : libs) {
                string resolvedFastq = PfMultiConfig::resolveFastqDir(lib.fastqs, P.pfMulti.crFastqRoot, fastqMap);
                string sampleName = lib.sample.empty() ? basenameOf(resolvedFastq) : lib.sample;
                string assignOut = assignBase + "/" + sanitizeDirName(sampleName);

                string mkAssignOutCmd = "mkdir -p \"" + assignOut + "\"";
                if (system(mkAssignOutCmd.c_str()) != 0) {
                    throw runtime_error("Failed to create assign output directory: " + assignOut);
                }

                string filteredRef = assignOut + "/feature_reference.filtered.csv";
                bool filtered = filterFeatureRefCsv(featureRef, spec.featureRefType, filteredRef);
                string refPath = filtered ? filteredRef : featureRef;
                if (!filtered) {
                    P.inOut->logMain << "WARNING: feature reference not filtered for " << spec.libraryType
                                     << "; using full reference\n";
                }

                int ret = PfMultiAssign::runAssignBarcodes(whitelist, refPath, resolvedFastq, assignOut, assignOpts);
                if (ret != 0) {
                    throw runtime_error("Failed to process feature library: " + spec.libraryType);
                }

                FeatureRun run;
                run.featureType = spec.featureRefType;
                run.assignOut = assignOut;
                run.featureRefPath = refPath;
                featureRuns.push_back(run);
            }
        }

        if (featureRuns.empty()) {
            throw runtime_error("No feature libraries found in multi config");
        }

        for (auto& run : featureRuns) {
            string realOut = findAssignOutputDir(run.assignOut);
            run.assignOut = realOut;
            int ret = PfMultiMexStub::processAssignOutput(run.assignOut, run.featureRefPath, run.featureType, false);
            if (ret != 0) {
                P.inOut->logMain << "WARNING: Failed to generate MEX stub for " << run.featureType << "\n";
            }
        }
        
        // Read GEX MEX from STARsolo output (both raw and filtered if available)
        string soloOut = outPrefix + "/Solo.out";
        string geneOut = soloOut + "/Gene";
        string geneFullOut = soloOut + "/GeneFull";
        string filteredOut = geneOut + "/filtered";
        string rawOut = geneOut + "/raw";
        string geneFullFiltered = geneFullOut + "/filtered";
        string geneFullRaw = geneFullOut + "/raw";

        if (P.pSolo.crGexFeature == ParametersSolo::CrGexGeneFull) {
            rawOut = geneFullRaw;
            filteredOut = geneFullFiltered;
            P.inOut->logMain << "NOTICE: --soloCrGexFeature=genefull (using GeneFull MEX for CR-compat merge)\n";
        } else if (P.pSolo.crGexFeature == ParametersSolo::CrGexGene) {
            rawOut = geneOut + "/raw";
            filteredOut = geneOut + "/filtered";
            P.inOut->logMain << "NOTICE: --soloCrGexFeature=gene (using Gene MEX for CR-compat merge)\n";
        }

        bool hasRaw = (P.pSolo.type != 0 && hasMexFiles(rawOut));
        bool hasFiltered = (P.pSolo.type != 0 && hasMexFiles(filteredOut));

        if (!hasRaw && !hasFiltered) {
            if (P.pSolo.crGexFeature == ParametersSolo::CrGexGeneFull) {
                P.inOut->logMain << "ERROR: GeneFull MEX directory not found for CR-compat merge\n";
                return 1;
            }
            if (P.pSolo.crGexFeature == ParametersSolo::CrGexGene) {
                P.inOut->logMain << "ERROR: Gene MEX directory not found for CR-compat merge\n";
                return 1;
            }
            // Fallback to GeneFull if Gene is absent (auto mode)
            bool hasGeneFullRaw = (P.pSolo.type != 0 && hasMexFiles(geneFullRaw));
            bool hasGeneFullFiltered = (P.pSolo.type != 0 && hasMexFiles(geneFullFiltered));
            if (hasGeneFullRaw || hasGeneFullFiltered) {
                rawOut = geneFullRaw;
                filteredOut = geneFullFiltered;
                hasRaw = hasGeneFullRaw;
                hasFiltered = hasGeneFullFiltered;
                P.inOut->logMain << "NOTICE: Using GeneFull MEX for CR-compat merge (Gene missing)\n";
            } else if (P.pSolo.type != 0 && hasMexFiles(geneOut)) {
                hasRaw = true;
                rawOut = geneOut;
            } else {
                P.inOut->logMain << "WARNING: GEX MEX directory not found, skipping merge\n";
                return 0;
            }
        }
        
        // Prefer in-memory filtered barcodes from Solo if available (avoids reading filtered MEX for barcode list)
        vector<string> filteredGexBarcodes;
        bool useFilteredGex = false;
        if (getFilteredBarcodesFromSolo(solo, P, filteredGexBarcodes)) {
            useFilteredGex = true;
            P.inOut->logMain << "NOTICE: Using GEX filtered barcodes from Solo (in-memory)\n";
        }

        // Read raw GEX MEX (required for raw_feature_bc_matrix)
        PfMultiMerge::MexData gexRawData;
        try {
            if (hasRaw) {
                gexRawData = PfMultiMerge::readMex(rawOut);
                // Filter to Gene Expression only if needed
                bool hasMultipleTypes = false;
                for (const auto& type : gexRawData.featureTypes) {
                    if (type != "Gene Expression" && !type.empty()) {
                        hasMultipleTypes = true;
                        break;
                    }
                }
                if (hasMultipleTypes) {
                    gexRawData = PfMultiMerge::filterByFeatureType(gexRawData, "Gene Expression");
                }
            }
        } catch (const exception& e) {
            P.inOut->logMain << "WARNING: Failed to read raw GEX MEX: " << e.what() << "\n";
            // If raw failed but filtered exists, we'll use filtered as fallback below
            if (!hasFiltered) {
                return 0;
            }
        }
        
        // Read filtered GEX MEX if needed (for counts fallback or barcode list)
        PfMultiMerge::MexData gexFilteredData;
        bool loadedFilteredMex = false;
        bool needFilteredMexForCounts = (!hasRaw || gexRawData.features.empty());
        bool needFilteredMexForBarcodes = (!useFilteredGex && hasFiltered);
        if (hasFiltered && (needFilteredMexForCounts || needFilteredMexForBarcodes)) {
            try {
                gexFilteredData = PfMultiMerge::readMex(filteredOut);
                // Filter to Gene Expression only if needed
                bool hasMultipleTypes = false;
                for (const auto& type : gexFilteredData.featureTypes) {
                    if (type != "Gene Expression" && !type.empty()) {
                        hasMultipleTypes = true;
                        break;
                    }
                }
                if (hasMultipleTypes) {
                    gexFilteredData = PfMultiMerge::filterByFeatureType(gexFilteredData, "Gene Expression");
                }
                loadedFilteredMex = true;
            } catch (const exception& e) {
                P.inOut->logMain << "WARNING: Failed to read filtered GEX MEX: " << e.what()
                                 << ", will use observed raw GEX barcodes for filtered output\n";
            }
        }
        
        // Determine primary GEX data for merging (use raw if available, otherwise filtered)
        PfMultiMerge::MexData gexData;
        if (hasRaw && !gexRawData.features.empty()) {
            gexData = gexRawData;
        } else if (hasFiltered && loadedFilteredMex && !gexFilteredData.features.empty()) {
            gexData = gexFilteredData;
            P.inOut->logMain << "WARNING: Raw GEX MEX not available, using filtered GEX MEX as merge base\n";
        } else {
            P.inOut->logMain << "ERROR: No valid GEX MEX data available for merging\n";
            return 1;
        }
        
        // Read feature MEX files
        vector<PfMultiMerge::MexData> featureDataVec;
        for (const auto& run : featureRuns) {
            try {
                PfMultiMerge::MexData featData = PfMultiMerge::readMex(run.assignOut);
                featureDataVec.push_back(featData);
            } catch (const exception& e) {
                P.inOut->logMain << "WARNING: Failed to read feature MEX for " << run.featureType
                                << ": " << e.what() << "\n";
            }
        }
        
        if (featureDataVec.empty()) {
            P.inOut->logMain << "WARNING: No feature MEX files found, skipping merge\n";
            return 0;
        }
        
        // Merge MEX files
        PfMultiMerge::MexData mergedData = PfMultiMerge::mergeMex(gexData, featureDataVec);
        
        // Extract GEM well from GEX library (with fallback logic)
        string gemWell = "1"; // Default
        vector<PfMultiConfig::LibraryEntry> gexLibs = config.getGexLibraries();
        if (!gexLibs.empty()) {
            gemWell = gexLibs[0].gem_well;
            // Check if multiple GEX entries disagree
            for (size_t i = 1; i < gexLibs.size(); ++i) {
                if (gexLibs[i].gem_well != gemWell) {
                    P.inOut->logMain << "WARNING: Multiple GEX libraries have different gem_well values ("
                                     << gemWell << " vs " << gexLibs[i].gem_well
                                     << "), using first: " << gemWell << "\n";
                    break;
                }
            }
        } else if (!config.libraries.empty()) {
            // Fallback to first library's gem_well
            gemWell = config.libraries[0].gem_well;
        }
        
        // Compute observed GEX barcodes from raw GEX triplets (for raw_feature_bc_matrix)
        vector<string> observedRawGexBarcodes;
        if (hasRaw && !gexRawData.features.empty()) {
            observedRawGexBarcodes = PfMultiMerge::computeObservedGexBarcodes(gexRawData);
        } else if (hasFiltered && !gexFilteredData.features.empty()) {
            // Fallback: use filtered GEX barcodes for raw output if raw is missing
            observedRawGexBarcodes = PfMultiMerge::computeObservedGexBarcodes(gexFilteredData);
            P.inOut->logMain << "WARNING: Raw GEX MEX not available, using filtered GEX barcodes for raw output\n";
        } else {
            P.inOut->logMain << "ERROR: Cannot compute observed GEX barcodes - no valid GEX data\n";
            return 1;
        }
        
        // Compute GEX barcodes for filtered output
        if (!useFilteredGex && hasFiltered && loadedFilteredMex && !gexFilteredData.features.empty()) {
            // Use filtered GEX barcodes (from filtered MEX)
            filteredGexBarcodes = PfMultiMerge::computeObservedGexBarcodes(gexFilteredData);
            useFilteredGex = true;
        }
        if (!useFilteredGex) {
            // Fallback to observed raw GEX barcodes (or filtered if raw missing)
            filteredGexBarcodes = observedRawGexBarcodes;
            P.inOut->logMain << "WARNING: Filtered GEX barcodes not available, using observed raw GEX barcodes for filtered output\n";
        }
        
        // Write raw_feature_bc_matrix (using observed raw GEX barcodes)
        string rawOutDir = outPrefix + "/outs/raw_feature_bc_matrix";
        int ret = PfMultiMerge::writeCombinedMex(rawOutDir,
                                                 mergedData,
                                                 gemWell,
                                                 P.inOut->logMain,
                                                 observedRawGexBarcodes,
                                                 effectiveChem,
                                                 outputChem);
        if (ret != 0) {
            throw runtime_error("Failed to write raw combined MEX");
        }
        P.inOut->logMain << "Raw MEX written to: " << rawOutDir << "\n";
        
        // Write filtered_feature_bc_matrix (using filtered GEX barcodes or fallback)
        string filteredOutDir = outPrefix + "/outs/filtered_feature_bc_matrix";
        ret = PfMultiMerge::writeCombinedMex(filteredOutDir,
                                             mergedData,
                                             gemWell,
                                             P.inOut->logMain,
                                             filteredGexBarcodes,
                                             effectiveChem,
                                             outputChem);
        if (ret != 0) {
            throw runtime_error("Failed to write filtered combined MEX");
        }
        P.inOut->logMain << "Filtered MEX written to: " << filteredOutDir << "\n";
        
        // Run CRISPR feature calling if CRISPR Guide Capture features were processed
        bool hasCrisprFeatures = false;
        for (const auto& run : featureRuns) {
            if (run.featureType == "CRISPR Guide Capture") {
                hasCrisprFeatures = true;
                break;
            }
        }
        
        if (hasCrisprFeatures) {
            string outsDir = outPrefix + "/outs";
            ret = runCrisprFeatureCalling(filteredOutDir, outsDir, P.pfMulti.crMinUmi, P.inOut->logMain);
            if (ret != 0) {
                P.inOut->logMain << "WARNING: CRISPR feature calling failed, continuing without crispr_analysis/\n";
            }
        }
        
        P.inOut->logMain << timeMonthDayTime() << " ..... finished pf-multi config processing\n";
        
    } catch (const exception& e) {
        P.inOut->logMain << "ERROR in pf-multi config processing: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
