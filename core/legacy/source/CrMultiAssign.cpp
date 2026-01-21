#include "CrMultiAssign.h"
#include "ErrorWarning.h"
#include "serviceFuns.cpp"
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <sys/stat.h>
#include <stdexcept>
using std::cerr;
using std::endl;

namespace CrMultiAssign {

int runAssignBarcodes(const string& assignBin, const string& whitelist,
                     const string& featureRef, const string& fastqDir,
                     const string& assignOut) {
    // Check inputs exist
    struct stat st;
    if (stat(assignBin.c_str(), &st) != 0) {
        ostringstream err;
        err << "assignBarcodes binary not found: " << assignBin;
        throw runtime_error(err.str());
    }
    if (stat(whitelist.c_str(), &st) != 0) {
        ostringstream err;
        err << "Whitelist file not found: " << whitelist;
        throw runtime_error(err.str());
    }
    if (stat(featureRef.c_str(), &st) != 0) {
        ostringstream err;
        err << "Feature reference file not found: " << featureRef;
        throw runtime_error(err.str());
    }
    if (stat(fastqDir.c_str(), &st) != 0 || !S_ISDIR(st.st_mode)) {
        ostringstream err;
        err << "FASTQ directory not found: " << fastqDir;
        throw runtime_error(err.str());
    }
    
    // Create output directory
    string cmd = "mkdir -p \"" + assignOut + "\"";
    int ret = system(cmd.c_str());
    if (ret != 0) {
        ostringstream err;
        err << "Failed to create output directory: " << assignOut;
        throw runtime_error(err.str());
    }
    
    // Build command
    ostringstream cmdStream;
    cmdStream << "\"" << assignBin << "\""
              << " --whitelist \"" << whitelist << "\""
              << " --featurelist \"" << featureRef << "\""
              << " --directory \"" << assignOut << "\""
              << " \"" << fastqDir << "\"";
    
    string fullCmd = cmdStream.str();
    
    // Execute
    ret = system(fullCmd.c_str());
    if (ret != 0) {
        ostringstream err;
        err << "assignBarcodes failed with exit code " << ret;
        throw runtime_error(err.str());
    }
    
    return 0;
}

int processFeatureLibraries(const CrMultiConfig::Config& config,
                          const string& assignBin, const string& whitelist,
                          const string& featureRef, const map<string, string>& fastqMap,
                          const string& fastqRoot, const string& outPrefix,
                          const vector<string>& featureTypes) {
    for (const auto& featureType : featureTypes) {
        vector<CrMultiConfig::LibraryEntry> libs = config.getFeatureLibraries(featureType);
        if (libs.empty()) {
            continue; // Skip if no libraries of this type
        }
        
        // Process first library of this type (could extend to handle multiple)
        const auto& lib = libs[0];
        string resolvedFastq = CrMultiConfig::resolveFastqDir(lib.fastqs, fastqRoot, fastqMap);
        
        // Create output directory name from feature type
        string featureTypeDir = featureType;
        // Replace spaces and special chars with underscores
        for (char& c : featureTypeDir) {
            if (c == ' ' || c == '/' || c == '\\') {
                c = '_';
            }
        }
        string assignOut = outPrefix + "/cr_assign/" + featureTypeDir;
        
        try {
            runAssignBarcodes(assignBin, whitelist, featureRef, resolvedFastq, assignOut);
        } catch (const exception& e) {
            cerr << "ERROR processing " << featureType << ": " << e.what() << endl;
            return 1;
        }
    }
    
    return 0;
}

} // namespace CrMultiAssign
