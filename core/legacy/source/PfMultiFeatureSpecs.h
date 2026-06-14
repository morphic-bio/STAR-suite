#ifndef PF_MULTI_FEATURE_SPECS_H
#define PF_MULTI_FEATURE_SPECS_H

#include "PfMultiConfig.h"
#include <map>
#include <set>
#include <sstream>

namespace PfMultiFeatureSpecs {

struct FeatureSpec {
    string libraryType;
    string featureRefType;
};

// Keys are normalized: lowercase, alphanumeric only (matching normalizedFeatureType output)
inline string normalizedFeatureRole(const string& value) {
    string normalized;
    normalized.reserve(value.size());
    for (unsigned char c : value) {
        if (std::isalnum(c)) {
            normalized.push_back(static_cast<char>(std::tolower(c)));
        }
    }
    return normalized;
}

inline bool isProteinAdtFeatureRole(const string& normalizedRole) {
    return normalizedRole == "antibodycapture"
        || normalizedRole == "adt"
        || normalizedRole == "protein";
}

inline bool shouldEmitAdtMexOutput(const string& featureRefType,
                                   const string& libraryType) {
    return isProteinAdtFeatureRole(normalizedFeatureRole(featureRefType))
        || isProteinAdtFeatureRole(normalizedFeatureRole(libraryType));
}

inline map<string, string> knownFeatureRefTypeMap() {
    return {
        {"crisprguidecapture", "CRISPR Guide Capture"},
        {"crisprguide", "CRISPR Guide Capture"},
        {"antibodycapture", "Antibody Capture"},
        {"adt", "Antibody Capture"},
        {"protein", "Antibody Capture"},
        {"cellplexcmo", "Multiplexing Capture"},
        {"multiplexingcapture", "Multiplexing Capture"}
    };
}

inline vector<FeatureSpec> buildFeatureSpecsFromConfig(
    const PfMultiConfig::Config& config, ostringstream& log)
{
    const map<string, string> knownMap = knownFeatureRefTypeMap();
    vector<FeatureSpec> specs;
    set<string> seen;

    for (const auto& lib : config.libraries) {
        string norm = lib.normalizedFeatureType();
        if (norm == "geneexpression" || norm == "gex") {
            continue;
        }
        if (norm.empty()) {
            throw runtime_error("Library row has empty feature_types"
                " (sample='" + lib.sample + "'"
                ", starLibraryId='" + lib.starLibraryId + "'"
                "); every non-GEX library must specify a feature_types value");
        }

        if (seen.count(norm)) {
            continue;
        }
        seen.insert(norm);

        FeatureSpec spec;
        spec.libraryType = lib.feature_types;

        auto it = knownMap.find(norm);
        if (it != knownMap.end()) {
            spec.featureRefType = it->second;
        } else {
            spec.featureRefType = lib.feature_types;
            log << "NOTICE: feature_types '" << lib.feature_types
                << "' is not a known 10x type; using verbatim as featureRefType\n";
        }
        specs.push_back(spec);
    }
    return specs;
}

} // namespace PfMultiFeatureSpecs

#endif // PF_MULTI_FEATURE_SPECS_H
