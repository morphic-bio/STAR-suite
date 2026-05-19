#include "OcmMultiConfig.h"
#include <cctype>
#include <set>
#include <sstream>
#include <stdexcept>

namespace OcmMultiConfig {

bool isValidOcmId(const string& id) {
    for (const char* valid : kValidOcmIds) {
        if (id == valid) {
            return true;
        }
    }
    return false;
}

bool isPathSafeSampleId(const string& sample_id) {
    if (sample_id.empty() || sample_id.size() > 128) {
        return false;
    }
    if (sample_id.find('/') != string::npos || sample_id.find('\\') != string::npos) {
        return false;
    }
    if (sample_id.find("..") != string::npos) {
        return false;
    }
    if (sample_id[0] == '.' || sample_id[0] == '-') {
        return false;
    }
    for (char c : sample_id) {
        const unsigned char uc = static_cast<unsigned char>(c);
        if (std::isalnum(uc) || c == '-' || c == '_' || c == '%') {
            continue;
        }
        return false;
    }
    return true;
}

PfMultiConfig::Config parseAndValidate(const string& configPath, std::ostream& logStream) {
    PfMultiConfig::Config config = PfMultiConfig::parseConfig(configPath);
    if (config.samples.empty()) {
        throw std::runtime_error("OCM multi config has no [samples] entries: " + configPath);
    }

    std::set<string> seenSampleIds;
    std::map<string, vector<string>> tagToSamples;

    for (const auto& sample : config.samples) {
        if (sample.sample_id.empty()) {
            throw std::runtime_error("OCM sample entry missing sample_id in " + configPath);
        }
        if (!isPathSafeSampleId(sample.sample_id)) {
            throw std::runtime_error("OCM sample_id is not path-safe '" + sample.sample_id + "' in " + configPath);
        }
        if (!seenSampleIds.insert(sample.sample_id).second) {
            throw std::runtime_error("Duplicate OCM sample_id '" + sample.sample_id + "' in " + configPath);
        }

        vector<string> ocmIds = sample.resolvedOcmIds();
        if (ocmIds.empty()) {
            throw std::runtime_error("OCM sample '" + sample.sample_id
                                   + "' has empty ocm_barcode_ids in " + configPath);
        }
        for (const auto& ocmId : ocmIds) {
            if (!isValidOcmId(ocmId)) {
                ostringstream err;
                err << "Invalid ocm_barcode_ids token '" << ocmId << "' for sample '"
                    << sample.sample_id << "'; expected one of OB1, OB2, OB3, OB4";
                throw std::runtime_error(err.str());
            }
            tagToSamples[ocmId].push_back(sample.sample_id);
        }
    }

    for (const auto& kv : tagToSamples) {
        if (kv.second.size() > 1) {
            ostringstream msg;
            msg << "NOTICE: OCM tag " << kv.first << " assigned to multiple samples: ";
            for (size_t i = 0; i < kv.second.size(); ++i) {
                if (i > 0) {
                    msg << ", ";
                }
                msg << kv.second[i];
            }
            msg << " (union sample design)\n";
            logStream << msg.str();
        }
    }

    return config;
}

} // namespace OcmMultiConfig
