#ifndef OCM_MULTI_CONFIG_H
#define OCM_MULTI_CONFIG_H

#include "PfMultiConfig.h"
#include <ostream>

namespace OcmMultiConfig {

static const char* kValidOcmIds[] = {"OB1", "OB2", "OB3", "OB4"};

bool isValidOcmId(const string& id);

/** Path-safe sample_id for outs/ and samples/ directory names. */
bool isPathSafeSampleId(const string& sample_id);

/** Parse and validate [samples] from a Cell Ranger multi config. */
PfMultiConfig::Config parseAndValidate(const string& configPath, std::ostream& logStream);

} // namespace OcmMultiConfig

#endif // OCM_MULTI_CONFIG_H
